import json
import subprocess
import shlex  # For safely formatting command arguments for display
import argparse
import os
import concurrent.futures  # For parallel processing


def execute_bcftools_query(vcf_file_path, region_string):
    """
    Executes a bcftools view query for a given region, excluding header lines.

    Args:
        vcf_file_path (str): Path to the VCF file.
        region_string (str): The region to query (e.g., "chr1:100-200").

    Returns:
        tuple: (list_of_variant_strings, bcftools_stderr_str)
               Returns (None, error_message) if bcftools command itself fails.
    """
    # Command to execute: bcftools view -H <file.vcf.gz> <chr:from-to>
    # The -H or --no-header option suppresses VCF header output
    command = ["bcftools", "view", "-H", vcf_file_path, region_string]
    # Worker PID is useful for debugging parallel processes
    print(f"  [Worker PID:{os.getpid()}] Executing: {' '.join(shlex.quote(arg) for arg in command)}")
    try:
        # Run the command
        process = subprocess.run(
            command,
            capture_output=True,  # Capture stdout and stderr
            text=True,  # Decode output as text
            check=False  # Do not raise an exception for non-zero exit codes automatically
            # We will check stderr manually
        )

        variants_list = []
        if process.stdout:
            # Split stdout into lines, removing any empty lines that might result from split
            # Since -H is used, no header lines should be present in stdout
            variants_list = [line for line in process.stdout.strip().split('\n') if line]

        return variants_list, process.stderr.strip()  # Return list of variants and any stderr output

    except FileNotFoundError:
        # This occurs if 'bcftools' executable is not found
        error_msg = f"Error: 'bcftools' command not found. Please ensure it is installed and in your PATH."
        print(error_msg)  # Print to console for immediate user feedback
        return None, error_msg  # Indicate critical failure
    except Exception as e:
        # Catch any other unexpected errors during subprocess execution
        error_msg = f"An unexpected error occurred while running bcftools for region {region_string}: {e}"
        print(error_msg)  # Print to console
        return None, error_msg


def process_item_for_bcftools(item, vcf_file_path):
    """
    Processes a single item from the input JSON, constructs the region,
    and prepares it for bcftools query.

    Args:
        item (dict): A dictionary from the input JSON (expected to have fields like
                     'top_level_node_id', 'haplotype_id', 'node_sequence_length',
                     'grch38_chromosome_region', 'grch38_start_position').
        vcf_file_path (str): Path to the VCF file.

    Returns:
        dict: A dictionary containing original item info, the query region string,
              a list of found variants, and any error messages from bcftools.
              Returns None if the item is fundamentally invalid for querying.
    """
    # Initialize the result entry with all expected fields from the input item
    base_result_entry = {
        "top_level_node_id": item.get("top_level_node_id"),
        "haplotype_id": item.get("haplotype_id"),
        "node_sequence_length": item.get("node_sequence_length"),
        "grch38_chromosome_region": item.get("grch38_chromosome_region"),
        "grch38_start_position": item.get("grch38_start_position"),
        "query_region": None,  # Will be populated if query proceeds
        "variants": [],  # Initialize as empty list for variants
        "bcftools_error": None  # Will store any errors
    }

    # Extract necessary fields from the input item
    start_pos_str = item.get("grch38_start_position")
    seq_len_val = item.get("node_sequence_length")
    chromosome = item.get("grch38_chromosome_region")
    item_hap_id_for_log = base_result_entry.get('haplotype_id', 'N/A')

    # Validate required fields for constructing the query region
    if not chromosome:
        base_result_entry["bcftools_error"] = "Skipped: Missing 'grch38_chromosome_region'."
        print(f"  Skipping item (HapID: {item_hap_id_for_log}): Missing 'grch38_chromosome_region'.")
        return base_result_entry  # Return the entry with error, it won't have variants.

    if start_pos_str is None:
        base_result_entry["bcftools_error"] = "Skipped: 'grch38_start_position' is null."
        print(f"  Skipping item (HapID: {item_hap_id_for_log}): 'grch38_start_position' is null.")
        return base_result_entry

    if seq_len_val is None:
        base_result_entry["bcftools_error"] = "Skipped: 'node_sequence_length' is null."
        print(f"  Skipping item (HapID: {item_hap_id_for_log}): 'node_sequence_length' is null.")
        return base_result_entry

    try:
        # Convert start position and sequence length to integers
        start_pos = int(start_pos_str)
        seq_len = int(seq_len_val)
        if seq_len < 0:
            base_result_entry[
                "bcftools_error"] = f"Skipped: Invalid 'node_sequence_length' ({seq_len}). Must be non-negative."
            print(f"  Skipping item (HapID: {item_hap_id_for_log}): Invalid 'node_sequence_length' ({seq_len}).")
            return base_result_entry
    except ValueError:
        base_result_entry[
            "bcftools_error"] = "Skipped: 'grch38_start_position' or 'node_sequence_length' is not a valid integer."
        print(
            f"  Skipping item (HapID: {item_hap_id_for_log}): 'grch38_start_position' or 'node_sequence_length' not a valid int.")
        return base_result_entry

    # Calculate end position based on the user's formula: grch38_start_position + node_sequence_length.
    # VCF coordinates are 1-based.
    # If node_sequence_length is the number of bases in a feature starting at start_pos,
    # the standard inclusive region is chr:start_pos-(start_pos + node_sequence_length - 1).
    # However, the user's formula implies the end coordinate is literally start_pos + node_sequence_length.
    # Example: start=100, length=1 -> user formula end=101 (region 100-101). Standard for 1bp feature: 100-100.
    # Example: start=100, length=0 -> user formula end=100 (region 100-100). This is fine for a point.
    end_pos = start_pos + seq_len

    # Sanity check for calculated end position
    if end_pos < start_pos:  # Should only happen if seq_len was negative (already checked)
        base_result_entry[
            "bcftools_error"] = f"Skipped: Calculated end position ({end_pos}) is less than start position ({start_pos}) with length ({seq_len})."
        print(
            f"  Skipping item (HapID: {item_hap_id_for_log}): Calculated end_pos {end_pos} < start_pos {start_pos} with length {seq_len}.")
        return base_result_entry

    # Construct the region string for bcftools
    query_region_str = f"{chromosome}:{start_pos}-{end_pos}"
    base_result_entry["query_region"] = query_region_str

    # Execute the bcftools query
    variants_list, bcftools_stderr_output = execute_bcftools_query(vcf_file_path, query_region_str)

    if variants_list is not None:  # Indicates bcftools command was attempted (not a FileNotFoundError for bcftools itself)
        base_result_entry["variants"] = variants_list  # Will be an empty list if no variants found

    if bcftools_stderr_output:  # Store any stderr output from bcftools (warnings or errors)
        if "bcftools command not found" in bcftools_stderr_output:
            # This error is critical and indicates a setup problem.
            # The FileNotFoundError will be raised by the main loop when future.result() is called.
            raise FileNotFoundError(bcftools_stderr_output)
        # Append to existing error or set if new, to preserve initial skip reasons.
        if base_result_entry["bcftools_error"]:
            base_result_entry["bcftools_error"] += f"; bcftools stderr: {bcftools_stderr_output}"
        else:
            base_result_entry["bcftools_error"] = f"bcftools stderr: {bcftools_stderr_output}"

    return base_result_entry


def main():
    """
    Main function to parse arguments and run the VCF querying process.
    """
    parser = argparse.ArgumentParser(
        description="Query a VCF file using bcftools view based on regions derived from an input JSON file. "
                    "The input JSON should be an array of objects, each with fields like "
                    "'grch38_chromosome_region', 'grch38_start_position', and 'node_sequence_length'. "
                    "Records with no variants found will not be saved to the output.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "input_json",
        help="Path to the input JSON file (output from the previous script containing regions)."
    )
    parser.add_argument(
        "vcf_file",
        help="Path to the gzipped and indexed VCF file (e.g., SMaHT_COLO829_SNV_truth_set_v1.0.vcf.gz)."
    )
    parser.add_argument(
        "-o", "--output_json",
        default="bcftools_query_results_with_variants.json",
        help="Path to save the output JSON file with bcftools query results (only includes records with variants)."
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=None,  # Will be set to CPU count or 1 if not provided
        help=(f"Number of parallel worker threads to use for bcftools queries. "
              f"(default: number of CPU cores, or 1 if CPU count cannot be determined. "
              f"Currently detected cores: {os.cpu_count() if os.cpu_count() else 'N/A'})")
    )

    args = parser.parse_args()

    # Set MAX_WORKERS based on args.threads or default to CPU count / 1
    max_workers = args.threads if args.threads is not None else (os.cpu_count() if os.cpu_count() else 1)
    if max_workers <= 0:  # Ensure at least one worker
        print("Warning: Number of threads must be positive. Defaulting to 1.")
        max_workers = 1

    # --- Print effective configuration ---
    print("--- Effective Configuration ---")
    print(f"Input JSON: {args.input_json}")
    print(f"VCF File: {args.vcf_file}")
    print(f"Output JSON: {args.output_json}")
    print(f"Max Workers (Threads): {max_workers}")
    print("-----------------------------")

    # --- Validate input files ---
    if not os.path.exists(args.input_json):
        print(f"Error: Input JSON file not found at '{args.input_json}'")
        return

    if not os.path.exists(args.vcf_file):
        print(f"Error: VCF file not found at '{args.vcf_file}'")
        return
    if not (os.path.exists(args.vcf_file + ".tbi") or os.path.exists(args.vcf_file + ".csi")):
        print(f"Error: Index file (.tbi or .csi) not found for '{args.vcf_file}'. "
              "Please index the VCF file using tabix or bcftools index.")
        return

    # --- Load input JSON data ---
    try:
        with open(args.input_json, 'r') as f:
            input_data = json.load(f)
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from '{args.input_json}'. Make sure it's a valid JSON file.")
        return

    if not isinstance(input_data, list):
        print(f"Error: Input JSON is not a list as expected. Found type: {type(input_data)}")
        return

    all_results_with_variants = []

    print(f"\nProcessing {len(input_data)} items from '{args.input_json}'...")
    print(f"Querying VCF: '{args.vcf_file}' with bcftools view (no header).")
    # Max workers already printed in effective configuration

    items_to_process = []
    skipped_due_to_input_validation = 0

    # Pre-filter items that can be processed
    for i, item in enumerate(input_data):
        item_hap_id_for_log = item.get('haplotype_id', f'index_{i}')
        start_pos_str = item.get("grch38_start_position")
        seq_len_val = item.get("node_sequence_length")
        chromosome = item.get("grch38_chromosome_region")

        if not chromosome or start_pos_str is None or seq_len_val is None:
            print(
                f"  [PRE-SKIP] Item (HapID: {item_hap_id_for_log}): Missing one or more required fields (chromosome, start_position, node_sequence_length).")
            skipped_due_to_input_validation += 1
            continue
        try:
            int(start_pos_str)
            if int(seq_len_val) < 0:
                print(
                    f"  [PRE-SKIP] Item (HapID: {item_hap_id_for_log}): Invalid 'node_sequence_length' ({seq_len_val}).")
                skipped_due_to_input_validation += 1
                continue
        except ValueError:
            print(
                f"  [PRE-SKIP] Item (HapID: {item_hap_id_for_log}): 'grch38_start_position' or 'node_sequence_length' not a valid int.")
            skipped_due_to_input_validation += 1
            continue
        items_to_process.append(item)

    if not items_to_process:
        print("No valid items to process after initial input validation.")
        if skipped_due_to_input_validation > 0:
            print(f"Skipped {skipped_due_to_input_validation} items due to missing/invalid fields in input JSON.")
        return

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_item_info = {
            executor.submit(process_item_for_bcftools, item, args.vcf_file):
                item.get('haplotype_id', f'original_index_{input_data.index(item) if item in input_data else "N/A"}')
            for item in items_to_process  # Only submit valid items
        }

        processed_count = 0
        skipped_due_to_no_variants = 0
        for future in concurrent.futures.as_completed(future_to_item_info):
            item_id_for_log = future_to_item_info[future]
            processed_count += 1
            try:
                result_entry = future.result()

                if result_entry.get("variants"):
                    all_results_with_variants.append(result_entry)
                else:
                    skipped_due_to_no_variants += 1
                    # Log if there was an error during processing, or just no variants
                    error_msg = result_entry.get("bcftools_error")
                    if error_msg and "Skipped:" not in error_msg:  # Don't re-log simple skips already printed
                        print(
                            f"  Item (HapID: {item_id_for_log}): No variants found or query failed. Error: {error_msg}")
                    elif not error_msg:  # No error, just no variants
                        print(
                            f"  Item (HapID: {item_id_for_log}): No variants found in region '{result_entry.get('query_region', 'N/A')}'.")

                if processed_count % 50 == 0 or processed_count == len(items_to_process):
                    print(f"  Processed {processed_count}/{len(items_to_process)} submitted tasks...")

            except FileNotFoundError as e:
                print(f"CRITICAL ERROR: {e}. Aborting further processing.")
                for f_cancel in future_to_item_info:
                    if not f_cancel.done() and not f_cancel.cancelled():
                        f_cancel.cancel()
                break
            except Exception as exc:
                print(f"An error occurred processing future for item associated with ID '{item_id_for_log}': {exc}")

    # --- Save all results to the output JSON file ---
    print(f"\n--- Summary ---")
    print(f"Total items in input JSON: {len(input_data)}")
    print(f"Items skipped due to initial validation: {skipped_due_to_input_validation}")
    print(f"Items submitted for bcftools query: {len(items_to_process)}")
    print(f"Items processed: {processed_count}")
    print(f"Records saved to output (with variants): {len(all_results_with_variants)}")
    print(f"Records skipped (no variants found or query skipped): {skipped_due_to_no_variants}")

    try:
        with open(args.output_json, 'w') as outfile:
            json.dump(all_results_with_variants, outfile, indent=4)
        print(f"\nSuccessfully saved {len(all_results_with_variants)} records to '{args.output_json}'")
    except IOError:
        print(f"Error: Could not write results to '{args.output_json}'.")
    except Exception as e:
        print(f"An unexpected error occurred while writing the output JSON: {e}")


if __name__ == "__main__":
    main()
