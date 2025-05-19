import json
import subprocess
import shlex
import argparse
import os
import concurrent.futures


def execute_bcftools_query(vcf_file_path, region_string):
    """
    Executes a bcftools view query for a given region.

    Args:
        vcf_file_path (str): Path to the VCF file.
        region_string (str): The region to query (e.g., "chr1:100-200").

    Returns:
        tuple: (list_of_variant_strings, bcftools_stderr_str)
               Returns (None, error_message) if bcftools command itself fails.
    """
    # bcftools view <file.vcf.gz> <chr:from-to>
    command = ["bcftools", "view", vcf_file_path, region_string]
    print(f"  Executing: {' '.join(shlex.quote(arg) for arg in command)}")
    try:
        process = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=False  # Don't raise exception for non-zero exit, we'll check stderr
        )

        variants_list = []
        if process.stdout:
            # Split stdout into lines, removing any empty lines
            variants_list = [line for line in process.stdout.strip().split('\n') if line]

        return variants_list, process.stderr.strip()

    except FileNotFoundError:
        error_msg = f"Error: 'bcftools' command not found. Please ensure it is installed and in your PATH."
        print(error_msg)
        return None, error_msg  # Special case for bcftools not found
    except Exception as e:
        error_msg = f"An unexpected error occurred while running bcftools for region {region_string}: {e}"
        print(error_msg)
        return None, error_msg


def process_item_for_bcftools(item, vcf_file_path):
    """
    Processes a single item from the input JSON, constructs the region,
    and prepares it for bcftools query.

    Args:
        item (dict): A dictionary from the input JSON.
        vcf_file_path (str): Path to the VCF file.

    Returns:
        dict: A dictionary containing original item info, query region,
              and bcftools variants/error.
    """
    base_result_entry = {
        "top_level_node_id": item.get("top_level_node_id"),
        "haplotype_id": item.get("haplotype_id"),
        "node_sequence_length": item.get("node_sequence_length"),
        "grch38_chromosome_region": item.get("grch38_chromosome_region"),
        "grch38_start_position": item.get("grch38_start_position"),
        "query_region": None,
        "variants": [],  # Initialize as empty list
        "bcftools_error": None
    }

    start_pos_str = item.get("grch38_start_position")
    seq_len_val = item.get("node_sequence_length")
    chromosome = item.get("grch38_chromosome_region")

    if not chromosome:
        base_result_entry["bcftools_error"] = "Skipped: Missing 'grch38_chromosome_region'."
        print(f"  Skipping item (HapID: {base_result_entry['haplotype_id']}): Missing 'grch38_chromosome_region'.")
        return base_result_entry

    if start_pos_str is None:
        base_result_entry["bcftools_error"] = "Skipped: 'grch38_start_position' is null."
        print(f"  Skipping item (HapID: {base_result_entry['haplotype_id']}): 'grch38_start_position' is null.")
        return base_result_entry

    if seq_len_val is None:
        base_result_entry["bcftools_error"] = "Skipped: 'node_sequence_length' is null."
        print(f"  Skipping item (HapID: {base_result_entry['haplotype_id']}): 'node_sequence_length' is null.")
        return base_result_entry

    try:
        start_pos = int(start_pos_str)
        seq_len = int(seq_len_val)
        if seq_len < 0:
            base_result_entry["bcftools_error"] = f"Skipped: Invalid 'node_sequence_length' ({seq_len})."
            print(
                f"  Skipping item (HapID: {base_result_entry['haplotype_id']}): Invalid 'node_sequence_length' ({seq_len}).")
            return base_result_entry
    except ValueError:
        base_result_entry[
            "bcftools_error"] = "Skipped: 'grch38_start_position' or 'node_sequence_length' is not a valid integer."
        print(
            f"  Skipping item (HapID: {base_result_entry['haplotype_id']}): 'grch38_start_position' or 'node_sequence_length' not a valid int.")
        return base_result_entry

    # Calculate end position: start_position + node_sequence_length
    # VCF coordinates are 1-based. If node_sequence_length is 1, the region is start_pos to start_pos.
    # The formula start_pos + node_sequence_length for the end is inclusive if we consider sequence length.
    # However, VCF region chr:start-end is typically inclusive.
    # If node_sequence_length is from a 0-indexed system or represents length,
    # end_pos = start_pos + seq_len -1 for a 1-based inclusive end.
    # But the user specified: grch38_start_position-grch38_start_position+node_sequence_length
    # This implies the end coordinate is literally start + length.
    # If length is 1, region is start-(start+1).
    # Let's assume the user's formula implies the end coordinate.
    # If node_sequence_length is, for example, 9 (from previous script), the region is start to start+9.
    end_pos = start_pos + seq_len

    # A common convention for VCF query is start_pos to (start_pos + length -1) if length is number of bases.
    # If user's `node_sequence_length` is truly the span, then `end_pos = start_pos + seq_len - 1`
    # However, sticking to the user's specified formula: `end_pos = start_pos + seq_len`
    # If `node_sequence_length` is 0, this means `start_pos` to `start_pos`.
    # If `node_sequence_length` is 1, this means `start_pos` to `start_pos+1`.
    # This might include one extra base than intended if node_sequence_length is a typical "length".
    # For now, strictly following: end = start + length

    if end_pos < start_pos:  # This should only happen if seq_len is negative, which is checked.
        base_result_entry[
            "bcftools_error"] = f"Skipped: Calculated end position ({end_pos}) is less than start position ({start_pos}) with length ({seq_len})."
        print(
            f"  Skipping item (HapID: {base_result_entry['haplotype_id']}): Calculated end_pos {end_pos} < start_pos {start_pos} with length {seq_len}.")
        return base_result_entry

    query_region_str = f"{chromosome}:{start_pos}-{end_pos}"
    base_result_entry["query_region"] = query_region_str

    variants, bcftools_stderr = execute_bcftools_query(vcf_file_path, query_region_str)

    if variants is not None:  # bcftools command was attempted
        base_result_entry["variants"] = variants  # This will be an empty list if no variants found
    if bcftools_stderr:
        base_result_entry["bcftools_error"] = bcftools_stderr
        if "bcftools command not found" in bcftools_stderr:
            raise FileNotFoundError(bcftools_stderr)

    return base_result_entry


def main():
    parser = argparse.ArgumentParser(
        description="Query a VCF file using bcftools view based on regions derived from an input JSON file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "input_json",
        help="Path to the input JSON file (output from the previous script)."
    )
    parser.add_argument(
        "vcf_file",
        help="Path to the gzipped and indexed VCF file (e.g., SMaHT_COLO829_SNV_truth_set_v1.0.vcf.gz)."
    )
    parser.add_argument(
        "-o", "--output_json",
        default="bcftools_query_results.json",
        help="Path to save the output JSON file with bcftools results."
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=min(4, os.cpu_count() if os.cpu_count() else 1),
        help="Number of parallel worker threads to use for bcftools queries."
    )

    args = parser.parse_args()

    # Validate VCF file existence and index
    # bcftools can use .tbi or .csi. Checking for .tbi for consistency,
    # but bcftools itself will error out if a suitable index is not found.
    if not os.path.exists(args.vcf_file):
        print(f"Error: VCF file not found at '{args.vcf_file}'")
        return
    if not (os.path.exists(args.vcf_file + ".tbi") or os.path.exists(args.vcf_file + ".csi")):
        print(
            f"Error: Index file (.tbi or .csi) not found for '{args.vcf_file}'. Please index the VCF file using tabix or bcftools index.")
        return

    try:
        with open(args.input_json, 'r') as f:
            input_data = json.load(f)
    except FileNotFoundError:
        print(f"Error: Input JSON file not found at '{args.input_json}'")
        return
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from '{args.input_json}'. Make sure it's a valid JSON file.")
        return

    if not isinstance(input_data, list):
        print(f"Error: Input JSON is not a list as expected. Found type: {type(input_data)}")
        return

    all_results = []

    print(f"Processing {len(input_data)} items from '{args.input_json}'...")
    print(f"Querying VCF: '{args.vcf_file}' with bcftools view.")
    print(f"Using up to {args.threads} parallel threads.")

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        future_to_item_info = {
            executor.submit(process_item_for_bcftools, item, args.vcf_file): item.get('haplotype_id', f'index_{i}')
            for i, item in enumerate(input_data)
        }

        processed_count = 0
        for future in concurrent.futures.as_completed(future_to_item_info):
            item_id_for_log = future_to_item_info[future]
            processed_count += 1
            try:
                result_entry = future.result()
                all_results.append(result_entry)
                if processed_count % 50 == 0 or processed_count == len(input_data):
                    print(f"  Processed {processed_count}/{len(input_data)} items...")
            except FileNotFoundError as e:
                print(f"CRITICAL ERROR: {e}. Aborting further processing.")
                # Cancel remaining futures
                for f_cancel in future_to_item_info:
                    if not f_cancel.done():
                        f_cancel.cancel()
                break
            except Exception as exc:
                print(f"An error occurred processing item associated with ID '{item_id_for_log}': {exc}")
                all_results.append({
                    "haplotype_id": item_id_for_log,  # Or more context from original item
                    "query_region": "N/A - processing error",
                    "variants": [],
                    "bcftools_error": str(exc)
                })

    try:
        with open(args.output_json, 'w') as outfile:
            json.dump(all_results, outfile, indent=4)
        print(f"\nSuccessfully saved {len(all_results)} results to '{args.output_json}'")
    except IOError:
        print(f"Error: Could not write results to '{args.output_json}'.")


if __name__ == "__main__":
    main()
