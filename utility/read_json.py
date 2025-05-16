import json
import subprocess
import shlex  # For safely formatting command arguments for display
import concurrent.futures
import os  # To get CPU count
import argparse  # For command-line argument parsing
import time

# --- Global Configuration (will be updated by CLI args) ---
GBZ_FILE_PATH = "hprc-v1.1-mc-grch38.d9.gbz"  # Default, overridden by CLI
OUTPUT_JSON_PATH = "grch38_node_positions_output.json"  # Default, overridden by CLI
MAX_WORKERS = os.cpu_count() if os.cpu_count() else 1  # Default, overridden by CLI
TASK_TIMEOUT_SECONDS = 300  # Default, overridden by CLI
ENABLE_COUNTING_DIAGNOSTICS = False  # Default, overridden by CLI


def load_json_to_dict(filepath):
    """
    Reads a JSON file and loads its content into a Python dictionary.
    """
    try:
        with open(filepath, 'r') as file_object:
            data_dictionary = json.load(file_object)
            return data_dictionary
    except FileNotFoundError:
        print(f"Error: The file '{filepath}' was not found.")
        return None
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from the file '{filepath}'. Make sure it's a valid JSON file.")
        return None
    except Exception as e:
        print(f"An unexpected error occurred while loading JSON: {e}")
        return None


def get_target_node_id_from_haplotype(haplotype_entry):
    """
    Extracts the target node ID from the thread_context of a haplotype entry.
    """
    if 'thread_context' in haplotype_entry and isinstance(haplotype_entry['thread_context'], list):
        for context_item in haplotype_entry['thread_context']:
            if isinstance(context_item, dict) and \
                    context_item.get('is_target') is True and \
                    'id' in context_item:
                return context_item['id']
    return None


def find_grch38_starting_position(node_id_to_find, chromosome_name_on_grch38, gbz_file_path_for_worker,
                                  worker_task_timeout):
    """
    Uses 'vg find' to get the starting position of a node on a specific GRCh38 path.
    This function is designed to be run in a separate process.
    Returns the position string, "TIMEOUT", or None.
    """
    path_name = f"GRCh38#0#{chromosome_name_on_grch38}"
    command = [
        "vg", "find",
        "-x", gbz_file_path_for_worker,
        "-n", str(node_id_to_find),
        "-P", path_name
    ]

    effective_timeout = worker_task_timeout if worker_task_timeout is not None else 300  # Default if None for some reason

    # This print will show up in the main console, potentially interleaved
    print(
        f"    [Worker PID:{os.getpid()}] Executing: {' '.join(shlex.quote(arg) for arg in command)} (Timeout: {effective_timeout}s)")

    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True, timeout=effective_timeout)
        output = result.stdout.strip()
        if output:
            parts = output.split()  # Splits by any whitespace

            if len(parts) >= 3:
                if str(parts[0]) == str(node_id_to_find):
                    return parts[2]
                else:
                    print(
                        f"    [Worker PID:{os.getpid()}] Warning: Output node ID {parts[0]} does not match queried node ID {node_id_to_find} for path {path_name}: '{output}'")
            elif len(parts) == 2:
                if str(parts[0]) == str(node_id_to_find):
                    print(
                        f"    [Worker PID:{os.getpid()}] Info: Handling abbreviated output for node {node_id_to_find} on {path_name}. Assuming format: node_id position. Output: '{output}'")
                    return parts[1]
                else:
                    print(
                        f"    [Worker PID:{os.getpid()}] Warning: Abbreviated output node ID {parts[0]} does not match queried node ID {node_id_to_find} for path {path_name}: '{output}'")
            else:
                print(
                    f"    [Worker PID:{os.getpid()}] Warning: Unexpected output format (not enough parts) for node {node_id_to_find} on {path_name}: '{output}'")
        else:
            print(
                f"    [Worker PID:{os.getpid()}] Warning: No output from vg find for node {node_id_to_find} on {path_name}.")
    except subprocess.TimeoutExpired:
        print(
            f"    [Worker PID:{os.getpid()}] Timeout: Command for node {node_id_to_find} on {path_name} exceeded {effective_timeout}s.")
        return "TIMEOUT"  # Special marker for timeout
    except subprocess.CalledProcessError as e:
        # This error means vg find ran but returned a non-zero exit code (e.g., node not found on path)
        # It's not necessarily a script error, but vg couldn't find what was asked.
        print(
            f"    [Worker PID:{os.getpid()}] 'vg find' failed for node {node_id_to_find} on {path_name}. Command: '{' '.join(shlex.quote(arg) for arg in e.cmd)}'. "
            f"Return code: {e.returncode}. Stdout: '{e.stdout.strip()}'. Stderr: '{e.stderr.strip()}'")
    except FileNotFoundError:
        # This is a critical error: 'vg' executable itself is not found
        print(
            f"    [Worker PID:{os.getpid()}] CRITICAL Error: 'vg' command not found. Ensure it is installed and in your system's PATH.")
        raise  # Re-raise to make the future fail explicitly, this indicates a setup problem
    except Exception as e:
        # Catch any other unexpected errors during the subprocess execution
        print(
            f"    [Worker PID:{os.getpid()}] Unexpected error during 'vg find' for node {node_id_to_find} on {path_name}: {e}")
    return None  # Indicates failure to get a position for any reason other than explicit timeout marker


def main():
    """
    Main function to parse arguments and run the processing.
    """
    global GBZ_FILE_PATH, OUTPUT_JSON_PATH, MAX_WORKERS, TASK_TIMEOUT_SECONDS, ENABLE_COUNTING_DIAGNOSTICS

    parser = argparse.ArgumentParser(
        description="Process node batch JSON data to find GRCh38 starting positions using 'vg find'.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter  # Shows default values in help
    )
    parser.add_argument(
        "json_file",  # Positional argument
        help="Path to the input JSON file (e.g., node_batch_info.json)."
    )
    parser.add_argument(
        "gbz_file",  # Positional argument
        help="Path to the GBZ graph file (e.g., hprc-v1.1-mc-grch38.d9.gbz)."
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        # Default is set after parsing if not provided by user
        help=f"Number of parallel worker threads/processes to use. (default: number of CPU cores, currently {os.cpu_count() if os.cpu_count() else 1})"
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=300,
        help="Timeout in seconds for each individual 'vg find' command."
    )
    parser.add_argument(
        "-o", "--output",
        default="grch38_node_positions_output.json",
        help="Path to save the output JSON file."
    )
    parser.add_argument(
        "--diag",
        action="store_true",
        default=False,
        help="Enable detailed diagnostics during the task counting phase."
    )

    args = parser.parse_args()

    # Set global configurations from parsed arguments
    json_file_path_cli = args.json_file
    GBZ_FILE_PATH = args.gbz_file
    MAX_WORKERS = args.threads if args.threads is not None else (os.cpu_count() if os.cpu_count() else 1)
    TASK_TIMEOUT_SECONDS = args.timeout
    OUTPUT_JSON_PATH = args.output
    ENABLE_COUNTING_DIAGNOSTICS = args.diag

    # --- Start of processing logic ---
    main_json_data = load_json_to_dict(json_file_path_cli)

    simplified_results_summary = []  # Store the simplified results
    tasks_to_submit_count = 0
    submitted_tasks_count = 0
    completed_tasks_count = 0
    successful_finds_count = 0
    timed_out_tasks_count = 0
    other_failed_tasks_count = 0

    if main_json_data:
        print("--- Task Identification Phase (Counting) ---")
        for top_level_node_id_str, node_specific_data in main_json_data.items():
            if not isinstance(node_specific_data, dict):
                if ENABLE_COUNTING_DIAGNOSTICS: print(
                    f"  [COUNT_SKIP] Top-level entry for '{top_level_node_id_str}' is not a dictionary.")
                continue
            if 'haplotypes' not in node_specific_data or not isinstance(node_specific_data['haplotypes'], list):
                if ENABLE_COUNTING_DIAGNOSTICS: print(
                    f"  [COUNT_SKIP] Node '{top_level_node_id_str}': Missing 'haplotypes' list or not a list.")
                continue

            haplotype_index = 0
            for haplotype_entry in node_specific_data['haplotypes']:
                haplotype_index += 1
                if not isinstance(haplotype_entry, dict):
                    if ENABLE_COUNTING_DIAGNOSTICS: print(
                        f"  [COUNT_SKIP] Node '{top_level_node_id_str}', Haplotype entry #{haplotype_index}: Not a dictionary.")
                    continue

                hap_id_for_log = haplotype_entry.get('haplotype_id', f'UnknownHapID_Index{haplotype_index}')
                is_grch38 = haplotype_entry.get('chromosome') == 'GRCh38'
                has_region = haplotype_entry.get('region') is not None and haplotype_entry.get('region') != ""
                target_node_id = get_target_node_id_from_haplotype(haplotype_entry)
                has_target_node = target_node_id is not None

                if is_grch38 and has_region and has_target_node:
                    tasks_to_submit_count += 1
                elif ENABLE_COUNTING_DIAGNOSTICS:
                    skip_reasons = []
                    if not is_grch38: skip_reasons.append(
                        f"Not GRCh38 (chromosome: {haplotype_entry.get('chromosome')})")
                    if not has_region: skip_reasons.append("Missing 'region'")
                    if not has_target_node: skip_reasons.append("No target node in thread_context")
                    print(
                        f"  [COUNT_SKIP] Node '{top_level_node_id_str}', Haplotype ID '{hap_id_for_log}': Skipped because ({'; '.join(skip_reasons)})")
        print(f"--- End Task Identification Phase ---")
        print(f"Identified {tasks_to_submit_count} potential 'vg find' tasks to execute.")

    if main_json_data and tasks_to_submit_count > 0:
        print(f"\nProcessing nodes from '{json_file_path_cli}' using GBZ file '{GBZ_FILE_PATH}'.")
        print(f"Using up to {MAX_WORKERS} parallel worker processes.")
        print(f"Individual task timeout set to: {TASK_TIMEOUT_SECONDS} seconds.")
        print(f"Output will be saved to: '{OUTPUT_JSON_PATH}'\n")

        with concurrent.futures.ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
            futures_to_metadata = {}  # Map future object to its descriptive metadata for context
            for top_level_node_id_str, node_specific_data in main_json_data.items():
                if not isinstance(node_specific_data, dict): continue
                if 'haplotypes' in node_specific_data and isinstance(node_specific_data['haplotypes'], list):
                    for haplotype_entry in node_specific_data['haplotypes']:
                        if not isinstance(haplotype_entry, dict): continue
                        if haplotype_entry.get('chromosome') == 'GRCh38':
                            grch38_region_name = haplotype_entry.get('region')
                            if not grch38_region_name: continue  # Ensure region is not None or empty

                            target_node_id_in_haplotype = get_target_node_id_from_haplotype(haplotype_entry)
                            if target_node_id_in_haplotype is not None:
                                future = executor.submit(
                                    find_grch38_starting_position,
                                    target_node_id_in_haplotype,
                                    grch38_region_name,
                                    GBZ_FILE_PATH,  # Pass the global GBZ_FILE_PATH
                                    TASK_TIMEOUT_SECONDS  # Pass the global TASK_TIMEOUT_SECONDS
                                )
                                # Store metadata needed for the simplified output, plus target for logging
                                metadata_for_future = {
                                    "top_level_node_id": top_level_node_id_str,
                                    "haplotype_id": haplotype_entry.get('haplotype_id', 'N/A'),
                                    "grch38_chromosome_region": grch38_region_name,
                                    "_target_node_id_for_logging": target_node_id_in_haplotype  # Internal for logging
                                }
                                futures_to_metadata[future] = metadata_for_future
                                submitted_tasks_count += 1
                                if submitted_tasks_count % 100 == 0 or submitted_tasks_count == tasks_to_submit_count:
                                    print(f"  Submitted {submitted_tasks_count}/{tasks_to_submit_count} tasks...")

            print(f"\nAll {submitted_tasks_count} tasks submitted. Waiting for completion (this may take a while)...\n")

            for future in concurrent.futures.as_completed(futures_to_metadata):
                completed_tasks_count += 1
                metadata = futures_to_metadata[future]  # Get the metadata associated with this future

                # Prepare data for the simplified output object
                output_entry = {
                    "top_level_node_id": metadata["top_level_node_id"],
                    "haplotype_id": metadata["haplotype_id"],
                    "grch38_chromosome_region": metadata["grch38_chromosome_region"],
                    "grch38_start_position": None  # Default to null
                }

                # For logging purposes
                log_prefix = (f"  Processed {completed_tasks_count}/{submitted_tasks_count} | "
                              f"TN {metadata['top_level_node_id']}, "
                              f"Hap {metadata['haplotype_id']}, "
                              f"Target {metadata['_target_node_id_for_logging']}: ")
                try:
                    # Add a small buffer to the future.result() timeout, as the worker has its own timeout
                    start_position_result = future.result(timeout=TASK_TIMEOUT_SECONDS + 20)  # Increased buffer

                    if start_position_result == "TIMEOUT":
                        print(f"{log_prefix}TIMED OUT (in worker)")
                        timed_out_tasks_count += 1
                        # output_entry["grch38_start_position"] remains None (null in JSON)
                    elif start_position_result is not None:  # Successfully got a position string
                        print(f"{log_prefix}SUCCESS -> Position: {start_position_result}")
                        output_entry["grch38_start_position"] = start_position_result
                        successful_finds_count += 1
                    else:  # Worker returned None (e.g. vg find failed, no output, unexpected format)
                        print(f"{log_prefix}FAILED (no position found or worker error, see worker logs)")
                        other_failed_tasks_count += 1
                        # output_entry["grch38_start_position"] remains None (null in JSON)

                except concurrent.futures.TimeoutError:  # future.result() timed out
                    print(f"{log_prefix}TIMED OUT (waiting for future result)")
                    timed_out_tasks_count += 1
                    # output_entry["grch38_start_position"] remains None (null in JSON)
                except Exception as exc:  # Other exceptions from task execution (e.g., vg not found if worker raised it)
                    print(f"{log_prefix}ERROR (task execution generated an exception: {exc})")
                    other_failed_tasks_count += 1
                    # output_entry["grch38_start_position"] remains None (null in JSON)

                simplified_results_summary.append(output_entry)

        print(f"\n--- Processing Complete ---")
        print(f"Total tasks identified by initial scan: {tasks_to_submit_count}")
        print(f"Total tasks submitted to executor: {submitted_tasks_count}")
        print(f"Total tasks completed (results processed): {completed_tasks_count}")
        print(f"Successfully found positions: {successful_finds_count}")
        print(f"Timed out tasks (worker or future): {timed_out_tasks_count}")
        print(f"Other failed tasks (no position/error): {other_failed_tasks_count}")

        try:
            with open(OUTPUT_JSON_PATH, 'w') as outfile:
                json.dump(simplified_results_summary, outfile, indent=4)
            print(f"\nResults successfully saved to '{OUTPUT_JSON_PATH}'")
        except IOError as e:
            print(f"\nError: Could not write results to '{OUTPUT_JSON_PATH}'. Reason: {e}")
        except Exception as e:
            print(f"\nAn unexpected error occurred while writing JSON output: {e}")

    elif tasks_to_submit_count == 0 and main_json_data:
        print(
            "No suitable GRCh38 haplotypes with target nodes found in the input JSON to process after diagnostic scan.")
    elif not main_json_data:
        print(f"Failed to load or parse JSON data from '{json_file_path_cli}'. Cannot proceed.")


# --- Entry point ---
if __name__ == "__main__":
    main()
