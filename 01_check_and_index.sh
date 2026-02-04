#!/usr/bin/env bash  # Shebang: Defines the script interpreter as Bash
set -euo pipefail    # Strict mode: exit on Error, Unset vars, or Pipe failures

# Usage: ./01_check_and_index.sh -i /path/to/bams -t 4
INPUT_DIR=""         # Initialises empty variable for the directory path
THREADS="${THREADS:-4}"  # Sets THREADS to 4 unless already set in the environment

usage() {            # Defines a function named 'usage' to show help text
  echo "Usage: $0 -i INPUT_DIR [-t THREADS]" >&2  # Prints help to 'stderr' (error output)
}                    # Closes the function definition

while [[ $# -gt 0 ]]; do  # Loop while number of arguments ($#) is greater than 0
  case "$1" in            # Start 'case' (multiple choice) block to check the current flag ($1)
    -i)                   # If the flag is '-i' (input)
        INPUT_DIR="$2"    # Assign the next argument ($2) to the variable
        shift 2           # Move the argument list forward by 2 to clear '-i' and its value
        ;;                # Double semicolon: ends the instructions for this specific case
    -t)                   # If the flag is '-t' (threads)
        THREADS="$2"      # Assign the next argument to the THREADS variable
        shift 2           # Move the argument list forward by 2
        ;;
    -h|--help)            # If the flag is '-h' or '--help'
        usage             # Run the 'usage' function defined above
        exit 0            # Exit the script with a 'success' code (0)
        ;;
    *)                    # The '*' is a wildcard (matches anything else/unknown flags)
        echo "Unknown option: $1" # Warn the user about the invalid input
        usage             # Show the correct usage
        exit 2            # Exit with an error code (2) for 'misuse of shell built-in'
        ;;
  esac                    # 'esac' is 'case' backwards: marks the end of the case block
done                      # Marks the end of the 'while' loop

[[ -z "${INPUT_DIR}" ]] && { usage; exit 1; }  # If INPUT_DIR is empty (-z), show help and exit
command -v samtools >/dev/null || { echo "samtools not found"; exit 1; } # Check if samtools command is available

# mapfile: reads lines into an array. -t: trims newlines.
# < <(...) is 'process substitution': it feeds the output of the find/sort commands into mapfile
mapfile -t BAMS < <(find "${INPUT_DIR}" -maxdepth 1 -type f -name "*.bam" | sort)

[[ ${#BAMS[@]} -gt 0 ]] || { echo "No BAMs in ${INPUT_DIR}"; exit 1; } # Exit if array length is 0

for BAM in "${BAMS[@]}"; do      # Start 'for' loop: process each file in the BAMS array
  echo ">> Checking ${BAM}"      # Print status message
  if ! samtools quickcheck -v "${BAM}"; then  # If 'samtools' returns an error (!)
    echo "!! ${BAM} failed samtools quickcheck — skip" # Warn user
    continue                     # Skip the rest of this loop and move to the next file
  fi                             # 'fi' is 'if' backwards: marks the end of the 'if' block

  if [[ -f "${BAM}.bai" ]]; then # If a file (-f) with the .bai extension already exists
    echo "== Index exists: ${BAM}.bai"  # Inform user and do nothing
  else                           # Otherwise (if the index file is missing)
    echo ">> Indexing ${BAM}"    # Inform user indexing is starting
    samtools index -@ "${THREADS}" "${BAM}" # Run indexing using specified threads (-@)
  fi                             # End of the 'if' block for indexing
done                             # End of the 'for' loop

echo "Done."                     # Print final completion message
