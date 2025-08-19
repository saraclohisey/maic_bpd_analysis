#!/bin/bash
#
# First, split your file into records. You will need to install coreutils on macos
# for this:
#   - brew install coreutils
#
# Then:
#   - gcsplit --elide-empty-files <input file> '/<record separator>/1' '{*}'
#
# Now use this script

input_dir="./" # Set the directory containing the small files
output_dir="./output" # Set the directory where the combined files will be saved
files_per_group=500
current_group=0
current_file_count=0

mkdir -p "$output_dir"

for file in "$input_dir"/xx*; do
  if [ "$current_file_count" -eq 0 ]; then
    combined_file="${output_dir}/combined_${current_group}.txt"
    touch "$combined_file"
  fi

  cat "$file" >> "$combined_file"
  current_file_count=$((current_file_count + 1))

  if [ "$current_file_count" -eq "$files_per_group" ]; then
    current_group=$((current_group + 1))
    current_file_count=0
  fi
done

echo "Combining complete. Files are saved in $output_dir."

