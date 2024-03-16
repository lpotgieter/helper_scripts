# This script aims to do count the proportion of scaffolds in an assembly and assess which proportion of the contigs are found in a local blast database
# It assumes that the files have the same basename and the blast results have _blastn at the end of the file name
# After counting, it determines the proportion of scaffolds/reads that have matches in the local blast db out of the total
# works with outfmt 6 from a local blast


# Output file
output_file=".log"
error_log="error_log.txt"

# Remove existing output file and error log
rm -f "$output_file" "$error_log"

# Specify the directory
directory=""

# Add headers to the output file
echo -e "File Name\tCount of >\tBlastn Line Count\tRatio (Column 3 / Column 2)" > "$output_file"

# Iterate through each fasta file in the specified directory
for fasta_file in "$directory"/*_metaspades.fasta; do
    # Extract file name without extension
    file_name=$(basename "$fasta_file" _metaspades.fasta)

    # Count the number of lines starting with ">" in the fasta file
    count=$(grep -c ">" "$fasta_file")

    # Check if the corresponding blastn file exists
    blastn_file="${file_name}_metaspades.fasta_blastn"
    if [ -e "$directory/$blastn_file" ]; then
        # Count the number of lines in the blastn file
        blastn_line_count=$(wc -l < "$directory/$blastn_file")

        # Calculate the ratio
        ratio=$(awk "BEGIN {printf \"%.2f\", $blastn_line_count/$count}")

        # Output to the result file with tab separation
        echo -e "$file_name\t$count\t$blastn_line_count\t$ratio" >> "$output_file"
    else
        echo "Error: $directory/$blastn_file does not exist for $fasta_file" >> "$error_log"
    fi
done
