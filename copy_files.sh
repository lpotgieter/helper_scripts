# This script pulls sequences with matching names from subdirectories in a directory and renames them to the subdirectory name
# This is particularly useful for pulling out all scaffolds.fasta files produced by SPAdes, for example
# Simply change the pattern you are looking for to match the pattern you need

# Set the source directory where the folders are located
source_dir=""

# Set the destination directory where you want to copy the files
destination_dir=""

# Iterate through each folder in the source directory
for folder in "$source_dir"/*; do
    # Check if the item in the source directory is a directory
    if [ -d "$folder" ]; then
        # Get the name of the parent folder
        parent_folder=$(basename "$folder")
        
        # Check if contigs.fasta file exists in the current folder
        if [ -e "$folder/scaffolds.fasta" ]; then
            # Copy contigs.fasta to the destination directory and rename it
            cp "$folder/scaffolds.fasta" "$destination_dir/$parent_folder.fasta"
            echo "Copied $folder/scaffolds.fasta to $destination_dir/$parent_folder.fasta"
        else
            echo "contigs.fasta not found in $folder"
        fi
    fi
done
