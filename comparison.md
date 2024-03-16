ISR2 Comparison
================
Lizel Potgieter
2024-03-16

## Libraries and organisation

``` r
#if you don't have pacman, uncomment next line
#install.packages("pacman")

pacman::p_load(here,
               tidyverse,
               UpSetR,
               reshape2)
```

``` r
here::i_am("comparison.md")
```

## Logic

To determine which amplicon sequencing method was best, I assembled the
filtered reads with mothur, megahit, and metaspades. I then created a
local blast database with target sequences. I performed a blastn, with
the maximum number of target sequences set to 1 on all 3 assembled
datasets. As a comparison, I converted the filtered fastq files to fasta
files and performed that same blast on them. Here I will show the
proportion of unique sequences retrieved by each method compared to the
unassembled sequences.

## File Loading

This assumes that there are 4 subdirectories in input called ‘raw’,
‘megahit’, ‘metaspades’, and ‘mothur’ and that files have the extension
’\_blastn_target1’.

### fastq

These results are on forward and reverse sequences. These outputs are
then combined to form a single dataframe.

``` r
# Set the directory path to your subfolder
subfolder <- "input/raw/"

# Get a list of all files in the subfolder
file_list <- list.files(subfolder, full.names = TRUE)

# Loop through the list of files
for (file in file_list) {
  # Extract the filename without extension
  filename <- tools::file_path_sans_ext(basename(file))
  
  # Load the file into a dataframe
  df <- read.delim(file, header = F)
  
  # Rename the dataframe with "_fastq" appended at the end
  df_name <- paste0(filename, "_fastq")
  
  # Remove the specific string from the dataframe name
  df_name <- gsub("\\.fq\\.gz\\.fasta_blastn_target1", "", df_name)
  
  # Assign the dataframe to the new name
  assign(df_name, df)
}

# Loop through all objects in the environment
for (obj in ls()) {
  # Check if the object is a dataframe and ends with "_1_fastq"
  if (is.data.frame(get(obj)) && grepl("_1_fastq$", obj)) {
    # Add a column called "orientation" filled with "forward"
    assign(obj, transform(get(obj), orientation = "forward"))
  }
}
for (obj in ls()) {
  # Check if the object is a dataframe and ends with "_1_fastq"
  if (is.data.frame(get(obj)) && grepl("_2_fastq$", obj)) {
    # Add a column called "orientation" filled with "forward"
    assign(obj, transform(get(obj), orientation = "reverse"))
  }
}



# Initialize a list to store merged dataframes
merged_dfs <- list()

# Extract unique base names from the list of dataframes
base_names <- unique(sub("_[12]_fastq$", "", ls(pattern = "_[12]_fastq")))

# Loop through each unique base name
for (base_name in base_names) {
  # Create a pattern for matching both _1_fastq and _2_fastq dataframes
  pattern <- paste0("^", base_name, "_[12]_fastq$")

  # Filter dataframes based on the pattern
  dfs <- mget(ls(pattern = pattern))
  
  # Check if both _1_fastq and _2_fastq dataframes exist
  if (length(dfs) == 2) {
    # Merge the two dataframes by row index
    merged_df <- do.call(rbind, dfs)
    
    # Remove the _1 and _2 suffixes from the base name
    new_base_name <- sub("_[12]_fastq$", "", base_name)
    
    # Add the merged dataframe to the list with the new base name
    merged_dfs[[new_base_name]] <- merged_df
    

    
    # Remove the individual dataframes used to make the combined one
    rm(list = ls(pattern = pattern))
  }
}

# Loop through the merged dataframes and assign them to the environment
for (base_name in names(merged_dfs)) {
  assign(paste0(base_name, "_fastq"), merged_dfs[[base_name]])
}
    # Remove the intermediate dataframes
    rm(dfs)
    rm(df)
    rm(merged_df)
    rm(merged_dfs)
```

### Megahit

``` r
# Set the directory path to your subfolder
subfolder <- "input/megahit/"

# Get a list of all files in the subfolder
file_list <- list.files(subfolder, full.names = TRUE)

excluded_files <- c()  # Initialize a vector to store excluded file names

# Loop through the list of files
for (file in file_list) {
  # Check if the file is empty
  if (file.size(file) == 0) {
    # If empty, print the filename and add it to excluded_files
    cat("Excluded empty file:", file, "\n")
    excluded_files <- c(excluded_files, file)
    next  # Move to the next iteration
  }
  
  # Extract the filename without extension
  filename <- tools::file_path_sans_ext(basename(file))
  
  # Load the file into a dataframe
  df <- read.delim(file, header = FALSE)
  
  # Rename the dataframe with "_fastq" appended at the end
  df_name <- paste0(filename, "_megahit")
  
  # Remove the specific string from the dataframe name
  df_name <- gsub("\\.fq\\.gz_megahit\\.fasta_blastn_target1", "", df_name)
  
  # Assign the dataframe to the new name
  assign(df_name, df)
}
```

    ## Excluded empty file: input/megahit/T2_1.fq.gz_megahit.fasta_blastn_target1 
    ## Excluded empty file: input/megahit/T27_1.fq.gz_megahit.fasta_blastn_target1 
    ## Excluded empty file: input/megahit/T37_1.fq.gz_megahit.fasta_blastn_target1

``` r
# Print the excluded files
cat("Excluded files:", excluded_files, "\n")
```

    ## Excluded files: input/megahit/T2_1.fq.gz_megahit.fasta_blastn_target1 input/megahit/T27_1.fq.gz_megahit.fasta_blastn_target1 input/megahit/T37_1.fq.gz_megahit.fasta_blastn_target1

### Metaspades

``` r
# Set the directory path to your subfolder
subfolder <- "input/metaspades/"

excluded_files <- c()  # Initialize a vector to store excluded file names

# Get a list of all files in the subfolder
file_list <- list.files(subfolder, full.names = TRUE)

# Loop through the list of files
for (file in file_list) {
  # Check if the file is empty
  if (file.size(file) == 0) {
    # If empty, print the filename and add it to excluded_files
    cat("Excluded empty file:", file, "\n")
    excluded_files <- c(excluded_files, file)
    next  # Move to the next iteration
  }
  
  # Extract the filename without extension
  filename <- tools::file_path_sans_ext(basename(file))
  
  # Load the file into a dataframe
  df <- read.delim(file, header = FALSE)
  
  # Rename the dataframe with "_fastq" appended at the end
  df_name <- paste0(filename, "_metaspades")
  
  # Remove the specific string from the dataframe name
  df_name <- gsub("\\.fq\\.gz_metaspades\\.fasta_blastn_target1", "", df_name)
  
  # Assign the dataframe to the new name
  assign(df_name, df)
}
```

    ## Excluded empty file: input/metaspades/T1_1.fq.gz_metaspades.fasta_blastn_target1 
    ## Excluded empty file: input/metaspades/T15_1.fq.gz_metaspades.fasta_blastn_target1 
    ## Excluded empty file: input/metaspades/T24_1.fq.gz_metaspades.fasta_blastn_target1 
    ## Excluded empty file: input/metaspades/T25_1.fq.gz_metaspades.fasta_blastn_target1 
    ## Excluded empty file: input/metaspades/T28_1.fq.gz_metaspades.fasta_blastn_target1 
    ## Excluded empty file: input/metaspades/T30_1.fq.gz_metaspades.fasta_blastn_target1 
    ## Excluded empty file: input/metaspades/T7_1.fq.gz_metaspades.fasta_blastn_target1 
    ## Excluded empty file: input/metaspades/T8_1.fq.gz_metaspades.fasta_blastn_target1 
    ## Excluded empty file: input/metaspades/T9_1.fq.gz_metaspades.fasta_blastn_target1

``` r
# Print the excluded files
cat("Excluded files:", excluded_files, "\n")
```

    ## Excluded files: input/metaspades/T1_1.fq.gz_metaspades.fasta_blastn_target1 input/metaspades/T15_1.fq.gz_metaspades.fasta_blastn_target1 input/metaspades/T24_1.fq.gz_metaspades.fasta_blastn_target1 input/metaspades/T25_1.fq.gz_metaspades.fasta_blastn_target1 input/metaspades/T28_1.fq.gz_metaspades.fasta_blastn_target1 input/metaspades/T30_1.fq.gz_metaspades.fasta_blastn_target1 input/metaspades/T7_1.fq.gz_metaspades.fasta_blastn_target1 input/metaspades/T8_1.fq.gz_metaspades.fasta_blastn_target1 input/metaspades/T9_1.fq.gz_metaspades.fasta_blastn_target1

### mothur

``` r
# Set the directory path to your subfolder
subfolder <- "input/mothur/"

# Get a list of all files in the subfolder
file_list <- list.files(subfolder, full.names = TRUE)

# Loop through the list of files
for (file in file_list) {
  # Extract the filename without extension
  filename <- tools::file_path_sans_ext(basename(file))
  
  # Load the file into a dataframe
  df <- read.delim(file, header = F)
  
  # Rename the dataframe with "_fastq" appended at the end
  df_name <- paste0(filename, "_mothur")
  
  # Remove the specific string from the dataframe name
  df_name <- gsub("\\.fq\\.trim\\.contigs\\.fasta_blastn_target1", "", df_name)
  
  # Assign the dataframe to the new name
  assign(df_name, df)
}

rm(df)
```

## Counting Accessions Recovered

### Fixing formating

``` r
# Loop through all dataframes in the environment
for (obj in ls()) {
  # Check if the object is a dataframe
  if (is.data.frame(get(obj))) {
    # Check if the column V2 exists in the dataframe
    if ("V2" %in% colnames(get(obj))) {
      # Remove everything after and including ":" in column V2
      assign(obj, transform(get(obj), V2 = sub(":.*", "", V2)))
    }
  }
}
```

### Counting

``` r
# Initialize an empty dataframe to store the counts
final_df <- data.frame(base_name = character(),
                        suffix = character(),
                        unique_value = character(),
                        count = numeric(),
                        stringsAsFactors = FALSE)

# Extract unique base names from the list of dataframes
base_names <- unique(sub("_.*$", "", ls()))

# Loop through each unique base name
for (base_name in base_names) {
  # Loop through each dataframe with the current base name
  for (df_name in ls(pattern = paste0("^", base_name, "_"))) {
    # Extract the suffix from the dataframe name
    suffix <- sub("^.*_", "", df_name)
    
    # Get the dataframe
    df <- get(df_name)
    
    # Check if column "V2" exists in the dataframe
    if ("V2" %in% colnames(df)) {
      # Count occurrences of each unique value in column V2
      counts <- table(df$V2)
      
      # Extract unique values and their counts
      unique_values <- names(counts)
      counts <- as.vector(counts)
      
      # Append the counts to the final dataframe
      final_df <- rbind(final_df, data.frame(base_name = base_name,
                                             suffix = suffix,
                                             unique_value = unique_values,
                                             count = counts,
                                             stringsAsFactors = FALSE))
    }
  }
}

# Rename columns
colnames(final_df) <- c("sample", "method", "accession", "count")

# Print the final dataframe
#Uncomment this next line if you want to see the dataframe
#print(final_df)
```

## Plotting

``` r
# Split the original dataframe by the 'sample' column
df_list <- split(final_df, final_df$sample)

# Loop through each dataframe in the list
for (sample_name in names(df_list)) {
  # Create a new dataframe name
  new_df_name <- paste0(sample_name, "_upset")
  
  # Extract the dataframe
  sample_df <- df_list[[sample_name]]
  
  # Get unique accessions
  unique_accessions <- unique(sample_df$accession)
  
  # Create an empty dataframe to store the UpSet data
  upset_data <- data.frame(accession = unique_accessions)
  
  # Loop through unique methods
  for (method in unique(sample_df$method)) {
    # Check if the method is present for each accession
    upset_data[method] <- ifelse(upset_data$accession %in% sample_df$accession[sample_df$method == method], 1, 0)
  }
  
  # Write the dataframe to a new variable
  assign(new_df_name, upset_data)
}
```

``` r
upset(Tx_upset, nsets =4)
```

## Conclusion

Through this analysis we can see that assemblies produced by mothur has
the best recovery rate between the 3 assembly methods compared to
sequences recovered from the unassembled reads. It is worth considering
that the target sequences in the database are highly similar, and the
higher abundance of accessions in the unassembled sets may not be true
and should be investigated more closely if the researcher chooses to
continue with these samples.
