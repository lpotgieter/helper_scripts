SASUF2024 Visualisation
================
Lizel Potgieter
2024-05-20

## Loading Packages and file

``` r
pacman::p_load(tidyverse, here, rtracklayer, treemap)
here::i_am("sasuf_visuals.Rmd")
```

    ## here() starts at C:/Users/llpo0001/Documents/sasuf/data

These are two really helpful packages to make your R journey more
pleasant! `Pacman` is a package that checks to see whether the packages
you are loading are installed. If they aren’t `pacman` installs them for
you before loading. `Here` is a package that sets your working directory
to the file you have set in the inverted commas.

``` r
#interpro table
beta <- read_tsv("beta_interpro.tsv", col_names = FALSE)
```

    ## Rows: 360928 Columns: 15
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (11): X1, X2, X4, X5, X6, X9, X11, X12, X13, X14, X15
    ## dbl  (3): X3, X7, X8
    ## lgl  (1): X10
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
#rename relevant columns
#you can rename as many columns as you'd like, but here we are just renaming the ones we will need for this visualisation
# Rename columns using indexing
colnames(beta)[1] <- "accession"
colnames(beta)[4] <- "tool"
colnames(beta)[5] <- "annotation"
colnames(beta)[6] <- "functional"
```

There are many ways of doing this. For today, we are only visualising
the Pfam accessions.

``` r
#interpro table
cerco <- read_tsv("cerco_interpro.tsv", col_names = FALSE)
```

    ## Rows: 111844 Columns: 15
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (11): X1, X2, X4, X5, X6, X9, X11, X12, X13, X14, X15
    ## dbl  (3): X3, X7, X8
    ## lgl  (1): X10
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# Rename columns using indexing
colnames(cerco)[1] <- "accession"
colnames(cerco)[4] <- "tool"
colnames(cerco)[5] <- "annotation"
colnames(cerco)[6] <- "functional"
```

We only care whether something has a signal peptide here, so we keep one
prediction per accession. If you care where the signal peptide is, you
can filter according to that.

``` r
#signalp <- read.delim("region_output.gff3", header = FALSE, stringsAsFactors = FALSE)

signalp <- readGFF("region_output.gff3")
signalp$seqid <- gsub("\\s.*", "", signalp$seqid)
colnames(signalp)[colnames(signalp) == "seqid"] <- "accession"
signalp <- unique(signalp)

# Extract unique accessions from signalp dataframe
unique_accessions <- unique(signalp$accession)

# Add signal_peptide column to cerco dataframe
cerco <- cerco %>%
  mutate(signal_peptide = ifelse(accession %in% unique_accessions, "yes", "no"))
```

We do the same for effectors

``` r
# Read the file
file_content <- readLines("cerco_effectorp")

# Find the index of the line containing "# Identifier"
identifier_index <- grep("^# Identifier", file_content)

# Extract lines from "# Identifier" to the end of the file
data_lines <- file_content[identifier_index:length(file_content)]

# Split each line into separate elements
data_elements <- lapply(strsplit(data_lines, "\t"), function(x) {
  if (length(x) < 13) {
    x <- c(x, rep(NA, 13 - length(x)))  # Pad with NA values if the number of elements is less than 13
  }
  return(x)
})

# Convert the list of elements into a dataframe
df <- as.data.frame(do.call(rbind, data_elements), stringsAsFactors = FALSE)

# Set column names based on the first row
colnames(df) <- df[1, ]

# Remove the first row (column names row)
df <- df[-1, ]

effectors <- df
colnames(effectors)[1] <- "accession"
effectors_filtered <- effectors[effectors$Prediction != "Non-effector", ]

effectors_filtered$accession <- gsub("\\s.*", "", effectors_filtered$accession)
unique_accessions <- unique(effectors_filtered$accession)


cerco <- cerco %>%
  mutate(effector = ifelse(accession %in% unique_accessions, "yes", "no"))
```

``` r
deeploc <- read.csv("cerco_deeploc.csv")
colnames(deeploc)[1] <- "accession"

# Assuming df1 and df2 are your data frames
merged_df <- merge(cerco, deeploc, by = "accession", all = TRUE)
```

Finally we have a dataframe that has effectors, signal peptides, and
localisation included. Now we can do some plotting!

Let us also have a look at the `merged_df`. Here you will see the effect
that prediction tools have. InterProScan did not annotated
XP_023448401.1 as having any particular features to annotate. Deeploc,
however, indicated that it is located in the cytoplasm. This is to show
that there are definitely benefits to using several prediction tools on
each dataset! We are merely making the computer look for patterns. If 1
out of 3 software notices there’s a pattern, it might be a false
positive. It just be a hypothetical protein.

``` r
head(merged_df)
```

    ##        accession                               X2  X3        tool annotation
    ## 1 XP_023448401.1                             <NA>  NA        <NA>       <NA>
    ## 2 XP_023448404.1 b7b28b12e4e17545ec70df5cc05d9d3a 443      PRINTS    PR00420
    ## 3 XP_023448404.1 b7b28b12e4e17545ec70df5cc05d9d3a 443      PRINTS    PR00420
    ## 4 XP_023448404.1 b7b28b12e4e17545ec70df5cc05d9d3a 443      PRINTS    PR00420
    ## 5 XP_023448404.1 b7b28b12e4e17545ec70df5cc05d9d3a 443      PRINTS    PR00420
    ## 6 XP_023448404.1 b7b28b12e4e17545ec70df5cc05d9d3a 443 SUPERFAMILY   SSF54373
    ##                                                         functional  X7  X8
    ## 1                                                             <NA>  NA  NA
    ## 2 Aromatic-ring hydroxylase (flavoprotein monooxygenase) signature   8  30
    ## 3 Aromatic-ring hydroxylase (flavoprotein monooxygenase) signature 162 177
    ## 4 Aromatic-ring hydroxylase (flavoprotein monooxygenase) signature 308 323
    ## 5 Aromatic-ring hydroxylase (flavoprotein monooxygenase) signature 323 339
    ## 6                         FAD-linked reductases, C-terminal domain 195 287
    ##        X9  X10        X11  X12  X13  X14  X15 signal_peptide effector
    ## 1    <NA>   NA       <NA> <NA> <NA> <NA> <NA>           <NA>     <NA>
    ## 2 4.5E-16 TRUE 08-05-2024    -    -    -    -             no       no
    ## 3 4.5E-16 TRUE 08-05-2024    -    -    -    -             no       no
    ## 4 4.5E-16 TRUE 08-05-2024    -    -    -    -             no       no
    ## 5 4.5E-16 TRUE 08-05-2024    -    -    -    -             no       no
    ## 6 1.73E-9 TRUE 08-05-2024    -    -    -    -             no       no
    ##   Localizations Signals Cytoplasm Nucleus Extracellular Cell.membrane
    ## 1     Cytoplasm            0.4871  0.4247        0.1228        0.0654
    ## 2          <NA>    <NA>        NA      NA            NA            NA
    ## 3          <NA>    <NA>        NA      NA            NA            NA
    ## 4          <NA>    <NA>        NA      NA            NA            NA
    ## 5          <NA>    <NA>        NA      NA            NA            NA
    ## 6          <NA>    <NA>        NA      NA            NA            NA
    ##   Mitochondrion Plastid Endoplasmic.reticulum Lysosome.Vacuole Golgi.apparatus
    ## 1        0.4422  0.0953                0.0501           0.2219          0.1762
    ## 2            NA      NA                    NA               NA              NA
    ## 3            NA      NA                    NA               NA              NA
    ## 4            NA      NA                    NA               NA              NA
    ## 5            NA      NA                    NA               NA              NA
    ## 6            NA      NA                    NA               NA              NA
    ##   Peroxisome
    ## 1     0.0242
    ## 2         NA
    ## 3         NA
    ## 4         NA
    ## 5         NA
    ## 6         NA

## Plotting effectors!

For simplicity, we will only consider the pfam annotations. This
excludes all of the proteins that InterProScan did not have annotations
predicted for. Depending on your project, you might want to change this
parameter. For this exercise, it is exactly what we are looking for.

``` r
cerco_pfam <- merged_df[grep("Pfam", merged_df$tool), ]
```

Let us see which proportion of all proteins have signal peptides

``` r
# Calculate the proportion of "yes" values in the signal_peptide column
proportion_signal_peptide_yes <- mean(cerco_pfam$signal_peptide == "yes")

# Calculate the proportion of "yes" values in the effector column
proportion_effector_yes <- mean(cerco_pfam$effector == "yes")

# Create a data frame for plotting
plot_data <- data.frame(
  Variable = c("Signal Peptide", "Effector"),
  Proportion = c(proportion_signal_peptide_yes, proportion_effector_yes)
)

# Create the combined bar plot
combined_plot <- ggplot(plot_data, aes(x = Variable, y = Proportion, fill = Variable)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(Proportion * 100, 2), "%")), vjust = -0.5) +
  labs(y = "Proportion", title = "Proportion of 'Yes' in Signal Peptide and Effector Columns") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())  # Remove x-axis labels and ticks

# Print the combined plot
print(combined_plot)
```

![](sasuf_visuals_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
# Subset dataframe to include only rows where the "effector" column has "yes" values
effector_yes_df <- cerco_pfam[cerco_pfam$effector == "yes", ]

# Tabulate frequencies of each annotation
annotation_freq <- table(effector_yes_df$annotation)

# Convert the table to a dataframe
annotation_freq_df <- as.data.frame(annotation_freq)
names(annotation_freq_df) <- c("Annotation", "Frequency")

# Create bar plot
ggplot(annotation_freq_df, aes(x = Annotation, y = Frequency, fill = Annotation)) +
  geom_bar(stat = "identity") +
  labs(y = "Frequency", title = "Frequency of Annotations for 'Yes' in Effector Column") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability
```

![](sasuf_visuals_files/figure-gfm/unnamed-chunk-5-1.png)<!-- --> This
is far too messy! Let us have a look at different types of
visualisations.

``` r
# Subset dataframe to include only rows where the "effector" column has "yes" values
effector_yes_df <- cerco_pfam[cerco_pfam$effector == "yes", ]

# Tabulate frequencies of each annotation
annotation_freq <- table(effector_yes_df$annotation)

# Convert the table to a dataframe
annotation_freq_df <- as.data.frame(annotation_freq)
names(annotation_freq_df) <- c("Annotation", "Frequency")

# Create treemap
treemap(annotation_freq_df,
        index = "Annotation",
        vSize = "Frequency",
        title = "Frequency of Annotations for 'Yes' in Effector Column")
```

![](sasuf_visuals_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
# Sort annotation_freq_df based on frequency in descending order
sorted_annotations <- annotation_freq_df[order(-annotation_freq_df$Frequency), ]

# Select top 10 annotations
top_10_annotations <- head(sorted_annotations$Annotation, 10)

# Subset cerco_pfam based on top 10 annotations
matching_rows <- cerco_pfam[cerco_pfam$annotation %in% top_10_annotations, ]

# Extract unique functions from the matching rows
unique_functions <- unique(matching_rows$functional)

# Print the unique top 10 functions
print(unique_functions)
```

    ##  [1] "Pathogen effector; putative necrosis-inducing factor"
    ##  [2] "Peroxidase, family 2"                                
    ##  [3] "Glycosyl hydrolases family 16"                       
    ##  [4] "Glycosyl hydrolases family 43"                       
    ##  [5] "LysM domain"                                         
    ##  [6] "Cysteine-rich secretory protein family"              
    ##  [7] "Hydrophobic surface binding protein A"               
    ##  [8] "Glycosyl hydrolases family 11"                       
    ##  [9] "Cutinase"                                            
    ## [10] "PAN domain"

``` r
# Sort annotation_freq_df based on frequency in descending order
sorted_annotations <- annotation_freq_df[order(-annotation_freq_df$Frequency), ]

# Select top 10 annotations
top_10_annotations <- head(sorted_annotations$Annotation, 10)

# Print top 10 annotations
print(top_10_annotations)
```

    ##  [1] PF01083 PF14295 PF01476 PF12296 PF00722 PF01328 PF04616 PF14856 PF00188
    ## [10] PF00457
    ## 103 Levels: PF00024 PF00026 PF00081 PF00082 PF00085 PF00089 PF00096 ... PF21203

``` r
# Subset cerco_pfam based on top 10 annotations
top_10_functions <- cerco_pfam[cerco_pfam$annotation %in% top_10_annotations, "function"]

# Check if any rows match the top 10 annotations in cerco_pfam
matching_rows <- cerco_pfam[cerco_pfam$annotation %in% top_10_annotations, ]

# Print matching rows to check
#print(matching_rows)
```

``` r
# Sort annotation_freq_df based on frequency in descending order
sorted_annotations <- annotation_freq_df[order(-annotation_freq_df$Frequency), ]

# Select top 10 annotations
top_10_annotations <- head(sorted_annotations$Annotation, 10)

# Print top 10 annotations
print(top_10_annotations)
```

    ##  [1] PF01083 PF14295 PF01476 PF12296 PF00722 PF01328 PF04616 PF14856 PF00188
    ## [10] PF00457
    ## 103 Levels: PF00024 PF00026 PF00081 PF00082 PF00085 PF00089 PF00096 ... PF21203

``` r
# Subset cerco_pfam based on top 10 annotations
top_10_functions <- cerco_pfam[cerco_pfam$annotation %in% top_10_annotations, "functional"]

# Check if any rows match the top 10 annotations in cerco_pfam
matching_rows <- cerco_pfam[cerco_pfam$annotation %in% top_10_annotations, ]

# Print matching rows with specific columns
print(matching_rows[, c("accession", "annotation", "functional", "signal_peptide", "effector", "Localizations")])
```

    ##             accession annotation
    ## 25     XP_023448410.1    PF14856
    ## 188    XP_023448431.1    PF01328
    ## 426    XP_023448463.1    PF01328
    ## 1612   XP_023448622.1    PF01328
    ## 2786   XP_023448754.1    PF00722
    ## 4560   XP_023448960.1    PF04616
    ## 4561   XP_023448960.1    PF04616
    ## 5890   XP_023449107.1    PF04616
    ## 6210   XP_023449133.1    PF04616
    ## 7307   XP_023449247.1    PF00722
    ## 9154   XP_023449460.1    PF00722
    ## 11835  XP_023449751.1    PF01328
    ## 12872  XP_023449859.1    PF01328
    ## 13767  XP_023449982.1    PF01476
    ## 13778  XP_023449982.1    PF01476
    ## 13783  XP_023449982.1    PF01476
    ## 14563  XP_023450069.1    PF01328
    ## 15908  XP_023450215.1    PF00188
    ## 16374  XP_023450266.1    PF12296
    ## 16451  XP_023450277.1    PF01328
    ## 16508  XP_023450285.1    PF12296
    ## 16846  XP_023450328.1    PF01328
    ## 18019  XP_023450473.1    PF01476
    ## 18020  XP_023450473.1    PF01476
    ## 18427  XP_023450522.1    PF12296
    ## 18714  XP_023450559.1    PF12296
    ## 18715  XP_023450560.1    PF12296
    ## 18721  XP_023450561.1    PF12296
    ## 18734  XP_023450563.1    PF12296
    ## 21034  XP_023450823.1    PF01476
    ## 21756  XP_023450896.1    PF01328
    ## 22782  XP_023451019.1    PF01328
    ## 23319  XP_023451084.1    PF00722
    ## 24951  XP_023451276.1    PF01328
    ## 25491  XP_023451335.1    PF01328
    ## 25744  XP_023451365.1    PF04616
    ## 27366  XP_023451554.1    PF01476
    ## 28209  XP_023451648.1    PF01328
    ## 32291  XP_023452110.1    PF12296
    ## 34120  XP_023452284.1    PF01328
    ## 34588  XP_023452340.1    PF00457
    ## 34750  XP_023452360.1    PF01083
    ## 35390  XP_023452437.1    PF01328
    ## 36069  XP_023452533.1    PF01328
    ## 37198  XP_023452679.1    PF01083
    ## 39676  XP_023452952.1    PF01083
    ## 41522  XP_023453188.1    PF01083
    ## 42059  XP_023453241.1    PF04616
    ## 42703  XP_023453334.1    PF12296
    ## 44659  XP_023453530.1    PF01083
    ## 45092  XP_023453572.1    PF00457
    ## 45931  XP_023453681.1    PF01083
    ## 46733  XP_023453770.1    PF04616
    ## 46959  XP_023453796.1    PF14295
    ## 47091  XP_023453817.1    PF04616
    ## 47335  XP_023453841.1    PF04616
    ## 50386  XP_023454160.1    PF14856
    ## 50964  XP_023454229.1    PF04616
    ## 51814  XP_023454327.1    PF00188
    ## 52251  XP_023454377.1    PF14295
    ## 52253  XP_023454377.1    PF14295
    ## 52714  XP_023454419.1    PF14856
    ## 52907  XP_023454436.1    PF01328
    ## 53088  XP_023454459.1    PF14295
    ## 53089  XP_023454459.1    PF14295
    ## 53099  XP_023454459.1    PF14295
    ## 53723  XP_023454516.1    PF12296
    ## 55683  XP_023454749.1    PF01476
    ## 55684  XP_023454749.1    PF01476
    ## 55685  XP_023454749.1    PF01476
    ## 55865  XP_023454771.1    PF14295
    ## 55866  XP_023454771.1    PF14295
    ## 55998  XP_023454787.1    PF12296
    ## 56943  XP_023454889.1    PF00722
    ## 59195  XP_023455148.1    PF12296
    ## 59863  XP_023455218.1    PF00457
    ## 61381  XP_023455429.1    PF01083
    ## 64718  XP_023455801.1    PF01328
    ## 64722  XP_023455802.1    PF04616
    ## 65617  XP_023455909.1    PF00722
    ## 65836  XP_023455938.1    PF01083
    ## 65971  XP_023455960.1    PF00722
    ## 66346  XP_023456003.1    PF00457
    ## 67393  XP_023456129.1    PF00722
    ## 67399  XP_023456130.1    PF00722
    ## 68093  XP_023456221.1    PF04616
    ## 70937  XP_023456545.1    PF14295
    ## 70938  XP_023456545.1    PF14295
    ## 70939  XP_023456545.1    PF14295
    ## 70940  XP_023456545.1    PF14295
    ## 70941  XP_023456545.1    PF14295
    ## 70942  XP_023456545.1    PF14295
    ## 71663  XP_023456632.1    PF01328
    ## 72225  XP_023456700.1    PF00722
    ## 72288  XP_023456710.1    PF01328
    ## 73021  XP_023456800.1    PF00188
    ## 73344  XP_023456825.1    PF01328
    ## 74368  XP_023456953.1    PF00722
    ## 74693  XP_023456984.1    PF04616
    ## 75341  XP_023457052.1    PF01083
    ## 79162  XP_023457465.1    PF01328
    ## 79922  XP_023457559.1    PF14856
    ## 80294  XP_023457598.1    PF12296
    ## 80581  XP_023457627.1    PF12296
    ## 80853  XP_023457665.1    PF01476
    ## 80855  XP_023457665.1    PF01476
    ## 80878  XP_023457666.1    PF01476
    ## 80892  XP_023457666.1    PF01476
    ## 82706  XP_023457835.1    PF00188
    ## 82769  XP_023457846.1    PF01328
    ## 84090  XP_023457997.1    PF01328
    ## 85156  XP_023458108.1    PF00722
    ## 91849  XP_023458810.1    PF14856
    ## 96061  XP_023459224.1    PF01083
    ## 96541  XP_023459282.1    PF14295
    ## 99898  XP_023459603.1    PF04616
    ## 100921 XP_023459710.1    PF01328
    ## 101220 XP_023459739.1    PF00722
    ## 105819 XP_023460210.1    PF01328
    ## 106000 XP_023460227.1    PF01328
    ## 106501 XP_023460274.1    PF12296
    ## 106803 XP_023460313.1    PF01328
    ## 106854 XP_023460320.1    PF12296
    ## 108684 XP_023460490.1    PF00722
    ## 108888 XP_023460511.1    PF12296
    ## 109761 XP_023460616.1    PF01328
    ## 109981 XP_023460642.1    PF14295
    ## 109982 XP_023460642.1    PF14295
    ## 111682 XP_023460829.1    PF01328
    ##                                                  functional signal_peptide
    ## 25     Pathogen effector; putative necrosis-inducing factor            yes
    ## 188                                    Peroxidase, family 2            yes
    ## 426                                    Peroxidase, family 2            yes
    ## 1612                                   Peroxidase, family 2            yes
    ## 2786                          Glycosyl hydrolases family 16            yes
    ## 4560                          Glycosyl hydrolases family 43            yes
    ## 4561                          Glycosyl hydrolases family 43            yes
    ## 5890                          Glycosyl hydrolases family 43            yes
    ## 6210                          Glycosyl hydrolases family 43             no
    ## 7307                          Glycosyl hydrolases family 16            yes
    ## 9154                          Glycosyl hydrolases family 16            yes
    ## 11835                                  Peroxidase, family 2            yes
    ## 12872                                  Peroxidase, family 2            yes
    ## 13767                                           LysM domain            yes
    ## 13778                                           LysM domain            yes
    ## 13783                                           LysM domain            yes
    ## 14563                                  Peroxidase, family 2            yes
    ## 15908                Cysteine-rich secretory protein family            yes
    ## 16374                 Hydrophobic surface binding protein A            yes
    ## 16451                                  Peroxidase, family 2            yes
    ## 16508                 Hydrophobic surface binding protein A            yes
    ## 16846                                  Peroxidase, family 2            yes
    ## 18019                                           LysM domain             no
    ## 18020                                           LysM domain             no
    ## 18427                 Hydrophobic surface binding protein A            yes
    ## 18714                 Hydrophobic surface binding protein A             no
    ## 18715                 Hydrophobic surface binding protein A            yes
    ## 18721                 Hydrophobic surface binding protein A            yes
    ## 18734                 Hydrophobic surface binding protein A            yes
    ## 21034                                           LysM domain             no
    ## 21756                                  Peroxidase, family 2            yes
    ## 22782                                  Peroxidase, family 2            yes
    ## 23319                         Glycosyl hydrolases family 16            yes
    ## 24951                                  Peroxidase, family 2            yes
    ## 25491                                  Peroxidase, family 2             no
    ## 25744                         Glycosyl hydrolases family 43            yes
    ## 27366                                           LysM domain            yes
    ## 28209                                  Peroxidase, family 2            yes
    ## 32291                 Hydrophobic surface binding protein A             no
    ## 34120                                  Peroxidase, family 2             no
    ## 34588                         Glycosyl hydrolases family 11            yes
    ## 34750                                              Cutinase            yes
    ## 35390                                  Peroxidase, family 2            yes
    ## 36069                                  Peroxidase, family 2             no
    ## 37198                                              Cutinase            yes
    ## 39676                                              Cutinase            yes
    ## 41522                                              Cutinase            yes
    ## 42059                         Glycosyl hydrolases family 43            yes
    ## 42703                 Hydrophobic surface binding protein A            yes
    ## 44659                                              Cutinase            yes
    ## 45092                         Glycosyl hydrolases family 11            yes
    ## 45931                                              Cutinase            yes
    ## 46733                         Glycosyl hydrolases family 43            yes
    ## 46959                                            PAN domain            yes
    ## 47091                         Glycosyl hydrolases family 43            yes
    ## 47335                         Glycosyl hydrolases family 43            yes
    ## 50386  Pathogen effector; putative necrosis-inducing factor            yes
    ## 50964                         Glycosyl hydrolases family 43            yes
    ## 51814                Cysteine-rich secretory protein family            yes
    ## 52251                                            PAN domain            yes
    ## 52253                                            PAN domain            yes
    ## 52714  Pathogen effector; putative necrosis-inducing factor            yes
    ## 52907                                  Peroxidase, family 2            yes
    ## 53088                                            PAN domain            yes
    ## 53089                                            PAN domain            yes
    ## 53099                                            PAN domain            yes
    ## 53723                 Hydrophobic surface binding protein A            yes
    ## 55683                                           LysM domain            yes
    ## 55684                                           LysM domain            yes
    ## 55685                                           LysM domain            yes
    ## 55865                                            PAN domain            yes
    ## 55866                                            PAN domain            yes
    ## 55998                 Hydrophobic surface binding protein A            yes
    ## 56943                         Glycosyl hydrolases family 16            yes
    ## 59195                 Hydrophobic surface binding protein A            yes
    ## 59863                         Glycosyl hydrolases family 11            yes
    ## 61381                                              Cutinase            yes
    ## 64718                                  Peroxidase, family 2            yes
    ## 64722                         Glycosyl hydrolases family 43            yes
    ## 65617                         Glycosyl hydrolases family 16            yes
    ## 65836                                              Cutinase            yes
    ## 65971                         Glycosyl hydrolases family 16             no
    ## 66346                         Glycosyl hydrolases family 11            yes
    ## 67393                         Glycosyl hydrolases family 16            yes
    ## 67399                         Glycosyl hydrolases family 16            yes
    ## 68093                         Glycosyl hydrolases family 43            yes
    ## 70937                                            PAN domain            yes
    ## 70938                                            PAN domain            yes
    ## 70939                                            PAN domain            yes
    ## 70940                                            PAN domain            yes
    ## 70941                                            PAN domain            yes
    ## 70942                                            PAN domain            yes
    ## 71663                                  Peroxidase, family 2            yes
    ## 72225                         Glycosyl hydrolases family 16            yes
    ## 72288                                  Peroxidase, family 2            yes
    ## 73021                Cysteine-rich secretory protein family            yes
    ## 73344                                  Peroxidase, family 2            yes
    ## 74368                         Glycosyl hydrolases family 16            yes
    ## 74693                         Glycosyl hydrolases family 43            yes
    ## 75341                                              Cutinase            yes
    ## 79162                                  Peroxidase, family 2            yes
    ## 79922  Pathogen effector; putative necrosis-inducing factor            yes
    ## 80294                 Hydrophobic surface binding protein A            yes
    ## 80581                 Hydrophobic surface binding protein A            yes
    ## 80853                                           LysM domain            yes
    ## 80855                                           LysM domain            yes
    ## 80878                                           LysM domain            yes
    ## 80892                                           LysM domain            yes
    ## 82706                Cysteine-rich secretory protein family            yes
    ## 82769                                  Peroxidase, family 2            yes
    ## 84090                                  Peroxidase, family 2             no
    ## 85156                         Glycosyl hydrolases family 16            yes
    ## 91849  Pathogen effector; putative necrosis-inducing factor            yes
    ## 96061                                              Cutinase            yes
    ## 96541                                            PAN domain            yes
    ## 99898                         Glycosyl hydrolases family 43             no
    ## 100921                                 Peroxidase, family 2             no
    ## 101220                        Glycosyl hydrolases family 16             no
    ## 105819                                 Peroxidase, family 2            yes
    ## 106000                                 Peroxidase, family 2            yes
    ## 106501                Hydrophobic surface binding protein A            yes
    ## 106803                                 Peroxidase, family 2            yes
    ## 106854                Hydrophobic surface binding protein A            yes
    ## 108684                        Glycosyl hydrolases family 16            yes
    ## 108888                Hydrophobic surface binding protein A             no
    ## 109761                                 Peroxidase, family 2            yes
    ## 109981                                           PAN domain            yes
    ## 109982                                           PAN domain            yes
    ## 111682                                 Peroxidase, family 2            yes
    ##        effector               Localizations
    ## 25          yes                   Cytoplasm
    ## 188          no                        <NA>
    ## 426          no                        <NA>
    ## 1612         no                        <NA>
    ## 2786        yes               Extracellular
    ## 4560         no                        <NA>
    ## 4561         no                        <NA>
    ## 5890         no                        <NA>
    ## 6210         no                        <NA>
    ## 7307         no                        <NA>
    ## 9154        yes               Extracellular
    ## 11835        no                        <NA>
    ## 12872       yes               Extracellular
    ## 13767       yes               Extracellular
    ## 13778       yes               Extracellular
    ## 13783       yes               Extracellular
    ## 14563        no                        <NA>
    ## 15908       yes               Extracellular
    ## 16374       yes                   Cytoplasm
    ## 16451        no                        <NA>
    ## 16508       yes               Extracellular
    ## 16846        no                        <NA>
    ## 18019        no                        <NA>
    ## 18020        no                        <NA>
    ## 18427        no                        <NA>
    ## 18714        no                        <NA>
    ## 18715        no                        <NA>
    ## 18721       yes     Cytoplasm|Extracellular
    ## 18734        no                        <NA>
    ## 21034        no                        <NA>
    ## 21756        no                        <NA>
    ## 22782       yes               Extracellular
    ## 23319        no                        <NA>
    ## 24951        no                        <NA>
    ## 25491        no                        <NA>
    ## 25744       yes               Extracellular
    ## 27366       yes               Extracellular
    ## 28209        no                        <NA>
    ## 32291        no                        <NA>
    ## 34120        no                        <NA>
    ## 34588       yes               Extracellular
    ## 34750       yes               Extracellular
    ## 35390        no                        <NA>
    ## 36069        no                        <NA>
    ## 37198       yes Extracellular|Cell membrane
    ## 39676       yes               Extracellular
    ## 41522       yes               Cell membrane
    ## 42059       yes               Extracellular
    ## 42703       yes               Extracellular
    ## 44659        no                        <NA>
    ## 45092       yes               Extracellular
    ## 45931       yes               Extracellular
    ## 46733        no                        <NA>
    ## 46959        no                        <NA>
    ## 47091        no                        <NA>
    ## 47335       yes               Extracellular
    ## 50386       yes               Extracellular
    ## 50964        no                        <NA>
    ## 51814       yes               Extracellular
    ## 52251       yes               Extracellular
    ## 52253       yes               Extracellular
    ## 52714       yes               Extracellular
    ## 52907        no                        <NA>
    ## 53088        no                        <NA>
    ## 53089        no                        <NA>
    ## 53099        no                        <NA>
    ## 53723        no                        <NA>
    ## 55683        no                        <NA>
    ## 55684        no                        <NA>
    ## 55685        no                        <NA>
    ## 55865       yes               Extracellular
    ## 55866       yes               Extracellular
    ## 55998       yes               Extracellular
    ## 56943        no                        <NA>
    ## 59195        no                        <NA>
    ## 59863       yes               Extracellular
    ## 61381       yes               Extracellular
    ## 64718        no                        <NA>
    ## 64722        no                        <NA>
    ## 65617        no                        <NA>
    ## 65836       yes               Extracellular
    ## 65971        no                        <NA>
    ## 66346       yes               Extracellular
    ## 67393       yes               Extracellular
    ## 67399       yes               Extracellular
    ## 68093       yes               Extracellular
    ## 70937        no                        <NA>
    ## 70938        no                        <NA>
    ## 70939        no                        <NA>
    ## 70940        no                        <NA>
    ## 70941        no                        <NA>
    ## 70942        no                        <NA>
    ## 71663        no                        <NA>
    ## 72225        no                        <NA>
    ## 72288        no                        <NA>
    ## 73021       yes               Extracellular
    ## 73344        no                        <NA>
    ## 74368        no                        <NA>
    ## 74693       yes               Extracellular
    ## 75341       yes               Extracellular
    ## 79162       yes               Extracellular
    ## 79922       yes               Extracellular
    ## 80294        no                        <NA>
    ## 80581       yes               Extracellular
    ## 80853        no                        <NA>
    ## 80855        no                        <NA>
    ## 80878       yes               Extracellular
    ## 80892       yes               Extracellular
    ## 82706       yes               Extracellular
    ## 82769       yes               Extracellular
    ## 84090        no                        <NA>
    ## 85156        no                        <NA>
    ## 91849       yes               Extracellular
    ## 96061       yes               Extracellular
    ## 96541       yes               Extracellular
    ## 99898        no                        <NA>
    ## 100921       no                        <NA>
    ## 101220       no                        <NA>
    ## 105819       no                        <NA>
    ## 106000       no                        <NA>
    ## 106501       no                        <NA>
    ## 106803       no                        <NA>
    ## 106854       no                        <NA>
    ## 108684      yes               Extracellular
    ## 108888       no                        <NA>
    ## 109761      yes               Extracellular
    ## 109981      yes               Extracellular
    ## 109982      yes               Extracellular
    ## 111682       no                        <NA>

From this table we can see that not all PFam families are effectors
(which is interesting, just as is). We can subsitute the `cerco_pfam`
dataframe with the `effector_yes_df` dataframe to have a look at only
the effectors and where they are localised.

``` r
# Sort annotation_freq_df based on frequency in descending order
sorted_annotations <- annotation_freq_df[order(-annotation_freq_df$Frequency), ]

# Select top 10 annotations
top_10_annotations <- head(sorted_annotations$Annotation, 10)

# Print top 10 annotations
print(top_10_annotations)
```

    ##  [1] PF01083 PF14295 PF01476 PF12296 PF00722 PF01328 PF04616 PF14856 PF00188
    ## [10] PF00457
    ## 103 Levels: PF00024 PF00026 PF00081 PF00082 PF00085 PF00089 PF00096 ... PF21203

``` r
# Subset cerco_pfam based on top 10 annotations
top_10_functions <- effector_yes_df[effector_yes_df$annotation %in% top_10_annotations, "functional"]

# Check if any rows match the top 10 annotations in cerco_pfam
matching_rows <- effector_yes_df[effector_yes_df$annotation %in% top_10_annotations, ]

# Print matching rows with specific columns
print(matching_rows[, c("accession", "annotation", "functional", "signal_peptide", "effector", "Localizations")])
```

    ##             accession annotation
    ## 25     XP_023448410.1    PF14856
    ## 2786   XP_023448754.1    PF00722
    ## 9154   XP_023449460.1    PF00722
    ## 12872  XP_023449859.1    PF01328
    ## 13767  XP_023449982.1    PF01476
    ## 13778  XP_023449982.1    PF01476
    ## 13783  XP_023449982.1    PF01476
    ## 15908  XP_023450215.1    PF00188
    ## 16374  XP_023450266.1    PF12296
    ## 16508  XP_023450285.1    PF12296
    ## 18721  XP_023450561.1    PF12296
    ## 22782  XP_023451019.1    PF01328
    ## 25744  XP_023451365.1    PF04616
    ## 27366  XP_023451554.1    PF01476
    ## 34588  XP_023452340.1    PF00457
    ## 34750  XP_023452360.1    PF01083
    ## 37198  XP_023452679.1    PF01083
    ## 39676  XP_023452952.1    PF01083
    ## 41522  XP_023453188.1    PF01083
    ## 42059  XP_023453241.1    PF04616
    ## 42703  XP_023453334.1    PF12296
    ## 45092  XP_023453572.1    PF00457
    ## 45931  XP_023453681.1    PF01083
    ## 47335  XP_023453841.1    PF04616
    ## 50386  XP_023454160.1    PF14856
    ## 51814  XP_023454327.1    PF00188
    ## 52251  XP_023454377.1    PF14295
    ## 52253  XP_023454377.1    PF14295
    ## 52714  XP_023454419.1    PF14856
    ## 55865  XP_023454771.1    PF14295
    ## 55866  XP_023454771.1    PF14295
    ## 55998  XP_023454787.1    PF12296
    ## 59863  XP_023455218.1    PF00457
    ## 61381  XP_023455429.1    PF01083
    ## 65836  XP_023455938.1    PF01083
    ## 66346  XP_023456003.1    PF00457
    ## 67393  XP_023456129.1    PF00722
    ## 67399  XP_023456130.1    PF00722
    ## 68093  XP_023456221.1    PF04616
    ## 73021  XP_023456800.1    PF00188
    ## 74693  XP_023456984.1    PF04616
    ## 75341  XP_023457052.1    PF01083
    ## 79162  XP_023457465.1    PF01328
    ## 79922  XP_023457559.1    PF14856
    ## 80581  XP_023457627.1    PF12296
    ## 80878  XP_023457666.1    PF01476
    ## 80892  XP_023457666.1    PF01476
    ## 82706  XP_023457835.1    PF00188
    ## 82769  XP_023457846.1    PF01328
    ## 91849  XP_023458810.1    PF14856
    ## 96061  XP_023459224.1    PF01083
    ## 96541  XP_023459282.1    PF14295
    ## 108684 XP_023460490.1    PF00722
    ## 109761 XP_023460616.1    PF01328
    ## 109981 XP_023460642.1    PF14295
    ## 109982 XP_023460642.1    PF14295
    ##                                                  functional signal_peptide
    ## 25     Pathogen effector; putative necrosis-inducing factor            yes
    ## 2786                          Glycosyl hydrolases family 16            yes
    ## 9154                          Glycosyl hydrolases family 16            yes
    ## 12872                                  Peroxidase, family 2            yes
    ## 13767                                           LysM domain            yes
    ## 13778                                           LysM domain            yes
    ## 13783                                           LysM domain            yes
    ## 15908                Cysteine-rich secretory protein family            yes
    ## 16374                 Hydrophobic surface binding protein A            yes
    ## 16508                 Hydrophobic surface binding protein A            yes
    ## 18721                 Hydrophobic surface binding protein A            yes
    ## 22782                                  Peroxidase, family 2            yes
    ## 25744                         Glycosyl hydrolases family 43            yes
    ## 27366                                           LysM domain            yes
    ## 34588                         Glycosyl hydrolases family 11            yes
    ## 34750                                              Cutinase            yes
    ## 37198                                              Cutinase            yes
    ## 39676                                              Cutinase            yes
    ## 41522                                              Cutinase            yes
    ## 42059                         Glycosyl hydrolases family 43            yes
    ## 42703                 Hydrophobic surface binding protein A            yes
    ## 45092                         Glycosyl hydrolases family 11            yes
    ## 45931                                              Cutinase            yes
    ## 47335                         Glycosyl hydrolases family 43            yes
    ## 50386  Pathogen effector; putative necrosis-inducing factor            yes
    ## 51814                Cysteine-rich secretory protein family            yes
    ## 52251                                            PAN domain            yes
    ## 52253                                            PAN domain            yes
    ## 52714  Pathogen effector; putative necrosis-inducing factor            yes
    ## 55865                                            PAN domain            yes
    ## 55866                                            PAN domain            yes
    ## 55998                 Hydrophobic surface binding protein A            yes
    ## 59863                         Glycosyl hydrolases family 11            yes
    ## 61381                                              Cutinase            yes
    ## 65836                                              Cutinase            yes
    ## 66346                         Glycosyl hydrolases family 11            yes
    ## 67393                         Glycosyl hydrolases family 16            yes
    ## 67399                         Glycosyl hydrolases family 16            yes
    ## 68093                         Glycosyl hydrolases family 43            yes
    ## 73021                Cysteine-rich secretory protein family            yes
    ## 74693                         Glycosyl hydrolases family 43            yes
    ## 75341                                              Cutinase            yes
    ## 79162                                  Peroxidase, family 2            yes
    ## 79922  Pathogen effector; putative necrosis-inducing factor            yes
    ## 80581                 Hydrophobic surface binding protein A            yes
    ## 80878                                           LysM domain            yes
    ## 80892                                           LysM domain            yes
    ## 82706                Cysteine-rich secretory protein family            yes
    ## 82769                                  Peroxidase, family 2            yes
    ## 91849  Pathogen effector; putative necrosis-inducing factor            yes
    ## 96061                                              Cutinase            yes
    ## 96541                                            PAN domain            yes
    ## 108684                        Glycosyl hydrolases family 16            yes
    ## 109761                                 Peroxidase, family 2            yes
    ## 109981                                           PAN domain            yes
    ## 109982                                           PAN domain            yes
    ##        effector               Localizations
    ## 25          yes                   Cytoplasm
    ## 2786        yes               Extracellular
    ## 9154        yes               Extracellular
    ## 12872       yes               Extracellular
    ## 13767       yes               Extracellular
    ## 13778       yes               Extracellular
    ## 13783       yes               Extracellular
    ## 15908       yes               Extracellular
    ## 16374       yes                   Cytoplasm
    ## 16508       yes               Extracellular
    ## 18721       yes     Cytoplasm|Extracellular
    ## 22782       yes               Extracellular
    ## 25744       yes               Extracellular
    ## 27366       yes               Extracellular
    ## 34588       yes               Extracellular
    ## 34750       yes               Extracellular
    ## 37198       yes Extracellular|Cell membrane
    ## 39676       yes               Extracellular
    ## 41522       yes               Cell membrane
    ## 42059       yes               Extracellular
    ## 42703       yes               Extracellular
    ## 45092       yes               Extracellular
    ## 45931       yes               Extracellular
    ## 47335       yes               Extracellular
    ## 50386       yes               Extracellular
    ## 51814       yes               Extracellular
    ## 52251       yes               Extracellular
    ## 52253       yes               Extracellular
    ## 52714       yes               Extracellular
    ## 55865       yes               Extracellular
    ## 55866       yes               Extracellular
    ## 55998       yes               Extracellular
    ## 59863       yes               Extracellular
    ## 61381       yes               Extracellular
    ## 65836       yes               Extracellular
    ## 66346       yes               Extracellular
    ## 67393       yes               Extracellular
    ## 67399       yes               Extracellular
    ## 68093       yes               Extracellular
    ## 73021       yes               Extracellular
    ## 74693       yes               Extracellular
    ## 75341       yes               Extracellular
    ## 79162       yes               Extracellular
    ## 79922       yes               Extracellular
    ## 80581       yes               Extracellular
    ## 80878       yes               Extracellular
    ## 80892       yes               Extracellular
    ## 82706       yes               Extracellular
    ## 82769       yes               Extracellular
    ## 91849       yes               Extracellular
    ## 96061       yes               Extracellular
    ## 96541       yes               Extracellular
    ## 108684      yes               Extracellular
    ## 109761      yes               Extracellular
    ## 109981      yes               Extracellular
    ## 109982      yes               Extracellular

We can also have a look at effectors that have one a single Pfam
annotation, indicating that these are rare in the genome

``` r
# Subset annotation_freq_df to include only rows where the frequency is 1
single_instance_annotations <- annotation_freq_df[annotation_freq_df$Frequency == 1, "Annotation"]

# Print single instance annotations
print(single_instance_annotations)
```

    ##  [1] PF00026 PF00081 PF00096 PF00160 PF00262 PF00313 PF00314 PF00383 PF00445
    ## [10] PF00544 PF00657 PF00684 PF00840 PF00967 PF01042 PF01095 PF01120 PF01423
    ## [19] PF01556 PF01607 PF01738 PF01764 PF01793 PF02221 PF02678 PF02777 PF03022
    ## [28] PF03067 PF03537 PF03664 PF03856 PF03928 PF04117 PF04420 PF05270 PF05390
    ## [37] PF05630 PF05726 PF05922 PF07249 PF07452 PF07749 PF07915 PF08212 PF08892
    ## [46] PF09206 PF09362 PF09430 PF09451 PF09769 PF09996 PF10528 PF10604 PF12454
    ## [55] PF12999 PF13015 PF13229 PF13577 PF13912 PF14040 PF14558 PF16541 PF19271
    ## [64] PF19535 PF21203
    ## 103 Levels: PF00024 PF00026 PF00081 PF00082 PF00085 PF00089 PF00096 ... PF21203

``` r
# Subset cerco_pfam based on single instance annotations
single_instance_functions <- effector_yes_df[effector_yes_df$annotation %in% single_instance_annotations, "functional"]

# Check if any rows match the single instance annotations in cerco_pfam
matching_rows_single <- effector_yes_df[effector_yes_df$annotation %in% single_instance_annotations, ]

# Print matching rows with specific columns
print(matching_rows_single[, c("accession", "annotation", "functional", "signal_peptide", "effector", "Localizations")])
```

    ##             accession annotation
    ## 529    XP_023448483.1    PF01120
    ## 556    XP_023448490.1    PF05270
    ## 561    XP_023448490.1    PF09206
    ## 2167   XP_023448691.1    PF00313
    ## 3211   XP_023448810.1    PF07249
    ## 3907   XP_023448891.1    PF01793
    ## 8035   XP_023449335.1    PF09996
    ## 9164   XP_023449462.1    PF13577
    ## 10384  XP_023449594.1    PF03928
    ## 11080  XP_023449678.1    PF19271
    ## 11626  XP_023449725.1    PF01764
    ## 12996  XP_023449879.1    PF00160
    ## 16829  XP_023450325.1    PF00840
    ## 17838  XP_023450452.1    PF04117
    ## 18472  XP_023450532.1    PF08892
    ## 21721  XP_023450893.1    PF00445
    ## 27287  XP_023451538.1    PF19535
    ## 30295  XP_023451869.1    PF00096
    ## 30304  XP_023451869.1    PF13912
    ## 32213  XP_023452102.1    PF00383
    ## 32744  XP_023452156.1    PF00967
    ## 34262  XP_023452297.1    PF00262
    ## 35625  XP_023452464.1    PF07452
    ## 36776  XP_023452630.1    PF01738
    ## 37149  XP_023452672.1    PF05390
    ## 37727  XP_023452730.1    PF05922
    ## 39038  XP_023452887.1    PF02221
    ## 39801  XP_023452969.1    PF12454
    ## 40269  XP_023453030.1    PF04420
    ## 42386  XP_023453281.1    PF09362
    ## 45693  XP_023453648.1    PF00544
    ## 46735  XP_023453771.1    PF08212
    ## 49764  XP_023454096.1    PF16541
    ## 51798  XP_023454325.1    PF01042
    ## 52191  XP_023454369.1    PF13229
    ## 53770  XP_023454524.1    PF00684
    ## 53772  XP_023454524.1    PF01556
    ## 53918  XP_023454539.1    PF00657
    ## 58379  XP_023455052.1    PF03067
    ## 59189  XP_023455146.1    PF12999
    ## 59191  XP_023455146.1    PF13015
    ## 61321  XP_023455419.1    PF09769
    ## 62032  XP_023455500.1    PF01607
    ## 65203  XP_023455859.1    PF05630
    ## 67579  XP_023456152.1    PF21203
    ## 71956  XP_023456671.1    PF00026
    ## 72800  XP_023456780.1    PF07749
    ## 74930  XP_023457011.1    PF09451
    ## 75068  XP_023457026.1    PF03664
    ## 79112  XP_023457460.1    PF01423
    ## 80756  XP_023457650.1    PF00314
    ## 80836  XP_023457663.1    PF14040
    ## 81237  XP_023457699.1    PF03022
    ## 84230  XP_023458014.1    PF09430
    ## 86220  XP_023458213.1    PF03856
    ## 87421  XP_023458349.1    PF01095
    ## 89011  XP_023458521.1    PF02777
    ## 89015  XP_023458521.1    PF00081
    ## 90695  XP_023458704.1    PF14558
    ## 94066  XP_023459027.1    PF10528
    ## 97730  XP_023459380.1    PF07915
    ## 99443  XP_023459552.1    PF02678
    ## 102010 XP_023459829.1    PF05726
    ## 104544 XP_023460079.1    PF03537
    ## 104549 XP_023460080.1    PF10604
    ##                                                         functional
    ## 529                                             Alpha-L-fucosidase
    ## 556                    Alpha-L-arabinofuranosidase B (ABFB) domain
    ## 561                       Alpha-L-arabinofuranosidase B, catalytic
    ## 2167                               'Cold-shock' DNA-binding domain
    ## 3211                                               Cerato-platanin
    ## 3907                        Glycolipid 2-alpha-mannosyltransferase
    ## 8035       Uncharacterized protein conserved in bacteria (DUF2237)
    ## 9164                                             SnoaL-like domain
    ## 10384                             Haem degrading protein HbpS-like
    ## 11080                                                  Nis1 family
    ## 11626                                             Lipase (class 3)
    ## 12996     Cyclophilin type peptidyl-prolyl cis-trans isomerase/CLD
    ## 16829                                  Glycosyl hydrolase family 7
    ## 17838                                         Mpv17 / PMP22 family
    ## 18472                                             YqcI/YcgG family
    ## 21721                                       Ribonuclease T2 family
    ## 27287                         Family of unknown function (DUF6060)
    ## 30295                                       Zinc finger, C2H2 type
    ## 30304                                        C2H2-type zinc finger
    ## 32213   Cytidine and deoxycytidylate deaminase zinc-binding region
    ## 32744                                                Barwin family
    ## 34262                                          Calreticulin family
    ## 35625                                                  CHRD domain
    ## 36776                                Dienelactone hydrolase family
    ## 37149       Yeast cell wall synthesis protein KRE9/KNH1 C-terminal
    ## 37727                                       Peptidase inhibitor I9
    ## 39038                                                    ML domain
    ## 39801                  GPI-anchored cell wall organization protein
    ## 40269                                            CHD5-like protein
    ## 42386                         Domain of unknown function (DUF1996)
    ## 45693                                                Pectate lyase
    ## 46735                                        Lipocalin-like domain
    ## 49764                              Alternaria alternata allergen 1
    ## 51798                                       Endoribonuclease L-PSP
    ## 52191                               Right handed beta helix region
    ## 53770                                          DnaJ central domain
    ## 53772                                       DnaJ C terminal domain
    ## 53918                               GDSL-like Lipase/Acylhydrolase
    ## 58379     Lytic polysaccharide mono-oxygenase, cellulose-degrading
    ## 59189                             Glucosidase II beta subunit-like
    ## 59191                     Glucosidase II beta subunit-like protein
    ## 61321                                             Apolipoprotein O
    ## 62032                          Chitin binding Peritrophin-A domain
    ## 65203                             Necrosis inducing protein (NPP1)
    ## 67579                       ER membrane protein complex subunit 10
    ## 71956                                 Eukaryotic aspartyl protease
    ## 72800       Endoplasmic reticulum protein ERp29, C-terminal domain
    ## 74930                                 Autophagy-related protein 27
    ## 75068                                 Glycosyl hydrolase family 62
    ## 79112                                                   LSM domain
    ## 80756                                             Thaumatin family
    ## 80836                                  Deoxyribonuclease NucA/NucB
    ## 81237                                    Major royal jelly protein
    ## 84230  ER membrane protein complex subunit 7, beta-sandwich domain
    ## 86220                                Beta-glucosidase (SUN family)
    ## 87421                                               Pectinesterase
    ## 89011      Iron/manganese superoxide dismutases, C-terminal domain
    ## 89015   Iron/manganese superoxide dismutases, alpha-hairpin domain
    ## 90695                                               ML-like domain
    ## 94066                                                 GLEYA domain
    ## 97730                     Glucosidase II beta subunit-like protein
    ## 99443                                                        Pirin
    ## 102010                               Pirin C-terminal cupin domain
    ## 104544                            Glycoside-hydrolase family GH114
    ## 104549          Polyketide cyclase / dehydrase and lipid transport
    ##        signal_peptide effector                           Localizations
    ## 529               yes      yes                           Extracellular
    ## 556               yes      yes                           Extracellular
    ## 561               yes      yes                           Extracellular
    ## 2167              yes      yes                               Cytoplasm
    ## 3211              yes      yes                           Extracellular
    ## 3907              yes      yes                           Extracellular
    ## 8035              yes      yes                               Cytoplasm
    ## 9164              yes      yes                               Cytoplasm
    ## 10384             yes      yes                               Cytoplasm
    ## 11080             yes      yes                           Extracellular
    ## 11626             yes      yes                           Extracellular
    ## 12996             yes      yes         Cytoplasm|Endoplasmic reticulum
    ## 16829             yes      yes                           Extracellular
    ## 17838             yes      yes                Mitochondrion|Peroxisome
    ## 18472             yes      yes                               Cytoplasm
    ## 21721             yes      yes                           Extracellular
    ## 27287             yes      yes                           Extracellular
    ## 30295             yes      yes                       Cytoplasm|Nucleus
    ## 30304             yes      yes                       Cytoplasm|Nucleus
    ## 32213             yes      yes                       Cytoplasm|Nucleus
    ## 32744             yes      yes                               Cytoplasm
    ## 34262             yes      yes                   Endoplasmic reticulum
    ## 35625             yes      yes                           Extracellular
    ## 36776             yes      yes                               Cytoplasm
    ## 37149             yes      yes                           Extracellular
    ## 37727             yes      yes                           Extracellular
    ## 39038             yes      yes                           Extracellular
    ## 39801             yes      yes             Extracellular|Cell membrane
    ## 40269             yes      yes                   Endoplasmic reticulum
    ## 42386             yes      yes                           Extracellular
    ## 45693             yes      yes                           Extracellular
    ## 46735             yes      yes                           Extracellular
    ## 49764             yes      yes                           Extracellular
    ## 51798             yes      yes                               Cytoplasm
    ## 52191             yes      yes   Cytoplasm|Extracellular|Cell membrane
    ## 53770             yes      yes Cytoplasm|Nucleus|Endoplasmic reticulum
    ## 53772             yes      yes Cytoplasm|Nucleus|Endoplasmic reticulum
    ## 53918             yes      yes                           Extracellular
    ## 58379             yes      yes                           Extracellular
    ## 59189             yes      yes                   Endoplasmic reticulum
    ## 59191             yes      yes                   Endoplasmic reticulum
    ## 61321             yes      yes                           Mitochondrion
    ## 62032             yes      yes                           Extracellular
    ## 65203             yes      yes                           Extracellular
    ## 67579             yes      yes                   Endoplasmic reticulum
    ## 71956             yes      yes                           Extracellular
    ## 72800             yes      yes                   Endoplasmic reticulum
    ## 74930             yes      yes  Endoplasmic reticulum|Lysosome/Vacuole
    ## 75068             yes      yes                           Extracellular
    ## 79112             yes      yes                       Cytoplasm|Nucleus
    ## 80756             yes      yes                           Extracellular
    ## 80836             yes      yes                 Cytoplasm|Extracellular
    ## 81237             yes      yes             Extracellular|Mitochondrion
    ## 84230             yes      yes                   Endoplasmic reticulum
    ## 86220             yes      yes                           Extracellular
    ## 87421             yes      yes                           Extracellular
    ## 89011             yes      yes                 Cytoplasm|Mitochondrion
    ## 89015             yes      yes                 Cytoplasm|Mitochondrion
    ## 90695             yes      yes                           Extracellular
    ## 94066             yes      yes                           Extracellular
    ## 97730             yes      yes                   Endoplasmic reticulum
    ## 99443             yes      yes                           Mitochondrion
    ## 102010            yes      yes                              Peroxisome
    ## 104544            yes      yes          Extracellular|Lysosome/Vacuole
    ## 104549            yes      yes                               Cytoplasm
