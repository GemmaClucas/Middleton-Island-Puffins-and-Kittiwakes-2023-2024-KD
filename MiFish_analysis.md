MiFish_analysis
================
Gemma Clucas
2025-08-21

## 1. Set up

The data is saved on my solid state hard drive.

    cd /Users/gc547/Dropbox/GitHub_copied/Fecal_metabarcoding/Middleton-Island_2023-2024_KD/MiFish/
    conda activate qiime2-amplicon-2024.10 

## 2. Import the data into Qiime2

The data is saved across multiple plates: Plate75_2, Plate76, Plate99,
Plate100, Plate101, Plate102, and Plate109.

Plate 76 had both Katelyn’s samples and some storm petrel samples on it,
which are saved in separate folders, so I’ll have to import that one
separately.

And the same for plate 109, which was also a mixed plate.

    for K in 75_2 99 100 101 102; do
      qiime tools import\
        --type 'SampleData[PairedEndSequencesWithQuality]'\
        --input-path /Volumes/Data_SS1/MiFish/Plate$K/reads/ \
        --input-format CasavaOneEightSingleLanePerSampleDirFmt\
        --output-path demux_Plate$K.qza
    done

    qiime tools import\
        --type 'SampleData[PairedEndSequencesWithQuality]'\
        --input-path /Volumes/Data_SS1/MiFish/Plate76/MID_reads/ \
        --input-format CasavaOneEightSingleLanePerSampleDirFmt\
        --output-path demux_Plate76.qza
        
    qiime tools import\
        --type 'SampleData[PairedEndSequencesWithQuality]'\
        --input-path /Volumes/Data_SS1/MiFish/Plate109/MID_reads/ \
        --input-format CasavaOneEightSingleLanePerSampleDirFmt\
        --output-path demux_Plate109.qza

Make visualisation files:

    for K in 75_2 76 99 100 101 102 109; do
      qiime demux summarize \
        --i-data demux_plate$K.qza \
        --o-visualization demux_Plate$K.qzv
    done

## 3. Trim primers using cutadapt

The MiFish sequences are:

F primer: GTCGGTAAAACTCGTGCCAGC (21 bp)  
R primer: CATAGTGGGGTATCTAATCCCAGTTTG (27 bp)

### Trim 3’ ends first

At the 3’ end of the read, the primer will have been read through after
reading the MiFish region. I need to be looking for the reverse
complement of the reverse primer in read 1 (—p-adapter-f) and the
reverse complement of the forward primer in R2 (—p-adapter-r).

F primer reverse complement: GCTGGCACGAGTTTTACCGAC  
R primer reverse complement: CAAACTGGGATTAGATACCCCACTATG

    for K in 75_2 76 99 100 101 102 109; do
      qiime cutadapt trim-paired \
        --i-demultiplexed-sequences demux_plate$K.qza \
        --p-adapter-f CAAACTGGGATTAGATACCCCACTATG \
        --p-adapter-r GCTGGCACGAGTTTTACCGAC \
        --o-trimmed-sequences trimd_Plate$K.qza \
        --verbose > cutadapt_out_Plate$K.txt
    done

To see how much data passed the filter for each sample:

    for K in 75_2 76 99 100 101 102 109; do
      grep "Total written (filtered):" cutadapt_out_Plate$K.txt 
    done

Mostly 77-80%, which is normal.

### Trim 5’ ends of reads

All R1 should begin with the forward primer: GTCGGTAAAACTCGTGCCAGC (21
bases). All R2 should begin with the reverse primer:
CATAGTGGGGTATCTAATCCCAGTTTG (27 bases).

Trim these with the following commands:

    for K in 75_2 76 99 100 101 102 109; do
      qiime cutadapt trim-paired \
        --i-demultiplexed-sequences trimd_Plate$K.qza \
        --p-front-f GTCGGTAAAACTCGTGCCAGC \
        --p-front-r CATAGTGGGGTATCTAATCCCAGTTTG \
        --o-trimmed-sequences trimd2_Plate$K.qza \
        --verbose > cutadapt_out2_Plate$K.txt
    done

To see how much data passed the filter for each sample:

    for K in 75_2 76 99 100 101 102 109; do
      grep "Total written (filtered):" cutadapt_out2_Plate$K.txt 
    done

88% - good.

## 4. Denoise with dada2

Same settings as usual.

Note, this step can be a little slow to run.

    for K in 75_2 76 99 100 101 102 109; do
      qiime dada2 denoise-paired \
        --i-demultiplexed-seqs trimd2_Plate$K.qza \
        --p-trunc-len-f 133 \
        --p-trunc-len-r 138 \
        --p-trim-left-f 0 \
        --p-trim-left-r 0 \
        --p-min-overlap 50 \
        --p-n-threads 8 \
        --o-representative-sequences rep-seqs_Plate$K \
        --o-table table_Plate$K \
        --o-denoising-stats denoise_Plate$K
    done

Create visualizations for the denoising stats.

    for K in 75_2 76 99 100 101 102 109; do  
      qiime metadata tabulate\
        --m-input-file denoise_Plate$K.qza\
        --o-visualization denoise_Plate$K.qzv
    done

Look good, blanks have low read numbers.

## 5. Merge across plates

    qiime feature-table merge \
      --i-tables table_Plate75_2.qza \
      --i-tables table_Plate76.qza \
      --i-tables table_Plate99.qza \
      --i-tables table_Plate100.qza \
      --i-tables table_Plate101.qza \
      --i-tables table_Plate102.qza \
      --i-tables table_Plate109.qza \
      --p-overlap-method sum \
      --o-merged-table merged-table.qza

    qiime feature-table summarize \
        --i-table merged-table.qza \
        --m-sample-metadata-file metadata.txt \
        --o-visualization merged-table
        
    qiime feature-table merge-seqs \
      --i-data rep-seqs_Plate75_2.qza \
      --i-data rep-seqs_Plate76.qza \
      --i-data rep-seqs_Plate99.qza \
      --i-data rep-seqs_Plate100.qza \
      --i-data rep-seqs_Plate101.qza \
      --i-data rep-seqs_Plate102.qza \
      --i-data rep-seqs_Plate109.qza \
      --o-merged-data merged_rep-seqs.qza

    qiime feature-table tabulate-seqs \
      --i-data merged_rep-seqs.qza \
      --o-visualization merged_rep-seqs.qzv

## 6. Assign taxonomy

I will use the newest version of the database here. I copied it into
this folder from the eider and guillemot folder where I made it.

    ./mktaxa_singlethreaded.py \
      ncbi-refseqs-withHuman.qza \
      ncbi-taxonomy-withHuman.qza \
      merged_rep-seqs.qza

    qiime metadata tabulate \
      --m-input-file superblast_taxonomy.qza \
      --o-visualization superblast_taxonomy

## 7. Make some barplots

First look at what’s in there.

    qiime taxa barplot \
      --i-table merged-table.qza \
      --i-taxonomy superblast_taxonomy.qza \
      --m-metadata-file metadata.txt \
      --o-visualization barplot_before_filtering.qzv

## 8. Remove non-food reads

Filter out any sequences from the bird, mammals (human), bacteria,
unnassigned etc. sequences since we’re not interested in these.

    qiime taxa filter-table \
      --i-table merged-table.qza \
      --i-taxonomy superblast_taxonomy.qza \
      --p-exclude Unassigned,Aves,Mammalia,Bacteria,Eukaryota,Actinomycetota \
      --o-filtered-table merged_table_noBirdsMammalsUnassigned.qza
      
    qiime feature-table summarize \
        --i-table merged_table_noBirdsMammalsUnassigned.qza \
        --m-sample-metadata-file metadata.txt \
        --o-visualization merged_table_noBirdsMammalsUnassigned

Barplot having removed bird/human/unassigned DNA:

    qiime taxa barplot \
      --i-table merged_table_noBirdsMammalsUnassigned.qza \
      --i-taxonomy superblast_taxonomy.qza \
      --m-metadata-file metadata.txt \
      --o-visualization barplot_noBirdsMammalsUnassigned.qzv

## 9. Read into R to calculate depth of samples and blanks

### Read in feature table, taxonomy, and metadata

Going to do this on the onlymetazoa version of the data as I should
include birds and mammals as potential contaminants. Make sure that file
paths do not overwrite the MiFish ones.

``` r
# Read in the QIIME 2 artifacts and metadata
# Replace file paths with your actual file paths
feature_table <- read_qza("MiFish/merged_table_noBirdsMammalsUnassigned.qza")
taxonomy_data <- read_qza("MiFish/superblast_taxonomy.qza")
metadata <- read_q2metadata("MiFish/metadata.txt")

# Convert feature table to data frame
feature_table_df <- as.data.frame(feature_table$data) %>%
  rownames_to_column("FeatureID")

# Clean up taxonomy data
taxonomy_clean <- taxonomy_data$data %>%
  parse_taxonomy() %>%  # This splits the taxonomy string into columns
  rownames_to_column("FeatureID")

# Reshape feature table to long format
feature_long <- feature_table_df %>%
  gather(key = "SampleID", value = "Abundance", -FeatureID)

# Merge all data together
complete_data <- feature_long %>%
  left_join(taxonomy_clean, by = "FeatureID") %>%
  left_join(metadata, by = "SampleID")
```

### Calculate read depth of samples vs blanks

``` r
# Calculate average read depths by Type and Plate
read_depth_comparison <- complete_data %>%
  group_by(SampleID, Type, Plate) %>%
  summarise(TotalReads = sum(Abundance)) %>%
  group_by(Type, Plate) %>%
  summarise(
    MeanReads = round(mean(TotalReads), 2),
    MedianReads = round(median(TotalReads), 2),
    SDReads = round(sd(TotalReads), 2),
    n = n()
  )
```

    ## `summarise()` has grouped output by 'SampleID', 'Type'. You can override using
    ## the `.groups` argument.
    ## `summarise()` has grouped output by 'Type'. You can override using the
    ## `.groups` argument.

``` r
# Create empty list to store results for each plate
plate_results <- list()

# Process each plate separately
for(current_plate in unique(complete_data$Plate)) {
  # Get all expected types from metadata for this plate
  expected_types <- metadata %>%
    filter(Plate == current_plate) %>%
    pull(Type) %>%
    unique()
  
  # Get actual total reads for blanks and mock (this only sums real data)
  ext_blank_total_reads <- complete_data %>%
    filter(Plate == current_plate, Type == "EXTBLANK") %>%
    group_by(SampleID) %>%
    summarise(TotalReads = sum(Abundance)) %>%
    pull(TotalReads) %>%
    sum()
    
  fld_blank_total_reads <- complete_data %>%
    filter(Plate == current_plate, Type == "FLDBLANK") %>%
    group_by(SampleID) %>%
    summarise(TotalReads = sum(Abundance)) %>%
    pull(TotalReads) %>%
    sum()
    
  pcr_blank_total_reads <- complete_data %>%
    filter(Plate == current_plate, Type == "PCRBLANK") %>%
    group_by(SampleID) %>%
    summarise(TotalReads = sum(Abundance)) %>%
    pull(TotalReads) %>%
    sum()
    
  # Get true numbers from metadata
  true_ext_blank_number <- metadata %>% 
    filter(Plate == current_plate, Type == "EXTBLANK") %>% 
    nrow()
    
  true_fld_blank_number <- metadata %>% 
    filter(Plate == current_plate, Type == "FLDBLANK") %>% 
    nrow()
    
  true_pcr_blank_number <- metadata %>% 
    filter(Plate == current_plate, Type == "PCRBLANK") %>% 
    nrow()
  
  # Filter data for current plate
  plate_data <- read_depth_comparison %>% filter(Plate == current_plate)
  
  # Calculate adjusted averages using true totals and true numbers
  adjusted_ext_blank_avg <- ext_blank_total_reads / true_ext_blank_number
  adjusted_fld_blank_avg <- fld_blank_total_reads / true_fld_blank_number
  adjusted_pcr_blank_avg <- pcr_blank_total_reads / true_pcr_blank_number
  sample_avg <- plate_data$MeanReads[plate_data$Type == "SAMPLE"]
  
  # Create base comparison dataframe with all expected types
  adjusted_comparison <- data.frame(
    Type = expected_types,
    Plate = current_plate,
    MeanReads = 0,
    MedianReads = 0,
    SDReads = 0,
    n = 0,
    PercentOfSample = 0
  )
  
  # Update with actual data where available
  for(type in unique(plate_data$Type)) {
    idx <- adjusted_comparison$Type == type
    if(any(idx)) {
      adjusted_comparison[idx, "MeanReads"] <- plate_data$MeanReads[plate_data$Type == type]
      adjusted_comparison[idx, "MedianReads"] <- plate_data$MedianReads[plate_data$Type == type]
      adjusted_comparison[idx, "SDReads"] <- plate_data$SDReads[plate_data$Type == type]
      adjusted_comparison[idx, "n"] <- plate_data$n[plate_data$Type == type]
    }
  }
  
  # Update blank rows with adjusted means and true n
  if("EXTBLANK" %in% expected_types) {
    adjusted_comparison$MeanReads[adjusted_comparison$Type == "EXTBLANK"] <- round(adjusted_ext_blank_avg, 2)
    adjusted_comparison$n[adjusted_comparison$Type == "EXTBLANK"] <- true_ext_blank_number
  }
  if("FLDBLANK" %in% expected_types) {
    adjusted_comparison$MeanReads[adjusted_comparison$Type == "FLDBLANK"] <- round(adjusted_fld_blank_avg, 2)
    adjusted_comparison$n[adjusted_comparison$Type == "FLDBLANK"] <- true_fld_blank_number
  }
  if("PCRBLANK" %in% expected_types) {
    adjusted_comparison$MeanReads[adjusted_comparison$Type == "PCRBLANK"] <- round(adjusted_pcr_blank_avg, 2)
    adjusted_comparison$n[adjusted_comparison$Type == "PCRBLANK"] <- true_pcr_blank_number
  }
  
  # Calculate percentages for all types
  adjusted_comparison$PercentOfSample <- round((adjusted_comparison$MeanReads / sample_avg) * 100, 2)
  adjusted_comparison$PercentOfSample[adjusted_comparison$Type == "SAMPLE"] <- 100
  
  # Store results for this plate
  plate_results[[as.character(current_plate)]] <- adjusted_comparison
}


# Calculate combined statistics across all plates
all_plates_summary <- complete_data %>%
  group_by(SampleID, Type) %>%
  summarise(TotalReads = sum(Abundance)) %>%
  group_by(Type) %>%
  summarise(
    MeanReads = round(mean(TotalReads), 2),
    MedianReads = round(median(TotalReads), 2),
    SDReads = round(sd(TotalReads), 2),
    n = n()
  )
```

    ## `summarise()` has grouped output by 'SampleID'. You can override using the
    ## `.groups` argument.

``` r
# Get true numbers from metadata across all plates
true_ext_blank_number <- metadata %>% 
  filter(Type == "EXTBLANK") %>% 
  nrow()
  
true_fld_blank_number <- metadata %>% 
  filter(Type == "FLDBLANK") %>% 
  nrow()
  
true_pcr_blank_number <- metadata %>% 
  filter(Type == "PCRBLANK") %>% 
  nrow()

# Calculate total reads for each type across all plates
ext_blank_total_reads <- complete_data %>%
  filter(Type == "EXTBLANK") %>%
  group_by(SampleID) %>%
  summarise(TotalReads = sum(Abundance)) %>%
  pull(TotalReads) %>%
  sum()
  
fld_blank_total_reads <- complete_data %>%
  filter(Type == "FLDBLANK") %>%
  group_by(SampleID) %>%
  summarise(TotalReads = sum(Abundance)) %>%
  pull(TotalReads) %>%
  sum()
  
pcr_blank_total_reads <- complete_data %>%
  filter(Type == "PCRBLANK") %>%
  group_by(SampleID) %>%
  summarise(TotalReads = sum(Abundance)) %>%
  pull(TotalReads) %>%
  sum()

# Calculate adjusted averages
adjusted_ext_blank_avg <- ext_blank_total_reads / true_ext_blank_number
adjusted_fld_blank_avg <- fld_blank_total_reads / true_fld_blank_number
adjusted_pcr_blank_avg <- pcr_blank_total_reads / true_pcr_blank_number
sample_avg <- all_plates_summary$MeanReads[all_plates_summary$Type == "SAMPLE"]

# Create base comparison dataframe with all types
all_plates_comparison <- data.frame(
  Type = unique(metadata$Type),
  MeanReads = 0,
  MedianReads = 0,
  SDReads = 0,
  n = 0,
  PercentOfSample = 0
)

# Update with actual data where available
for(type in unique(all_plates_summary$Type)) {
  idx <- all_plates_comparison$Type == type
  if(any(idx)) {
    all_plates_comparison[idx, "MeanReads"] <- all_plates_summary$MeanReads[all_plates_summary$Type == type]
    all_plates_comparison[idx, "MedianReads"] <- all_plates_summary$MedianReads[all_plates_summary$Type == type]
    all_plates_comparison[idx, "SDReads"] <- all_plates_summary$SDReads[all_plates_summary$Type == type]
    all_plates_comparison[idx, "n"] <- all_plates_summary$n[all_plates_summary$Type == type]
  }
}

# Update blank rows with adjusted means and true n
all_plates_comparison$MeanReads[all_plates_comparison$Type == "EXTBLANK"] <- round(adjusted_ext_blank_avg, 2)
all_plates_comparison$n[all_plates_comparison$Type == "EXTBLANK"] <- true_ext_blank_number

all_plates_comparison$MeanReads[all_plates_comparison$Type == "FLDBLANK"] <- round(adjusted_fld_blank_avg, 2)
all_plates_comparison$n[all_plates_comparison$Type == "FLDBLANK"] <- true_fld_blank_number

all_plates_comparison$MeanReads[all_plates_comparison$Type == "PCRBLANK"] <- round(adjusted_pcr_blank_avg, 2)
all_plates_comparison$n[all_plates_comparison$Type == "PCRBLANK"] <- true_pcr_blank_number

# Calculate percentages for all types
all_plates_comparison$PercentOfSample <- round((all_plates_comparison$MeanReads / sample_avg) * 100, 2)
all_plates_comparison$PercentOfSample[all_plates_comparison$Type == "SAMPLE"] <- 100

# Display plate-specific results first
for(plate in names(plate_results)) {
  cat(sprintf("\n### Read Depth Summary - Plate %s\n", plate))
  print(knitr::kable(plate_results[[plate]], 
                     caption = sprintf("Summary of read depths by sample type - Plate %s", plate),
                     col.names = c("Type", "Plate", "Mean Reads", "Median Reads", "SD Reads", "n", "% of Sample Reads"),
                     align = c('l', 'l', 'r', 'r', 'r', 'r', 'r')))
  cat("\n")
}
```

    ## 
    ## ### Read Depth Summary - Plate 101
    ## 
    ## 
    ## Table: Summary of read depths by sample type - Plate 101
    ## 
    ## |Type     |Plate | Mean Reads| Median Reads| SD Reads|  n| % of Sample Reads|
    ## |:--------|:-----|----------:|------------:|--------:|--:|-----------------:|
    ## |MOCK     |101   |  120662.00|       120662|       NA|  1|            119.03|
    ## |PCRBLANK |101   |       0.00|            0|     0.00|  2|              0.00|
    ## |SAMPLE   |101   |  101372.71|        89519| 66468.92| 85|            100.00|
    ## |EXTBLANK |101   |       0.33|            2|       NA|  6|              0.00|
    ## |FLDBLANK |101   |   30901.00|        30901|       NA|  1|             30.48|
    ## 
    ## 
    ## ### Read Depth Summary - Plate 75
    ## 
    ## 
    ## Table: Summary of read depths by sample type - Plate 75
    ## 
    ## |Type     |Plate | Mean Reads| Median Reads|  SD Reads|  n| % of Sample Reads|
    ## |:--------|:-----|----------:|------------:|---------:|--:|-----------------:|
    ## |MOCK     |75    |  339103.00|       339103|        NA|  1|            256.68|
    ## |PCRBLANK |75    |       0.00|            0|      0.00|  2|              0.00|
    ## |SAMPLE   |75    |  132110.25|       112035| 104761.16| 68|            100.00|
    ## |EXTBLANK |75    |     104.83|           25|    327.68|  6|              0.08|
    ## |FLDBLANK |75    |    1884.33|         1713|    455.83|  3|              1.43|
    ## 
    ## 
    ## ### Read Depth Summary - Plate 99
    ## 
    ## 
    ## Table: Summary of read depths by sample type - Plate 99
    ## 
    ## |Type     |Plate | Mean Reads| Median Reads| SD Reads|  n| % of Sample Reads|
    ## |:--------|:-----|----------:|------------:|--------:|--:|-----------------:|
    ## |MOCK     |99    |   66816.00|      66816.0|       NA|  1|            135.75|
    ## |PCRBLANK |99    |       0.00|          0.0|     0.00|  3|              0.00|
    ## |SAMPLE   |99    |   49220.06|       8918.5| 87641.52| 72|            100.00|
    ## |EXTBLANK |99    |       0.00|          0.0|     0.00|  6|              0.00|
    ## |FLDBLANK |99    |    2143.00|       3214.5|  2345.47|  3|              4.35|
    ## 
    ## 
    ## ### Read Depth Summary - Plate 100
    ## 
    ## 
    ## Table: Summary of read depths by sample type - Plate 100
    ## 
    ## |Type     |Plate | Mean Reads| Median Reads| SD Reads|  n| % of Sample Reads|
    ## |:--------|:-----|----------:|------------:|--------:|--:|-----------------:|
    ## |MOCK     |100   |   130655.0|       130655|       NA|  1|            109.26|
    ## |PCRBLANK |100   |        0.0|            0|      0.0|  2|              0.00|
    ## |SAMPLE   |100   |   119578.2|       106767| 100882.1| 73|            100.00|
    ## |EXTBLANK |100   |        0.0|            0|      0.0|  6|              0.00|
    ## |FLDBLANK |100   |     1092.5|         2185|       NA|  2|              0.91|
    ## 
    ## 
    ## ### Read Depth Summary - Plate 102
    ## 
    ## 
    ## Table: Summary of read depths by sample type - Plate 102
    ## 
    ## |Type     |Plate | Mean Reads| Median Reads| SD Reads|  n| % of Sample Reads|
    ## |:--------|:-----|----------:|------------:|--------:|--:|-----------------:|
    ## |MOCK     |102   |   386488.0|       386488|       NA|  1|            199.67|
    ## |PCRBLANK |102   |        0.0|            0|      0.0|  1|              0.00|
    ## |SAMPLE   |102   |   193564.1|       141453| 174110.9| 34|            100.00|
    ## |EXTBLANK |102   |        0.0|            0|      0.0|  2|              0.00|
    ## 
    ## 
    ## ### Read Depth Summary - Plate 76
    ## 
    ## 
    ## Table: Summary of read depths by sample type - Plate 76
    ## 
    ## |Type     |Plate | Mean Reads| Median Reads| SD Reads|  n| % of Sample Reads|
    ## |:--------|:-----|----------:|------------:|--------:|--:|-----------------:|
    ## |MOCK     |76    |   190818.0|     190818.0|       NA|  1|            130.94|
    ## |PCRBLANK |76    |        0.0|          0.0|      0.0|  2|              0.00|
    ## |SAMPLE   |76    |   145724.5|     143198.5| 102663.7| 28|            100.00|
    ## |EXTBLANK |76    |        0.0|          0.0|      0.0|  2|              0.00|
    ## 
    ## 
    ## ### Read Depth Summary - Plate 109
    ## 
    ## 
    ## Table: Summary of read depths by sample type - Plate 109
    ## 
    ## |Type   |Plate | Mean Reads| Median Reads| SD Reads|  n| % of Sample Reads|
    ## |:------|:-----|----------:|------------:|--------:|--:|-----------------:|
    ## |SAMPLE |109   |   12443.38|        13300|  4615.72| 13|               100|

``` r
# Display combined results
cat("\n### Read Depth Summary - All Plates Combined\n")
```

    ## 
    ## ### Read Depth Summary - All Plates Combined

``` r
print(knitr::kable(all_plates_comparison,
                   caption = "Summary of read depths by sample type - All Plates Combined",
                   col.names = c("Type", "Mean Reads", "Median Reads", "SD Reads", "n", "% of Sample Reads"),
                   align = c('l', 'r', 'r', 'r', 'r', 'r')))
```

    ## 
    ## 
    ## Table: Summary of read depths by sample type - All Plates Combined
    ## 
    ## |Type     | Mean Reads| Median Reads|  SD Reads|   n| % of Sample Reads|
    ## |:--------|----------:|------------:|---------:|---:|-----------------:|
    ## |MOCK     |  205757.00|     160736.5| 128727.69|   6|            188.58|
    ## |PCRBLANK |       0.00|          0.0|      0.00|  12|              0.00|
    ## |SAMPLE   |  109105.79|      84203.0| 108491.56| 373|            100.00|
    ## |EXTBLANK |      21.03|         20.5|    286.99|  30|              0.02|
    ## |FLDBLANK |    5018.67|       2185.0|  10843.06|   9|              4.60|

``` r
# Save all results if requested
if(knitr::opts_chunk$get("save_output")) {
  for(plate in names(plate_results)) {
    write.csv(plate_results[[plate]],
              sprintf("MiFish/adjusted_read_depth_comparison_plate_%s.csv", plate),
              row.names = FALSE)
  }
  write.csv(all_plates_comparison,
            "MiFish/adjusted_read_depth_comparison_all_plates.csv",
            row.names = FALSE)
}
```

The field blank on plate 101 has quite a few reads, but overall the
field blank depth is 4.6% so it’s fine.

## 13. Calculate alpha rarefaction curves

Must be done on collapsed table (because we don’t care about ASV
diversity, just species diversity).

    qiime taxa collapse \
      --i-table merged_table_noBirdsMammalsUnassigned.qza \
      --i-taxonomy superblast_taxonomy.qza \
      --p-level 7 \
      --o-collapsed-table merged_table_noBirdsMammalsUnassigned_collapsed.qza

    qiime diversity alpha-rarefaction \
      --i-table merged_table_noBirdsMammalsUnassigned_collapsed.qza \
      --m-metadata-file metadata.txt \
      --p-min-depth 100 \
      --p-max-depth 5000 \
      --o-visualization alpha-rarefaction-100-5000

Repeat up to 1000

    qiime diversity alpha-rarefaction \
      --i-table merged_table_noBirdsMammalsUnassigned_collapsed.qza \
      --m-metadata-file metadata.txt \
      --p-min-depth 100 \
      --p-max-depth 1000 \
      --o-visualization alpha-rarefaction-100-1000

Species richness plateaus at 300 reads, so we should drop samples with
fewer than 300.

## 14. Sample filtering

#### Drop mock community and blanks and samples with fewer than 300 reads

    qiime feature-table filter-samples \
      --i-table merged_table_noBirdsMammalsUnassigned.qza \
      --p-min-frequency 300 \
      --m-metadata-file metadata.txt \
      --p-where "Type='SAMPLE'" \
      --o-filtered-table merged_table_noBirdsMammalsUnassigned_minfreq300
      
    qiime taxa barplot \
      --i-table merged_table_noBirdsMammalsUnassigned_minfreq300.qza  \
      --m-metadata-file metadata.txt \
      --i-taxonomy superblast_taxonomy.qza \
      --o-visualization barplot_noBirdsMammalsUnassigned_minfreq300

## 16. Abundance filtering

``` r
library(qiime2R)
library(tidyverse)
library(biomformat)

# Read in the QIIME 2 artifacts and metadata
feature_table <- read_qza("MiFish/merged_table_noBirdsMammalsUnassigned_minfreq300.qza")
taxonomy_data <- read_qza("MiFish/superblast_taxonomy.qza")
metadata <- read_q2metadata("MiFish/metadata.txt")

# Convert feature table to data frame
feature_table_df <- as.data.frame(feature_table$data) %>%
  rownames_to_column("FeatureID")

# Extract taxonomy information
taxonomy_df <- data.frame(
  FeatureID = taxonomy_data$data[,1],
  Taxonomy = taxonomy_data$data[,2],
  stringsAsFactors = FALSE
)


# Reshape feature table to long format
feature_long <- feature_table_df %>%
  gather(key = "SampleID", value = "Abundance", -FeatureID) %>%
  # Join with taxonomy
  left_join(taxonomy_df, by = "FeatureID")

# Aggregate abundances by taxonomy and sample
species_level_abundances <- feature_long %>%
  group_by(SampleID, Taxonomy) %>%
  summarise(Abundance = sum(Abundance)) %>%
  ungroup()
```

    ## `summarise()` has grouped output by 'SampleID'. You can override using the
    ## `.groups` argument.

``` r
# Apply 1% filter per sample at species level
filtered_species <- species_level_abundances %>%
  group_by(SampleID) %>%
  mutate(
    TotalReads = sum(Abundance),
    RelativeAbundance = Abundance / TotalReads,
    # Set abundance to 0 if relative abundance is < 1%
    Abundance = if_else(RelativeAbundance < 0.01, 0, Abundance)
  ) %>%
  ungroup()

# Get the list of taxonomy strings that passed the filter in any sample
retained_taxa <- filtered_species %>%
  filter(Abundance > 0) %>%
  pull(Taxonomy) %>%
  unique()

# Filter the original feature-level data based on taxonomy
filtered_features <- feature_long %>%
  filter(Taxonomy %in% retained_taxa) %>%
  select(FeatureID, SampleID, Abundance, Taxonomy)

# Convert back to wide format for QIIME2
filtered_feature_table <- filtered_features %>%
  select(FeatureID, SampleID, Abundance) %>%
  pivot_wider(
    names_from = SampleID,
    values_from = Abundance,
    values_fill = 0
  )

# Create output files for QIIME2
if(knitr::opts_chunk$get("save_output")) {
  # Create feature table in BIOM format
  filtered_feature_matrix <- as.matrix(filtered_feature_table[,-1])
  rownames(filtered_feature_matrix) <- filtered_feature_table$FeatureID
  
  # Create BIOM object
  biom_obj <- make_biom(data = filtered_feature_matrix)
  
  # Write BIOM file
  write_biom(biom_obj, "MiFish/filtered_feature_table.biom")
  
  # Write taxonomy file with correct header
  filtered_taxonomy <- taxonomy_df %>%
    filter(FeatureID %in% filtered_feature_table$FeatureID)
  
  # Create the data without headers first
  filtered_taxonomy_output <- data.frame(
    FeatureID = filtered_taxonomy$FeatureID,
    Taxonomy = filtered_taxonomy$Taxonomy
  )
  
  # Write to file with explicit column names
  write.table(filtered_taxonomy_output, 
              "MiFish/filtered_taxonomy.tsv", 
              sep="\t", 
              quote=FALSE, 
              row.names=FALSE,
              col.names = c("Feature ID", "Taxon"))
}

# Print summary statistics
cat("Original number of features:", nrow(feature_table_df), "\n")
```

    ## Original number of features: 2248

``` r
cat("Number of features after filtering:", nrow(filtered_feature_table), "\n")
```

    ## Number of features after filtering: 2134

``` r
cat("Number of features removed:", nrow(feature_table_df) - nrow(filtered_feature_table), "\n")
```

    ## Number of features removed: 114

``` r
# Calculate and print the number of unique taxa before and after filtering
n_taxa_before <- length(unique(taxonomy_df$Taxonomy))
n_taxa_after <- length(retained_taxa)
cat("\nNumber of unique taxa before filtering:", n_taxa_before, "\n")
```

    ## 
    ## Number of unique taxa before filtering: 140

``` r
cat("Number of unique taxa after filtering:", n_taxa_after, "\n")
```

    ## Number of unique taxa after filtering: 48

``` r
cat("Number of taxa removed:", n_taxa_before - n_taxa_after, "\n")
```

    ## Number of taxa removed: 92

``` r
# Calculate and print the number of reads before and after filtering
total_reads_before <- sum(feature_long$Abundance)
total_reads_after <- sum(filtered_features$Abundance)
cat("\nTotal reads before filtering:", total_reads_before, "\n")
```

    ## 
    ## Total reads before filtering: 40694710

``` r
cat("Total reads after filtering:", total_reads_after, "\n")
```

    ## Total reads after filtering: 40675573

``` r
cat("Percentage of reads retained:", round(total_reads_after/total_reads_before * 100, 2), "%\n")
```

    ## Percentage of reads retained: 99.95 %

### Import back into qiime

    # Import feature table
    qiime tools import \
      --input-path filtered_feature_table.biom \
      --type 'FeatureTable[Frequency]' \
      --input-format BIOMV100Format \
      --output-path filtered_table_minfreq300_minabund1.qza

    # Import taxonomy
    qiime tools import \
      --input-path filtered_taxonomy.tsv \
      --type 'FeatureData[Taxonomy]' \
      --input-format TSVTaxonomyFormat \
      --output-path filtered_taxonomy_minfreq300_minabund1.qza

Make a barplot to see how it looks now

    qiime taxa barplot \
      --i-table filtered_table_minfreq300_minabund1.qza \
      --i-taxonomy filtered_taxonomy_minfreq300_minabund1.qza \
      --m-metadata-file metadata.txt \
      --o-visualization barplot_minfreq300_minabund1.qzv

## 15. Taxonomy edits

First, filter down the rep-seqs to just those that are found in our
remaining samples after abundance filtering.

    qiime feature-table filter-seqs \
      --i-data merged_rep-seqs.qza \
      --i-table filtered_table_minfreq300_minabund1.qza \
      --o-filtered-data filtered_rep-seqs_minfreq300_minabund1.qza

Next make a table with sequences and taxonomy strings:

    qiime metadata tabulate \
      --m-input-file filtered_rep-seqs_minfreq300_minabund1.qza \
      --m-input-file filtered_taxonomy_minfreq300_minabund1.qza \
      --o-visualization sequence_taxonomy_minfreq300_minabund1.qzv

Make the edits needed after reading the taxonomy artifact into R.

``` r
tax <- read_qza("MiFish/superblast_taxonomy.qza")
```

Changes that I made for the first set of samples: \* Clupea harengus to
Clupea pallasii.  
\* Sprattus sprattus to Clupea pallasii.  
\* Gadus chalcogrammus to Gadus sp. (but most likely walleye pollock).  
\* Microgradus proximus to Gadidae sp. (but most likely Pacific
tomcod).  
\* Stenobrachius leucopsarus to Stenobrachius sp.  
\* Stenobrachius nannochir to Stenobrachius sp.  
\* Spirinchus lanceolatus to Thaleichthys pacificus.  
\* Hexagrammos agrammus to Hexagrammos sp.  
\* Hexagrammos octogrammus to Hexagrammos sp.  
\* Sebastes babcocki to Sebastes sp.  
\* Sebastes crameri to Sebastes sp.  
\* Sebastes mentella to Sebastes sp.  
\* Anoplarchus insignis to Anoplarchus sp.  
\* Xiphister mucosus to Stichaeidae sp.  
\* Oncorhynchus gorbuscha to Oncorhynchus sp.  
\* Oncorhynchus keta to Oncorhynchus sp.  
\* Oncorhynchus kisutch to Oncorhynchus sp.  
\* Oncorhynchus mykiss to Oncorhynchus sp.  
\* Oncorhynchus nerka to Oncorhynchus sp.  
\* Lycodapus mandibularis to Lycodapus sp.  
\* Lycodapus microchir to Lycodapus sp.  
\* Ammodytes personatus to A. dubius.

Additional changes for the second set of samples:  
\* Poromitra cristiceps to Poromitra crassiceps based on range. \* Gadus
macrocephalus to Gadus sp.  
\* Gadus ogac to Gadus sp.  
\* Lampanyctus tenuiformis to Lampanyctus sp. (recent genus name
changes). \* Nannobrachium regale to Lampanyctus sp. (recent genus name
changes). \* Ruscarius creaseri to Cottidae (no close match). \*
Cryptacanthodes aleutensis to Cryptacanthodes sp.  
\* Polypera greeni to Liparis greeni (name change).  
\* Sebastes miniatus to Sebeastes sp.

Katelyn finished checking these and removed those that we don’t need
anymore (after applying the 1% filter first this time around).

``` r
tax$data$Taxon <- tax$data$Taxon %>%
  str_replace_all("g__Poromitra;s__cristiceps", "g__Poromitra;s__crassiceps") %>%
  str_replace_all("g__Clupea;s__harengus", "g__Clupea;s__pallasii") %>%
  str_replace_all("g__Sprattus;s__sprattus", "g__Clupea;s__pallasii") %>%
  str_replace_all("g__Gadus;s__chalcogrammus", "g__Gadus;s__") %>%
  str_replace_all("g__Gadus;s__macrocephalus", "g__Gadus;s__") %>%
  str_replace_all("g__Gadus;s__ogac", "g__Gadus;s__") %>%
  str_replace_all("g__Microgadus;s__proximus", "g__;s__") %>%
  str_replace_all("g__Lampanyctus;s__tenuiformis", "g__Lampanyctus;s__") %>%
  str_replace_all("g__Nannobrachium;s__regale", "g__Lampanyctus;s__") %>%
  str_replace_all("g__Stenobrachius;s__leucopsarus", "g__Stenobrachius;s__") %>%
  str_replace_all("g__Stenobrachius;s__nannochir", "g__Stenobrachius;s__") %>%
  str_replace_all("g__Ruscarius;s__creaseri", "g__;s__") %>%
  str_replace_all("g__Cryptacanthodes;s__aleutensis", "g__Cryptacanthodes;s__") %>%
  str_replace_all("g__Hexagrammos;s__agrammus", "g__Hexagrammos;s__") %>%
  str_replace_all("g__Hexagrammos;s__decagrammus", "g__Hexagrammos;s__") %>%
  str_replace_all("g__Polypera;s__greeni", "g__Liparis;s__greeni") %>%
  str_replace_all("g__Sebastes;s__miniatus", "g__Sebastes;s__") %>%
  str_replace_all("g__Anoplarchus;s__insignis", "g__Anoplarchus;s__") %>%
  str_replace_all("g__Lycodapus;s__mandibularis", "g__Lycodapus;s__") %>%
  str_replace_all("g__Lycodapus;s__microchir", "g__Lycodapus;s__") %>%
  str_replace_all("g__Atheresthes;s__evermanni", "g__Atheresthes;s__stomias") %>%
  str_replace_all("g__Glyptocephalus;s__stelleri", "g__;s__") %>%
  str_replace_all("g__Oncorhynchus;s__gorbuscha", "g__Oncorhynchus;s__") %>%
  str_replace_all("g__Oncorhynchus;s__kisutch", "g__Oncorhynchus;s__") %>%
  str_replace_all("g__Oncorhynchus;s__nerka", "g__Oncorhynchus;s__") %>%
  str_replace_all("g__Ammodytes;s__japonicus", "g__Ammodytes;s__personatus") %>%
  str_replace_all("g__Hyperoplus;s__immaculatus", "g__Ammodytes;s__personatus")
```

To check that it worked, you can search for the old or new names using
grepl:

``` r
tax$data %>% filter(grepl("Ammodytes;s__dubius", Taxon)) %>% count()
```

    ##   n
    ## 1 1

``` r
tax$data %>% filter(grepl("Ammodytes;s__personatus", Taxon)) %>% count()
```

    ##    n
    ## 1 72

``` r
tax$data %>% filter(grepl("g__Oncorhynchus;s__gorbuscha|g__Oncorhynchus; s__keta|g__Oncorhynchus; s__kisutch", Taxon)) %>% count()
```

    ##   n
    ## 1 0

``` r
tax$data %>% filter(grepl("g__Oncorhynchus;s__", Taxon)) %>%  count()
```

    ##    n
    ## 1 24

``` r
tax$data %>% filter(grepl("g__Clupea;s__pallasii", Taxon)) %>%  count()
```

    ##     n
    ## 1 557

Export to check (commented out for knitting).

``` r
# write.table(tax$data,
#           quote = FALSE,
#           row.names = FALSE,
#           file = "MiFish/superblast_taxonomy_edited.tsv",
#           sep = "\t")
```

Remove the period and reload into qiime (run in terminal):

    sed -i.bak 's/Feature.ID/Feature ID/g' superblast_taxonomy_edited.tsv

    conda activate qiime2-amplicon-2024.10 

    qiime tools import \
      --input-path superblast_taxonomy_edited.tsv \
      --output-path superblast_taxonomy_edited.qza \
      --type 'FeatureData[Taxonomy]'
      
    qiime metadata tabulate \
      --m-input-file superblast_taxonomy_edited.qza \
      --o-visualization superblast_taxonomy_edited

## Final barplot

    qiime taxa barplot \
      --i-table filtered_table_minfreq300_minabund1.qza \
      --i-taxonomy superblast_taxonomy_edited.qza \
      --m-metadata-file metadata.txt \
      --o-visualization barplot_minfreq300_minabund1_taxedited.qzv

Downloaded this final dataset as a csv from the qiime viewer and this is
what Katelyn is working with.
