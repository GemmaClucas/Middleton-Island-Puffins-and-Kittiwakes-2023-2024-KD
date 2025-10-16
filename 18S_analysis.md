18S_commands_2023&2024
================
Gemma Clucas
2025-08-29

## 1. Import the data

Just plate 75 for now. Will need to update to include all plates.

    cd /Users/gc547/Dropbox/GitHub_copied/Fecal_metabarcoding/Middleton-Island_2023-2024_KD/18S
    conda activate qiime2-amplicon-2024.10 

    for K in 75 76 99 100 101; do
    qiime tools import\
      --type 'SampleData[PairedEndSequencesWithQuality]'\
      --input-path /Volumes/Data_SS1/18S/Plate$K/reads/ \
      --input-format CasavaOneEightSingleLanePerSampleDirFmt\
      --output-path demux_Plate$K.qza
      
    for K in 75 76 99 100 101; do
      qiime demux summarize \
        --i-data demux_Plate$K.qza \
        --o-visualization demux_Plate$K.qzv
    done

## 2. Trim primers using cutadapt

The 18S sequences are:

F primer: GGTCTGTGATGCCCTTAGATG (21 bp)  
R primer: GGTGTGTACAAAGGGCAGGG (20 bp)

### Trim 3’ ends first

At the 3’ end of the read, the primer will have been read through after
reading the 18S region. I need to be looking for the reverse complement
of the reverse primer in read 1 (—p-adapter-f) and the reverse
complement of the forward primer in R2 (—p-adapter-r).

F primer reverse complement: CATCTAAGGGCATCACAGACC  
R primer reverse complement: CCCTGCCCTTTGTACACACC

    for K in 75 76 99 100 101; do
      qiime cutadapt trim-paired \
        --i-demultiplexed-sequences demux_Plate$K.qza \
        --p-adapter-f CCCTGCCCTTTGTACACACC \
        --p-adapter-r CATCTAAGGGCATCACAGACC \
        --o-trimmed-sequences trimd_Plate$K.qza \
        --verbose > cutadapt_out_Plate$K.txt
    done

To see how much passed the filters:

    for K in 75 76 99 100 101; do
      grep "Total written (filtered):" cutadapt_out_Plate$K.txt 
    done

This is very consistent with 74 - 77% passing the filters, which is
typical for this marker.

### Trim 5’ ends of reads

All R1 should begin with the forward primer: GGTCTGTGATGCCCTTAGATG (21
bases). All R2 should begin with the reverse primer:
GGTGTGTACAAAGGGCAGGG (20 bases).

Trim these with the following commands:

    for K in 75 76 99 100 101; do
      qiime cutadapt trim-paired \
        --i-demultiplexed-sequences trimd_Plate$K.qza \
        --p-front-f GGTCTGTGATGCCCTTAGATG \
        --p-front-r GGTGTGTACAAAGGGCAGGG \
        --o-trimmed-sequences trimd2_Plate$K.qza \
        --verbose > cutadapt_out2_Plate$K.txt
    done

To see how much passed the filters:

    for K in 75 76 99 100 101; do
      grep "Total written (filtered):" cutadapt_out2_Plate$K.txt 
    done

89% passes here, also normal.

## 3. Denoise with Dada2

    for K in 75 76 99 100 101; do
      qiime dada2 denoise-paired \
        --i-demultiplexed-seqs trimd2_Plate$K.qza \
        --p-trunc-len-f 150 \
        --p-trunc-len-r 150 \
        --p-trim-left-f 0 \
        --p-trim-left-r 0 \
        --p-min-overlap 50 \
        --p-n-threads 8 \
        --o-representative-sequences rep-seqs_Plate$K \
        --o-table table_Plate$K \
        --o-denoising-stats denoise_Plate$K
    done

Create visualizations for the denoising stats. 76 99 100 101

    for K in 75; do  
      qiime metadata tabulate\
        --m-input-file denoise_Plate$K.qza\
        --o-visualization denoise_Plate$K.qzv
    done

## 4. Merge across plates

    qiime feature-table merge \
      --i-tables table_Plate75.qza \
      --i-tables table_Plate76.qza \
      --i-tables table_Plate99.qza \
      --i-tables table_Plate100.qza \
      --i-tables table_Plate101.qza \
      --p-overlap-method sum \
      --o-merged-table merged-table.qza

    qiime feature-table summarize \
        --i-table merged-table.qza \
        --m-sample-metadata-file metadata.txt \
        --o-visualization merged-table
        
    qiime feature-table merge-seqs \
      --i-data rep-seqs_Plate75.qza \
      --i-data rep-seqs_Plate76.qza \
      --i-data rep-seqs_Plate99.qza \
      --i-data rep-seqs_Plate100.qza \
      --i-data rep-seqs_Plate101.qza \
      --o-merged-data merged_rep-seqs.qza

    qiime feature-table tabulate-seqs \
      --i-data merged_rep-seqs.qza \
      --o-visualization merged_rep-seqs.qzv

## 5. Assign taxonomy using Naive Bayes classifier

I am going to use the same classifier that I trained for Will’s puffin
samples. I just copied it from the storm petrel folder into here. Notes
on how I trained it are in the repo for his analysis.

I copied the metadata file over from the MiFish folder and changed
underscores to hyphens.

    qiime feature-classifier classify-sklearn \
      --i-classifier classifier.qza \
      --i-reads rep-seqs_Plate75.qza \
      --o-classification sklearn_taxonomy.qza
      
    qiime taxa barplot\
       --i-table table_Plate75.qza\
       --i-taxonomy sklearn_taxonomy.qza\
       --m-metadata-file metadata.txt\
       --o-visualization barplot_sklearn-taxa.qzv
       
       
    # merged version
    qiime feature-classifier classify-sklearn \
      --i-classifier classifier.qza \
      --i-reads merged_rep-seqs.qza \
      --o-classification sklearn_taxonomy.qza
      
    qiime taxa barplot\
       --i-table merged_table.qza\
       --i-taxonomy sklearn_taxonomy.qza\
       --m-metadata-file metadata.txt\
       --o-visualization barplot_sklearn-taxa.qzv

The blanks seems to mostly have fungal and a load of other
micro-organismal DNA in them, so that’s not too bad. One has some fish.

## 5. Filter out non-prey sequences

    qiime taxa filter-table \
      --i-table table_Plate75.qza \
      --i-taxonomy sklearn_taxonomy.qza \
      --p-include Metazoa \
      --o-filtered-table merged-table_onlymetazoa.qza
      
    qiime taxa barplot\
          --i-table merged-table_onlymetazoa.qza\
          --i-taxonomy sklearn_taxonomy.qza\
          --m-metadata-file metadata.txt\
          --o-visualization barplot_sklearn-taxa_onlymetazoa

Need to filter for just prey next with Katelyn. I downloaded the level 6
version for her to work on. This is what we decided to keep:

    qiime taxa filter-table \
      --i-table merged-table_onlymetazoa.qza \
      --i-taxonomy sklearn_taxonomy.qza \
      --p-include Neopterygii,Hydroidolina,Phyllodocida,Copepoda,Decapodiformes,Eumalacostraca,Thecostraca,Gastropoda \
      --o-filtered-table merged-table_onlyprey.qza
      
    qiime taxa barplot\
          --i-table merged-table_onlyprey.qza \
          --i-taxonomy sklearn_taxonomy.qza\
          --m-metadata-file metadata.txt\
          --o-visualization barplot_onlyprey 

Need to download and check before sending to Katelyn.

Out of all 11 blanks on plate 75 (8 lab blanks, 3 field) there were 4678
reads (found in 2 lab and 1 field blanks). Average depth of blanks was
425 reads. The average depth of the samples was 26,214 reads. Therefore
the depths of the blanks was 1.62% the depth of the samples.

## 6. Alpha rarefaction

    qiime taxa collapse \
      --i-table merged-table_onlyprey.qza \
      --i-taxonomy sklearn_taxonomy.qza \
      --p-level 6 \
      --o-collapsed-table merged-table_onlyprey_collapsed.qza

    qiime diversity alpha-rarefaction \
      --i-table merged-table_onlyprey_collapsed.qza \
      --m-metadata-file metadata.txt \
      --p-min-depth 100 \
      --p-max-depth 10000 \
      --o-visualization alpha-rarefaction-100-10000


    qiime diversity alpha-rarefaction \
      --i-table merged-table_onlyprey_collapsed.qza \
      --m-metadata-file metadata.txt \
      --p-min-depth 50 \
      --p-max-depth 500 \
      --o-visualization alpha-rarefaction-50-500

Even at 50 reads we are capturing the diversity in these samples, so I
will drop two that had fewer than 50 reads (P53 and P59).

## 7. Final filtering

Drop the mock and blanks.

    qiime feature-table filter-samples \
      --i-table merged-table_onlyprey.qza \
      --m-metadata-file metadata.txt \
      --p-where "Species='BLKI'" \
      --o-filtered-table BLKI_table_onlyprey.qza 

Retain only features with at least 1% abundance in at least 1 sample (we
have 92 samples so 1/74 = 0.0135).

    qiime feature-table filter-features-conditionally \
      --i-table BLKI_table_onlyprey.qza \
      --p-abundance 0.01 \
      --p-prevalence 0.01 \
      --o-filtered-table BLKI_table_onlyprey_filtered.qza
      
    qiime taxa barplot \
      --i-table BLKI_table_onlyprey_filtered.qza  \
      --m-metadata-file metadata.txt \
      --i-taxonomy sklearn_taxonomy.qza \
      --o-visualization barplot_BLKI_onlyprey_filtered

So the file I am sending to Katelyn with the final data for 18S is
called `level6_sklearn-taxa_onlyprey_filtered.csv`. I deleted P53 and
P59 from this as it was easiest to do this by hand rather than create
another file in qiime.
