## Create Merged Sample Table ##
# Modeled on code from the Therkildsen Lab at Cornell

# Merge samples in sample table

library(tidyverse)

## Define base directory and reference name
#basedir <- "/projects/f_mlp195/kyraf/croaker/lcwgs/"
refname <- "hifi"

## Read in unmerged sample tables and combine pe and se
sample_table<-read_tsv("Sample_Table_Combined.tsv") %>%
  mutate(sample_seq_id=paste(sample_id,seq_id,lane_number, sep = "_"))

## Add a sample_id_corrected column, just in case that the same sample got assigned different IDs, or a few of the IDs are wrong
# This is not the case for the Greenland cod project
# When this happends, just edit "wrong_id" and "correct_id"
sample_table <- mutate(sample_table, sample_id_corrected=ifelse(sample_id!="wrong_id", sample_id, "correct_id"))

## Create a merged table by keeping only one row for each unique sample
# seq_id, lane_number, and data_type are all replaced with "merged" for duplicated samples
sample_table_merged <- group_by(sample_table, sample_id_corrected) %>%
  summarise(population=unique(population), seq_id=ifelse(n()==1,seq_id, "merged"), lane_number=ifelse(length(unique(lane_number))==1,unique(lane_number), "merged"), data_type=paste0(unique(data_type), collapse = "")) %>%
  mutate(sample_seq_id=paste(sample_id_corrected, seq_id, lane_number, data_type, sep = "_")) %>%
  select(sample_seq_id, lane_number, seq_id, sample_id_corrected, population, data_type)

## Write the merged table
write_tsv(sample_table_merged, "sample_table_merged.tsv")

## Create bam lists as inputs for future steps
bam_list_merged <- paste0(sample_table_merged$sample_seq_id, "_bt2_", refname, "_minq0_sorted", ".bam")

bam_list_dedup_overlapclipped <- transmute(sample_table_merged, suffix=ifelse(data_type=="se", paste0("_bt2_", refname, "_minq0_sorted_dedup.bam"), paste0("_bt2_", refname, "_minq0_sorted_dedup_overlapclipped.bam"))) %>%
                         .$suffix %>%
                         paste0(sample_table_merged$sample_seq_id, .)

bam_list_realigned <- transmute(sample_table_merged, suffix=ifelse(data_type=="se", paste0("_bt2_", refname, "_minq0_sorted_dedup_realigned.bam"), paste0("_bt2_", refname, "_minq0_sorted_dedup_overlapclipped_realigned.bam"))) %>%
  .$suffix %>%
  paste0(sample_table_merged$sample_seq_id, .)

write_lines(bam_list_merged, "bam_list_merged.tsv")
write_lines(bam_list_dedup_overlapclipped, "bam_list_dedup_overlapclipped.tsv")
write_lines(bam_list_realigned, "bam_list_realigned.tsv")

## Create Merging Script

## Find all duplicated samples
duplicated_samples <- (sample_table$sample_id_corrected)[duplicated(sample_table$sample_id_corrected)]
duplicated_samples_seq_ids <- sample_table_merged[match(duplicated_samples,sample_table_merged$sample_id_corrected),] %>%
  .$sample_seq_id
final_duplicated_samples_seq_ids <- unique(duplicated_samples_seq_ids)
merging_script<-NULL

## Loop through all duplicated samples
for (i in 1:length(duplicated_samples)){
  duplicated_sample <- duplicated_samples[i]
  duplicated_samples_seq_id <- final_duplicated_samples_seq_ids[i]
  ## Extract the bam file names from the unmerged sample table
  input <- filter(sample_table, sample_id_corrected==duplicated_sample) %>%
    mutate(unmerged_bam=paste(sample_id, seq_id, lane_number, data_type, "bt2", refname, "minq0_sorted.bam", sep = "_")) %>%
    # Note that sample_id is used in here instead of sample_id_corrected, since the unmerged bam file still uses sample_id as part of its name, not the corrected one.
    .$unmerged_bam %>%
    paste0( .) %>%
    paste(collapse = " ")
  ## Paste together the command line
  merging_script[i] <- paste0("$SAMTOOLS merge -@ $THREADS ",duplicated_samples_seq_id, "_bt2_", refname, "_minq0_sorted", ".bam ", input)
}

## Write the script
write_lines(merging_script, "merge_bam.sh")
