# Atlantic_croaker_range
## Repository for bioinformatics and genomic analysis of Atlantic croaker (Micropogonias undulatus) collected off the Eastern United States. 

This repository provides the code for the third chapter of the doctoral dissertation Conservation Genomics of Marine Fish Populations: https://doi.org/doi:10.7282/t3-wgyk-5s02. This research was supported by a NOAA Margaret A. Davidson fellowship (Grant # NA22NOS4200056). Samples were obtained by the NEAMAP and SEAMAP-SA surveys under VIMS IACUC-2020-02-24-14108-jxgart and Rutgers IACUC Protocol PROTO999900001.

As genomic data is too large to be stored on GitHub, sequence data were uploaded to the NCBI SRA and metadata to GEOME-DB.

Metadata is available at: https://n2t.net/ark:/21547/GjB2.

PacBio Hifi Sequencing data for one individual (used to generate a de novo assembly): BioProject PRJNA1460334 https://www.ncbi.nlm.nih.gov/sra/PRJNA1460334.

Illumina reads (low-coverage whole genome sequencing) for 400 individuals: BioProject PRJNA1457870 https://www.ncbi.nlm.nih.gov/bioproject/1457870.

Genome assembly, pre-processing, mapping, and downstream analysis were run using Amarel, Rutgers University's high performance computing system.
Statistical analysis was run using R v.4.3.2 on a MacBook Pro (Apple M3 Pro Chip with 36 GB memory).

Contact: kyrasfitz at gmail dot com

# De Novo Genome Assembly
## Assembly using Hifiasm

Hifiasm v. 0.19.3-r572  
sbatch hifiasm_code.sh

### Evaluating Hifiasm assemblies
QUAST v. 5.2.0  
sbatch quast_code.sh  

Selected primary contig assembly to move forward based on highest N50, lowest L50, and a total length similar to other croaker species. N50 18,479,188; L50 17; GC% 42.08; # contigs 3921; total length 782,013,506.

BUSCO v. 5.4.6 using lineage actinopterygii on primary contig Hifiasm assembly (hifi4523.asm.bp.p_ctg.fa)  
sbatch busco_code.sh  
Complete:98.9%[Single:97.4%,Duplicate:1.5%],Fragemented:0.2%,Missing:0.9%,n:3640

BlobToolKit: identifying and filtering out contamination

Running Blastn as first step to identifying contamination (following https://blobtoolkit.genomehubs.org/install/#databases):  
Blast v. 2.6.0  
sbatch blastn_code.sh  

sbatch blobtools_create.sh 
Plots show reads with GC content < 40% are not vertebrates. Will remove these reads. 

Filtering out contamination with Blobtools:  
sbatch blobtools_create_filter.sh  

Filtered assembly: N50 19,238,659; L50 17, GC% 42.54; # of contigs 2860; total length 731,051,241  

Moving forward with filtered assembly for mapping of Illumina reads: hifi4523.asm.bp.p_ctg.filtered.fa

### Annotating genome assembly

Using Liftoff v. 1.6.3 and the Collicthys lucidus assembly GCA_004119915.2 (Gan et al. 2021), the closet relative with an annotated
genome, as a reference. https://github.com/agshumate/Liftoff

Produces annotation hifi.gff

# Pre-processing and Mapping of Illumina reads

Scripts are adapted from Nina Therkildsen's lab's data-processing repository (https://github.com/therkildsen-lab/data-processing).  

Reads were run through pre-processing steps in two batches (Run1 and Run2) based on sequencing runs and merged after mapping.
See Therkildsen repository above for a description of the sample tables and lists used in the scripts.

## Pre-processing

### Adapter trimming
trimmomatic v. 0.39  
Sample Lists/Tables: Sample_Table_Run1.tsv and Sample_List_Run1.tsv; Sample_Table_Run2.tsv and Sample_List_run2.tsv
sbatch adapter_clipping.sh

### Quality filtering
fastp v. 0.23.4
multiqc v. 1.14  

Trimming to quality scores > 30  
sbatch quality_filtering_q30.sh  

One additional trimming step to remove index barcodes at the beginning of reads  
sbatch trim2.sh

## Mapping

bowtie2 v. 2.5.1  
samtools v. 1.17

### Build Bowtie Reference Index

sbatch build_bowtie_ref_index.sh  
Outputs .bt2, .dict, and .fai files  

### Mapping with Bowtie2

Sample_List_Combined.tsv and Sample_Table_Combined.tsv  
sbatch low_coverage_mapping.sh  

Running samtools conversion to bam file  
sbatch samtools_mapping.sh

### Merging bam files
bam_list_merged.tsv  
sbatch merge_bam.sh

### Overlap clipping of merged bam files
bamutil v. 1.0.15  
sbatch clipoverlap.sh

## In-del realignment
GATK v. 3.7.0  
sbatch realign_indels.sh  

The resulting merged, overlap clipped, and realigned bam files are used as the inputs for ANGSD.

## Downstream analysis using ANGSD

ANGSD v. 0.940  
PCAngsd v. 0.98.2  

### SNP Calling with ANGSD

sbatch angsd_global_snp_calling.sh  
Generates list of SNPs to be used for subsequent analyses.

### Genotype-likelihood estimation 

sbatch get_beagle.sh  
Generates Beagle file: bam_list_realigned_mindp132_maxdp4000_minind0.beagle.gz

### PCAngsd: PCA, Admixture, Selection Scans

Outputs from PCAngsd are used in the Statistical Analysis in R section  

sbatch run_pcangsd.sh  
Ouputs covariance matrix: pcangsd_bam_list_realigned_mindp132_maxdp4000_minind0.cov 

Analyze and plot data in R using individual_pca_functions.R script  
Uses pcangsd_bam_list_realigned_mindp132_maxdp4000_minind0.cov and pop_labels.csv as input files  

Selection Scan  
run_pcangsd_selection.sh  
Outputs npy file: pcangsd_bam_list_realigned_mindp132_maxdp4000_minind0.selection.npy  
Analyze in R with selection_scans.R script

Admixture with PCAngsd  
sbatch run_pcangsd_admix.sh  
Selects K=2 based on MAP analysis  
Outputs: pcangsd_bam_list_realigned_mindp132_maxdp4000_minind0_e.admix.Q.npy, pcangsd_bam_list_realigned_mindp132_maxdp4000_minind0_e.cov  

Analyze and plot in R with admixture.R script

### SFS generation, FST calculations, and Isolation-by-Distance

Generate site frequency spectrum (SFS)  
sbatch sfs_angsd.sh  
Produces .saf file, run get_sfs.sh script to get sfs from saf file  
sbatch get_sfs.sh 

Plotting sfs in R using sfs.R script (uses .sfs files as inputs)  

Calculate pairwise Fst between populations  
sbatch get_fst.sh script

Use pairwise Fst results to test for isolation-by-distance using the ibd_analysis.R script  
Uses Lat_Lon_9pop.csv and Pairwise_Fst_9pop.csv files as inputs  

### Genotype-environnment association (GEA) analysis

Uses individual allele frequencies generated by PCAngsd  
Include -indf_save command in PCAngsd code to save this frequency matrix as an output  
sbatch run_pcangsd_freqmatrix.sh  
Using *.indf.npy output for subsequent analysis in R  

Loading in environmental data from Copernicus Marine Services for the GEA analysis in env_data.R script  
Environmental variables are latitude, depth, bottom temperature, bottom salinity, and bottom dissolved oxygen  
Uses mercatorbiomer4v2r1_global_mean_bio_202210.nc and glo12_rg_1m-m_202210-202210_2D_hcst.nc files as inputs  

Running GEA analyses (RDA and LFMM) in rda_analysis.R and lfmm_analysis.R scripts  
rda_analysis.R script uses rdadapt.R script to identify candidate SNPs  
Uses pcangsd_minmaf0.05_bam_list_realigned_mindp132_maxdp4000_minind0_minmaf0.indf.npy, bam_list_realigned_mindp132_maxdp4000_minind0_minmaf0.05.mafs.txt, and bam_list_realigned_mindp132_maxdp4000_minind0.50_minq20.pos.gz files as inputs  

Annotating candidate SNPS identified by RDA with SnpEff v.4_3s, using the hifi.gff file produced in the Genome Annotation section.  
sbatch run_snpeff.sh

Outputs 2320snps.ann.vcf and snpEff_genes.txt which give the SNP annotations.


