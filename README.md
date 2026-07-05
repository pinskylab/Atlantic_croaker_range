# Atlantic_croaker_range
## Repository for bioinformatics and genomic analysis of Atlantic croaker (Micropogonias undulatus) collected off the Eastern United States. 

This repository provides the code for the third chapter of the doctoral dissertation Conservation Genomics of Marine Fish Populations: https://doi.org/doi:10.7282/t3-wgyk-5s02. This research was supported by a NOAA Margaret A. Davidson fellowship. 

As genomic data is too large to be stored on GitHub, sequence data were uploaded to the NCBI SRA and metadata to GEOME-DB.

Metadata is available at: https://n2t.net/ark:/21547/GjB2.

PacBio Hifi Sequencing data for one individual (used to generate a de novo assembly): BioProject PRJNA1460334 https://www.ncbi.nlm.nih.gov/sra/PRJNA1460334.

Illumina reads (low-coverage whole genome sequencing) for 400 individuals: BioProject PRJNA1457870 https://www.ncbi.nlm.nih.gov/bioproject/1457870.

Genome assembly, pre-processing, mapping, and downstream analysis were run using Amarel, Rutgers University's high performance computing system.
Statistical analysis was run using R v.4.4.1 on a MacBook Pro (Apple M3 Pro Chip with 36 GB memory).

# De Novo Genome Assembly
## Assembly using Hifiasm

Hifiasm v. 0.19.3-r572  
sbatch hifiasm_code.sh

### Evaluating Hifiasm assemblies
QUAST v. 5.2.0  
sbatch quast_code.sh  

Selected primary contig assembly to move forward based on highest N50, lowest L50, and a total length similar to other croaker species. N50 18,479,188; L50 17; # contigs 3921; total length 782,013,506.

BUSCO v. 5.4.6 using lineage actinopterygii on primary contig Hifiasm assembly (hifi4523.asm.bp.p_ctg.fa)  
sbatch busco_code.sh  
Complete:98.9%[Single:97.4%,Duplicate:1.5%],Fragemented:0.2%,Missing:0.9%,n:3640

