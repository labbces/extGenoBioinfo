wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/101/765/GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0/GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.fna.gz
gunzip GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.fna.gz

#RNASeq READS
wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR292/067/SRR29288667/SRR29288667_1.fastq.gz"
wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR292/067/SRR29288667/SRR29288667_2.fastq.gz"

#PacBio Sequel Reads
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR120/097/SRR12072097/SRR12072097_subreads.fastq.gz

conda activate transcriptomics
SAMPLE=SRR29288667
bbduk.sh in=${SAMPLE}_1.fastq.gz in2=${SAMPLE}_1.fastq.gz out=${SAMPLE}_clean_1.fastq.gz out2=${SAMPLE}_clean_2.fastq.gz ref=adapters refstats=${SAMPLE}_clean_adapters_refstats ktrim=r threads=10
hisat2-build -p 10 GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.fna GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic
hisat2 --threads 10 --fast -x GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic -1 ${SAMPLE}_clean_1.fastq.gz -2 ${SAMPLE}_clean_2.fastq.gz | samtools view -bS - | samtools sort -o ${SAMPLE}_sorted.bam 
samtools index ${SAMPLE}_sorted.bam
samtools faidx GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.fna
conda deactivate

conda activate jupiterplot
minimap2 -a -x map-pb GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.fna SRR12072097_subreads.fastq.gz | samtools view -bS - | samtools sort -o SRR12072097_PacBio_sorted.bam

samtools index SRR12072097_PacBio_sorted.bam
conda deactivate
