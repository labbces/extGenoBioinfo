#Download Spodoptera frugiperda genome and annotation
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/101/765/GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0/GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/101/765/GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0/GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.gff.gz
gunzip GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.fna.gz
gunzip GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.gff.gz
sort -k1,1 -k4,4n  GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.gff > GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.sorted.gff
bgzip GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.sorted.gff
tabix -p gff GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.sorted.gff.gz

#Download Reads
#RNASeq READS Illumina
wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR292/067/SRR29288667/SRR29288667_1.fastq.gz"
wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR292/067/SRR29288667/SRR29288667_2.fastq.gz"
#GENOMIC READS Illumina
wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR291/066/SRR29141966/SRR29141966_1.fastq.gz"
wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR291/066/SRR29141966/SRR29141966_2.fastq.gz"

#GENOMIC PacBio Sequel Reads
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR120/097/SRR12072097/SRR12072097_subreads.fastq.gz

#Generate mapping files
conda activate transcriptomics
SAMPLE=SRR29288667
bbduk.sh in=${SAMPLE}_1.fastq.gz in2=${SAMPLE}_2.fastq.gz out=${SAMPLE}_clean_1.fastq.gz out2=${SAMPLE}_clean_2.fastq.gz ref=adapters refstats=${SAMPLE}_clean_adapters_refstats ktrim=r threads=10
hisat2-build -p 10 GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.fna GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic
hisat2 --threads 10 --fast -x GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic -1 ${SAMPLE}_clean_1.fastq.gz -2 ${SAMPLE}_clean_2.fastq.gz | samtools view -bS - | samtools sort -o ${SAMPLE}_RNASeq_sorted.bam 
samtools index ${SAMPLE}_RNASeq_sorted.bam
samtools faidx GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.fna

SAMPLE2=SRR29141966
bbduk.sh in=${SAMPLE2}_1.fastq.gz in2=${SAMPLE2}_2.fastq.gz out=${SAMPLE2}_clean_1.fastq.gz out2=${SAMPLE2}_clean_2.fastq.gz ref=adapters refstats=${SAMPLE2}_clean_adapters_refstats ktrim=r threads=10
bwa index GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.fna
bwa mem -t 10 GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.fna ${SAMPLE2}_clean_1.fastq.gz ${SAMPLE2}_clean_2.fastq.gz | samtools view -@ 10 -bS - | samtools sort -@ 10 -o ${SAMPLE2}_Illumina_sorted.bam
samtools index ${SAMPLE2}_Illumina_sorted.bam

conda deactivate

conda activate jupiterplot
minimap2 -a -x map-pb GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.fna SRR12072097_subreads.fastq.gz | samtools view -bS - | samtools sort -o SRR12072097_PacBio_sorted.bam

samtools index SRR12072097_PacBio_sorted.bam
conda deactivate


#extract mappings to chromosome 2 only, to fasta display in IGV
conda activate emboss
seqret GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.fna:NC_064213.1 NC_064213.1.fasta
conda deactivate

conda activate jupiterplot
bamtools filter -in SRR12072097_PacBio_sorted.bam -out SRR12072097_PacBio_NC_064213.1_sorted.bam -region NC_064213.1
samtools index SRR12072097_PacBio_NC_064213.1_sorted.bam
bamtools filter -in SRR29288667_RNASeq_sorted.bam -out SRR29288667_RNASeq_NC_064213.1_sorted.bam -region NC_064213.1
samtools index SRR29288667_RNASeq_NC_064213.1_sorted.bam
bamtools filter -in SRR29141966_Illumina_sorted.bam -out SRR29141966_Illumina_NC_064213.1_sorted.bam -region NC_064213.1
samtools index SRR29141966_Illumina_NC_064213.1_sorted.bam

conda deactivate
