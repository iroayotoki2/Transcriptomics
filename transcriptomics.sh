#Assignment 2 Code script
#Generating the fastq files from the SRR files using the SRA toolkit
mkdir -p tmp
for SRR in SRR*
do
fasterq-dump "$SRR" --temp tmp 
rm -rf tmp/*
done
#Fastqc and multiqc reports
fastqc SRR10551658.fastq SRR10551660.fastq SRR10551662.fastq SRR10551664.fastq SRR10551657.fastq SRR10551659.fastq  SRR10551661.fastq  SRR10551663.fastq  SRR10551665.fastq
multiqc .
#STAR Alignment
gzip *.fastq #converting to .gz to conserve storage 
 STAR --runMode genomeGenerate --genomeDir ~/binf6110/Assignment2 --genomeFastaFiles Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa --sjdbGTFfile Saccharomyces_cerevisiae.R64-1-1.113.gtf --sjdbOverhang ReadLength-1 

for FILE in SRR*.fastq.gz
do
    STAR \
      --genomeDir ~/binf6110/Assignment2 \
      --readFilesIn "$FILE" \
      --readFilesCommand zcat \
      --runThreadN 2 \
      --quantMode TranscriptomeSAM \
	  --outSAMtype None \
      --outFileNamePrefix "${FILE%.gz}_"
done
#RSEM Quantification
#Prepare reference
rsem-prepare-reference \
  --star \
  --gtf Saccharomyces_cerevisiae.R64-1-1.113.gtf \
  Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
  Saccharomyces_ref
#Quantification loop
for BAM in *_Aligned.toTranscriptome.out.bam
do
  SAMPLE=${BAM%_Aligned.toTranscriptome.out.bam}

  rsem-calculate-expression \
    --bam \
    --no-bam-output \
    --estimate-rspd \
    -p 2 \
    "$BAM" \
    Saccharomyces_ref \
    "$SAMPLE"
done
