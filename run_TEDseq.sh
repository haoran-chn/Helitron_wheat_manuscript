#!/bin/bash
set -euo pipefail

MAX_JOBS=12
THREADS_PER_SAMPLE=20

REF=TAIR10.fa
TES=(X1 F8)

SAMPLES=(
  WTH_S1
  X1_L9_S2
  F8_L2_S3
  F8_L6_S4
  X1F8_L3_S5
  X1F8_L12_S6
)

mkdir -p work
echo -e "Sample\tTE\tInput_reads\tPrimer_pass\tPrimer_%\tExtremity_pass\tExtremity_%" > primer_stats.tsv

# prepare reference
bwa index $REF
samtools faidx $REF

run_sample() {

  S=$1
  TE=$2

  PRIMER=${TE}_primer.fa
  EXT=${TE}_extremity.fa

  PRIMER_LEN=$(grep -v ">" $PRIMER | tr -d '\n' | wc -c)
  EXT_LEN=$(grep -v ">" $EXT | tr -d '\n' | wc -c)

  mkdir -p work/$TE/{01_primer,02_ext,03_R1map,04_R2map,05_validated,06_collapsed}

  R1=${S}_R1_001.fastq.gz
  R2=${S}_R2_001.fastq.gz

  echo ">>> $S / $TE"

  # =============================
  # Step1: primer anchoring
  # =============================
  cutadapt \
    -g file:$PRIMER \
    --error-rate 0 \
    --overlap $PRIMER_LEN \
    --discard-untrimmed \
    -o work/$TE/01_primer/${S}_R1_primer.fastq.gz \
    -p work/$TE/01_primer/${S}_R2_primer.fastq.gz \
    $R1 $R2 \
    > work/$TE/01_primer/${S}.log

  IN=$(zcat $R1 | echo $((`wc -l`/4)))
  OUT=$(zcat work/$TE/01_primer/${S}_R1_primer.fastq.gz | echo $((`wc -l`/4)))
  PCT=$(awk -v o=$OUT -v i=$IN 'BEGIN{printf "%.4f",o/i*100}')

  # =============================
  # Step2: extremity anchoring + trim
  # =============================
  cutadapt \
    -g file:$EXT \
    --error-rate 0 \
    --overlap $EXT_LEN \
    --minimum-length 15 \
    --discard-untrimmed \
    -o work/$TE/02_ext/${S}_R1_ext.fastq.gz \
    -p work/$TE/02_ext/${S}_R2_ext.fastq.gz \
    work/$TE/01_primer/${S}_R1_primer.fastq.gz \
    work/$TE/01_primer/${S}_R2_primer.fastq.gz

  EXT_OUT=$(zcat work/$TE/02_ext/${S}_R1_ext.fastq.gz | echo $((`wc -l`/4)))
  EXT_PCT=$(awk -v o=$EXT_OUT -v i=$OUT 'BEGIN{printf "%.4f",o/i*100}')

  echo -e "$S\t$TE\t$IN\t$OUT\t$PCT\t$EXT_OUT\t$EXT_PCT" >> primer_stats.tsv

  # =============================
  # Step3: R1 mapping (primary insertion sites, NO clipped reads)
  # =============================
  bwa mem -M -t $THREADS_PER_SAMPLE $REF work/$TE/02_ext/${S}_R1_ext.fastq.gz | \
    samtools sort -@ $THREADS_PER_SAMPLE -o work/$TE/03_R1map/${S}_R1.bam
  samtools index work/$TE/03_R1map/${S}_R1.bam

  samtools view -F 4 work/$TE/03_R1map/${S}_R1.bam | \
    awk '$6!~/[SH]/{
      strand=(and($2,16)?"-":"+");
      print $3"\t"$4"\t"$4+1"\t"strand
    }' | sort | uniq -c | \
    awk '{print $2"\t"$3"\t"$3+1"\t"$1"\t"$5}' \
    > work/$TE/03_R1map/${S}_R1_sites.bed
    # chr start end R1_support strand

  # =============================
  # Step4: R2 mapping (paired validation, NO clipped reads)
  # =============================
  bwa mem -M -t $THREADS_PER_SAMPLE $REF work/$TE/02_ext/${S}_R2_ext.fastq.gz | \
    samtools sort -@ $THREADS_PER_SAMPLE -o work/$TE/04_R2map/${S}_R2.bam
  samtools index work/$TE/04_R2map/${S}_R2.bam

  samtools view -F 4 work/$TE/04_R2map/${S}_R2.bam | \
    awk '$6!~/[SH]/{print $3"\t"$4"\t"$4+1}' | sort | uniq -c | \
    awk '{print $2"\t"$3"\t"$3+1"\t"$1}' \
    > work/$TE/04_R2map/${S}_R2_sites.bed
    # chr start end R2_support

  # =============================
  # Step5: R1 filtered by R2 ±400 bp
  # =============================
  bedtools window -w 400 \
    -a work/$TE/03_R1map/${S}_R1_sites.bed \
    -b work/$TE/04_R2map/${S}_R2_sites.bed | \
    awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' | sort -u \
    > work/$TE/05_validated/${S}_validated.bed

  # =============================
  # Step6: ±50 bp collapse (strand majority)
  # =============================
  sort -k1,1 -k2,2n work/$TE/05_validated/${S}_validated.bed | \
    bedtools cluster -d 50 | \
    awk '
    {
      k=$6;
      chr[k]=$1;
      sum[k]+=$2;
      n[k]++;
      sup[k]+=$4;
      strand[k][$5]++;
    }
    END{
      for(k in chr){
        pos=int(sum[k]/n[k]);
        max=0;
        for(s in strand[k]) if(strand[k][s]>max){max=strand[k][s]; ms=s}
        print chr[k]"\t"pos"\t"pos+1"\t"sup[k]"\t"ms;
      }
    }' | sort -k1,1 -k2,2n \
    > work/$TE/06_collapsed/${S}_final_insertions.bed
}

export -f run_sample
export REF THREADS_PER_SAMPLE

parallel -j $MAX_JOBS run_sample ::: "${SAMPLES[@]}" ::: "${TES[@]}"

