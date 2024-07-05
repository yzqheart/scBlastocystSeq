fq1=$1 # trimmed r1
fq2=$2 # trimmed r2
outdir=$3 # output dir
sample=$4 # sample name

PICARD="java -Xmx120g -jar /mnt/data/yanzhiqiang/142/software/picard-tools-1.141/picard.jar"
index=/mnt/data/yanzhiqiang/138/DB/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa

ID=$sample
PL=ILLUMINA
LB=$sample
SM=$sample
rg="@RG\tID:$ID\tPL:$PL\tLB:$LB\tSM:$SM"

bwa mem -v 1 -R $rg -M -t 12 $index $fq1 $fq2 > ${outdir}/${sample}/${sample}.sam

samtools view -h -b -q 10 ${outdir}/${sample}/${sample}.sam > ${outdir}/${sample}/${sample}.uniq.bam
samtools sort -o ${outdir}/${sample}/${sample}.uniq.sort.bam ${outdir}/${sample}/${sample}.uniq.bam
samtools index ${outdir}/${sample}/${sample}.uniq.sort.bam
samtools rmdup ${outdir}/${sample}/${sample}.uniq.sort.bam ${outdir}/${sample}/${sample}.rmdup.bam
samtools index ${outdir}/${sample}/${sample}.rmdup.bam

rm ${outdir}/${sample}/${sample}.uniq.bam ${outdir}/${sample}/${sample}.uniq.sort.bam* ${outdir}/${sample}/${sample}.sam
