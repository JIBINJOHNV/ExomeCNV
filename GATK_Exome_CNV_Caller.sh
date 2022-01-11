

#Reference; https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-common-and-rare-germline-copy-number-variants

docker run -v /mnt/Data/NGS_Analysis_Database_Tools/Human_Exome_hg19:/gatk/data \
-v /mnt/Data/NGS_Analysis_Working/Sandeep/Burden_Test_CH/CNV_TEST:/gatk/bam  -it broadinstitute/gatk:latest


Fasta="/gatk/data/hg19.fa"
Bed="/gatk/data/CReV2_Regions.bed"
BamFolder="/gatk/bam/BamFiles/"

N=5

ls ${BamFolder}*.bam >bamfile.txt
rm Sample_Names.txt temp.txt
cat temp.txt

while read Fastq; do echo $Fastq | awk -F"/" '{print $NF}' | cut -d "." -f 1 >> temp.txt ;done < bamfile.txt
cat temp.txt |sort| uniq |sed '/^$/d' > Sample_Names.txt
cat -n Sample_Names.txt




#Preprocesses Intervels ; For exome data, pad target regions, e.g. with 250 bases.
echo -e "1. Padd target region with 250 bases:\n"

mkdir -p GATK_CNV
gatk PreprocessIntervals \
        -R ${Fasta} \
        -L ${Bed} \
        --bin-length 0 \
        -imr OVERLAPPING_ONLY \
        -O GATK_CNV/targets.preprocessed.interval_list


#Count reads per bin using CollectReadCounts
echo -e "2. Count reads per bin using CollectReadCounts"

mkdir -p GATK_CNV/2.GatkCoverage
while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
    gatk CollectReadCounts \
    -L GATK_CNV/targets.preprocessed.interval_list -R ${Fasta} \
    -imr OVERLAPPING_ONLY \
    -I ${BamFolder}${Samplename}.bam \
    --format TSV \
    -O GATK_CNV/2.GatkCoverage/${Samplename}.tsv & 

done < Sample_Names.txt ; wait


#AnnotateIntervals with GC content
echo -e "3. AnnotateIntervals with GC content"
gatk AnnotateIntervals \
    -L GATK_CNV/targets.preprocessed.interval_list \
    -R ${Fasta} \
    -imr OVERLAPPING_ONLY \
    -O GATK_CNV/targets.preprocessed.interval_list.annotated.tsv


#FilterIntervals based on GC-content and cohort extreme counts
echo -e "4. FilterIntervals based on GC-content and cohort extreme counts"
gatk FilterIntervals \
    -L GATK_CNV/targets.preprocessed.interval_list \
    $echo$(ls GATK_CNV/2.GatkCoverage/ | tr "\t" "\n" | awk -v pwd="$PWD" '{print "-I","GATK_CNV/2.GatkCoverage/"$0}' | tr "\n" " ") \
    -imr OVERLAPPING_ONLY \
    -O GATK_CNV/cohort.gc.filtered.interval_list


#Call autosomal and allosomal contig ploidy with DetermineGermlineContigPloidy

echo -e "5. DetermineGermlineContigPloidy in COHORT MODE"

gatk DetermineGermlineContigPloidy \
-L GATK_CNV/cohort.gc.filtered.interval_list \
--interval-merging-rule OVERLAPPING_ONLY \
$echo$(ls GATK_CNV/2.GatkCoverage/ | tr "\t" "\n" | awk '{print "-I","GATK_CNV/2.GatkCoverage/"$0}' | tr "\n" " ") \
    --contig-ploidy-priors bam/priors.tsv \
    --output GATK_CNV/ \
    --output-prefix ploidy \
    --verbosity DEBUG




gatk GermlineCNVCaller \
    --run-mode COHORT \
    -L GATK_CNV/cohort.gc.filtered.interval_list \
$echo$(ls GATK_CNV/2.GatkCoverage/ | tr "\t" "\n" | awk '{print "-I","GATK_CNV/2.GatkCoverage/"$0}' | tr "\n" " ") \
    --contig-ploidy-calls GATK_CNV/ploidy-calls \
    --annotated-intervals GATK_CNV/targets.preprocessed.interval_list.annotated.tsv \
    --interval-merging-rule OVERLAPPING_ONLY \
    --output GATK_CNV/GATK_CNV \
    --output-prefix GATK_CNV \
    --verbosity DEBUG



##Post processes

INDEX=0
N=5
while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
gatk PostprocessGermlineCNVCalls \
    --model-shard-path GATK_CNV/GATK_CNV/GATK_CNV-model \
    --calls-shard-path  GATK_CNV/GATK_CNV/GATK_CNV-calls \
    --allosomal-contig chrX --allosomal-contig chrY \
    --output-genotyped-intervals GATK_CNV/genotyped-intervals-cohort_${Samplename}.vcf.gz \
    --output-genotyped-segments  GATK_CNV/genotyped-segments-cohort_${Samplename}.vcf.gz \
    --sequence-dictionary /gatk/data/hg19.dict \
    --contig-ploidy-calls GATK_CNV/ploidy-calls \
    --sample-index $INDEX \
    --output-denoised-copy-ratios GATK_CNV/Denoise-cohort_${Samplename} &
    INDEX=$(($INDEX+1))
done < Sample_Names.txt ; wait


java -jar -Xmx4g BuildBamIndex.jar I=foo.bam O=foo.bam.bai
