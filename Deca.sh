

#####deca
DECA_JAR='/mnt/Data/NGS_Analysis_Working/Sandeep/Burden_Test_CH/CNV_TEST/Software/deca/deca-cli/target/deca-cli_2.11-0.2.1-SNAPSHOT.jar'
pseq='/mnt/Data/NGS_Analysis_Working/Sandeep/Burden_Test_CH/CNV_TEST/Software/plinkseq/pseq'
fasta='/mnt/Data/NGS_Analysis_Database_Tools/Human_Exome_hg19/hg19.fa'
TargetBed='Twist_Exome_Target_hg38.bed'
xhmm='/mnt/Data/NGS_Analysis_Working/Sandeep/Burden_Test_CH/CNV_TEST/Software/statgen-xhmm-998f7c405974/xhmm'
gatk='java -jar /mnt/Data/NGS_Analysis_Working/Sandeep/Burden_Test_CH/CNV_TEST/Software/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar'
deca='/mnt/Data/NGS_Analysis_Working/Sandeep/Burden_Test_CH/CNV_TEST/Software/deca/bin/deca-submit'

export PATH=$PATH:/mnt/Data/NGS_Analysis_Working/Sandeep/Burden_Test_CH/CNV_TEST/Software/deca/spark-2.1.0-bin-hadoop2.7/bin


awk -F"\t" '{print $1":"$2"-"$3}' ${TargetBed} > CReV2_Regions.target_gatk.bed
mkdir -p Deca_CNV


${deca} \
--master local[16] \
--driver-class-path $DECA_JAR \
--conf spark.local.dir=. \
--conf spark.driver.maxResultSize=0 \
--conf spark.kryo.registrationRequired=true \
--driver-memory 16G \
-- coverage \
-L ${TargetBed} \
-I *.bam \
-o Deca_CNV/DECA.RD.txt



#-------------------------------------------------------------------------------------------------------------------------------
#Optional: run Plink/Seq to calculate sequence complexity of targets ; since seqdb for hg38 is not available this step not running 
#sed "1i#CHR\tPOS1\tPOS2\tID" CReV2_Regions.target.bed > CReV2_Regions.target.PlinkSeq.reg
#pseq='/mnt/Data/NGS_Analysis_Working/Sandeep/Burden_Test_CH/CNV_TEST/Software/plinkseq/pseq'

#${pseq} . loc-load --locdb EXOME.targets.LOCDB --file CReV2_Regions.target.PlinkSeq.reg --group targets

#${pseq} . loc-stats --locdb EXOME.targets.LOCDB --group targets --seqdb seqdb | awk "{if (NR > 1) print $_}" |
#sort -k1 -g | awk "{print $10}" | paste CReV2_Regions.target.bed - |
#awk "{print $1"t"$2}" > DATA.locus_complexity.txt



# run GATK to calculate GC content of targets

#${gatk} -T GCContentByInterval -L ${TargetBed} -R ${fasta} -o TargetRegion_GC.txt
    # Create a list of all exome targets to be excluded, based on their GC content: 
#cat TargetRegion_GC.txt | awk "{if ($2 < 0.1 || $2 > 0.9) print $1}" > extreme_gc_targets.txt

#cat TargetRegion_GC.txt | awk '{if ($2< 0.1 || $2> 0.9) print "blah" }' > extreme_gc_targets.txt

#----------------------------------------------------------------------------------------------------------------------------------------------



${deca} \
--master local[16] \
--driver-class-path $DECA_JAR \
--conf spark.local.dir=. \
--conf spark.driver.maxResultSize=0 \
--conf spark.kryo.registrationRequired=true \
--driver-memory 16G \
-- normalize_and_discover \
-min_some_quality 29.5 \
-I Deca_CNV/DECA.RD.txt \
-o Deca_CNV/DECA.gff3

#-exclude_targets exclude_targets.txt \


mkdir -p Deca_CNV/AnnotatedCNVs

cat Deca_CNV/DECA.gff3 | awk -F"\t" 'BEGIN{OFS="\t"};NR>1{print $1, $4-1,$5,
$3,$2,$6}' | sed '1 i#Chrom\tStart\tEnd\tSV_type\tSamples_ID\tScore' > Deca_CNV/AnnotatedCNVs/DECA_CNV.bed


export ANNOTSV=/mnt/Data/NGS_Analysis_Database_Tools/AnnotSV

$ANNOTSV/bin/AnnotSV -SVinputFile Deca_CNV/AnnotatedCNVs/DECA_CNV.bed \
        -outputFile Deca_CNV/AnnotatedCNVs/DECA_CNV.annotated_75.tsv -svtBEDcol 4 \
        -genomeBuild GRCh38 -overlap 75 -REreport yes




#-overlap           Minimum overlap (%) between the user features and the annotated SV to be reported 
                    #Range values: [0-100], default = 100

#-minTotalNumber    Minimum number of individuals tested to consider a benign SV for the ranking
                    #Range values: [100-1000], default = 500

#-promoterSize:     Number of bases upstream from the transcription start site
                    #Default = 500

#-promoterSize:     Number of bases upstream from the transcription start site
                    #Default = 500
#-REreport:         Create a report to link the annotated SV and the overlapped regulatory elements (coordinates and sources)
                    #alues: no (default) or yes

#-SVminSize:        #SV minimum size (in bp)   Default = 50
#-svtBEDcol:        #Number of the column describing the SV type (DEL, DUP) if the input SV file is a BED
                    #Range values: [4-[, default = -1 (value not given)





conda activate DECoNT_linux
