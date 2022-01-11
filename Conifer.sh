
#Reference 

bed="Twist_Exome_Target_hg38.bed"
mkdir -p CONIFER

docker run \
-v /mnt/Data/NGS_Analysis_Database_Tools/Hg38/GATK_hg38_Resources:/home/bio/Ref \
-v /mnt/Data/bam_all/CNV_Analysis:/home/bio/Data \
-v /mnt/Data/bam_all/:/mnt/Data/bam_all \
-it molecular/conifer

bed="Twist_Exome_Target_hg38.bed"

#1. Calculate RPKM values for all samples
for i in $(ls Data/*bam) ; do
sample=${i::-10}
sample=$(echo ${sample} |awk -F"/" '{print $2}')

python /home/bio/conifer_v0.2.2/conifer.py  rpkm --probes Data/${bed} \
--input Data/${sample}_recal.bam --output ${sample}_rpkm.txt
done

mkdir -p CONIFER/rpkm
mv *rpkm.txt CONIFER/rpkm


#2 Analyze all RPKM values for all samples and create SVD-ZRPKM values. This steps create a single 
#unified conifer data file which contains all data, probes and samples for downstream analysis.

python /home/bio/conifer_v0.2.2/conifer.py analyze --probes /home/bio/Data/${bed} \
--rpkm_dir /home/bio/CONIFER/rpkm/ \
--output /home/bio/CONIFER/CONIFER_analysis.hdf5 \
--svd 0 \
--write_svals /home/bio/CONIFER/CONIFER_singular_values.txt \
--plot_scree /home/bio/CONIFER/CONIFER_screeplot.png \
--write_sd /home/bio/CONIFER/CONIFER_sd_values.txt

#Making calls:
python /home/bio/conifer_v0.2.2/conifer.py call \
  --input /home/bio/CONIFER/CONIFER_analysis.hdf5 \
  --output /home/bio/CONIFER/CONIFER_calls.txt  \
   --threshold 0.95


#From a different terminal
docker ps
#cOPY REQUIRED FOLDER TO THE SYSTEM
# 51fbffcbafea IS CONTAINER ID

sudo docker cp 51fbffcbafea:/home/bio/CONIFER .


cat CONIFER/CONIFER_calls.txt | awk -F"\t" 'BEGIN{OFS="\t"};NR>1{print $2,$3,$4,$5,$1,25}' |  \
sed '1 i#Chrom\tStart\tEnd\tSV_type\tSamples_ID\tScore' > CONIFER/CONIFER_calls.bed

#Annotate the CNV using ANNOTSV
export ANNOTSV=/mnt/Data/NGS_Analysis_Database_Tools/AnnotSV
$ANNOTSV/bin/AnnotSV -SVinputFile CONIFER_calls.bed \
        -outputFile CONIFER_calls.annotated_75.tsv -svtBEDcol 4 \
        -genomeBuild GRCh38 -overlap 75 -REreport yes

