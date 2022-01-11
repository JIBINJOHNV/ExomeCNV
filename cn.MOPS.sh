
to kill docker sudo systemctl restart docker.socket docker.service

DataFolder='/mnt/Data/bam_all/CNV_Analysis/'
mkdir -p cn.mops

docker run \
-v /mnt/Data/bam_all/CNV_Test/:/home/bio/Data \
-it dceoy/cn.mops:latest

library(cn.mops)

setwd("/home/bio/Data/")
BAMFiles <- list.files(,pattern=".bam$")

segments <- read.table("Twist_Exome_Target_hg38.bed",
                       sep="\t",as.is=TRUE)
gr<-GRanges(segments[,1],IRanges(segments[,2],segments[,3]))
X<-getSegmentReadCountsFromBAM(BAMFiles,GR=gr)

resCNMOPS<-exomecn.mops(X)
resCNMOPS<-calcIntegerCopyNumbers(resCNMOPS)
CNVs<-as.data.frame(cnvs(resCNMOPS))
write.csv(CNVs,file="cnmops.csv")


docker ps
sudo docker cp 07a8d574af8e:/home/bio/Data/cnmops.csv .

cat cn.mops/cnmops.csv | awk -F"," 'BEGIN{OFS="\t"};NR>1{print $2,$3,$4,$10,\
$7}' | sed 's/"//g' |  sed '1 i#Chrom\tStart\tEnd\tSV_type\tSamples_ID' > cn.mops/cnmops_CNV.bed



mkdir -p cn.mops/AnnotatedCNVs/
export ANNOTSV=/mnt/Data/NGS_Analysis_Database_Tools/AnnotSV

$ANNOTSV/bin/AnnotSV -SVinputFile cn.mops/cnmops_CNV.bed \
        -outputFile cn.mops/AnnotatedCNVs/cnmops_CNV.annotated.tsv -svtBEDcol 4 \
        -genomeBuild GRCh38 



sudo systemctl restart docker.socket docker.service # run from different  terminal
mkdir cnmops
mv cnmops.csv cnmops
