#### 
##https://api.gdc.cancer.gov/data/9bd7cbce-80f9-449e-8007-ddc9b1e89dfb下载snp6.na35.remap.hg38.subset.txt.gz
##https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files下载snp6.na35.remap.hg38.subset.txt.gz
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6801 download probe info. of snp6.0:GPL6801-4019.txt(hg19)

### step1 prepare Marker File
snp <- read.table('/mdshare/node9/yanzijun/CRU/ped_M/data/CNV/annoProbe/snp6.na35.remap.hg38.subset.txt',header = T)
snp <- snp[,1:3]
colnames(snp) <- c('Marker name','chromosome','Marker position')
write.table(snp,'/mdshare/node9/yanzijun/CRU/ped_M/data/CNV/annoProbe/Marker.txt',sep='\t',quote = F,col.names = F,row.names = F)

### step2 make bed format from GPL6801-4019.txt
cd /mdshare/node9/yanzijun/CRU/ped_M/data/CNV/annoProbe
grep -E '^SNP|AFFX-SNP' GPL6801-4019.txt | cut -f 2,7 | grep -v "^---" | sed -e 's/^/chr/' > GPL6801_SNP_left.bed
cut -f 2 GPL6801_SNP_left.bed > GPL6801_SNP_right.bed
paste GPL6801_SNP_left.bed GPL6801_SNP_right.bed > GPL6801_SNP.bed

grep '^CN' GPL6801-4019.txt | cut -f 2,4,5 | grep -v "^---" | sed -e 's/^/chr/' > GPL6801_CN.bed


cat GPL6801_SNP.bed GPL6801_CN.bed > GPL6801_hg19.bed

## liftover hg19 to hg38
/local/yzj/software/liftOver GPL6801_hg19.bed /local/zy/tools/files_liftOver/hg19ToHg38.over.chain.gz GPL6801_hg38.bed GPL6801_hg38_unmap.bed
