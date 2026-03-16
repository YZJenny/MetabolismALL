cd /mdshare/node9/yanzijun/software/GISTIC2/
nohup ./gistic2 \
  -b /mdshare/node10/yzj/CRU/ped_M/res/GISTIC2/ALL/ \
  -seg /mdshare/node10/yzj/CRU/ped_M/res/GISTIC2/segmentationfile_ALL.txt \
  -refgene /mdshare/node9/yanzijun/software/GISTIC2/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat \
  -mk /mdshare/node9/yanzijun/CRU/ped_M/data/CNV/annoProbe/Marker.txt \
  -alf /mdshare/node10/yzj/CRU/ped_M/res/GISTIC2/arraylistfile_ALL.txt \
  -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 \
  -armpeel 1 -savegene 1 -gcm extreme > /mdshare/node10/yzj/CRU/ped_M/res/GISTIC2/ALL.log 2>&1 &
  