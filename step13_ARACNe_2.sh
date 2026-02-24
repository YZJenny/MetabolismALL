#### 2.run ARACNe ####
###node10上运行
#Step1 Calculate threshold with a fixed seed
# java -Xmx5G -jar /local/yzj/software/ARACNe-AP/ARACNe-AP-master/dist/aracne.jar -e /mdshare/node10/yzj/CRU/ped_M/res_ALL/ARACNe/mat.txt -o /mdshare/node10/yzj/CRU/ped_M/res_ALL/ARACNe/Network --tfs /mdshare/node10/yzj/CRU/ped_M/data/Human_TF/TF_names_v_1.01.txt --pvalue 1E-8 --seed 1 --calculateThreshold &
#   
#Step2: Run 100 reproducible bootstraps
for i in {1..100}
do
java -Xmx5G -jar /local/yzj/software/ARACNe-AP/ARACNe-AP-master/dist/aracne.jar -e /mdshare/node10/yzj/CRU/ped_M/res_ALL/ARACNe/mat.txt -o /mdshare/node10/yzj/CRU/ped_M/res_ALL/ARACNe/Network --tfs /mdshare/node10/yzj/CRU/ped_M/data/Human_TF/TF_names_v_1.01.txt --pvalue 1E-8 --seed $i --threads 60
done

#Step3: Consolidate bootstraps in the output folder
java -Xmx5G -jar /local/yzj/software/ARACNe-AP/ARACNe-AP-master/dist/aracne.jar -o /mdshare/node10/yzj/CRU/ped_M/res_ALL/ARACNe/Network --consolidate --threads 60

