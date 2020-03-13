for i in {1..22} X
do
	bgzip -c RUMap_chr${i}.txt > RUMap_chr${i}.txt.gz
	tabix  -s1 -b7 -e7 -c# RUMap_chr${i}.txt.gz
done 
