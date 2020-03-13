for i in {1..22} X 
do 
	awk -F'\t' -v chromosome="$i" 'BEGIN {OFS="\t"} {if (NR==1) {print "#chr",$0} else {if ($2=="SNP") {print chromosome,$1,$2,$3,$4,$5,$6,$7,$8,$9,$6,$7,$8,$9,$14,$15} else {print chromosome,$0}}}' RUMapv3_B137_chr${i}.txt > RUMap_chr${i}.txt 
done
