help:
	@echo "This is test program for SEQLinkage package. Use auto-complete to view options"
regular:
	@seqlink --fam MissingData/real100quads.ped --vcf MissingData/rep1.vcf.gz -f MERLIN --output output/DISC1nppnomiss --prephased   #family:100 nuclear families; vcf: simulated from ped file with missing_call rate set to 0
missing:
	@seqlink --fam MissingData/real100quads.ped --vcf MissingData/rep3.vcf.gz -f MERLIN --output output/DISC1nppmiss --prephased    #same family structure; vcf: simulated from ped file with missing_call rate set to 0.05
clean:
	rm -rf cache output
