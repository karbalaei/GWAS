Imputation will be done using the Michigan Imputation Server (BIG THANKS TO THEM FOR MAKING OUR LIVES SO MUCH EASIER!!)

https://imputationserver.sph.umich.edu/index.html

Using Reference panel: HRC 1.1 2016 and Eagle v2.3 Phasing

### More info on reference panel here :

https://imputationserver.readthedocs.io/en/latest/reference-panels/ and http://www.haplotype-reference-consortium.org/

### Output of Imputation server

Three files per input chromosomes e.g. for chromosome 21:

```
chr21.dose.vcf.gz.tbi
chr21.dose.vcf.gz
chr21.info.gz
```
The .dose.vcf.gz are your imputed genotypes with dosage information
The .dose.vcf.gz.tbi is the index file of your .vcf.gz file (which some programs require)
The .info.gz file is a file that provided some information on the included variants such as quality and frequency

Example of the .info.gz file (Note we copied and pasted some random variants here)
```
SNP	REF(0)	ALT(1)	ALT_Frq	MAF	AvgCall	Rsq	Genotyped	LooRsq	EmpR	EmpRsq	Dose0	Dose1
9:11739	G	A	0.00004	0.00004	0.99996	0.00023	Imputed	-	-	-	-	-
9:14665	G	A	0.20555	0.20555	0.82608	0.20321	Imputed	-	-	-	-	-
9:62179	G	A	0.22655	0.22655	0.99994	0.99976	Genotyped	0.998	0.999	0.99717	0.99920	0.00078
```
```
Most important columns are:
SNP = variant name and in this case Chromosome:Basepair
REF(0) = Reference allele
ALT(1) = Alternative allele
ALT_Frq = Alternative allele frequency
MAF = Minor Allele frequency
Rsq = This is the imputation quality, depending on who you ask, >0.3, >0.6 or >0.8 is good enough for GWAS.
Genotyped = This says either Imputed or Genotyped. Genotyped means that it was included in the input data.
```

### References:

Imputation server: https://www.ncbi.nlm.nih.gov/pubmed/27571263

HRC panel: https://www.ncbi.nlm.nih.gov/pubmed/27548312

Eagle Imputation: https://www.ncbi.nlm.nih.gov/pubmed/27694958

