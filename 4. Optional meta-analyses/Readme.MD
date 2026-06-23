# 4. Optional meta-analyses
Meta-analyses are used when you have multiple datasets.

## Prep of your GWAS results (based on RVTEST file)

merge all GWAS files per chromosome and merge all .info files
```
cat *.SingleWald.assoc | grep -v 'N_INFORMATIVE' > allChrs_FILE.assoc 
cat maf001rsq03minimums_chr*.info | grep -v 'Rsq' > allChrs_FILE.Info
```
Then in reformat in R
```
R
# note1 you can speed this up with the data.table package and changing read.table and write.table to fread and fwrite
# note2 we set here beta filtering at <5 and >-5 since typically these are unrealistic beta's coming from GWAS
infos <- read.table(paste("allChrs_FILE.Info"))
colnames(infos) <- c("SNP","ALT_Frq","Rsq")
assoc <- read.table(paste("allChrs_FILE.assoc"))
colnames(assoc) <- c("CHROM","POS","REF","ALT","N_INFORMATIVE","Test","Beta","SE","Pvalue")
data <- merge(infos, assoc, by.x = "SNP", by.y = "Test", all.y = T)
dat <- subset(data, Beta < 5 & Beta > -5 & !is.na(data$Pvalue))
dat$chr <- paste("chr",dat$CHROM, sep = "")
dat$markerID <- paste(dat$chr,dat$POS, sep = ":")
dat$minorAllele <- ifelse(dat$ALT_Frq <= 0.5, as.character(dat$ALT), as.character(dat$REF))
dat$majorAllele <- ifelse(dat$ALT_Frq <= 0.5, as.character(dat$REF), as.character(dat$ALT))
dat$beta <- ifelse(dat$ALT_Frq <= 0.5, dat$Beta, dat$Beta*-1)
dat$se <- dat$SE
dat$maf <- ifelse(dat$ALT_Frq <= 0.5, dat$ALT_Frq, 1 - dat$ALT_Frq)
dat$P <- dat$Pvalue
dat0 <- dat[,c("markerID","minorAllele","majorAllele","beta","se","maf","P")]
write.table(dat0, file=paste("toMeta.FILE.tab"), quote = F, sep = "\t", row.names = F)
```

## Create a metal file 
Create a metal file that looks like below and save it as my_metal_analysis.txt

```
#../generic-metal/metal metalAll.txt
#THIS SCRIPT EXECUTES AN ANALYSIS OF EIGHT STUDIES
#THE RESULTS FOR EACH STUDY ARE STORED IN FILES Inputfile1.txt THROUGH Inputfile8.txt
SCHEME  STDERR
AVERAGEFREQ ON
MINMAXFREQ ON
LABEL TotalSampleSize as N # If input files have a column for the sample size labeled as 'N'
# LOAD THE FIRST SEVEN INPUT FILES

# UNCOMMENT THE NEXT LINE TO ENABLE GenomicControl CORRECTION
# GENOMICCONTROL ON

# === DESCRIBE AND PROCESS THE FIRST INPUT FILE ===
MARKER markerID
ALLELE minorAllele majorAllele
FREQ   maf
EFFECT beta
STDERR se
PVALUE P
WEIGHT N 
PROCESS toMeta.GWAS1.tab

# === DESCRIBE AND PROCESS THE SECOND INPUT FILE ===
MARKER markerID
ALLELE minorAllele majorAllele
FREQ   maf
EFFECT beta
STDERR se
PVALUE P
WEIGHT N
PROCESS toMeta.GWAS2.tab


# === DESCRIBE AND PROCESS THE SIXTH INPUT FILE ===
MARKER markerID
ALLELE minorAllele majorAllele
FREQ   maf
EFFECT beta
STDERR se
PVALUE P
WEIGHT N
PROCESS toMeta.GWAS3.tab


OUTFILE MY_THREE_GWAS_FILES_META .tbl
ANALYZE HETEROGENEITY

QUIT
``` 

## Run metal

Then run metal like this:

```
metal metal.txt
```

## Output will look something like this:
These are a couple of lines from one a recent GWAS, including 17 datasets. 

```
MarkerName      Allele1 Allele2 Freq1   FreqSE  MinFreq MaxFreq Effect  StdErr  P-value Direction       HetISq  HetChiSq        HetDf   HetPVal TotalSampleSize
chr4:90666041   t       c       0.6036  0.0229  0.5505  0.6295  0.6980  0.1169  2.348e-09       +-++++++-++++++++       40.4    26.827  16      0.04344 138030
chr3:8379719	a	g	0.2688	0.0101	0.2534	0.2884	0.0006	0.1582	0.997	-+-?+?-+??--++--+	0.0	8.536	12	0.742 140900

```
A couple of notes:

For the Direction column three options are possible +,- and ?. ? means the variant was not present in that dataset, + means positive beta and - means negative beta value.
The HetDf value is the number of included datasets for the variant -1 
Note the it is also recommend to filter for HetISq, anything higher than 80 is not reliable.

Also optional sorting and filtering:
```
# column 10 is the p-value column
sort -gk 10 META_ALL_DATA1.tbl > SORTED_META_ALL_DATA1.tbl
```

### References:

METAL: https://www.ncbi.nlm.nih.gov/pubmed/20616382
