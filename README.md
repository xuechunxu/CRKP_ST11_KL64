# CRKP_ST11_KL64
# Genome-centric Analysis Workflow for the CRKP_ST11_KL64 paper

This repository contains a collection of code and scripts used in the paper **A point mutation in recC promotes subclonal replacement of carbapenem-resistant Klebsiella pneumoniae ST11 in China** (DOI: XXX) by Kai Zhou, Chun-Xu Xue et al..

Adaptation to selective pressures is crucial for clinically important pathogens to establish epidemics, but the underlying evolutionary drivers remain poorly understood. The current epidemic of carbapenem-resistant Klebsiella pneumoniae (CRKP) poses a significant threat to public health. In this study we analyzed the genome sequences of 794 CRKP bloodstream isolates collected in 40 hospitals in China between 2014 and 2019. We uncovered a subclonal replacement in the predominant clone ST11, where the previously prevalent subclone OL101:KL47 was replaced by O2v1:KL64 over time in a stepwise manner. O2v1:KL64 carried a higher load of mobile genetic elements, and a canonical mutation exclusively detected in the recC of O2v1:KL64 significantly promotes recombination proficiency. The epidemic success of O2v1:KL64 was further bolstered by a selective hypervirulence sublineage with enhanced resistance to phagocytosis, sulfamethoxazole-trimethoprim, and tetracycline. The phenotypic alterations were linked to the overrepresentation of hypervirulence determinants and antibiotic genes conferred by the acquisition of an rmpA-positive pLVPK-like virulence plasmid and an IncFII-type multidrug-resistant plasmid. The dissemination of the sublineage was further promoted by more frequent inter-hospital transmission. The results collectively demonstrate that the expansion of O2v1:KL64 is driven by a repertoire of genomic alterations convergent in a selective population with evolutionary advantages.

The links below the sub-headings lead to the scripts needed for the corresponding steps. Most of the scripts were developed for running on the Huawei FusionServer Pro 5885H V5 server. You may download and adapt the scripts to suit your own requirements.

## 1. Software used in this workflow

- [Perl](https://www.perl.org/)
- [Python3](https://www.python.org/)
- [Trimmomatic](https://github.com/timflutre/trimmomatic)
- [SPAdes](https://github.com/ablab/spades)
- [Prokka](https://github.com/tseemann/prokka)
- [Prodigal](https://github.com/hyattpd/Prodigal)
- [HMMER](http://hmmer.org/)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [RAxML](https://evomics.org/learning/phylogenetics/raxml/)
- [Pandas](https://pandas.pydata.org/)
- [Numpy](https://numpy.org/)
- [CheckM](https://ecogenomics.github.io/CheckM/)
- [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk)

>Take the KP16932 isolate as an example.

## 2. Dataset
All assembled Illumina sequence data have been deposited in GenBank under the BioProject accession number [PRJNA778807](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA778807).

## 3. read trimming
### Trimmomatic
```bash
java -jar trimmomatic-0.36.jar PE -threads 5 KP16932_raw_1.fq.gz KP16932_raw_2.fq.gz KP16932_clean_1.fq.gz KP16932__unpaired_1.fq.gz KP16932_clean_2.fq.gz KP16932__unpaired_2.fq.gz
```

## 4. Assembly
### SPAdes
```bash
spades.py -1 KP16932_clean_1.fq.gz -2 KP16932_clean_2.fq.gz --isolate --cov-cutoff auto -o KP16932.fasta
```
### Unicycle
```bash
unicycler -1 KP16932_1.clean_1.fq.gz -2 KP16932_2.clean_1.fq.gz -l KP16932.nanopore.fq.gz -o KP16932.unicycle.fasta
```

## 5. Taxonomy assignment
### GTDB
```bash
nohup gtdbtk classify_wf --genome_dir fasta_dir/ --out_dir fasta_dir.GTDB.out --extension fasta &
# fasta_dir, the input directory containing a set of genomic assembly sequences.
# fasta_dir.GTDB.out, output directory
```
## 6. Amino acid identity (ANI) calculation
### fastANI
```bash
fastANI --ql quer_genome.list --rl ref_genome.list -o FastANI.out -t 40
```

## 7. Genome annotation
### Prokka
```bash
prokka KP16932.fasta --prefix KP16932 --outdir KP16932.prokka.out/KP16932 --compliant
```

## 8. ST assignment
### Kleborate
```bash
kleborate --all -o kleborate.results.txt -a fasta_dir/*.fasta
# fasta_dir, the input directory containing a set of genomic assembly sequences.
```

## 9. Identification of ARGs,
### Abricate
```bash
mkdir ARG_dir
for f in `ls fasta_dir`; do abricate -db resfinder --nopath --minid 50 --mincov 70 --quiet fasta_dir/${f} > ARG_dir/${f%%.fasta}.tab; done
abricate --nopath --summary ARG_dir/*tab > ARG.tab
# fasta_dir, the input directory containing a set of genomic assembly sequences.
```

## 10. Core SNP Phylogenetic analysis
### Snippy, Gubbins, RAxML
```bash
# call SNPs for multiple isolates from the same reference KP16932.
snippy-multi input.tab --ref KP16932.fa  --cpu 24  > runme.sh
# input.tab, a tab separated input file as follows
# input.tab = ID assembly.fasta
# Isolate	/path/to/contigs.fasta
less runme.sh   # check the script makes sense
sh ./runme.sh   # leave it running over lunch

# remove all the "weird" characters and replace them with N
snippy-clean_full_aln core.full.aln > clean.full.aln 

### Gubbins
# detect recombination region
run_gubbins.py -f 50 -p gubbins clean.full.aln

# remove recombination region
snp-sites -c gubbins.filtered_polymorphic_sites.fasta > clean.core.aln
# -c only output columns containing exclusively ACGT

### RAxML
# build core SNP tree
raxmlHPC -f a -x 12345 -p 12345 -# 100 -m GTRGAMMAX -s clean.core.aln -n tree
```

## 11. Ancestral state reconstruction of mobile genetic element (MGE) number
### Phytools, R
```R
# core_SNP.tre and mge_count.csv can be found in `MGE_ancestral_state` dictionary,
# selected replicon as example
setwd("path/work_dictionary")
library(phytools)

tree <- read.tree("core_SNP.tre")
#the phylogenetic tree built in 10 above.

mge <- read.csv("mge_count.csv",row.names=1) # input MGE number of each isolates.
mge<-as.matrix(mge)[,1] # selected replicon as example

# estimate ancestral states and compute variances & 95% confidence intervals for each node:
fit<-fastAnc(tree,mge,vars=TRUE,CI=TRUE)
fit

# projection of the reconstruction onto the edges of the tree
obj<-contMap(tree,mge,plot=FALSE)
plot(obj,legend=0.7*max(nodeHeights(tree)),
     fsize=c(0.1,0.9), lwd=1, outline = F, leg.txt="Replicons",ftype="off")
# fsize, set font size
# outline, logical value indicated whether or not to outline the plotted color bar with a 1 pt
line.
# leg.txt, set title of legend.
# ftype="off", don't show leave's name.

# OR set colors manually
obj<-setMap(obj,c("red", "#fffc00", "green", "purple", "blue", "#d7ff00", "black"))
plot(obj,legend=0.7*max(nodeHeights(tree)),
     fsize=c(0.1,0.9), lwd=1, outline = F, leg.txt="replicons",ftype="off")
```

## 12. Plasmid coverage across isolates
### Blastn, Seqkit
```bash
# Blastn each genome to plasmid sequence
makeblastdb -in plasmid.fasta -dbtype nucl -parse_seqids -out plasmid_db
mkdir blastn_result
for i in fasta_dir/*.fasta; do blastn -query $i -db plasmid_db -out blastn_result/${i##*/}.blastn.out -outfmt 6; done

# calculate coverage
seqkit fx2tab --length --name --header-line plasmid.fasta # calculate the length of plasmid
cd blastn_result
for i in *.out; do python /home/xuechunxu/pipeline/coverage.py -i $i -l <plasmid length>; done > ../coverage_result.tab
# -l, length of plasmid
```
