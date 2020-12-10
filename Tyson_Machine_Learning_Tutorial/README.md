# Steps for Tutorial </br>
Note: For easy installation of packages, this tutorial assumes you have conda installed. Visit here: https://docs.conda.io/en/latest/miniconda.html first if this is not the case </br> </br>
1. Download the fasta file in this github directory </br>
2. Create a Conda environment </br>
```conda create -n machineLearningTutorial -c bioconda snp-sites vcftools mafft```</br>
3. Run MAFFT to align the genomes</br>
```mafft AandBGenomes.fasta > aligned.fasta```</br>
4. Run snp-sites to generate a vcf file of your SNP locations</br>
```snp-sites -v -o vcf.tsv aligned.fasta```</br>
5. Run vcftools to filter your vcf by minor allele frequency>0.05</br>
```vcftools --vcf vcf.tsv --maf 0.00005 --recode --out vcf_0.00005maf.tsv```</br>
6. Rename the output for simplicity</br>
```mv vcf_0.00005maf.tsv.recode.vcf ncbivcf0.05MAF.tsv```</br>
7. Open ```ncbivcf0.05MAF.tsv``` in your favorite text editor (ie vim, nano) and remove the first three lines as well as the ```#``` at the beginning of the fourth line. Save.</br>
8. Run the R script in this github directory </br>
### Congrats! You just did a simple machine learning exercise to predict clade with SNP data!
