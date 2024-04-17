# MLPheWAS (Multi-level PheWAS)

A flexiable tool running phenome-wide association study at gene and pathway level.  
MLPheWAS simply runs as one line command of R in linux environment.

## Introduction
The accumulation of large-scale genome sequencing datasets and electronic medical records (EMR) information has enabled an unprecedented potential to discover associations between genetics and disease by phenome-wide association studies (PheWAS). Current PheWAS approaches perform association tests at the variant level, which is powered mostly for using common variants. Therefore, the automated aggregation and association of rare and common high impact variants into genes and pathways with EMR phenotypes is likely to significantly enhance novel PheWAS discoveries identify etiologies of human diseases.  
We present Multi-level PheWAS (MLPheWAS), a user-friendly and highly effective framework automating variant-, gene-, and pathway-level burden  tests , PheWAS, and principal component analysis (PCA) for population-specific associations   using common and rare variants as well as user-defined covariates. 

### Flowchart
<img src="/img/gitlab-flowchart.png"  width="1000" height="500">

MLPheWAS reads genotypes, phenotypes, gene-variant pairs as input for conducting variant-, gene- and pathway level phenome-wide association analysis by simply leveraging the variant set file. It runs PCA prior to the PheWAS to generate principal components as covariates in running phewas to control the population stratification. Moreover，it supports any other covariates defined by users as well. A high resolution plot in PDF format will be generated for quick displaying the top association other than the raw results table. Meanwhile, it's flexible for extending usages by senior users as SKAT and PheWAS functions are fully retained in MLPheWAS.  
## Installation for dependencies
- Installing PheWAS in R
```
 install.packages(c("dplyr","tidyr","ggplot2","MASS","meta","ggrepel","DT"))
 devtools::install_github("PheWAS/PheWAS")
```
- Installing SKAT in R
```
install.packages("SKAT")
```
- Installing getopt in R
```
install.packages("getopt")
```
- Installing PLINK 1.9 in Linux
```
sudo apt-get install plink/1.90
```

## Usage

> **MLPheWAS**  

MLPheWAS running simply as one-line commmand in linux:
```
Rscript MLPheWAS.R -i genotypes -p phenotypes -s variant-sets -m methods -c covariates -o output-filename
-i | --input            : This is the prefix of the plink files
-p | --pheno            : This is the phenotype code table
-s | --set              : This is the file name of the prepared setID
-m | --method           : Algorithm for test, 'SKAT','Burden' or 'SKATO',SKATO is defult setting
-c | --covariates       : Specifying the '-c covariates.txt' as the user defined covariates, or '-c suppress' to not include any covariates.   
                          defult setting is to run PCA prior to conducting phewas.
-o | --output           : Specifying a name for output files

-m and -c are optional parameters, if they are not specified, the defaul setting will be used.
```
A running template:  
`Rscript MLPheWAS.R -i test-inputs -p simulated-phenotypes.csv -s simulated-variant-sets.setID -m SKAT -c suppress -o MyOutputs`  

'test-inputs' files in the repository of master branch. 'test-inputs' are plinks of 1000 genome samples (hg38) with filtered sites for reducing data size.   
'simulated-phenotypes.csv' and 'simulated-variant-sets.setID' can be found in repository of main branch. Both files were created for users for testing. 

### Input files
  - --input is the prefix of plink binary files,for example "Input1", the files include Input1.bed, Input1.bim and Input1.fam. The plink files should include as many   variants as possible for the purpose to population structure inference, no need to prune or filter out variants.
  - --pheno Specifying a table of phenotypes for samples. It is the diagnosis information (encoded as ICD10,ICD9,phecodes) of individuals which was summarized from Biobank, it needs to be prepared as a csv file including three columns with headers, as shown below:
```
---------------------------
id,icd10,count
Sample1_ID,D369,5
Sample1_ID,M7070,9
Sample1_ID,K5909,1
Sample1_ID,J40,1
Sample1_ID,R609,8
Sample1_ID,Z794,42
Sample2_ID,I739,41
Sample2_ID,I25118,9
.....
---------------------------
The diagnosis information should be in comma-delimited format. It should include a header for IDs of samples, phenotype related ICD10 and the count of ICD10 summarized from Biobank.
The sample ID should be consistent with the individual ID (IID) in Input1.fam.
The ICD10 codes were used for sorting out cases and controls for each phenotype.
The count of ICD10 can be employed for collecting reliable cases, the default cutoff was set to 1. Changing values of 'min.code.count = 1' in line 47 of MLPheWAS.R to adjust the cutoff, more detail can be found from 'createPhewasTable' function in PheWAS package.

```
  
  - --set is the variant set, users can leverage this set file for the purpose of conducting gene-level or pathway-level phenome wide tests. 
      For gene level PheWAS, the variant set file inlcudes two columns with no header, the first column is the gene names, the second column is the filtered variants (prepared by users) for each gene.

```
-------gene level----------
ABHD17A 19:1881248:C:T
ABHD17A 19:1881278:T:C
ABHD17A 19:1881371:A:T
ABHD17A 19:1881377:C:T
ACAT2   6:159762979:A:G
ACAT2   6:159767081:G:T
ACAT2   6:159768515:C:G
ACAT2   6:159768613:C:T
.....
---------------------------
```
For pathway level test, the first column is the name of pathway, the second column is variants belonging to genes in pathways, but the gene info is omitted here.
```
------Pathway level--------
Pathway1 19:1881248:C:T
Pathway1 19:1881278:T:C
Pathway1 19:1881371:A:T
Pathway1 19:1881377:C:T
Pathway1  6:159762979:A:G
Pathway1  6:159767081:G:T
Pathway1  6:159768515:C:G
Pathway1  6:159768613:C:T
...
PathwayN  ...
---------------------------
```
  - --output Specifying a prefix of the output file which can be defined as any name by user, a csv file named as 'Input5-MLPhewas-Output.csv' will be created after successfully running MLPheWAS.

- Other necessay input files for running MLPheWAS.

"high-LD-region-hg38.txt", the high LD regions based on genome build of hg38, this file should be located the same directory with MLPheWAS.R.  
"high-LD-region-hg19.txt", the high LD regions based on genome build of hg19, this file should be located the same directory with MLPheWAS.R.  
"high-LD-region-hg38.txt" is the default input in line 37 of MLPheWAS, users can choose the a proper one by modifying the file name in line 37.  
"icd10cm_codes_2019.csv", this file inldudes interpretations of each ICD10 code, which help to add detailed phenotype information for each association. It should be located at the same directory of MLPheWAS.R
### Outputs
- --output Specifying a name for output files.

For example, using "My-outputs", then we will get:
- "My-outputs-MLPhewas-Output.csv",which includes results of gene-level or pathway level association tests, will be created for all results after successfully running MLPheWAS.  
Interpretation of headers in 'My-outputs-MLPhewas-Output.csv':  
```
codeinfo: Descriptions of the each phenotype.
codes: The ICD10 for each phenotype.
numcases: The number of cases for the phenotype summarized from 'Code-table.csv'
numcontrols: The number of controls used for association test.
SetID: The names of the variant sets, they are euqal to the gene name or pathway ID in MLPheWAS.
P.value: P value in SKATO test.
N.Marker.All: a number of SNPs in the genotype matrix.
N.Marker.Test: a number of SNPs used for the test. It can be different from N.Marker.All when some markers are monomorphic or have higher missing rates than the missing_cutoff.
MAC: total minor allele count.
m: the number of individuals with minor alleles.
Method.bin: a type of method to compute a p-value (default="Hybrid"),see more detail in manuels of SKAT.
MAP: the minimum possible p-values. It is available when the method.bin="ER" and m is sufficiently small, ,see more detail in manuels of SKAT.
Bonferroni: "Ture" or "FALSE" indicating if a assoication passed the Bonferroni-adjusted significance (Phenomoe wide)
```

- "My-outputs-tops-in-PheWAS.pdf"  
A pdf file for displaying results created by successfully running MLPheWAS, it displays the top 50 associations in results accoring to p-values. The number of displayed items can be specified in line 120 of MLPheWAS.R.
The red dash line in plot denotes the phenome wide significance.

### Optional parameters
Users have some options in running MLPheWAS other than the functions in defaul setting. 

- --m  Specifying what test is going to be used: 'SKAT', 'Burden' or 'SKATO'.  
    User can choose the algorithms of SKAT, Burden or optimal SKAT(SKATO). 
    SKATO will be used if not specified.

- --covariates Running with user defined covariates  
    User can sepcify their own covariates using this parameter. The covariate file should be named as 'covariates.txt'.
    The format of the covariate.txt should be tab-delimited txt, the first two columns should be the same with --input *fam file. From the  3rd column are values of covariates.

  

> ****Geting ORs of variant sets****  

We provide function for calculating the gene-level and pathway level OR for any set of interested variants.  
```
Rscript Get-ORs-pairs.r input1 input input3 pairs-interested 
```
- Inputs  
'Input1' and 'input3' are the same with above input files for --input and --set respectively, 'input' is extract the word "input", which is the name of temp file for counting variant carriers using plink. 'pairs-interested' is a tab delimited txt file including ICD10 and genes as two columns. For each pair, we will get the OR of the gene for the phenotype (ICD10).  
'plinks-for-OR.bash' should be placed at the same directory with Get-ORs-pairs.r
- Output  
'all-ors-output', a tab delimited file gives ORs for each phenotype-gene pair.


## Computing environment
The MLPheWAS running in linux platforms.  
MLPheWAS was developed and tested within following computing environment:  
Ubuntu 18.04  
RAM: 512GB.  >128GB RAM is recommended.  
CPU: Intel(R) Xeon(R) Gold 6148 CPU @ 2.40GHz  
R 4.1.0  
SKAT 2.0.1 in R  
PheWAS 0.99.5.5 in R  
PLINK 1.90  
# contract
Email: yiming.wu@mssm.edu  
Principal investigator: Yuval Itan (yuval.itan@mssm.edu)  
Institude: The Charles Bronfman Institute for Personalized Medicine, Icahn School of Medicine at Mount Sinai, New York, NY, USA
