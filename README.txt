0. Installing required packages

0.1 Installing PheWAS in R
 install.packages(c("dplyr","tidyr","ggplot2","MASS","meta","ggrepel","DT"))
 devtools::install_github("PheWAS/PheWAS")

0.2 Installing SKAT in R
 install.packages("SKAT")

0.3 Installing ggplot2 in R
 install.packages("ggplot2")

0.4 Installing PLINK 1.9 in Linux
 sudo apt-get install plink/1.90

1. Runing of MLPheWAS
MLPheWAS running simply as one-line commmand in linux:

Rscript MLPheWAS.R Input1 Input2 Input3 Input4

2. Input files

2.1 Input1 is the prefix of plink binary files including Input1.bed, Input1.bim and Input1.fam. The plink files should include as many variants as possible for the purpose to population structure inference, no need to prune or filter out variants.

2.2 Input2 can be any name as the intermediate plink files generated for each cases versus controls.

2.3 Input3 is the variant set, users can leverage this set file for the purpose of conducting gene-level or pathway-level phenome wide tests.
For gene level PheWAS, the variant set file inlcudes two columns with no header, the first column is the gene names, the second column is the filtered variants (prepared by users) for each gene.
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

For pathway level test, the first column is the name of pathway, the second column is variants belonging to genes in pathways, but the gene info is omitted here.
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

2.4 Input4 is the prefix of the output file which can be defined as any name by user, a csv file named as 'Input5-MLPhewas-Output.csv' will be created after successfully running MLPheWAS.

2.5 "Code-table.csv"
It is the diagnosis information (encoded as ICD10) of individuals which was summarized from Biobank, it needs to be prepared as a csv file including three columns with headers, as shown below:

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


3. The Outputs.

3.1 "Input5-MLPhewas-Output.csv",which includes results of gene-level or pathway level association tests, will be created for all results after successfully running MLPheWAS.

Interpretation of headers in 'Input5-MLPhewas-Output.csv',
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

3.2 "input4-tops-in-PheWAS.pdf"

A pdf file will be created after successfully running MLPheWAS, it displays the top 50 associations in results accoring to p-values. The number of displayed items can be specified in line 94 of MLPheWAS.R. 

4. Other necessay input files working with MLPheWAS.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"high-LD-region-hg38.txt", the high LD regions based on genome build of hg38, this file should be located the same directory with MLPheWAS.R.
"high-LD-region-hg19.txt", the high LD regions based on genome build of hg19, this file should be located the same directory with MLPheWAS.R.
"high-LD-region-hg38.txt" is the default input in line 37 of MLPheWAS, users can choose the a proper one by modifying the file name in line 37.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"icd10cm_codes_2019.csv", this file inldudes interpretations of each ICD10 code, which help to add detailed phenotype information for each association. It should be located at the same directory of MLPheWAS.R. Moreover, the ICD10 in this table should be in same format with the ICD10 in input3. 

5. Some options in MLPheWAS.
Users have some options in running MLPheWAS other than the functions in defaul setting. Some important parts have been mentioned in source code as comments.
5.1 SKAT or Burden test. Users can choose the algorithms of SKAT, Burden or optimal SKAT(SKATO) in line 76 of MLPheWAS.
5.2 Running without any covariates
Comment out line 31-36,72,74,meanwhile uncomment line 75.
5.3 Running with user defined covariates.
Comment out line 31-36, replace the 'covariates-PCA.txt' with the user defined covariates, modifying the number of 'FAM_Cov$COV' in line 74 to mathch with the number of covariates.
5.4 More usages can be digged up by checking manuals of SKAT and PheWAS packages.

6. Computing environment
The MLPheWAS running in linux platforms.
MLPheWAS was developed and tested within following computing environment:
Ubuntu 18.04
RAM: 512GB.  >128GB RAM is recommended.
CPU: Intel(R) Xeon(R) Gold 6148 CPU @ 2.40GHz
R 4.1.0
SKAT 2.0.1 in R
PheWAS 0.99.5.5 in R
PLINK 1.90



