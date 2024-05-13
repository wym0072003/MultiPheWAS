############################
# Install required packages#
############################
#-------------------------------------------#
##1. Installing PheWAS in R##
# install.packages(c("dplyr","tidyr","ggplot2","MASS","meta","ggrepel","DT"))
# devtools::install_github("PheWAS/PheWAS")
#-------------------------------------------#

#-------------------------------------------#
##2. Installing SKAT in R##
# install.packages("SKAT")
#-------------------------------------------#

#-------------------------------------------#
##3. Installing getopt in R##
# install.packages("getopt")
#-------------------------------------------#

#-------------------------------------------#
##3. Installing ggplot2 in R##
# install.packages("ggplot2")
#-------------------------------------------#

#-------------------------------------------#
##4. Installing PLINK 1.9 in Linux##
# sudo apt-get install plink/1.90
#-------------------------------------------#



library(getopt)

spec <- matrix(

  c("input", "i", 2, "character", "This is the prefix of the plink files",
    "pheno", "p", 2, "character",  "This is the phenotype code table",
    "set", "s", 2, "character",  "This is the file name of the prepared setID",
	"method", "m", 1, "character",  "Algorithm for test, 'SKAT','Burden' or 'SKATO',SKATO is defult setting",
	"covariates", "c", 1, "character",  "Specifying the 'covariates.txt' as the user defined covariates, or 'suppress' to not include any covariates",
    "output", "o", 2, "character",  "Specifying a name for output files",
    "help",  "h", 0, "logical",  "This is Help"),
  byrow=TRUE, ncol=5)

opt <- getopt(spec=spec)
if( is.null(opt$input) || is.null(opt$set) || is.null(opt$pheno) ){
  cat("Missing core input files,please check names following -i -p -s -o\n")
  cat(paste(getopt(spec=spec, usage = T), "\n"))
  quit()
}
input1=opt$input;input2=opt$pheno;setid=opt$set;output=opt$output;covs=opt$covariates;test=opt$method;
command<-paste(
  "plink1.9 --bfile", paste(input1,collapse=""), "--exclude high-LD-region-hg38.txt --range --indep-pairwise 50 5 0.2 --out temp.in",";",
  "plink1.9 --bfile", paste(input1,collapse=""), "--extract temp.in.prune.in --make-bed --out temp.indep",";",
  "plink1.9 --bfile temp.indep --geno 0.05 --maf 0.02 --make-bed --out pca-input",";",
  "plink1.9 --bfile pca-input --pca --out pca-output",";","mv pca-output.eigenvec covariates-PCA.txt")
if(is.null(covs)){system(command)}else{}
if(covs == "covariates.txt"){} else if (covs == "suppress"){}else{system(command)}
#####################################
#           Running PheWAS          #
#####################################

library ("PheWAS")
library("SKAT")
library("ggplot2")
csv.phenotypes=read.csv(paste(input2,collapse=""),colClasses=c("character","character","integer"))
codedict<-read.csv(file="./icd10cm_codes_2019.csv",header = F)
phenotypes=createPhewasTable(csv.phenotypes,min.code.count = 1,translate = F,add.exclusions = F)
allsamples<-read.table(file=paste(input1,".fam",sep=""))
outneed<-data.frame();
for (i in 2:ncol(phenotypes)){ 
  needcases<-phenotypes[phenotypes[,i]==TRUE,c(1,i)]
  controls<-phenotypes[phenotypes[,i]==FALSE,c(1,i)]
  set.seed(99)
  if(nrow(controls)<(nrow(needcases)*5)){needcontrols=controls} #!Ratio of cases to controls.
  else{ needcontrols<-controls[controls$id%in%sample(controls$id,size=nrow(needcases)*5,replace = F),]} #!Ratio of cases to controls.
  #Agjusting the number in "(needcases)*5" in above 2 lines for the ratio of cases to controls. Eg. here 5 means the ratio of cases to controls is 1:5.
  phenocase<-cbind(allsamples[allsamples$V2 %in% needcases$id,c(1,2)],rep(2,nrow(needcases)));colnames(phenocase)=c("fid","iid","pheno");
  phenocontrol<-cbind(allsamples[allsamples$V2 %in% needcontrols$id,c(1,2)],rep(1,nrow(needcontrols)));colnames(phenocontrol)=c("fid","iid","pheno");
  phenos<-rbind(phenocase,phenocontrol)
  write.table(phenos,paste(input1,"pheno-keep-update.txt",sep = ""),row.names = F,col.names = F,quote = F,sep = "\t")
  command<-paste("plink1.9 --bfile", paste(input1,collapse=""), "--keep", paste(input1,"pheno-keep-update.txt",sep = ""), "--pheno", paste(input1,"pheno-keep-update.txt",sep=""), "--make-bed --out intermediate")
  system(command)
  ##skat/burden/skatO
  File.Bed<-"intermediate.bed";
  File.Bim<-"intermediate.bim";
  File.Fam<-"intermediate.fam";
  File.SetID<-setid;
  File.SSD<-"intermediate.SSD";
  File.Info<-"intermediate.SSD.info";
  Generate_SSD_SetID(File.Bed,File.Bim,File.Fam,File.SetID,File.SSD,File.Info);
  FAM<-Read_Plink_FAM(File.Fam,Is.binary = TRUE)
  SSD.INFO<-Open_SSD(File.SSD,File.Info)
  y<-FAM$Phenotype;
  if(covs == "suppress"){obj<-SKAT_Null_Model(y~1,out_type = "D")
  }else if(covs == "covariates.txt"){  File.Cov<-"./covariates.txt";
  FAM_Cov<-Read_Plink_FAM_Cov(File.Fam,File.Cov,Is.binary = T,cov_header = FALSE);ncov=ncol(FAM_Cov)-6;
  obj<-SKAT_Null_Model(y~noquote(paste("FAM_Cov$COV",1:ncov,sep = "",collapse = "+")),out_type = "D")
  }else{File.Cov<-"./covariates-PCA.txt";
  FAM_Cov<-Read_Plink_FAM_Cov(File.Fam,File.Cov,Is.binary = T,cov_header = FALSE);
  obj<-SKAT_Null_Model(y~FAM_Cov$COV1+FAM_Cov$COV2+FAM_Cov$COV3+FAM_Cov$COV2+FAM_Cov$COV4+FAM_Cov$COV5,out_type = "D")}; 
  
  if(is.null(test)){out.skato<-SKATBinary.SSD.All(SSD.INFO,obj,method="SKATO")
  }else{out.skato<-SKATBinary.SSD.All(SSD.INFO,obj,method=paste(test,collapse=""))}   
  codes<-colnames(phenotypes[,i]);
  numcases<-sum(y==1);numcontrols<-sum(y==0);
  need=cbind(codes,numcases,numcontrols,out.skato$results);
  outneed<-rbind(outneed,need)
  rm(needcases,phenocontrol)
}
Bonferroni<-outneed$P.value<(0.05/ncol(codedict));
codeinfo<-codedict[match(outneed$codes,codedict$V1),2];
finalneed<-cbind(codeinfo,outneed,Bonferroni);
write.csv(finalneed,paste(output,"-MLPhewas-Output.csv",sep = ""),quote = FALSE)

#####################################
##Generating plot for PheWAS results#
#####################################
data<-read.csv(file=paste(output,"-MLPhewas-Output.csv",sep = ""),header = T)
data<-data[data$numcases>=20,] 
dataInTop<-head(data[order(data$P.value),],n=50);V3<-matrix(nrow(dataInTop):1,byrow = 1); #"n=50" shows the top 50 associations. Modifying the number to show the items user like.
combinedcodeinfo=paste(dataInTop$codeinfo,dataInTop$SetID,sep = " -- ")
dataInTop<-data.frame(dataInTop,combinedcodeinfo)
log10=-log10(dataInTop$P.value)
phenowideSig=-log10(0.05/length(unique(data$codes)))
finaldata<-data.frame(dataInTop$combinedcodeinfo,log10,V3)
colnames(finaldata)<-c("V1","V2","V3")
library("ggplot2")
pdf(file =paste(output,"-tops-in-PheWAS.pdf",sep = ""),height = 8,width = 11)
ggplot(finaldata,aes(x=V2,y=reorder(V1,V3))) +geom_point(size=2,aes(colour=V1))+theme_bw()+
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey60",linetype = "dashed",size = 0.2),
        axis.text.y = element_text(size=11,colour = "black",face="italic")
  )+xlab("-log10(P-value)")+ylab("Phenotypes (ICD10) -- Gene")+geom_vline(xintercept =phenowideSig,linetype="dashed",color="red",size=0.2)
dev.off()
