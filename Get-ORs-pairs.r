args<-commandArgs(TRUE);
library ("PheWAS")
library("SKAT")
input1=args[1];input2=args[2];wholesets=args[3];anypairs=args[4]; #input1 is plinks, input2 is 'input', input3 is top icd10, setid is the gene sets
#anypairs is for generating setID of a single gene
pair<-read.table(file=paste(anypairs,sep=""),header = F)
sets<-read.table(file=paste(wholesets,sep=""),header=F)
csv.phenotypes=read.csv("id.icd10-check.clean.nopoint.csv",colClasses=c("integer","character","integer"))
#codetable<-read.csv(file="/home/yiming/Documents/icd10cm_codes_2019.csv",header = F)
#library(RODBC)
phenotypes=createPhewasTable(csv.phenotypes,min.code.count = 1,translate = F,add.exclusions = F)
allors<-data.frame();
for (i in 1:nrow(pair)){
onecode=as.character(pair[i,1]);
onegene=as.character(pair[i,2]);
oneset=sets[sets$V1==onegene,]; #checking this
snpsinOneSet=oneset$V2;
write.table(snpsinOneSet,file="temp.snp.list",row.names = F,col.names = F,quote = F);
needcases<-phenotypes[phenotypes[,onecode]==TRUE,c("id",onecode)];
controls<-phenotypes[phenotypes[,onecode]==FALSE,c("id",onecode)];

if(nrow(controls)<(nrow(needcases)*5)){needcontrols=controls}
else{ needcontrols<-controls[controls$id%in%sample(controls$id,size=nrow(needcases)*5,replace = F),]} 

phenocase<-cbind(needcases$id,needcases$id,rep(2,ncol(needcases)))
phenocontrol<-cbind(needcontrols$id,needcontrols$id,rep(1,ncol(needcontrols)))
phenos<-rbind(phenocase,phenocontrol)
write.table(phenos,paste(input1,"pheno-keep-update.txt",sep = ""),row.names = F,col.names = F,quote = F,sep = "\t")

command<-paste("plink1.9 --bfile", paste(input1,collapse=""), "--keep", paste(input1,"pheno-keep-update.txt",sep = ""), "--pheno", paste(input1,"pheno-keep-update.txt",sep=""), "--extract temp.snp.list --make-bed --out", paste(input2,collapse=""))
system(command)

system("bash plinks-for-or.bash")
or=read.table(file="or.txt");
allors=rbind(allors,or);
}
results<-cbind(pair,allors);
write.table(results,file="all-ors-output",row.names=F,col.names=F,quote=F)
