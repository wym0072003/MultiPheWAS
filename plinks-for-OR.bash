#!usr/bin/bash
rm temp.*;
plink1.9 --bfile input --filter-controls --recode rlist --out controls;
cat controls.rlist | grep -v "NIL"|awk '{for(i=6;i<=NF;i+=2){print $i"\t"$i}}' >> temp.control.txt;
controlcarrier=`cat temp.control.txt| sort|uniq| wc -l`;
plink1.9 --bfile input --filter-cases --recode rlist --out cases;
cat cases.rlist | grep -v "NIL"|awk '{for(i=6;i<=NF;i+=2){print $i"\t"$i}}' >> temp.cases.txt;
casescarrier=`cat temp.cases.txt|sort|uniq|wc -l`;
allcases=`grep "2$" input.fam| wc -l`;
allcontrols=`grep "1$" input.fam| wc -l`;
caseodds=`awk 'BEGIN{printf "%.8f\n",(('$casescarrier'+0.5)/(('$allcases'-'$casescarrier')+0.5))}'`;
controlodds=`awk 'BEGIN{printf "%.8f\n",(('$controlcarrier'+0.5)/(('$allcontrols'-'$controlcarrier')+0.5))}'`;
awk 'BEGIN{printf "%.4f\n",('$caseodds'/'$controlodds')}' > or.txt;