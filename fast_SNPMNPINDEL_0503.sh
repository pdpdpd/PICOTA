#!/bin/bash

#lets make it simple! (maybe slow)
#because I cant remove insertions, right now only allows SNPs 

chrarr=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
indarr=(HG02810 HG03833 HG00109 HG01885 HG01556) #HG02810 HG00109 HG01885 HG01556
cd /home/yp11/Desktop/genomes/1000genome/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/SNPMNPINDEL
#finally changed address and the path
#bcftools norm -m-any testbcf_chr10_HG02810.vcf.gz | bcftools norm -Ov --check-ref w -f ~/Downloads/hg19/chr10.fa > chr10_HG2810_OUT.VCF
#for filtering
#ignoring MT at this point

for ind in ${indarr[@]}
do
  mkdir hap1_${ind}_all
  mkdir hap2_${ind}_all

  for chrs in ${chrarr[@]}
  do
      if [ ! -f ./chr${chrs}_${ind}_all.vcf.gz ]; then
        echo bcftools view -c1 -Oz -s ${ind} -v snps,mnps,indels -o chr${chrs}_${ind}_all.vcf.gz chr${chrs}.vcf.gz
        bcftools view -c1 -Oz -s ${ind} -v snps,mnps,indels -o chr${chrs}_${ind}_all.vcf.gz chr${chrs}.vcf.gz
      fi 
         tabix -f chr${chrs}_${ind}_all.vcf.gz
         cat /home/yp11/Desktop/genomes/hg19/chromFa_numbers/chr${chrs}.fa  | bcftools consensus -e 'ALT~"<.*>"' -H 1pIu -s ${ind} chr${chrs}_${ind}_all.vcf.gz > hap1_${ind}_all/chr${chrs}.fa 
         cat /home/yp11/Desktop/genomes/hg19/chromFa_numbers/chr${chrs}.fa  | bcftools consensus -e 'ALT~"<.*>"' -H 2pIu -s ${ind} chr${chrs}_${ind}_all.vcf.gz > hap2_${ind}_all/chr${chrs}.fa
	 #-e 'ALT~"<.*>"'exclude all symbolic alleles.
	 # 1pIu   2pIu :http://samtools.github.io/bcftools/bcftools.html#consensus

         #cat hap1_${ind}_all_chr${chrs}.fa >> hap1_${ind}_all.fa
         #cat q_hap2_${ind}_all_chr${chrs}.fa >> q_hap2_${ind}_all.fa
done

#Then need to determine if this individual is male or female, but ours are all male so doesn't matter
bcftools view -c1 -Oz -s ${ind} -v snps,mnps,indels -o chrX_${ind}_all.vcf.gz chrX.vcf.gz
tabix chrX_${ind}_all.vcf.gz
bcftools view -c1 -Oz -s ${ind} -v snps,mnps,indels -o chrY_${ind}_all.vcf.gz chrY.vcf.gz
tabix chrY_${ind}_all.vcf.gz
#bcftools view -c1 -Oz -s ${ind} -v snps,mnps,indels -o chrM_${ind}_all.vcf.gz chrM.vcf.gz
#tabix chrM_${ind}_all.vcf.gz

cat /home/yp11/Desktop/genomes/hg19/chromFa_numbers/chrX.fa | bcftools consensus -e 'ALT~"<.*>"' -s ${ind} chrX_${ind}_all.vcf.gz > hap1_${ind}_all/chrX.fa
cat /home/yp11/Desktop/genomes/hg19/chromFa_numbers/chrY.fa | bcftools consensus -e 'ALT~"<.*>"' -s ${ind} chrY_${ind}_all.vcf.gz > hap1_${ind}_all/chrY.fa
#cat /home/yp11/Desktop/genomes/hg19/chromFa_numbers/chrM.fa | bcftools consensus -e 'ALT~"<.*>"' -s ${ind} chrM_${ind}_all.vcf.gz > hap1_${ind}_all/chrM.fa
cp hap1_${ind}_all/chrX.fa hap2_${ind}_all/chrX.fa
cp hap1_${ind}_all/chrY.fa hap2_${ind}_all/chrY.fa
#cp hap1_${ind}_all/chrM.fa hap2_${ind}_all/chrM.fa

#rename and index
bash /home/yp11/Documents/scripts/chrtoPERSONAL_totagger_0504.sh /home/yp11/Desktop/genomes/1000genome/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/SNPMNPINDEL/hap1_${ind}_all ${ind} 1 hg19${ind}hap1 Y
bash /home/yp11/Documents/scripts/chrtoPERSONAL_totagger_0504.sh /home/yp11/Desktop/genomes/1000genome/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/SNPMNPINDEL/hap2_${ind}_all ${ind} 2 hg19${ind}hap2 Y


done

#a quick version with qs. just to make a quick test for circle-seq pipeline. :3
