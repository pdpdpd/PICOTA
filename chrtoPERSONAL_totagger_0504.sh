#!/bin/bash



#for our data that doesn't have Ref or NC..
#input 0: dic of target site
#input 1: #sample
#input 2: #hap
#input 3: species(default:Homo sapiens)
#input 4: has Y? (default:Y)
#PWD="/home/yp11/Desktop/genomes/mm10/mm10_tagger/"
PWD=${1:-/home/yp11/Desktop/genomes/mm10/mm10_tagger/}
samplen=${2:-default}
   hapn=${3:-1}
species=${4:-Homosapiens}
  haveY=${5:-Y}

#determine if the input has chromosome Y
if [ "$haveY" == 'Y' ]; then
 Filename=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
#should be MTfor chrMT..
else
 Filename=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
fi

cd $PWD

name=${samplen}'_'${hapn}'_all_chr.fa'
  
echo -n >  $name

i=0

for x in "${Filename[@]}"
do
 i=$((i + 1))

#the file to be change
 ofile='chr'$x'.fa'

#change the header and combine into a huge file with the listed order
 ref_position='>gi|0000000'$i'|chr'$x':NC_000000.'$i'|'$species

 sed -i "1s/.*/$ref_position/" $ofile
 
#gi|000000000|chrI:NC_000000.1|c elegans ce10 WS220

 cat $ofile >> $name

done

pwd
#tagscan:
 echo " ~/Desktop/useful_tools/tagger-1.3/genwin -a ${name/_all_chr.fa/_idx.txt} -o ${name/_all_chr.fa/_out.bin} $name"
 ~/Desktop/useful_tools/tagger-1.3/genwin -a ${name/_all_chr.fa/_idx.txt} -o ${name/_all_chr.fa/_out.bin} $name
 echo " ~/Desktop/useful_tools/tagger-1.3/sortGWI -o ${name/_all_chr.fa/_idx.bin} -u ${name/_all_chr.fa/_out.bin}"
 ~/Desktop/useful_tools/tagger-1.3/sortGWI -o ${name/_all_chr.fa/_idx.bin} -u ${name/_all_chr.fa/_out.bin}
 echo " ~/Desktop/useful_tools/tagger-1.3/fetchGWI -a ${name/_all_chr.fa/_idx.txt} -i ${name/_all_chr.fa/_idx.bin} -m 1 CCTGTCGTCCTTGAAGAA"
 ~/Desktop/useful_tools/tagger-1.3/fetchGWI -a ${name/_all_chr.fa/_idx.txt} -i ${name/_all_chr.fa/_idx.bin} -m 1 CCTGTCGTCCTTGAAGAA

#${f/Mus_musculus.GRCm38.75.dna.chromosome./chr}



#useit  % compile the hg19 genome using tagscan
#system('..\genwin -a ../hg19_idx.txt -o ../hg19_out.bin hg19_all_chr.fa');
#system('..\sortGWI -o ../hg19_idx.bin -u ../hg19_out.bin');

#% test use of it
#system('..\fetchGWI -a ../hg19_idx.txt -i ../hg19_idx.bin -m 2 CCTGTCGTCCTTGAAGAA');

#return
