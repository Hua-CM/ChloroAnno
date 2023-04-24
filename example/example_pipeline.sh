#!/bin/bash

# iso
# # Two assembly need check the isomer
# |inpath1	| refpath |
# |/mnt/distributed/assembly/chloroplast/202211/result/G430225LY0081_1.fasta |	/mnt/distributed/assembly/chloroplast/annotation/ref/NC_056970.1.fasta |
# |/mnt/distributed/assembly/chloroplast/202211/result/G620802LY0210_1.fasta |	/mnt/distributed/assembly/chloroplast/annotation/ref/NC_035719.1.fasta |


# curate
realpath result/* > curate.info
sed -i '1i\inpath1\' curate.info
mkdir curated
python /home/hzy/software/ChloroAnno/ChloroAnno -t curate -i curate.info -o curated

# correct
# #PGA annotation
mkdir anno_PGA anno_tmp
realpath curated/* > seq.lst
while read line
do
seqname=$(echo "$line"| awk -F "/" '{print $NF}')
perl ~/tools/PGA/PGA.pl \
    -r ~/tools/PGA/plastomes \
    -p 85 \
    -i 100000 \
    -t ${line} \
    -o anno_tmp
mv anno_tmp/*.gb anno_PGA
done<seq.lst

# #CPGAVAS2 annotation
mkdir anno_cpg anno_cpg2
singularity shell -B /mnt:/mnt /home/assembly/tools/CPGAVAS2.sif
while read line
do
    seqname=$(echo "$line"| awk -F "/" '{print $NF}' | cut -f1 -d.)
    singularity exec -B /mnt:/mnt /home/assembly/tools/CPGAVAS2.sif cpgavas-cluster \
        -pid ${seqname} \
        -in ${line} \
        -db 1 \
        -outdir /mnt/distributed/assembly/chloroplast/annotation/anno_cpg
    cp anno_cpg/dir_{seqname}/${seqname}.gbf anno_cpg2/${seqname}.gb
done<seq.lst
while read line
do
    seqname=$(echo "$line"| awk -F "/" '{print $NF}' | cut -f1 -d.)
    cp anno_cpg/dir_${seqname}/${seqname}.gbf anno_cpg2/${seqname}.gb
done<seq.lst

# #correct
# #Output tbl format to make it easier for upload.
mkdir corrected
python /home/hzy/software/ChloroAnno/ChloroAnno -t correct -i demo.info --ouffmt tbl -o corrected

