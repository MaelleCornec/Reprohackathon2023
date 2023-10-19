cd data_echantillons

cd perister_replicate_1
trim_galore -q 20 --phred33 --length 25 SRR10379721_1.fastq
cd -

cd perister_replicate_2
trim_galore -q 20 --phred33 --length 25 SRR10379722_1.fastq
cd -

cd perister_replicate_3
trim_galore -q 20 --phred33 --length 25 SRR10379723_1.fastq
cd -

cd control_1
trim_galore -q 20 --phred33 --length 25 SRR10379724_1.fastq
cd -

cd control_2
trim_galore -q 20 --phred33 --length 25 SRR10379725_1.fastq
cd -

cd control_3
trim_galore -q 20 --phred33 --length 25 SRR10379726_1.fastq
cd -

