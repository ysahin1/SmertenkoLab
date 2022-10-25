for f in $(ls *.fastq | rev | cut -c 9- | rev | uniq);
	do
	 echo "$f file is under progress"
	 hisat2 -t --new-summary --summary ${f}_hisat_summary.txt -x sbv3_referenceGenome -1 ${f}_1.fastq -2 ${f}_2.fastq -S ${f}.sam
done
for s in *.sam;
	do 
		echo "$s is processing"
		samtools sort -@ 4 \
		-o ./BAM/$s.bam \
		$s
done
cd ./BAM
bash featureCounts
