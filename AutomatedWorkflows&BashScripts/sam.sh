for s in *.sam;
	do 
		echo "$s is processing"
		samtools sort -@ 4 \
		-o ./BAM/$s.bam \
		$s
done
cd ./BAM
bash featureCounts
