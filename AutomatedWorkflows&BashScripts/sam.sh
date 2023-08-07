for s in *.sam;
	do 
		echo "$s is processing"
		samtools sort -@ 4 \
		-o ./BAM/$s.bam \
		$s
done
cd ./BAM
bash featureCounts

for f in $(ls ../StressInSpikes/SAM/*.sam);
 do
  var1=${f}
  output=${var1##*/}
  samtools view -Sb $f > ../StressInSpikes/Bam/${output}.bam
  samtools view -h -b -f 4 ../StressInSpikes/Bam/${output}.bam > ../StressInSpikes/UnmappedBam/unmapped_${output}.bam
  samtools view -h -b -q 20 ../StressInSpikes/Bam/${output}.bam > ../StressInSpikes/FilteredBam/${output}.bam
done

