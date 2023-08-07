for fs in $(find ../StressInSpikes/StressInspikesRawFastq/ -type f -name "*.fq.gz");
	do
	 echo "$fs file is under progress"
	 for f in $fs; do
	  var1=${f}
          output=${var1##*/}
	  echo "${output}"
	  hisat2 -p 20 -t --new-summary --summary ../StressInSpikes/HisatSummary/${output}.txt -x ../wheat_reference_genome_2.1v/ref/treasv2 -U ${f} -S ../StressInSpikes/SAM/${output}.sam
         done
done
