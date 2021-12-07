#!/bin/bash

# read from command line the unfiltered and unsortde bam file

MODE=$1
SAMPLE=$2
p1=$3
BLACK=$4


out1prefix=$(echo $p1 | sed 's/\.bam$//')
out1="${out1prefix}.sorted.bam"

samtools view -u -q 30 "${p1}" | sambamba sort  --out ${out1} /dev/stdin
samtools index $out1
#samtools rmdup

out2=$(echo ${out1} | sed 's/\.bam$/.nodup.bam/')


picard MarkDuplicates VERBOSITY=WARNING \
	INPUT=${out1} OUTPUT=${out2} \
	METRICS_FILE="${SAMPLE}_dups.log" \
	REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
# index
samtools index $out2

samtools idxstats $out2 > ${SAMPLE}_idxstats.log

out2m=$(echo $out1 | sed 's/\.bam$/.nodup.noM.temp.bam/')
out3=$(echo $out1 | sed 's/\.bam$/.nodup.noM.bam/')
export CHROMOSOMES=$(samtools view -H $out2 | grep '^@SQ' | cut -f 2 | grep -v -e _ -e chrM -e 'VN:' | sed 's/SN://' | xargs echo)

if [[ $MODE = "PE" ]]
then 
  samtools view -b -h -f 3 -F 4 -F 8 -F 256 -F 1024 -F 2048 -q 30 $out2 $CHROMOSOMES > ${out3}	
else
  samtools view -b -h -F 4 -F 256 -F 1024 -F 2048 -q 30 $out2 $CHROMOSOMES > ${out3}
fi
samtools index $out3
out4=$(echo $out1 | sed 's/\.bam$/.nodup.noM.black.bam/')
bedtools subtract -A -a ${out3} -b $BLACK > ${out4}	
samtools index ${out4}

finalPrefix=$(basename "$p1" .bam)
sambamba sort --out ${SAMPLE}_final.bam ${out4}

if [[ $MODE = "PE" ]]
then
	picard CollectInsertSizeMetrics \
		I=${SAMPLE}_final.bam \
		O=${SAMPLE}_insertSizes.log \
		H=${SAMPLE}_insertHist.pdf
fi

rm $out1 $out2 $out3 $out4 *bai 

