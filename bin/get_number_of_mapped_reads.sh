function getMappedReads {
	bamfile=$1
	sample=$2
    outFile=$3
	echo "Sample;num_reads;num_mapped_reads" > $outFile
    numReads=$(samtools view -c $bamfile)
    numMappedReads=$(samtools view -c -F 260 $bamfile)
	echo "$sample;$numReads;$numMappedReads" >> $outFile 
}
getMappedReads $1 $2 $3