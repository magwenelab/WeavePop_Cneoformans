#/usr/bin/env xonsh

# Bash code to combine fastqs into one. Needs to be translated to Xonsh:
mkdir fastqs_combined
cat PRJNA382844_SRStoSRR.txt | while read line
do 
    SRS=$(echo $line | cut -d' ' -f1)
    SRR=$(echo $line | cut -d' ' -f2)
    echo "Processing $SRS"
    cat fastqs/${SRR}_1.fastq >> fastqs_combined/${SRS}_1.fastq
    cat fastqs/${SRR}_2.fastq >> fastqs_combined/${SRS}_2.fastq
done

# Count lines in combined fastqs and compare them to the sum of the combined files line count