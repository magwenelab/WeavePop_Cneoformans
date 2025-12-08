sample=$1
reference=$2
bin_size="500"
bam_path="/FastData/czirion/persvade_test/data/bam/${sample}.bam"
reference_fasta="/FastData/czirion/persvade_test/data/references/${reference}.fasta"
output="/FastData/czirion/persvade_test/results_cnvpytor"

cores="1"
date_start=$(date +%s)

if [ -f $bam_path ]; then
    if [ -f ${output}/${sample}_${bin_size}_calls.tsv ]; then
        echo "CNVpytor already done. Skipping this step."
    else
        echo "Analysing ${sample}"
        cnvpytor -conf ${output}/${reference}_conf.py -j ${cores} -root ${output}/${sample}.pytor -rd $bam_path
        cnvpytor -conf ${output}/${reference}_conf.py -j ${cores} -root ${output}/${sample}.pytor -his $bin_size
        cnvpytor -conf ${output}/${reference}_conf.py -j ${cores} -root ${output}/${sample}.pytor -partition $bin_size
        cnvpytor -conf ${output}/${reference}_conf.py -j ${cores} -root ${output}/${sample}.pytor -call $bin_size 1> ${output}/${sample}_${bin_size}_calls.tsv
        date_end=$(date +%s)
        cnvpytor -conf ${output}/${reference}_conf.py -j ${cores} -root ${output}/${sample}.pytor -plot rd $bin_size -o ${output}/${sample}_${bin_size}_rdplot.svg
        echo "Finished with ${sample}"

        seconds=$(($date_end - $date_start))
        echo "time_elapsed",$seconds
    fi
else
    echo "Bam file not found"
    exit
fi


# Other CNVpytor commands
# # cnvpytor -conf ${output}/${reference}_conf.py -root ${output}/${sample}.pytor -plot manhatan $bin_size -o ${output}/${sample}_${bin_size}_manahatanplot.svg
# # cnvpytor -conf ${output}/${reference}_conf.py -root ${sample}.pytor -view $bin_size
# # Combined call using RD and SNPs heterozygosity as CNV signal, not useful for homozygous samples (ploidy 1, GT 1/1)
# # cnvpytor -conf ${output}/${reference}_conf.py -j 1 -root ${output}/${sample}.pytor -snp /FastData/czirion/WeavePop_Cneoformans/Crypto_Desjardins/results/01.Samples/snippy/${sample}/snps.vcf -noAD
# # cnvpytor -conf ${output}/${reference}_conf.py -j 1 -root ${output}/${sample}.pytor -pileup /FastData/czirion/WeavePop_Cneoformans/Crypto_Desjardins/results/01.Samples/snippy/${sample}/snps.bam
# # cnvpytor -conf ${output}/${reference}_conf.py -j 1 -root ${output}/${sample}.pytor -baf $bin_size -nomask
# # cnvpytor -conf ${output}/${reference}_conf.py -j 1 -root ${output}/${sample}.pytor -call combined $bin_size