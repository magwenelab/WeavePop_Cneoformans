
sample=$1
reference=$2
algorithms="HMMcopy,AneuFinder"
window_size=500
ploidy=1
mem=0.1
threads=1

ref="/FastData/czirion/WeavePop_Cneoformans/analyses/data/cnv_software/references/${reference}.fasta"
bam="/FastData/czirion/WeavePop_Cneoformans/analyses/data/cnv_software/bam/${sample}.bam"
wd="/FastData/czirion/WeavePop_Cneoformans/analyses/results/cnv_software/results_HMMcopy_AneuFinder"
output_dir="$wd/${sample}"


date_start=$(date +%s)

cd $wd
mkdir -p $output_dir

# call_CNVs
# if call_CNVs/perSVade_finished_file.txt exists, skip this step
echo "Starting CNV calling module..."
if [ -f "$output_dir/call_CNVs/perSVade_finished_file.txt" ]; then
    echo "CNV calling already done. Skipping this step."
else
    echo "Calling CNVs"
    python /perSVade/scripts/perSVade call_CNVs \
    --fraction_available_mem $mem \
    --threads $threads \
    -o $output_dir/call_CNVs \
    -r $ref \
    --mitochondrial_chromosome no_mitochondria \
    -sbam $bam \
    --window_size_CNVcalling $window_size \
    -p $ploidy \
    --cnv_calling_algs $algorithms
fi

date_end=$(date +%s)
seconds=$(($date_end - $date_start))
echo "time_elapsed",$seconds
