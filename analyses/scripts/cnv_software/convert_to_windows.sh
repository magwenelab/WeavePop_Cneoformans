# CNVpytor
cut -f1 /FastData/czirion/WeavePop_Cneoformans/analyses/data/processed/haploid_samples_with_cnv.tsv  | while read sample
#sample="SRS885916"
do
    tsv='/FastData/czirion/persvade_test/data/depth_by_windows/'${sample}'.tsv'
    cnv='/FastData/czirion/persvade_test/results_cnvpytor/'${sample}'_500_calls.tsv'
    cnv_bed='/FastData/czirion/persvade_test/results_cnvpytor/'${sample}'_500_calls.bed'
    windows_bed='/FastData/czirion/persvade_test/results_cnvpytor/'${sample}'_500_windows.bed'

    if [ -f $windows_bed ]; then
        echo "File present. Skipping CNVpytor for $sample."
    else
        # separate second column(accession:start-end) into accession start end and print columns accession start end first fourth
        awk '{split($2,a,":|-"); print a[1], a[2], a[3], $4, $1}' OFS='\t' $cnv > $cnv_bed
        # intersect a base bed file of windows with the CNV regions to get the calls per window
        awk 'NR > 1 {start = $2 + 1; print $1, start, $3}' OFS='\t' $tsv | bedtools intersect -a stdin -b $cnv_bed -loj | cut -f1,2,3,7,8 > $windows_bed
    fi
done

# AneuFinder
cut -f1 /FastData/czirion/WeavePop_Cneoformans/analyses/data/processed/haploid_samples_with_cnv.tsv | while read sample
# sample="ERS542306"
do
    tsv='/FastData/czirion/persvade_test/data/depth_by_windows/'${sample}'.tsv'
    cnv='/FastData/czirion/persvade_test/results_AneuFinder/'${sample}'/call_CNVs/final_CNVcalling.tab'
    cnv_bed='/FastData/czirion/persvade_test/results_AneuFinder/'${sample}'_final_CNVcalling.bed'
    windows_bed='/FastData/czirion/persvade_test/results_AneuFinder/'${sample}'_500_windows.bed'

    if [ -f $windows_bed ]; then
        echo "File present. Skipping AneuFinder for $sample"
    else
        # select columns with coordinates, depth and CNV call
        awk 'NR > 1 {print $1, $3, $4, $2, $9}' OFS='\t' $cnv > $cnv_bed
        # intersect a base bed file of windows with the CNV regions to get the calls per window
        awk 'NR > 1 {start = $2 + 1; print $1, start, $3}' OFS='\t' $tsv | bedtools intersect -a stdin -b $cnv_bed -loj | cut -f1,2,3,7,8 | sed 's/DEL/deletion/g' | sed 's/DUP/duplication/g' > $windows_bed
    fi
done

# HMMcopy
cut -f1 /FastData/czirion/WeavePop_Cneoformans/analyses/data/processed/haploid_samples_with_cnv.tsv | while read sample
# sample="ERS542306"
do
    tsv='/FastData/czirion/persvade_test/data/depth_by_windows/'${sample}'.tsv'
    cnv='/FastData/czirion/persvade_test/results_HMMcopy/'${sample}'/call_CNVs/final_CNVcalling.tab'
    cnv_bed='/FastData/czirion/persvade_test/results_HMMcopy/'${sample}'_final_CNVcalling.bed'
    windows_bed='/FastData/czirion/persvade_test/results_HMMcopy/'${sample}'_500_windows.bed'

    if [ -f $windows_bed ]; then
        echo "File present. Skipping HMMcopy for $sample."
    else
        # select columns with coordinates, depth and CNV call
        awk 'NR > 1 {print $1, $3, $4, $2, $9}' OFS='\t' $cnv > $cnv_bed
        # intersect a base bed file of windows with the CNV regions to get the calls per window
        awk 'NR > 1 {start = $2 + 1; print $1, start, $3}' OFS='\t' $tsv | bedtools intersect -a stdin -b $cnv_bed -loj | cut -f1,2,3,7,8 | sed 's/DEL/deletion/g' | sed 's/DUP/duplication/g' > $windows_bed
    fi
done

# HMMcopy_AneuFinder
cut -f1 /FastData/czirion/WeavePop_Cneoformans/analyses/data/processed/haploid_samples_with_cnv.tsv | while read sample
# sample="ERS542306"
do
    tsv='/FastData/czirion/persvade_test/data/depth_by_windows/'${sample}'.tsv'
    cnv='/FastData/czirion/persvade_test/results_HMMcopy_AneuFinder/'${sample}'/call_CNVs/final_CNVcalling.tab'
    cnv_bed='/FastData/czirion/persvade_test/results_HMMcopy_AneuFinder/'${sample}'_final_CNVcalling.bed'
    windows_bed='/FastData/czirion/persvade_test/results_HMMcopy_AneuFinder/'${sample}'_500_windows.bed'

    if [ -f $windows_bed ]; then
        echo "File present. Skipping HMMcopy_AneuFinder for $sample."
    else
        # select columns with coordinates, depth and CNV call
        awk 'NR > 1 {print $1, $3, $4, $2, $10}' OFS='\t' $cnv > $cnv_bed
        # intersect a base bed file of windows with the CNV regions to get the calls per window
        awk 'NR > 1 {start = $2 + 1; print $1, start, $3}' OFS='\t' $tsv | bedtools intersect -a stdin -b $cnv_bed -loj | cut -f1,2,3,7,8 | sed 's/DEL/deletion/g' | sed 's/DUP/duplication/g' > $windows_bed
    fi
done