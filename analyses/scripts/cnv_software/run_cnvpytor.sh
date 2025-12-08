# Generate results_cnvpytor/{lineage}_conf.py file before running this
# conda activate cnvpytor
cd /FastData/czirion/persvade_test/CNVpytor
mkdir -p /FastData/czirion/WeavePop_Cneoformans/analyses/results/cnv_software/results_cnvpytor/logs_cnvpytor

for reference in VNI VNII VNBI VNBII
do
    if [ -f "/FastData/czirion/WeavePop_Cneoformans/analyses/results/cnv_software/results_cnvpytor/$reference.pytor" ]; then
        echo "Reference pytor file present. Skipping this step."
    else
        echo "Creating pytor file for $reference"
        cnvpytor --root /FastData/czirion/WeavePop_Cneoformans/analyses/results/cnv_software/results_cnvpytor/$reference.pytor \
            -gc /FastData/czirion/WeavePop_Cneoformans/analyses/data/cnv_software/references/$reference.fasta -make_gc_file \
            &> "/FastData/czirion/WeavePop_Cneoformans/analyses/results/cnv_software/results_cnvpytor/logs_cnvpytor/$reference.log"
    fi
done

parallel -j 60 --colsep '\t' -a /FastData/czirion/WeavePop_Cneoformans/analyses/results/tables/haploid_samples_with_cnv.tsv \
'bash /FastData/czirion/WeavePop_Cneoformans/analyses/scripts/cnv_software/cnvpytor_commands.sh {1} {3} &> \
/FastData/czirion/WeavePop_Cneoformans/analyses/results/cnv_software/results_cnvpytor/logs_cnvpytor/{1}.log'