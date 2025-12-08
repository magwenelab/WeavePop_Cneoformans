# CNV software benchmarking

1) Get list of samples of Desjardins and Ashton with significant CNV content
2) Link references
```
mkdir data scripts notebooks
ln -s /FastData/czirion/WeavePop_Cneoformans/Crypto_Desjardins/data/references data/references
```
3) Link tables depth_by_windows.tsv
```
grep "Desjardins" /FastData/czirion/WeavePop_Cneoformans/analyses/results/tables/haploid_samples_with_cnv.tsv | cut -f1 | bash /FastData/shared/WeavePop_global/symlinks2weavepop.sh -d /FastData/czirion/WeavePop_Cneoformans/Crypto_Desjardins/results/01.Samples/cnv -f depth_by_windows.tsv -m f -r s -o analyses/data/cnv_software/depth_by_windows
grep "Ashton" /FastData/czirion/WeavePop_Cneoformans/analyses/results/tables/haploid_samples_with_cnv.tsv | cut -f1 | bash /FastData/shared/WeavePop_global/symlinks2weavepop.sh -d /FastData/czirion/WeavePop_Cneoformans/Crypto_Ashton/results/01.Samples/cnv -f depth_by_windows.tsv -m f -r s -o analyses/data/cnv_software/depth_by_windows
```
3) Link snps.bam
```
grep "Desjardins" /FastData/czirion/WeavePop_Cneoformans/analyses/results/tables/haploid_samples_with_cnv.tsv | cut -f1 | bash /FastData/shared/WeavePop_global/symlinks2weavepop.sh -d /FastData/czirion/WeavePop_Cneoformans/Crypto_Desjardins/results/01.Samples/mapping_and_variants -f snps.bam -m f -r s -o analyses/data/cnv_software/bam
grep "Desjardins" /FastData/czirion/WeavePop_Cneoformans/analyses/results/tables/haploid_samples_with_cnv.tsv | cut -f1 | bash /FastData/shared/WeavePop_global/symlinks2weavepop.sh -d /FastData/czirion/WeavePop_Cneoformans/Crypto_Desjardins/results/01.Samples/mapping_and_variants -f snps.bam.bai -m f -r s -o analyses/data/cnv_software/bam
grep "Ashton" /FastData/czirion/WeavePop_Cneoformans/analyses/results/tables/haploid_samples_with_cnv.tsv | cut -f1 | bash /FastData/shared/WeavePop_global/symlinks2weavepop.sh -d /FastData/czirion/WeavePop_Cneoformans/Crypto_Ashton/results/01.Samples/mapping_and_variants -f snps.bam -m f -r s -o analyses/data/cnv_software/bam
grep "Ashton" /FastData/czirion/WeavePop_Cneoformans/analyses/results/tables/haploid_samples_with_cnv.tsv | cut -f1 | bash /FastData/shared/WeavePop_global/symlinks2weavepop.sh -d /FastData/czirion/WeavePop_Cneoformans/Crypto_Ashton/results/01.Samples/mapping_and_variants -f snps.bam.bai -m f -r s -o analyses/data/cnv_software/bam
```
# CNVpytor

4) Create config files for each lineage for CNVpytor: `results/{lineage}_confg.py`  

5) Run `bash analyses/scripts/cnv_software/run_cnvpytor.sh`  


# PerSVade
## Installation and test
### Installation with Singularity
```
mamba create -n singularity singularity
conda activate singularity
singularity build --docker-login ./mikischikora_persvade_v1.02.6.sif docker://mikischikora/persvade:v1.02.6
```
### Example commands
List contents of the /perSVAde directory inside the container:
```
singularity exec -e mikischikora_persvade_v1.02.6.sif ls /perSVade
```
Run perSVade help command:
```
singularity exec -e mikischikora_persvade_v1.02.6.sif bash -c 'source /opt/conda/etc/profile.d/conda.sh && conda activate perSVade_env && python /perSVade/scripts/perSVade --help'
```

### Run perSVade test installation

```
mkdir perSVade_testing_outputs #0
singularity exec \ #1
    -B ./perSVade_testing_outputs:/perSVade/installation/test_installation/testing_outputs \ #2
    -e mikischikora_persvade_v1.02.6.sif \ #3
    bash -c \ #4
    'source /opt/conda/etc/profile.d/conda.sh && \ #5
    conda activate perSVade_env && \ #6
    /perSVade/installation/test_installation/test_installation_modules.py' \ #7
    &> test.log
```
0: Make the directory.   
1: Execute the singularity command.  
2: Bind the output directory to the container.  
3: Specify the singularity image file.  
4: Use bash to run the command.  
5: Source the conda profile script.  
6: Activate the perSVade environment.  
7: Run the test installation script.  

## Run for Desjardins and Ashton samples with heterogeneous depth

6) Execute `bash run_persvade.sh` in the singularity environment.

## Run notebook

7) Render `quarto render scripts/cnv_software/agreement_comparison.qmd`