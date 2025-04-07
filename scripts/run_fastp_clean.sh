# Clean Desjardins dataset
cd /FastData/czirion/Crypto_Diversity_Pipeline/
tail -n +2 Crypto_Desjardins/config/metadata.csv | cut -d',' -f2 > Crypto_Desjardins/config/samples.txt
conda activate fastp
bash scripts/fastp_clean.sh Crypto_Desjardins/config/samples.txt _1.fq.gz _2.fq.gz Crypto_Desjardins/data/samples Crypto_Desjardins/data/samples_clean Crypto_Desjardins/data/fastp_reports 16 &> Crypto_Desjardins/data/fastp_clean.log

# Clean Ashton dataset
tail -n +2 Crypto_Ashton/config/metadata.csv | cut -d',' -f1 > Crypto_Ashton/config/samples.txt
conda activate fastp
bash scripts/fastp_clean.sh Crypto_Ashton/config/samples.txt _1.fq.gz _2.fq.gz Crypto_Ashton/data/samples Crypto_Ashton/data/samples_clean Crypto_Ashton/data/fastp_reports 16 &> Crypto_Ashton/data/fastp_clean.log