# From: https://github.com/magwenelab/Crypto_Lineage_Assignment_Pipeline/blob/main/generate_sample_names.sh

# Script to clean reads using fastp

# conda activate fastp

# Usage: bash fastp_clean.sh <sample_names_file> <suffix_1> <suffix_2> <out> <reports>
# <sample_names_file> is a file containing the sample names
# <suffix_1> is the suffix of the forward reads
# <suffix_2> is the suffix of the reverse reads
# <out> is the output directory
# <reports> is the directory to save the reports
# <jobs> is the number of jobs to run in parallel

# Check if the correct number of arguments is provided
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <sample_names_file> <suffix_1> <suffix_2> <out> <reports> <jobs>"
    exit 1
fi

# Assign the command line argument to a variable
sample_names_file=$1
suffix_1=$2
suffix_2=$3
in=$4
out=$5
reports=$6
jobs=$7

# Check if the provided file exists
if [ ! -f $sample_names_file ]; then
    echo "File not found!"
    exit 1
fi

# Check if the output directory exists
if [ ! -d $out ]; then
    mkdir -p $out
fi

# Check if the reports directory exists
if [ ! -d $reports ]; then
    mkdir -p $reports
fi


parallel -j $jobs  "fastp \
--dont_overwrite \
--dedup \
--in1 $in/{}$suffix_1 \
--in2 $in/{}$suffix_2 \
--out1 $out/{}_cleaned$suffix_1 \
--out2 $out/{}_cleaned$suffix_2 \
-h $reports/{}_fastp.html \
-j $reports/{}_fastp.json"  \
:::: $sample_names_file