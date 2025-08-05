
#!/bin/bash

metadata_file=$1 #"/FastData/czirion/WeavePop/test/config/metadata.csv"
input_dir=$2 #"test/results"
output_dir=$3 #"test/results_samples"


cut -d',' -f3 $metadata_file | tail -n +2 | while read sample; do
    if [ ! -d "${output_dir}${sample}" ]; then
        mkdir -p ${output_dir}/${sample}
    fi

    ls ${input_dir}/01.Samples/annotation/${sample}/* | while read file; do
        filename=$(basename $file)
        ln -sr $file ${output_dir}/${sample}/${filename}
    done

    ls ${input_dir}/01.Samples/cnv/${sample}/* | while read file; do
        filename=$(basename $file)
        ln -sr $file ${output_dir}/${sample}/${filename}
    done

    ls ${input_dir}/01.Samples/plots/${sample}/* | while read file; do
        filename=$(basename $file)
        ln -sr $file ${output_dir}/${sample}/${filename}
    done

    ls ${input_dir}/01.Samples/depth_quality/${sample}/* | while read file; do
        filename=$(basename $file)
        ln -sr $file ${output_dir}/${sample}/${filename}
    done


    ln -sr ${input_dir}/01.Samples/snippy/${sample} ${output_dir}/${sample}/snippy

done
