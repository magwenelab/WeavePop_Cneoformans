
# Command to get the start time and end time of all jobs of a rule

time_to_seconds () {
h=$(echo $1 | cut -d":" -f1)
m=$(echo $1 | cut -d":" -f2)
s=$(echo $1 | cut -d":" -f3)
echo $(( "${h#0}" * 3600 + "${m#0}" * 60 + "${s#0}" ))
}

# Mosdepth
log="Crypto_Desjardins/logs/weavepop.log"

grep "localrule mosdepth" -A6 -B1  $log \
| grep -E 'jobid|\[' | sed -z 's/]\n/,/g' | cut -d" " -f4,10 | while read line; 
do jobid=$(echo $line | cut -d" " -f2); 
start_1=$(echo $line |cut -d" " -f1); 
end_1=$(grep -B1 "Finished jobid: $jobid" $log | \
grep -v "Finished" | cut -d" " -f4); 
start=$(time_to_seconds $start_1)
end=$(time_to_seconds $end_1)
echo $(( $end - $start )) ; 
done > analyses/results/tables/mosdepth_time.csv

# CNV calling
grep "localrule cnv_calling" -A6 -B1 $log | \
grep -E 'jobid|\[' | sed -z 's/]\n/,/g' | cut -d" " -f4,10 | while read line; 
do jobid=$(echo $line | cut -d" " -f2); 
start_1=$(echo $line |cut -d" " -f1); 
end_1=$(grep -B1 "Finished jobid: $jobid" $log | \
grep -v "Finished" | cut -d" " -f4); 
start=$(time_to_seconds $start_1)
end=$(time_to_seconds $end_1)
echo $(( $end - $start )) ; 
done > analyses/results/tables/cnv_time.csv

# PerSVade and CNVpytor
software_results="" #FIXME
for folder $software_results
do
grep "time_elapsed" $folder/logs_*/*.log
done | while read line
do
software=$(echo $line | cut -d"/" -f1)
sample=$(echo $line | cut -d"/" -f3 | cut -d"." -f1)
time=$(echo $line | cut -d"," -f2)
echo $software,$sample,$time 
done > analyses/results/tables/software_time.csv

