SAMPLE_NAMES=$(cat milich2021-ftp-urls.txt | rev | cut -d "/" -f 1 | rev)
BAM_DIR=data/bams/

echo $SAMPLE_NAMES
for sample in "${SAMPLE_NAMES}"
do
  echo ${sample[2]}
done


# for file in "${BAM_DIR}*"
# do
#   # echo ${file}
#   TEXT={$(echo $file | tr "/" "\n" | head )}
#   for i in "${!TEXT[@]}"
#   do
#     echo "${TEXT[$i]}"
#     if [[ "${TEXT[$i]}" == "bam" ]]
#     then
#       echo "${i}"
#     fi
#   done
# done
