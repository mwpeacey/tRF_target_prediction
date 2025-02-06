#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=4G
#$ -N miranda_wrapper
#$ -o miranda_wrapper_output.txt
#$ -e miranda_wrapper_output.txt

SMALL_RNA_FASTA=$1
TRANSCRIPTOME_FASTA=$2
RUN_ID=$3
RUN_MODE=$4

cd /grid/schorn/home/mpeacey/tRF_targets/data/miranda_output

mkdir ${RUN_ID}
cd ${RUN_ID}

cat ${SMALL_RNA_FASTA} | grep '>' > temp_sRNA_list.txt

i=1
for line in $(cat temp_sRNA_list.txt); do

        result=${line#*>}

        if [ $i = 1 ]; then

                echo $result > sRNA_list.txt

        else

            	echo $result >> sRNA_list.txt

        fi

	((i=i+1))

done

rm temp_sRNA_list.txt

cat ${TRANSCRIPTOME_FASTA} | grep '>' > temp_lines.txt

i=1
for line in $(cat temp_lines.txt); do

        result=${line#*>}

        if [ $i = 1 ]; then

               	echo $result > transcript_list.txt

        else

            	echo $result >> transcript_list.txt

        fi

        ((i=i+1))

done

rm temp_lines.txt

for sRNA in $(cat sRNA_list.txt); do

	cd /grid/schorn/home/mpeacey/tRF_targets/data/miranda_output/${RUN_ID}

	echo "Submitting job for sRNA ${sRNA}..."

	mkdir ${sRNA}
	cd ${sRNA}

	cat ${SMALL_RNA_FASTA} | sed -n "/${sRNA}/,/>/p" | head -2 > temp_${sRNA}.fasta

	qsub -N ${RUN_ID}_${sRNA}_miranda \
	-o ${RUN_ID}_${sRNA}_miranda_output.txt -e ${RUN_ID}_${sRNA}_miranda_output.txt \
	/grid/schorn/home/mpeacey/scripts/miranda/miranda_v3/miranda_v3.sh \
	${SMALL_RNA_FASTA} ${TRANSCRIPTOME_FASTA} ${RUN_ID} \
	${sRNA} ${RUN_MODE}

done
