#!/bin/bash

#$1 is fasta file for chrM (in one sequence)
#$2 is the read length
#$3 is the path to first read fastq
#$4 if the path to second read fastq

#This script is used to generate transposition frequence bigwig files for ATACseq data
#Orignal data required: fastq file(s) from ATAcseq, fasta file for chrM in a single sequence
fasta=0
read_1=0
read_2=0
read_length=0
dnase=0
threads=0
output=0
scale=0

usage="CircAlign\n---------\n\tRequired Arguments:\n\t\t-f [FASTA_FILE] Fasta file with Chromosome M in one sequence\n\t\t-1 [FASTQ_FILE] Fastq file for read (Read 1 if Paired End)\n\t\t-r [READ_LENGTH] Length of reads found in [FASTQ_FILE]\n\tOptional Arguments:\n\t\t-2 [FASTQ_FILE] Fastq file for read 2 (Paired end data)\n\t\t-o [OUTPUT_DIRECTORY] Directory to place output results into (./ is default)\n\t\t-t [THREAD_#] Number of threads to use\n\t\t-d This will run the data as DNase (ATACseq is default)\n\t\t-s This will force the BigWigs generated to not be scaled (Default is to scale by read number aligned to [FASTA_FILE])\n"

printf "\n"
while getopts ":hf:1:2:r:o:t:ds" opt; do
    case $opt in
        h)
            printf "$usage\n"
            exit 0
            ;;
        f)
            printf -- "-Fasta file for Chromosome M: $OPTARG\n"
            fasta=$OPTARG
            ;;
        1)
            printf -- "-Fastq file (Read 1): $OPTARG\n"
            read_1=$OPTARG
            ;;
        2)
            printf -- "-Fastq file (Read 2): $OPTARG\n"
            read_2=$OPTARG
            ;;
        o)
            printf -- "-Output Directory is: $OPTARG\n"
            output=${OPTARG%*/}
            ;;
        r)
            printf -- "-Read length is: $OPTARG\n"
            read_length=$OPTARG
            ;;
        t)
            printf -- "-Threads to use: $OPTARG\n"
            threads=$OPTARG
            ;;
        d)
            printf -- "-Data type set to DNase\n"
            dnase=1
            ;;
        s)
            printf -- "-BigWig files will not be scaled\n"
            scale=1
            ;;
        \?)
            printf "\nError: Invalid option: -$OPTARG. \n$usage"
            exit 1
            ;;
        :)
            printf "\nError: Option -$OPTARG requires an argument.\n"
            exit 1
            ;;
    esac
done
shift $((OPTIND-1))

if [[ ! -f $fasta ]]; then
    printf "Error: fasta file for chrM does not exist. Please check path and run again.\n"
    exit 1
fi
if [[ ! -f $read_1 ]]; then
    printf "Error: Sequence file for read 1 does not exist. Please check path and run again.\n"
    exit 1
fi
if [[ $read_2 != 0 && ! -f $read_2 ]]; then
    printf "Error: Sequence file provided for read 2 does not exist. Please check path and run again or run without read 2.\n"
    exit 1
fi
if [[ $read_length == 0 ]]; then
    printf "Error: Read length was not provided. Please specify read length with '-r' and run again.\n"
    exit 1
fi
if [[ $output == 0 ]]; then
    output=$(pwd)
    printf -- "-Output will be located in ${output}\n"
fi
if [[ $threads == 0 ]]; then
    printf -- "-Threads to use: 1\n"
    threads=1
fi

#step 0: setting enviroment
mkdir -p "$output"
source /home/segil_lab/.profile_litao
fasta_file_name="${fasta##*/}"

#step 1: circulating chrM by shifting 2x(read length) bp of chrM fasta from start to the end

temp="${output}/temp.fasta"
original="${output}/${fasta_file_name%.*}_original.fasta"
shifted="${output}/${fasta_file_name%.*}_shifted.fasta"
chrM_genome="${output}/chrM.genome.txt"

sed -e 's/\(^>.*$\)/#\1#/' "$fasta" | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' | sed 1d > "${temp}"
echo ">chrM" | cat - "$temp" > "$original"
seq_length=$(tail -1 "$original"| wc -c)
seq_length=$((${seq_length} - 1))
bwa index "$original"

shiftlength=$((${read_length}+${read_length}))
fragment1=$(cut -c1-$shiftlength "${temp}")
fragment2=$(cut -c $shiftlength- "${temp}" | cut -c 2-)
echo -e ">chrM\n${fragment2}${fragment1}" > "$shifted"
bwa index "$shifted"

#make a genome file just for chrM
genomeLength=$(head -2 "$original" | tail -1 | awk '{print length($0)}')
echo -e "chrM\t${genomeLength}" > "$chrM_genome"
rm "${temp}"

#step2: align reads to chrM and shifted-chrM using bwa
original_align="${output}/originalAlign"
shifted_align="${output}/shiftedAlign"
mkdir -p "$original_align"
mkdir -p "$shifted_align"

original_read1_sai="${original_align}/read_1_align.sai"
original_read2_sai="${original_align}/read_2_align.sai"
original_sample_sam="${original_align}/sample_align.sam"
original_sample_bed="${original_align}/sample_align.bed"
original_forward_bgr="${original_align}/sample_cov_forward.bgr"
original_reverse_bgr="${original_align}/sample_cov_reverse.bgr"

shifted_read1_sai="${shifted_align}/read_1_align.sai"
shifted_read2_sai="${shifted_align}/read_2_align.sai"
shifted_sample_sam="${shifted_align}/sample_align.sam"
shifted_sample_bed="${shifted_align}/sample_align.bed"
shifted_forward_bgr="${shifted_align}/sample_cov_forward.bgr"
shifted_reverse_bgr="${shifted_align}/sample_cov_reverse.bgr"

# Due to small genome, bwa can only use about 4 cores anyway so giving 8 is useless and would be better used split
if [[ $threads -gt 7 ]]; then
    half_threads=$(echo "$threads/2" | bc)
    bwa aln -t "$half_threads" "$original" "$read_1" > "$original_read1_sai" &
    bwa aln -t "$half_threads" "$shifted" "$read_1" > "$shifted_read1_sai" &
    wait
else
    bwa aln -t "$threads" "$original" "$read_1" > "$original_read1_sai"
    bwa aln -t "$threads" "$shifted" "$read_1" > "$shifted_read1_sai"
fi


if [[ $read_2 != 0 ]]; then
    # Due to small genome, bwa can only use about 4 cores anyway so giving 8 is useless and would be better used split
    if [[ $threads -gt 7 ]]; then
        half_threads=$(echo "$threads/2" | bc)
        bwa aln -t "$half_threads" "$original" "$read_2" > "$original_read2_sai" &
        bwa aln -t "$half_threads" "$shifted" "$read_2" > "$shifted_read2_sai" &
        wait
    else
        bwa aln -t "$threads" "$original" "$read_2" > "$original_read2_sai"
        bwa aln -t "$threads" "$shifted" "$read_2" > "$shifted_read2_sai"
    fi

    if [[ $threads > 1 ]]; then
        bwa sampe "$original" "$original_read1_sai" "$original_read2_sai" "$read_1" "$read_2" > "$original_sample_sam" &
        bwa sampe "$shifted" "$shifted_read1_sai" "$shifted_read2_sai" "$read_1" "$read_2" > "$shifted_sample_sam" &
        wait
    else
        bwa sampe "$original" "$original_read1_sai" "$original_read2_sai" "$read_1" "$read_2" > "$original_sample_sam"
        bwa sampe "$shifted" "$shifted_read1_sai" "$shifted_read2_sai" "$read_1" "$read_2" > "$shifted_sample_sam"
    fi
else
    if [[ $threads > 1 ]]; then
        bwa samse "$original" "$original_read1_sai" "$read_1" > "$original_sample_sam" &
        bwa samse "$shifted" "$shifted_read1_sai" "$read_1" > "$shifted_sample_sam" &
        wait
    else
        bwa samse "$original" "$original_read1_sai" "$read_1" > "$original_sample_sam"
        bwa samse "$shifted" "$shifted_read1_sai" "$read_1" > "$shifted_sample_sam"
    fi
fi

rm -f "${orignal}*" "${shifted}*"

#step3: make transposition frequence bigwig files original alignment
if [[ $dnase == 1 ]]; then
    if [[ $threads -gt 7 ]]; then
        half_threads=$(echo "$threads/2" | bc)
        samtools view -@ "$half_threads" -b "$original_sample_sam" | bamToBed -i stdin > "$original_sample_bed" &
        samtools view -@ "$half_threads" -b "$shifted_sample_sam" | bamToBed -i stdin > "$shifted_sample_bed" &
        wait
    else
        samtools view -@ "$threads" -b "$original_sample_sam" | bamToBed -i stdin > "$original_sample_bed"
        samtools view -@ "$threads" -b "$shifted_sample_sam" | bamToBed -i stdin > "$shifted_sample_bed"
    fi
else
    if [[ $threads -gt 7 ]]; then
        half_threads=$(echo "$threads/2" | bc)
        samtools view -@ "$half_threads" -b "$original_sample_sam" | bamToBed -i stdin | awk '{if($6=="+") print $1"\t"$2+4"\t"$3"\t"$4"\t"$5"\t"$6; else print $1"\t"$2"\t"$3-4"\t"$4"\t"$5"\t"$6}' > "$original_sample_bed" &
        samtools view -@ "$half_threads" -b "$shifted_sample_sam" | bamToBed -i stdin | awk '{if($6=="+") print $1"\t"$2+4"\t"$3"\t"$4"\t"$5"\t"$6; else print $1"\t"$2"\t"$3-4"\t"$4"\t"$5"\t"$6}' > "$shifted_sample_bed" &
        wait
    else
        samtools view -@ "$threads" -b "$original_sample_sam" | bamToBed -i stdin | awk '{if($6=="+") print $1"\t"$2+4"\t"$3"\t"$4"\t"$5"\t"$6; else print $1"\t"$2"\t"$3-4"\t"$4"\t"$5"\t"$6}' > "$original_sample_bed"
        samtools view -@ "$threads" -b "$shifted_sample_sam" | bamToBed -i stdin | awk '{if($6=="+") print $1"\t"$2+4"\t"$3"\t"$4"\t"$5"\t"$6; else print $1"\t"$2"\t"$3-4"\t"$4"\t"$5"\t"$6}' > "$shifted_sample_bed"
    fi

fi

junction_bed="${output}/junction_align.bed"
split_junc_bed="${output}/split_junction_align.bed"
junction_id="${output}/junction_ids.txt"
complete_bed="${output}/complete_sample_align.bed"
complete_bgr="${output}/complete_sample_align.bgr"
complete_bw="${output}/complete_sample_align.bw"

# bedtools intersect -v -a "$shifted_sample_bed" -b "$original_sample_bed" > "$junction_bed"
junction_posistion=$((${seq_length} - ${shiftlength}))
awk -F "\t" -v var="$junction_posistion" '{if($2 < var && $3 > var) print}' "$shifted_sample_bed" > "$junction_bed"
awk -F "\t" '{print $4}' "$junction_bed" > "$junction_id"
awk -F "\t" -v var="$seq_length" -v var2="$shiftlength" 'BEGIN{OFS = "\t"}{if($3 > (var - var2)) print $1, $2 + var2, var, $4, $5, $6"\n"$1, "1", ($3 + var2) - var, $4, $5, $6; print $0}' "$junction_bed" > "$split_junc_bed"
grep -v -F -f "$junction_id" "$original_sample_bed" | cat - "$split_junc_bed" > "$complete_bed"

rm -f "$junction_id" "$junction_bed" "$split_junc_bed"

if [[ $scale == 1 ]]; then
    scaleFactor=1
else
    readNumber=$(wc -l "$original_sample_bed" | awk '{print $1}')
    scaleFactor=$(echo "1000000/${readNumber}" | bc -l)
fi

echo "Scale factor for bigwig is ${scaleFactor}"
bedtools genomecov -bg -scale "$scaleFactor" -strand + -5 -i "$original_sample_bed" -g "$chrM_genome" > "$original_forward_bgr"
bedtools genomecov -bg -scale "$scaleFactor" -strand - -5 -i "$original_sample_bed" -g "$chrM_genome" > "$original_reverse_bgr"

bedtools genomecov -bg -scale "$scaleFactor" -strand + -5 -i "$shifted_sample_bed" -g "$chrM_genome" > "$shifted_forward_bgr"
bedtools genomecov -bg -scale "$scaleFactor" -strand - -5 -i "$shifted_sample_bed" -g "$chrM_genome" > "$shifted_reverse_bgr"

bedtools genomecov -bg -scale "$scaleFactor" -i "$complete_bed" -g "$chrM_genome" > "$complete_bgr"

#step4: fill the gap of original alignment by data from shifted alignment

head="${output}/head.bed"
head_shifted="${output}/head_shifted.bed"
tail="${output}/tail.bed"
tail_shifted="${output}/tail_shifted.bed"

forward_bgr="${output}/forward.bgr"
forward_bw="${output}/forward.bw"
forward_head="${output}/forward_head.bgr"
forward_middle="${output}/forward_middle.bgr"
forward_tail="${output}/forward_tail.bgr"

reverse_bgr="${output}/reverse.bgr"
reverse_bw="${output}/reverse.bw"
reverse_head="${output}/reverse_head.bgr"
reverse_middle="${output}/reverse_middle.bgr"
reverse_tail="${output}/reverse_tail.bgr"

a=$(echo "${genomeLength}-${read_length}+1" | bc)
b=$(echo "${a}-${read_length}-1" | bc)
c=$(echo "${a}-1" | bc)
d=$(echo "${b}-${read_length}+1" | bc)
echo -e "chrM\t0\t${read_length}" > "$head"
echo -e "chrM\t${a}\t${genomeLength}" > "$tail"
echo -e "chrM\t${b}\t${c}" > "$head_shifted"
echo -e "chrM\t${d}\t${b}" > "$tail_shifted"

bedtools intersect -v -a "$original_forward_bgr" -b "$head" | bedtools intersect -v -a stdin -b "$tail" > "$forward_middle"
bedtools intersect -a "$shifted_forward_bgr" -b "$head_shifted" | awk -v var="$b" '{print $1"\t"$2-var"\t"$3-var"\t"$4}' > "$forward_head"
bedtools intersect -a "$shifted_forward_bgr" -b "$tail_shifted" | awk -v var="$shiftlength" '{print $1"\t"$2+var"\t"$3+var"\t"$4}' > "$forward_tail"
cat "$forward_head" "$forward_middle" "$forward_tail" > "$forward_bgr"
rm -f "$forward_head" "$forward_middle" "$forward_tail"

bedtools intersect -v -a "$original_reverse_bgr" -b "$head" | bedtools intersect -v -a stdin -b "$tail" > "$reverse_middle"
bedtools intersect -a "$shifted_reverse_bgr" -b "$head_shifted" | awk -v var="$b" '{print $1"\t"$2-var"\t"$3-var"\t"$4}' > "$reverse_head"
bedtools intersect -a "$shifted_reverse_bgr" -b "$tail_shifted" | awk -v var="$shiftlength" '{print $1"\t"$2+var"\t"$3+var"\t"$4}' > "$reverse_tail"
cat "$reverse_head" "$reverse_middle" "$reverse_tail" > "$reverse_bgr"
rm -f "$reverse_head" "$reverse_middle" "$reverse_tail"
rm -f "$head" "$tail" "$head_shifted" "$tail_shifted"

#step5: generate bigwig file from bedGraph
bedGraphToBigWig "$forward_bgr" "$chrM_genome" "$forward_bw"
bedGraphToBigWig "$reverse_bgr" "$chrM_genome" "$reverse_bw"
bedGraphToBigWig "$complete_bgr" "$chrM_genome" "$complete_bw"

rm -f "$chrM_genome"
