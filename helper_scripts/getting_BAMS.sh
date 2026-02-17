#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# SRA to BAM pipeline using BWA-MEM2
#
# This script depends on bashbone. Activate bashbone (More: https://github.com/Hoffmann-Lab/bashbone) before running, 
# for example: source /opt/bashbone/latest/bashbone/activate.sh -c true
#
# Requirements:
#	 - idx: BWA index for GRCh38_v106. (Not provided)
###############################################################################

source <path/of/installation>/latest/bashbone/activate.sh -c true

# -----------------------
# User configuration
# -----------------------

threads=56
get_cmd="wget -c -q --timeout=60 --waitretry=10 --tries=10 --retry-connrefused -O /dev/stdout"
idx=../Data/GRCh38.fa.bwa.idx/bwa
do_rmdup=true
outdir="../Data/ATACseq/"

mkdir -p "$outdir"

# -----------------------
# Build commands
# -----------------------

declare -a cmds=()
while read -r id; do
	commander::makecmd -a cmds -v id -v threads -v get_cmd -v idx -v do_rmdup -v outdir -c <<-'CMD'
		echo "working on $id"		
		read -r l n < <(fastq-dump --split-3 --defline-seq '@$ac.$si[.$sn] \$ri' --defline-qual '+' -X 1 --stdout "$id" 2>/dev/null | parallel --pipe --tee {} ::: "head -2 | tail -1 | wc -c" "wc -l" | xargs echo)		
		minscore=$(echo $l | awk '{print $1-sprintf("%.0d",5/100*$1+1)*5}')		
		url=$([[ $(echo -n $id | wc -c) -lt 10 ]] && echo ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${id:0:6}/$id || echo ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${id:0:6}/$(printf '%03i' $(echo ${id:9} | sed 's/^0*//'))/$id)		
		if [[ $n -gt 4 ]]; then		
			#PE
			rmdup_cmd=$(${do_rmdup:-false} && echo "samtools rmdup - -" || echo "cat")
			awk -v f=<($get_cmd $url/${id}_2.fastq.gz | pigz -p 1 -kdc | sed '/^@/{s/\/[12]$//}' | paste - - - -) '{print $0; getline < f; print $0}' <($get_cmd $url/${id}_1.fastq.gz | pigz -p 1 -kdc | sed '/^@/{s/\/[12]$//}' | paste - - - -) | tr '\t' '\n' | cutadapt -j $threads -q 20 -m 18 --trim-n -a CTGTCTCTTATACACATCT -a AGATCGGAAGAGC -A CTGTCTCTTATACACATCT -A AGATCGGAAGAGC --interleaved -o - - | mbuffer -l /dev/null -q -Q -m 10M | stdbuf -o0 -e0 bwa-mem2 mem -p -T $minscore -R '@RG\tID:A1\tSM:sample1\tLB:library1\tPU:unit1\tPL:illumina' -a -Y -t $threads $idx - | samtools view -@ $threads -q 1 -F 256 -F 2048 -u | samtools sort -u -@ $threads -T /data/tmp/samtools_tmp.$id | $rmdup_cmd | samtools view -@ $threads -b --write-index -o "$outdir/$id.bam##idx##$outdir/$id.bai"
		else		
			#SE 
			rmdup_cmd=$(${do_rmdup:-false} && echo "samtools rmdup -s - -" || echo "cat")
			$get_cmd $url/${id}.fastq.gz | pigz -p 1 -kdc | cutadapt -j $threads -q 20 -m 18 --trim-n -a CTGTCTCTTATACACATCT -a AGATCGGAAGAGC -o - - | mbuffer -l /dev/null -q -Q -m 10M | stdbuf -o0 -e0 bwa-mem2 mem -T $minscore -R '@RG\tID:A1\tSM:sample1\tLB:library1\tPU:unit1\tPL:illumina' -a -Y -t $threads $idx - | samtools view -@ $threads -q 1 -F 256 -F 2048 -u | samtools sort -u -@ $threads -T /data/tmp/samtools_tmp.$id | $rmdup_cmd | samtools view -@ $threads -b > $outdir/$id.bam && samtools index -@ $threads $outdir/$id.bam $outdir/$id.bai			
		fi
	CMD
done < ./srr.list

# -----------------------
# Submit jobs via bashbone
# -----------------------

commander::printcmd -a cmds
commander::qsubcmd -c bwa -v -n sra2Bam -o sge -r -a cmds -i 10 -p threads -t $threads