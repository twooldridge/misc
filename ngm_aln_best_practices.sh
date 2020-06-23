#!/bin/bash 

usage(){ echo "                                             
Usage: quick.sh [-m|--memory] [-t|--time] [-n|--threads] [-i|--jobid] "command.here"
        
        -m: job memory
        -t: job time
        -n: threads
        -i: sample id
        -b: input bam file
        -r: reference genome                                                           
"                                                           
1>&2;exit 1;}

die() {
    printf '%s\n' "$1" >&2
    exit 1
}

# Initialize all the option variables.
# This ensures we are not contaminated by variables from the environment.
###Defaults###
memory=80000
time=72:00:00
threads=16
##############


while :; do
        case $1 in
                -h|-\?|--help)
                        usage    # Display a usage synopsis.
                        exit1
                        ;;
                -m|--memory)
                        if [ "$2" ];then
                                memory=$2
                                shift
                        fi
                        ;;
                -t|--time)
                        if [ "$2" ];then
                                time=$2
                                shift
                        fi
                        ;;
                -n|--threads)
                        if [ "$2" ];then
                                threads=$2
                                shift
                        fi
                        ;;
                -i|--sampleid)
                        if [ "$2" ];then
                                outid=$2
                                shift
                        else
                                echo "output ID needed";exit1
                        fi
                        ;;
                -b|--bam)
			if [ "$2" ];then
                                ubamfile=$2
                                shift
                        else
                                echo "ubamfile needed";exit1
                        fi
                        ;;
                -r|--ref)
                        if [ "$2" ];then
                                refgenome=$2
                                shift
                        else
                                echo "reference genome needed";exit1
                        fi
                        ;;
                --)
                        shift
                        break
                        ;;
                -?*)
                        printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
                        ;;
                *)
                        break
                        ;;
        esac
        shift
done


cat <<EOF > ${outid}.aln.slurm
#!/bin/bash
#SBATCH -p hoekstra,commons,shared,general
#SBATCH -t ${time}
#SBATCH --mem=${memory}
#SBATCH -n ${threads}
#SBATCH -e ./logs/${outid}.aln.e
#SBATCH -o ./logs/${outid}.aln.o
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -J ${outid}.aln

source activate py36
module load centos6/0.0.1-fasrc01; module load java
module load samtools

if [ -f ${outid}.sorted.ubam.md5 ]; then echo "${outid}.sorted.ubam succesfully created previously,skipping" ; else java -Xmx8G -XX:ParallelGCThreads=1 -Djava.io.tmpdir=`pwd`/tmp -jar ~/Software/picard/2.18.4/picard.jar SortSam CREATE_MD5_FILE=true SORT_ORDER=queryname I=${ubamfile} O=${outid}.sorted.ubam;fi

if [ -f ${outid}.markedadapters.ubam.md5 ]; then echo "${outid}.markedadapters.ubam succesfully created previously,skipping" ; else java -Xmx8G -XX:ParallelGCThreads=1 -Djava.io.tmpdir=`pwd`/tmp -jar ~/Software/picard/2.18.4/picard.jar MarkIlluminaAdapters CREATE_MD5_FILE=true I=${outid}.sorted.ubam O=${outid}.markedadapters.ubam M=MZ11012.metrics.txt;fi

if [ -f ${outid}.samtofastq_interleaved.fq.md5 ]; then echo "${outid}.samtofastq_interleaved.fq succesfully created previously,skipping" ; else java -XX:ParallelGCThreads=1 -Djava.io.tmpdir=`pwd`/tmp -Dsamjdk.buffer_size=131072 -Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4G -jar ~/Software/picard/2.18.4/picard.jar SamToFastq CREATE_INDEX=true CREATE_MD5_FILE=true I=${outid}.markedadapters.ubam FASTQ=${outid}.samtofastq_interleaved.fq CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true;fi

if [ -f ${outid}.ngm_mem.sam ]; then echo "${outid}.ngm_mem.sam created previously, skipping mapping"; else ngm --sensitive -t ${threads} -r ${refgenome} -p -q ${outid}.samtofastq_interleaved.fq -o ${outid}.ngm_mem.sam;fi

samtools quickcheck ${outid}.ngm_mem.sam && ( echo "nothing wrong with ${outid}.ngm_mem.sam"; ) || ( echo "${outid}.ngm_mem.sam corrupted, exiting"; rm ${outid}.ngm_mem.sam; exit 1; )

if [ -f ${outid}.ngm_mem.sorted.bam.md5 ]; then echo "${outid}.sorted.ubam succesfully created previously,skipping" ; else java -Xmx8G -XX:ParallelGCThreads=1 -Djava.io.tmpdir=`pwd`/tmp -jar ~/Software/picard/2.18.4/picard.jar SortSam CREATE_MD5_FILE=true SORT_ORDER=queryname I=${outid}.ngm_mem.sam O=${outid}.ngm_mem.sorted.bam;fi

if [ -f ${outid}.mdup2500.bam.md5 ]; then echo "${outid}.mdup2500.bam.md5 succesfully created previously,skipping"; else java -Dsamjdk.buffer_size=131072 -Dsamjdk.use_async_io=true -Dsamjdk.compression_level=1 -XX:+UseStringCache -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:ParallelGCThreads=1 -Djava.io.tmpdir=`pwd`/tmp -Xmx32G -jar ~/Software/picard/2.18.4/picard.jar MarkDuplicates INPUT=${outid}.ngm_mem.sorted.bam OUTPUT=${outid}.mdup2500.bam METRICS_FILE=MZ11012_markduplicates_metrics.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true CREATE_MD5_FILE=true;fi

EOF

mkdir -p logs
#sbatch ${outid}.aln.slurm
