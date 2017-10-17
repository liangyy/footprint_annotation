srr=$1
dirname=$2
log=$3
numLines=$(fastq-dump -X 1 -Z --split-spot $srr | wc -l)
if [ $numLines -eq 4 ]
then
  echo "$srr is single-end"
  fastq-dump -O $dirname --gzip $srr > $log 2>&1
else
  echo "$srr is paired-end"
  fastq-dump -O $dirname --gzip $srr --split-3 > $log 2>$1
fi
