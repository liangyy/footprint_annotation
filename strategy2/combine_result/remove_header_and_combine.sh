args=("$@")
o=${args[0]}
cmd='cat'
for ((i=1; i < $#; i++))
{
    tail -n +2 <(zcat ${args[$i]}) > ${args[$i]}.temp
    cmd="$cmd ${args[$i]}.temp"
}
$cmd | gzip > $o
for ((i=1; i < $#; i++))
{
    rm ${args[$i]}.temp
}


