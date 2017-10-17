trim=$1
input1=$2
input2=$3
output1=$4
output2=$5
trimp=$6

if [ $trim -eq 'No' ]
then
  echo "No trim need to be done"
  cp $input1 $output1
  cp $input2 $output2
else
  echo "Do trim" $trim
  $trim $input1 $output1 $input2 $output2 $trimp
fi
