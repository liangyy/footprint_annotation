trim=$1
input1=$2
input2=$3
output1=$4
output2=$5
trimp=$6
echo $1 $2 
echo $3 $4 $5 $6
if [ "$trim" == 'No' ]
then
  echo "No trim need to be done"
  cp $input1 $output1
  cp $input2 $output2
else
  echo "Do trim" $trim
  $trim $input1 $input2 $output1 $output1.unpaired $output2 $output2.unpaired $trimp
fi
