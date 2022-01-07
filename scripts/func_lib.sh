function cal_time()
{
 rtime1=$1
 id=0
 ih=0
 im=0
 ih=$(($rtime1/60))
 im=`expr $rtime1 % 60`
 if [ "$ih" -ge "24" ]; then
  id=$(($ih/24))
  ih=`expr $ih % 24`
 fi
 if [ "$im" -lt "10" ]; then
  mm="0$im"
 elif [ "$im" -lt "60" ]; then
  mm="$im"
 else
  echo " minute index calculation is incorrect ...exit 1"
  exit 1
 fi
 if [ "$ih" -lt 10 ]; then
  hh="0$ih"
 elif [ "$ih" -lt 24 ]; then
  hh="$ih"
 else
  echo " hour index calculation is incorrect ...exit 1"
  exit 1
 fi
 if [ "$id" -lt 10 ]; then
  dd="0$id"
 elif [ "$id" -lt 100 ]; then
  dd="$id"
 else
  echo " day index calculation is out of range 100 days ...exit 2"
  exit 2
 fi
 get_file="$dd:$hh:$mm"
}

function cal_index()
{
 input=$1
 if [ "$input" -lt 10 ]; then
  get_index="00$input"
 elif [ "$input" -lt 100 ]; then
  get_index="0$input"
 elif [ "$input" -lt 1000 ]; then
  get_index="$input"
 else
  echo " TO MANY MEMBER...exit 10"
  exit 10
 fi
}
