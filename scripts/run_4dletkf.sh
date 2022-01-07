#!/bin/sh
echo ""
echo "################################################################"
echo "##                                                            ##"
echo "## NOTE:                                                      ##"
echo "##      This is to run an assimilation system that performs   ##"
echo "##      the 4D Local Ensemble Transform Kalman Filter for     ##"
echo "##      the primitive equation 1-layer model                  ##"
echo "##                                                            ##"
echo "##      This system works under the following assumptions:    ##"
echo "##      1. Observational errors are random with the same err  ##"
echo "##         cov at all time slices;                            ##"
echo "##      2. Obs are given at all model grid points such that   ##"
echo "##         H = I for all time slices. Need to change the core ##"
echo "##         if obs distribution is change (module_4dletkf.f90; ##"
echo "##      3. The adaptive inflation/localization is not used;   ##"
echo "##      4. Model is perfect.                                  ##"
echo "##      5. The total integration < 100 days. This is because  ##"
echo "##         all day index is limited at 2 digits, which is     ##"
echo "##         merely because that is enough for testing. Need to ##"
echo "##         modify code and script as 0000:00:00 if longer     ##"
echo "##         integration is needed.                             ##"
echo "##                                                            ##"
echo "## CONTACT:                                                   ##"
echo "##      Author: Chanh Q. Kieu, Lecturer                       ##"
echo "##      Lab for Climate and Weather Research                  ##"
echo "##      Vietnam National University, Hanoi                    ##"
echo "##      email: chanhkq@vnu.edu.vn                             ##"
echo "##                                                            ##"
echo "################################################################"
echo ""
ninterval=6
nensemble=30
restart="Y"
nt=72
ncycle=$(($nt/$ninterval+1))
cold_start="T"
history="L"
t_window=3600
dt_window=3600
mpi_run="Y"
nslices=$(($t_window/$dt_window+1))
#nslices1=$nslices
nslices1=1
main_dir="/state/partition1/home/kieuc/model/da_4dletkf_swe"
echo "SUMMARY OF 4D-LETKF CONFIGURATION"
echo "==========================================================="
echo " Forecast time is                                : $nt"
echo " Number of cycles that the script will run is    : $ncycle"
echo " Number of the ensemble members for this test is : $nensemble"
echo " Interval for doing the analysis is              : $ninterval"
echo " Cold start option is                            : $cold_start"
echo " Re-start option is                              : $restart"
echo " MPI-RUN option is                               : $mpi_run"
echo " Time window for 4D-LETKF is                     : $t_window"
echo " Time window interval for 4D-LETKF is            : $dt_window"
echo " History option is (L-long; S-short)             : $history"       
echo "==========================================================="
if [ "$t_window" -ge "$(($ninterval*3600))" ]; then
 echo "Assimilation window has to be smaller than assimilation cycle...exit"
 echo "$nensemble < $t_window"
 exit 1
fi
#
# create a shared namelist
#
mv -f namelist.swe namelist.swe.bk
cat > namelist.swe << EOF
debug        = 0               ! debuging level
restart      = ${ninterval}               ! reatart inteval for assimilation [h] 
dy           = 208e+3          ! grid distance in the y-direction [m](dy=1.875 degree)
zamp         = 110             ! amplitude of z-anomaly           [m]
zscale       = 500e+3          ! scale of z-anomaly               [m]
icen         = 63              ! i-center of the anomaly          [] 
jcen         = 6               ! j-center of the anomaly          []
obs_err_u    = 1.0             ! obs error for u                  [m/s]
obs_err_v    = 1.0             ! obs error for v                  [m/s]
obs_err_z    = 3.0             ! obs error magnitude              [m] 
bgd_err_u    = 5.0             ! background error magnitude       [m]
bgd_err_v    = 5.0             ! background error for v           [m/s]
bgd_err_z    = 18.             ! background error magnitude       [m]
model_flag   = 0               ! flag for model: 0-perfect, 1-imperfect
ini_flag     = 1               ! flag for initial condition: 0-perfect, 1-imperfect
rscale       = 3               ! scale of localization
ifactor      = 0.1             ! inflation factor to increase model error 
nme          = 0               ! number of ensemble member for model error calculation
timeout      = 100             ! output interval for the model (steps)
tlm_flag     = 1               ! option for the TLM model: 1-first order,2-second order
no           = 1938            ! number of local observation
ne           = ${nensemble}              ! number of ensemble members for LETKF
nxl          = 5               ! size of the local patch
slat         = 10.             ! start latitude
nx           = 114             ! grid point in x-direction
ny           = 17              ! grid point in y-direction
tfcst        = ${nt}            ! length of forecast               [h]
dt           = 300.            ! model timestep                    [s]
obs_flag     = 0               ! option for obs output: 0:full, 1:vortex 
oscale       = 700e+03         ! radius of obs influence
da_flag      = 4               ! 0 (no da); 1(u); 2(v); 3(z); 4(all)               
t_window     = ${t_window}           ! 4d assimilation window          [s]
dt_window    = ${dt_window}            ! 4d assimilation window interval [s]
nt           = ${nslices1}               ! number of time slices for 4d-letkf
para5        = 0               ! extra slot for later new para
para6        = 0               ! extra slot for later new para
para7        = 0               ! extra slot for later new para
para8        = 0               ! extra slot for later new para
para9        = 0               ! extra slot for later new para
para10       = 0               ! extra slot for later new para
para11       = 0               ! extra slot for later new para
para12       = 0               ! extra slot for later new para
para13       = 0               ! extra slot for later new para
para14       = 0               ! extra slot for later new para
para15       = 0               ! extra slot for later new para
para16       = 0               ! extra slot for later new para
para17       = 0               ! extra slot for later new para
para18       = 0               ! extra slot for later new para
para19       = 0               ! extra slot for later new para
EOF
#
# set up first some base run
#
echo "Cleaning up previous outputs..."
rm -rf log*
rm -rf ${main_dir}/ana/*
rm -rf ${main_dir}/fsc/*
rm -rf ${main_dir}/bgd/*
rm -rf ${main_dir}/dig/*:*:*.dat
if [ "$restart" = "Y" ]; then
#
# start running the truth first to serve as a base
#
 echo "Generating a truth base ..."
 cd ${main_dir}/truth
 rm -rf tru_*:*.dat truth.dat bvortex_truth.dat bvortex.dat
 ./bvortex.exe > ${main_dir}/run/log.truth
 ./truth.exe >>  ${main_dir}/run/log.truth
#
# now create observational data from the given truth obtained
# above
#
 echo "Generating observation data ..."
 cd ${main_dir}/obs/
 rm -rf tru_*:*.dat obs_*:*.dat obs.dat
 ln -sf ${main_dir}/truth/tru_*:*.dat ./
 ./obs.exe > ${main_dir}/run/log.obs
fi 
#
# copy some initial startup for the cold start mode of KF
# cycle
#
id=0
ih=0
mm="00"
if [ "$ih" -lt 10 ]; then
 hh="0$ih"
elif [ "$ih" -lt 24 ]; then
 hh="$ih"
else 
 echo "-> hour index calculation is incorrect ...exit 1"
 exit 1
fi
if [ "$id" -lt 10 ]; then
 dd="0$id"
elif [ "$id" -lt 100 ]; then
 dd="$id"
else
 echo "-> day index calculation is out of range 100 days ...exit 2"
 exit 2
fi
ofile="$dd:$hh:$mm"
if [ "$cold_start" = "T" ]; then
 echo "Creating an ensemble of initial cold start ..."
 cd ${main_dir}/ini/
 rm -rf bgd_*:* ini_*:*.dat
 ./ini.exe > ${main_dir}/run/log.ini
 ./bgd.exe >> ${main_dir}/run/log.ini
 cd ${main_dir}/bgd
 if [ -d ${ofile} ]; then
  rm -rf ${ofile} 
 fi
 mkdir ${ofile}
 mv ${main_dir}/ini/bgd_*.dat ./${ofile}/
fi
#
# running a control with no assimilation for comparsion
#
cd ${main_dir}/ctl/
ln -sf ${main_dir}/bgd/${ofile}/bgd_001_${ofile}.dat ./bgd.dat
./ctl.exe > ${main_dir}/run/log.ctl
#
# start running the model with assimilation now. Begin with
# looping over the entire cycles.
#
icycle=1
rtime=0
echo "Perfoming the LETKF cycle run now ..."
. ${main_dir}/run/func_lib.sh
while [ "$icycle" -le "$ncycle" ]
do
 echo "====== CYCLE: ${icycle}th @ $ofile ========"
#
# checking the background data first
#
 echo "checking background data ..."
 rm -rf ${main_dir}/letkf/bgd_*.dat
 rm -rf ${main_dir}/letkf/ana_*.dat
 rm -rf ${main_dir}/letkf/obs_*.dat 
 rm -rf ${main_dir}/utils/mem_*.dat
 ie=1
 while [ "$ie" -le "$nensemble" ]
 do
  cal_index $ie
  imember="_${get_index}"
  it=1
  timeslice=0
  while [ "$it" -le "$nslices" ]; do
   cal_time $timeslice
   cal_index $it
   islice=$get_index
   bgdfile="bgd$imember"_"${get_file}.dat"
   bgdletkf="bgd${imember}_t${islice}.dat"
   if [ -f "${main_dir}/bgd/$ofile/$bgdfile" ]; then
    # echo " Background file $bgdfile exists...linking $bgdletkf"
    ln -sf ${main_dir}/bgd/$ofile/$bgdfile ${main_dir}/letkf/$bgdletkf
    if [ "$it" -eq 1 ]; then
     ln -sf ${main_dir}/bgd/$ofile/$bgdfile ${main_dir}/utils/mem$imember.dat
    fi
   else
    echo " -> background file $bgdfile does not exist...exit 2"
    exit 2
   fi
   it=$(($it + 1))
   timeslice=$(($timeslice+$dt_window/60))
  done
  ie=$(($ie + 1))
 done
#
# checking the observation data now
#
 echo "linking observational file ..."
 it=1
 timeslice=0
 while [ "$it" -le "$nslices" ]; do
  echo "time slice = $timeslice"
  cal_time $timeslice
  cal_index $it
  islice=$get_index 
  obsfile="obs_${get_file}.dat"
  if [ -f "${main_dir}/obs/$obsfile" ]; then
   echo "obs $obsfile file exists ...linking"
   rm -rf ${main_dir}/letkf/obs_${islice}.dat
   ln -sf ${main_dir}/obs/$obsfile ${main_dir}/letkf/obs_${islice}.dat
  else
   echo " -> obs file $obsfile does not exist. exit 3"
   exit 3
  fi
  it=$(($it + 1))
  timeslice=$(($timeslice+$dt_window/60))
 done
#
# perform the LETKF assimilation now and compute the ensemble
# mean for the background now.
#
 echo "LETKF is being called ..."
 cd ${main_dir}/letkf/
 rm -rf ana_*.dat
 if [ "$mpi_run" == "Y" ]; then
  time mpirun -np 2 ./4dletkf_mpi.exe >& ${main_dir}/run/log.letkf
 else
  time ./4dletkf_mpi.exe >& ${main_dir}/run/log.letkf
 fi
 cd ${main_dir}/utils
 ./mean.exe > ${main_dir}/run/log.utils
 mv mean.dat ${main_dir}/dig/bgd_$ofile.dat
 mv spread.dat ${main_dir}/dig/spb_$ofile.dat
#
# running the model now
#
 ie=1
 if [ -d ${main_dir}/fsc/${ofile} ]; then
  rm -rf ${main_dir}/fsc/${ofile}
 fi
 mkdir ${main_dir}/fsc/${ofile}
 if [ -d ${main_dir}/ana/${ofile} ]; then
  rm -rf ${main_dir}/ana/${ofile}
 fi
 mkdir ${main_dir}/ana/${ofile}
 cd ${main_dir}/model/
 rm -rf ./ana.dat ./fsc_*:*:*.dat ./bgd_*.dat
 while  [ "$ie" -le "$nensemble" ]
 do
  echo "running swe.exe for member ${ie} ..."
  cal_index $ie
  imember="_${get_index}"
  ln -sf ${main_dir}/letkf/ana$imember.dat ./ana.dat
  ./swe.exe > ${main_dir}/run/log.swe
#
# add the member extention to the output from the model runs. Also
# backup the analysis output at each each
#
  fscfile="fsc$imember"_"$ofile.dat"
  anafile="ana$imember"_"$ofile.dat"
  mkdir -p ${main_dir}/fsc/${ofile}/mem${imember}
  mv -f swe.dat ${main_dir}/fsc/${ofile}/mem${imember}/
  mv -f fsc_*:*.dat ${main_dir}/fsc/${ofile}/mem${imember}/
  mv -f ${main_dir}/letkf/ana$imember.dat ${main_dir}/ana/${ofile}/$anafile 
  ln -sf ${main_dir}/ana/${ofile}/$anafile ${main_dir}/utils/mem$imember.dat
  it=1
  timeslice=0
  while [ "$it" -le "$nslices" ]; do
   cal_time $timeslice
   cal_index $it
   islice=$get_index
   bgdfile="bgd${imember}_${get_file}.dat"
   mv -f bgd_t${islice}.dat ./$bgdfile
   it=$(($it + 1))
   timeslice=$(($timeslice+$dt_window/60))
  done
  ie=$(($ie + 1))
 done
 rm -rf ana.dat
#
# compute the ensemble mean for the analysis now
#
 cd ${main_dir}/utils/
 ./mean.exe >> ${main_dir}/run/log.utils
 mv mean.dat ${main_dir}/dig/ana_$ofile.dat
 mv spread.dat ${main_dir}/dig/spa_$ofile.dat
 rm mem_*.dat
#
# update the next cycle now
#
 icycle=$(($icycle+1))
 rtime=$(($rtime+$ninterval*60))
#
# update the background at the next cycle now
#
 cal_time $rtime
 ofile=$get_file
 if [ -d ${main_dir}/bgd/${ofile} ]; then
  rm -rf ${main_dir}/bgd/${ofile}
 fi
 mkdir ${main_dir}/bgd/${ofile}
 echo "new processing time is $ofile"
 mv ${main_dir}/model/bgd_*:*:*.dat ${main_dir}/bgd/$ofile/
 echo ""
done
#
# doing some analysis now
#
echo "doing analysis now ..."
cd ${main_dir}/dig
ln -sf ${main_dir}/truth/tru_*:*.dat ./
ln -sf ${main_dir}/obs/obs_*:*.dat ./
ln -sf ${main_dir}/ctl/ctl_*:*.dat ./
./ana.exe > ${main_dir}/run/log.dig
cd ${main_dir}/run/
#
# cleaning history now
#
if [ "$history" = "S" ]; then
 ./clean.sh
fi
echo "LETKF DONE"
exit 0
