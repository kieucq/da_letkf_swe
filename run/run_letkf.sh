#!/bin/sh
#####################################################################
#                                                            
# NOTE: This is to run an assimilation system that performs  
#       the Local Ensemble Transform Kalman Filter for the  
#       primitive equation 1-layer model.
#
# HOW TO RUN:
#       - Step 1: change main_dir to your local system
#       - Step 2: make sure the entire system is compilied with PGI
#       - Step 3: set main options for the DA system as well as LETKF 
#                 namelist. Note that the LETKF options are controled
#                 by the namelist.swe, while the DA cycling is set
#                 by the options at the begining of this script.
#       - Step 4: run sh ./run_letkf.sh
#       - Step 5: check diagnostic analysis in the dig dir
#
# HIST: - Sep 01, 08: created by CK at UMD from the SWE model
#                     of Krishnamurthy and Bounoua (Chapter 8). 
#       - Jan 30, 10: added cycling DA for demonstration of
#                     LETKF with vortex assimiation at VNU
#       - Sep 09, 10: simplied the workflow for DA teaching 
#                     at VNU                     
#       - Feb 18, 11: developed options for 4D-LETKF to compare
#                     the system with 4DVAR. 
#       - Mar 05, 11: freezed this version and moved the entire 
#                     4DLETKF development to a new code 
#                     structure (see da_4dletkf_swe). 
#       - Jan 05, 22: revised and added more options for OSSE for DA
#                     a DA course at IU
#                                                        
# AUTH: Chanh Kieu, Email: ckieu@iu.edu. Tel: 812-856-5704.                     
#                                                   
#####################################################################
#
# set up namelist options
#
main_dir="/N/u/ckieu/Karst/model/da_letkf_swe"
ninterval=3                     # DA window/restart interval [h]  
nensemble=30                    # number of ensemble members 
restart="Y"                     # fresh run opt from begining with truth ref
nt=48                           # forecast time for model [h]
ncycle=$(($nt/$ninterval+1))    # number of DA cycles 
cold_start="T"                  # cold start option
history="L"                     # long or short history
t_window=0                      # 4D DA window [s] - must be equal 0 for LETKF
dt_window=1                     # obs frequency in each 4D DA cycle [s] - must be 1 for LETKF 
mpi_run="N"                     # MPI option for LETKF
osse=1                          # option for obs osse design: 0-all grid point (no=nx*ny)
                                #                             1-vortex obs around (icen,jcen)
                                #                             2-a single obs at (icen,jcen)
echo "SUMMARY OF LETKF CONFIGURATION"
echo "==========================================================="
echo " Forecast time is                                : $nt"
echo " Number of cycles that the script will run is    : $ncycle"
echo " Number of the ensemble members for this test is : $nensemble"
echo " Interval for doing the analysis is              : $ninterval"
echo " Cold start option is                            : $cold_start"
echo " Re-start option is                              : $restart"
echo " MPI-RUN option is                               : $mpi_run"
echo " OSSE observational option is                    : $osse"
echo " History option is (L-long; S-short)             : $history"       
echo "==========================================================="
#if [ $history == "L" ]; then
#   set -x
#fi
#
# create a shared namelist
#
mv -f namelist.swe namelist.swe.bk
cat > namelist.swe << EOF
debug        = 0               ! debuging level
restart      = ${ninterval}               ! reatart inteval for assimilation [h] 
dy           = 208e+3          ! grid distance in the y-direction [m] (dy =1.875 degree from tape23)
zamp         = 110             ! amplitude of z-anomaly           [m]
zscale       = 500e+3          ! scale of z-anomaly               [m]
icen         = 72              ! i-center of the anomaly          [] 
jcen         = 7               ! j-center of the anomaly          []
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
tlm_flag     = 1               ! option for the TLM model: 1-first order, 2-second order
no           = 100             ! number of local observation
ne           = ${nensemble}              ! number of ensemble members for LETKF
nxl          = 5               ! size of the local patch
slat         = 10.             ! start latitude
nx           = 114             ! grid point in x-direction
ny           = 17              ! grid point in y-direction
tfcst        = ${nt}             ! length of forecast               [h]
dt           = 300.            ! model timestep
obs_flag     = ${osse}               ! osse type: 0-full grid, 1-vortex, 2-single obs 
oscale       = 700e+03         ! radius of obs influence
da_flag      = 4               ! 0 (no da); 1(u); 2(v); 3(z); 4(all)               
t_window     = ${t_window}             ! 4d assimilation window          [s]
dt_window    = ${dt_window}             ! 4d assimilation obs frequency [s]
nt           = $(($t_window/$dt_window+1))               ! number of time slices for 4d-letkf
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
rm -rf ${main_dir}/dig/*.dat
if [ "$restart" = "Y" ]; then
#
# start running the truth first to serve as a base. Here the 
# truth assumes a vortex generated at the center of the domain.
#
   echo "Generating a truth base ..."
   cd ${main_dir}/truth
   rm -rf tru_*:*.dat truth.dat bvortex_truth.dat bvortex.dat
   cp ${main_dir}/run/namelist.swe ./
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
   cp ${main_dir}/run/namelist.swe ./
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
 cp ${main_dir}/run/namelist.swe ./
 rm -rf bgd_*:* ini.dat
 ./ini_letkf.exe > ${main_dir}/run/log.ini
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
cp ${main_dir}/run/namelist.swe ./
./ctl.exe > ${main_dir}/run/log.ctl
#
# start running the model with assimilation now. Begin with
# looping over the entire cycles.
#
icycle=1
rtime=0
echo "Perfoming the LETKF cycle run now ..."
while [ "$icycle" -le "$ncycle" ]
do
 echo "====== CYCLE: ${icycle}th @ $ofile ========"
#
# checking the background data first
#
 echo "checking background data ..."
 rm -rf ${main_dir}/letkf/bgd_*.dat ${main_dir}/utils/mem_*.dat
 ie=1
 while [ "$ie" -le "$nensemble" ]
 do
  if [ "$ie" -lt 10 ]; then
   imember="_00$ie"
  elif [ "$ie" -lt 100 ]; then
   imember="_0$ie"
  elif [ "$ie" -lt 1000 ]; then
   imember="_$ie"
  else
   echo " TO MANY MEMBER...exit 10"
   exit 10
  fi
  bgdfile="bgd$imember"_"$ofile.dat"
  if [ -f "${main_dir}/bgd/$ofile/$bgdfile" ]; then
   #  echo " Background file $bgdfile exists...linking"
   ln -sf ${main_dir}/bgd/$ofile/$bgdfile ${main_dir}/letkf/bgd$imember.dat
   ln -sf ${main_dir}/bgd/$ofile/$bgdfile ${main_dir}/utils/mem$imember.dat
  else
   echo " -> background file $bgdfile does not exist...exit 2"
   exit 2
  fi
  ie=$(($ie + 1))
 done
#
# checking the observation data now
#
 echo "linking observational file ..."
 obsfile="obs_$ofile.dat"
 if [ -f "${main_dir}/obs/$obsfile" ]; then
  echo "obs $obsfile file exists ...linking"
  rm -rf ${main_dir}/letkf/obs.dat
  ln -sf ${main_dir}/obs/$obsfile ${main_dir}/letkf/obs.dat
 else
  echo " -> obs file $obsfile does not exist. exit 3"
  exit 3
 fi
#
# perform the LETKF assimilation now and compute the ensemble
# mean for the background now.
#
 echo "LETKF is being called ..."
 cd ${main_dir}/letkf/
 cp ${main_dir}/run/namelist.swe ./
 rm -rf ana_*.dat
 if [ "$mpi_run" == "Y" ]; then
  time mpirun -np 2 ./letkf_mpi.exe >& ${main_dir}/run/log.letkf
 else
  time ./letkf_serial.exe >& ${main_dir}/run/log.letkf
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
 cp ${main_dir}/run/namelist.swe ${main_dir}/model/
 while  [ "$ie" -le "$nensemble" ]
 do
  echo "running swe.exe for member ${ie} ..."
#
# creat an extension for the member first
#
  if [ "$ie" -lt 10 ]; then
   imember="_00$ie"
  elif [ "$ie" -lt 100 ]; then
   imember="_0$ie"
  else
   imember="_$ie"
  fi
  cd ${main_dir}/model/
  rm -rf ./ana.dat
  ln -sf ${main_dir}/letkf/ana$imember.dat ./ana.dat
  cd ${main_dir}/model/
  ./swe.exe > ${main_dir}/run/log.swe
#
# add the member extention to the output from the model runs. Also
# backup the analysis output at each each
#
  fscfile="fsc$imember"_"$ofile.dat"
  bgdfile="bgd$imember.dat"
  anafile="ana$imember"_"$ofile.dat"
  #echo " swe.dat will be renamed as $fscfile ..."
  #echo " bgd.dat will be renamed as $bgdfile ..."
  mkdir -p ${main_dir}/fsc/${ofile}/mem${imember}
  mv -f swe.dat ${main_dir}/fsc/${ofile}/mem${imember}/
  mv -f fsc_*:*.dat ${main_dir}/fsc/${ofile}/mem${imember}/
  mv -f bgd_t001.dat ./$bgdfile
  mv -f ${main_dir}/letkf/ana$imember.dat ${main_dir}/ana/${ofile}/$anafile 
  ln -sf ${main_dir}/ana/${ofile}/$anafile ${main_dir}/utils/mem$imember.dat
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
 rtime=$(($rtime+$ninterval))
#
# update the background at the next cycle now
#
 if [ "$rtime" -lt 24 ]; then
  ih=$rtime
 else
  ih=$(($rtime-24))
  id=$(($id+1))
  rtime=$(($rtime-24))
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
 ofile="$dd:$hh:$mm"
 if [ -d ${main_dir}/bgd/${ofile} ]; then
  rm -rf ${main_dir}/bgd/${ofile}
 fi
 mkdir ${main_dir}/bgd/${ofile}
 echo "new processing time is $ofile"
 ie=1
 while  [ "$ie" -le "$nensemble" ]
 do
   if [ "$ie" -lt 10 ]; then
    imember="_00$ie"
   elif [ "$ie" -lt 100 ]; then
    imember="_0$ie"
   else
    imember="_$ie"
   fi
   bgdfile="bgd$imember"_"$ofile.dat"
   #  echo " background for the next analysis is $bgdfile"
   mv ${main_dir}/model/bgd$imember.dat ${main_dir}/bgd/$ofile/$bgdfile
   ie=$(($ie+1))
 done
 echo ""
done
#
# doing some analysis now
#
echo "doing analysis now ..."
cd ${main_dir}/dig
cp ${main_dir}/run/namelist.swe ./
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
