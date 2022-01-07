!
! This is the common interface for all codes that use read namelist. 
! If namelist.swe is changed, this include file has to be modified
! properly.
!
      IMPLICIT NONE              
      REAL                    :: dy                   ! grid distance in the y-direction [m]
      INTEGER                 :: debug                ! debuging
      REAL                    :: restart              ! restart inverval
      INTEGER                 :: icen,jcen            ! x and y-location of the anomaly
      CHARACTER*100           :: temc                 ! tem var
      REAL                    :: zamp                 ! amplitude of z-anomaly           [m]
      REAL                    :: zscale               ! scale of z-anomaly               [m]
      REAL                    :: obs_err_u            ! obs error for u                  [m/s]
      REAL                    :: obs_err_v            ! obs error for u                  [m/s]
      REAL                    :: obs_err_z            ! obs error for u                  [m/s]
      REAL                    :: bgd_err_u            ! background error for u           [m/s]
      REAL                    :: bgd_err_v            ! background error for u           [m/s]
      REAL                    :: bgd_err_z            ! background error for u           [m/s]
      REAL                    :: ifactor              ! inflation factor
      INTEGER                 :: model_flag           ! flag for model: 0-perfect, 1-imperfect
      INTEGER                 :: ini_flag             ! flag for initial condition: 0-perfect, 1-imperfect
      INTEGER                 :: nme                  ! number of ensemble for model error calculation
      INTEGER                 :: timeout              ! output interval for the model (steps)
      REAL                    :: rscale               ! localization scale around a point
      REAL                    :: slat                 ! start latitude
      INTEGER                 :: tlm_flag             ! option for the TLM model
      INTEGER                 :: no                   ! number of global obs
      INTEGER                 :: ne                   ! number of ensemble members for LETKF
      INTEGER                 :: nxl                  ! size of local patch
      INTEGER                 :: nx,ny                ! model domain in x and y directions
      REAL                    :: dt                   ! model time step
      REAL                    :: tfcst                ! length of forecast                [h]
      INTEGER                 :: obs_flag             ! option for obs output: 0:full, 1:vortex
      REAL                    :: oscale               ! radius of obs influence
      INTEGER                 :: da_flag              ! 0(no da); 1(u); 2(v); 3(z); 4(all)
      INTEGER                 :: i,j
      REAL                    :: t_window             ! 4d assimilation window
      REAL                    :: dt_window            ! 4d assimilation window interval
      INTEGER                 :: nt                   ! number of time slices for 4d-letkf
      OPEN(14,file='namelist.swe',status='old')
      READ(14,*)temc,temc,debug
      READ(14,*)temc,temc,restart
      READ(14,*)temc,temc,dy
      READ(14,*)temc,temc,zamp
      READ(14,*)temc,temc,zscale
      READ(14,*)temc,temc,icen
      READ(14,*)temc,temc,jcen
      READ(14,*)temc,temc,obs_err_u
      READ(14,*)temc,temc,obs_err_v
      READ(14,*)temc,temc,obs_err_z
      READ(14,*)temc,temc,bgd_err_u
      READ(14,*)temc,temc,bgd_err_v
      READ(14,*)temc,temc,bgd_err_z
      READ(14,*)temc,temc,model_flag
      READ(14,*)temc,temc,ini_flag
      READ(14,*)temc,temc,rscale
      READ(14,*)temc,temc,ifactor
      READ(14,*)temc,temc,nme
      READ(14,*)temc,temc,timeout
      READ(14,*)temc,temc,tlm_flag
      READ(14,*)temc,temc,no
      READ(14,*)temc,temc,ne
      READ(14,*)temc,temc,nxl
      READ(14,*)temc,temc,slat
      READ(14,*)temc,temc,nx
      READ(14,*)temc,temc,ny
      READ(14,*)temc,temc,tfcst
      READ(14,*)temc,temc,dt
      READ(14,*)temc,temc,obs_flag
      READ(14,*)temc,temc,oscale
      READ(14,*)temc,temc,da_flag
      READ(14,*)temc,temc,t_window
      READ(14,*)temc,temc,dt_window
      READ(14,*)temc,temc,nt
      REWIND(14)
