MODULE module_interface
!=======================================================================
!
! NOTE:    MPI interface for 4d-letkf communication
!
! HISTORY: 08 Sep 2011: created by CK
!         
!
!=======================================================================
!
! model domain
!
  integer                          :: nx                  ! model domain size in x direction
  integer                          :: ny                  ! model domain size in y direction
  integer                          :: ng                  ! model model global vector length     
  integer                          :: ne                  ! number of ensemble members
  integer                          :: nt                  ! assimilation window time intervals
  integer                          :: nxl                 ! local patch dimension (NEED TO BE ODD)
  integer                          :: no                  ! number of observations
  integer                          :: nol                 ! local patch of obs
!
! 2D state arrays
!
  real, allocatable                :: ub(:,:,:,:)         ! background u
  real, allocatable                :: vb(:,:,:,:)         ! background v
  real, allocatable                :: zb(:,:,:,:)         ! background z
  real, allocatable                :: ua(:,:,:)           ! analysis u
  real, allocatable                :: va(:,:,:)           ! analysis v
  real, allocatable                :: za(:,:,:)           ! analysis z
  real, allocatable                :: uo(:,:)             ! obs u
  real, allocatable                :: vo(:,:)             ! obs v
  real, allocatable                :: zo(:,:)             ! obs z
!
! vector state arrays
!      
  real, allocatable                :: ub_g(:,:,:)         ! global u background (forecast)
  real, allocatable                :: vb_g(:,:,:)         ! global v background (forecast)
  real, allocatable                :: zb_g(:,:,:)         ! global z background (forecast)
  real, allocatable                :: ub_gm(:,:)          ! global ensemble mean of background of u
  real, allocatable                :: vb_gm(:,:)          ! global ensemble mean of background of v
  real, allocatable                :: zb_gm(:,:)          ! global ensemble mean of background of z
  real, allocatable                :: ua_g(:,:)           ! global analysis for u
  real, allocatable                :: va_g(:,:)           ! global analysis for v
  real, allocatable                :: za_g(:,:)           ! global analysis for z
  real, allocatable                :: ua_gm(:)            ! global ensemble mean of analysis for u
  real, allocatable                :: va_gm(:)            ! global ensemble mean of analysis for v
  real, allocatable                :: za_gm(:)            ! global ensemble mean of analysis for z
  real, allocatable                :: uo_g(:,:,:)         ! global background obs ensemble for u
  real, allocatable                :: vo_g(:,:,:)         ! global background obs ensemble for v
  real, allocatable                :: zo_g(:,:,:)         ! global background obs ensemble for z
  real, allocatable                :: uo_gm(:,:)          ! global ensemble mean of ensemble obs for u
  real, allocatable                :: vo_gm(:,:)          ! global ensemble mean of ensemble obs for v
  real, allocatable                :: zo_gm(:,:)          ! global ensemble mean of ensemble obs for z
  real, allocatable                :: Xa(:,:)             ! local analysis
  real, allocatable                :: Xa_m(:)             ! local ensemble mean of analysis
  real, allocatable                :: Xb(:,:)             ! local background (forecast) 
  real, allocatable                :: Xb_m(:)             ! local ensemble mean of background
  real, allocatable                :: Yb(:,:,:)           ! local background obs ensemble
  real, allocatable                :: Yb_m(:,:)           ! local ensemble mean of ensemble obs
  real, allocatable                :: Yo(:,:)             ! local observation data
  real, allocatable                :: Ro(:,:)             ! local obs error cov mtx
  real, allocatable                :: rho(:)              ! localization matrix
  real, allocatable                :: olon(:,:),olat(:,:) ! observation location
  real, allocatable,dimension(:,:) :: tem2d1,tem2d2       ! tem var
!
! system options
!
  integer                          :: model_flag          ! flag for model: 0-perfect,1-imperfect
  integer                          :: ini_flag            ! flag for initial condition: 0-perfect, 1-imperfect
  integer                          :: np                  ! total local patch dimension = nxl*nxl*num_var
  integer                          :: nme                 ! number of model ensemble
  real                             :: ifactor             ! inflation factor
  character*100                    :: ifile,ofile,temc    ! I/O files
  integer                          :: debug               ! debuging
  integer                          :: irec                ! output record
  real                             :: obs_err_u           ! observation variance for u
  real                             :: obs_err_v           ! observation variance for v
  real                             :: obs_err_z           ! observation variance for z
  real                             :: rscale              ! scale of background variance matrix
  integer                          :: da_flag             ! flag for choosing assimilation variables
  real                             :: mis                 ! missing value for ploting
CONTAINS
END MODULE module_interface
