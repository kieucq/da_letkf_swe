!
! NOTE:     
!      This program performs a local ensemble transform Kalman filter
!      for the shallow water equation. Note particular that this program
!      is designed only for a model with 3 varirables only. 
!
! HISTORY:  
!    - 02 Feb 2010: Updated from the LETKF for L40 model.
!    - 20 Feb 2010: Code runs well with local domain = 1 and ne < 9.
!                   Somehow the allocation does not allow for bigger
!                   ensemble size. Will add localization now so that
!                   local domain > 1.
!    - 23 Feb 2010: Add an option that allows for assimilating only u
!                   or v or z via localization factor.
!
! REFERENCE:
!    - Hunt et al., 2006: arxiv:physics/0511236
!
! AUTHOR: 
!      Chanh Q. Kieu, Research Associate
!      Dept. of Atmospheric and Oceanic science
!      Univ. of Maryland, College Park, MD
!      email: kieucq@atmos.umd.edu
!
! COPYRIGHT: (C) 2009
!
!=================================================================
!
  program LETKF_3variables
  use common
  use common_mtx
  implicit none
  integer                          :: nx                  ! model domain size in x direction
  integer                          :: ny                  ! model domain size in y direction
  integer                          :: ng                  ! model model global vector length     
  integer                          :: ne                  ! number of ensemble members
  integer                          :: nxl                 ! local patch dimension of model state (need to be odd)
  integer                          :: no                  ! number of observations
  integer                          :: nol                 ! local patch of obs
!
! 2D state arrays
!
  real, allocatable                :: ub(:,:,:)           ! background u
  real, allocatable                :: vb(:,:,:)           ! background v
  real, allocatable                :: zb(:,:,:)           ! background z
  real, allocatable                :: ua(:,:,:)           ! analysis u
  real, allocatable                :: va(:,:,:)           ! analysis v
  real, allocatable                :: za(:,:,:)           ! analysis z
  real, allocatable                :: uo(:),vo(:),zo(:)   ! obs u,v,z
!
! vector state arrays
!      
  real, allocatable                :: ub_g(:,:)           ! global u background (forecast)
  real, allocatable                :: vb_g(:,:)           ! global v background (forecast)
  real, allocatable                :: zb_g(:,:)           ! global z background (forecast)
  real, allocatable                :: ub_gm(:)            ! global ensemble mean of background of u
  real, allocatable                :: vb_gm(:)            ! global ensemble mean of background of v
  real, allocatable                :: zb_gm(:)            ! global ensemble mean of background of z
  real, allocatable                :: ua_g(:,:)           ! global analysis for u
  real, allocatable                :: va_g(:,:)           ! global analysis for v
  real, allocatable                :: za_g(:,:)           ! global analysis for z
  real, allocatable                :: ua_gm(:)            ! global ensemble mean of analysis for u
  real, allocatable                :: va_gm(:)            ! global ensemble mean of analysis for v
  real, allocatable                :: za_gm(:)            ! global ensemble mean of analysis for z
  real, allocatable                :: uo_g(:,:)           ! global background obs ensemble for u
  real, allocatable                :: vo_g(:,:)           ! global background obs ensemble for v
  real, allocatable                :: zo_g(:,:)           ! global background obs ensemble for z
  real, allocatable                :: uo_gm(:)            ! global ensemble mean of ensemble obs for u
  real, allocatable                :: vo_gm(:)            ! global ensemble mean of ensemble obs for v
  real, allocatable                :: zo_gm(:)            ! global ensemble mean of ensemble obs for z
  real, allocatable                :: Xa(:,:)             ! local analysis
  real, allocatable                :: Xa_m(:)             ! local ensemble mean of analysis
  real, allocatable                :: Xb(:,:)             ! local background (forecast) 
  real, allocatable                :: Xb_m(:)             ! local ensemble mean of background
  real, allocatable                :: Yb(:,:)             ! local background obs ensemble
  real, allocatable                :: Yb_m(:)             ! local ensemble mean of ensemble obs
  real, allocatable                :: Yo(:)               ! local observation data
  real, allocatable                :: Ro(:,:)             ! local obs error cov mtx
  real, allocatable                :: rho(:)              ! localization matrix
  real, allocatable                :: olon(:),olat(:)     ! observation location
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
  integer                          :: i,j,m,n,k,ie        ! indexing
  integer                          :: debug               ! debuging
  integer                          :: irec                ! output record
  real                             :: obs_err_u           ! observation variance for u
  real                             :: obs_err_v           ! observation variance for v
  real                             :: obs_err_z           ! observation variance for z
  real                             :: rscale              ! scale of background variance matrix
  integer                          :: da_flag             ! flag for choosing assimilation variables
  real                             :: mis                 ! missing value for ploting
  mis           = -99999
  irec          = 1
!
! reading namelist now
!
  call input_namelist(debug,obs_err_u,obs_err_v,obs_err_z,model_flag, &
                      ini_flag,ifactor,rscale,nme,no,ne,nx,ny,nxl,da_flag)
  print*,'letkf.exe: input namelist is'
  print*,'letkf.exe: debug         =',debug
  print*,'letkf.exe: obs_err_u     =',obs_err_u
  print*,'letkf.exe: obs_err_v     =',obs_err_v
  print*,'letkf.exe: obs_err_z     =',obs_err_z
  print*,'letkf.exe: model_flag    =',model_flag
  print*,'letkf.exe: ini_flag      =',ini_flag
  print*,'letkf.exe: ifactor       =',ifactor
  print*,'letkf.exe: rscale        =',rscale  
  print*,'letkf.exe: da_flag       =',da_flag
  print*,'letkf.exe: nme           =',nme
  print*,'letkf.exe: no            =',no
  print*,'letkf.exe: ne            =',ne
  print*,'letkf.exe: nxl           =',nxl
  print*,'letkf.exe: nx            =',nx
  print*,'letkf.exe: ny            =',ny
  if (debug.ge.1) read*
  if (ini_flag.eq.0) then
   print*,'WARNING: INPUT IS PERFECT? WHY NEED ASSIMILATION...STOP'
   stop
  endif 
  if (nme.gt.1.and.ifactor.gt.1e-10) then
   print*,'WARNING: CANNOT HAVE BOTH INFLATION AND PERTURBED FORCING...STOP'
   stop
  endif
  open(92,file='letkf.dat',access='direct',form='unformatted',recl=nx*ny*4)
  ng            = nx*ny  
  np            = nxl*nxl*3
!
! allocate arrays now
!
  allocate(ub(nx,ny,ne),vb(nx,ny,ne),zb(nx,ny,ne))
  allocate(ua(nx,ny,ne),va(nx,ny,ne),za(nx,ny,ne))
  allocate(uo(no),vo(no),zo(no),olat(no),olon(no))
  allocate(ub_g(ng,ne),vb_g(ng,ne),zb_g(ng,ne))
  allocate(ub_gm(ng),vb_gm(ng),zb_gm(ng))
  allocate(ua_g(ng,ne),va_g(ng,ne),za_g(ng,ne)) 
  allocate(ua_gm(ng),va_gm(ng),za_gm(ng))   
  allocate(uo_g(no,ne),vo_g(no,ne),zo_g(no,ne))
  allocate(uo_gm(no),vo_gm(no),zo_gm(no))
  allocate(Xb(np,ne),Xb_m(np),Xa(np,ne),Xa_m(np))
!
! reading observation data
!
  open(90,file='obs.dat')
  i              = 1
2 continue
  read(90,*,end=3)olon(i),olat(i),uo(i),vo(i),zo(i)
  if (debug.ge.2) then
   if (i.eq.1) print*,'letkf.exe: Checking the obs'
   print*,i,olon(i),olat(i),uo(i),vo(i),zo(i)
  endif
  i             = i + 1
  goto 2
3 print*,'letkf.exe: reading obs returns',i,' data points'
  close(90)
  if (i.ne.no) then
     print*,'letkf.exe: actual # of data points differ from no. Reset no'
     no = i-1
     print*,'letkf.exe: # of data points now is no = ',no
  endif
!
! quick check of obs data
!
  allocate(tem2d1(nx,ny))
  tem2d1        = mis
  do i          = 1,no
   tem2d1(nint(olon(i)),nint(olat(i))) = uo(i)
  enddo
  write(92,rec=irec)((tem2d1(i,j),i=1,nx),j=1,ny)
  irec          = irec + 1
  tem2d1        = mis
  do i          = 1,no
   tem2d1(nint(olon(i)),nint(olat(i))) = vo(i)
  enddo
  write(92,rec=irec)((tem2d1(i,j),i=1,nx),j=1,ny)
  irec          = irec + 1
  tem2d1        = mis
  do i          = 1,no
   tem2d1(nint(olon(i)),nint(olat(i))) = zo(i)
  enddo
  write(92,rec=irec)((tem2d1(i,j),i=1,nx),j=1,ny)
  irec          = irec + 1
  deallocate(tem2d1)
!
! reading forecast files
!
  ifile='bgd_000.dat'
  ie            = 0
8 continue
  ie            = ie + 1
  if (ie.le.9) then
   write(ifile(7:7),'(1I1)')ie
  elseif (ie.le.99) then
   write(ifile(6:7),'(1I2)')ie 
  elseif (ie.le.999) then
   write(ifile(5:7),'(1I3)')ie
  else
   print*,'letkf.exe: Stop as too many ensemble members'
   stop
  endif
  if (debug.ge.1) print*,'letkf.exe: input forecast file: ',ifile(1:20)
  open(91,file=ifile)
  read(91,'(6E13.6)')((ub(i,j,ie),i=1,nx),j=1,ny)
  read(91,'(6E13.6)')((vb(i,j,ie),i=1,nx),j=1,ny)
  read(91,'(6E13.6)')((zb(i,j,ie),i=1,nx),j=1,ny)
  close(91)
  if (ie.lt.ne) goto 8
  do ie          = 1,ne
   write(92,rec=irec)((ub(i,j,ie),i=1,nx),j=1,ny)
   irec          = irec + 1
  enddo
  do ie          = 1,ne 
   write(92,rec=irec)((vb(i,j,ie),i=1,nx),j=1,ny)
   irec          = irec + 1
  enddo
  do ie          = 1,ne
   write(92,rec=irec)((zb(i,j,ie),i=1,nx),j=1,ny)
   irec          = irec + 1
  enddo
!
! convert from 2d arrays to 1d column vectors
!
  call convert_array_vector(ub,vb,zb,nx,ny,ne,ub_g,vb_g,zb_g,ng)
!
! create the background first for the analysis
!
  ua_g           = ub_g
  va_g           = vb_g
  za_g           = zb_g
! 
! compute the global pertrubation background observation ensemble 
! from the global background. (step 1 in H04)
! 
  call global_observation(uo_g,vo_g,zo_g,uo_gm,vo_gm,zo_gm,  &
                          ub,vb,zb,olat,olon,ne,nx,ny,ng,no)
!
! compute the global mean and ensebmle of pertubation background.
! Note that the global background will not be needed from now on,
! so the Xb_g array will now store perturbation, not the total.
! This is step 2 in H04.
!
  ub_gm(:)       = 0
  vb_gm(:)       = 0
  zb_gm(:)       = 0
  do i           = 1,ne
   ub_gm(:)      = ub_gm(:) + ub_g(:,i)
   vb_gm(:)      = vb_gm(:) + vb_g(:,i)
   zb_gm(:)      = zb_gm(:) + zb_g(:,i)       
  enddo
  ub_gm(:)       = ub_gm(:)/ne
  vb_gm(:)       = vb_gm(:)/ne
  zb_gm(:)       = zb_gm(:)/ne
  do i           = 1,ne
   ub_g(:,i)     = ub_g(:,i) - ub_gm(:)
   vb_g(:,i)     = vb_g(:,i) - vb_gm(:)
   zb_g(:,i)     = zb_g(:,i) - zb_gm(:)
  enddo
!
! Now loop over all grid point with the corresponding
! local patch.
!
  grid_loop: do i = 1,ng
   print*,'global grid loop @ i =',i
!
! compute first the grid location and var.
!
   n             = i/nx + 1
   m             = mod(i,nx)
   if (m.eq.0) then 
    n            = n - 1
    m            = nx
   endif             
   if (debug.ge.1) print*,'Working at the local patch around',i,m,n
!
! finding the number of local observation within each patch. Note that
! the size of the local patch in 2D array will be (m-nxl/2,m+nxl/2) x
! (n-nxl/2,n+nxl/2). This step is merely for defining the local size 
!
   call local_obs(m,n,nx,ny,nxl,olon,olat,no,nol,debug)
   if (debug.ge.1) print*,'letkf.exe: Number of local obs is',nol
   if (nol.lt.1) goto 10
!
! now we will project from global to local patch at each grid point (m,n)
! The cross correlation can be handled separately to create a new set of
! e.g., un from z. (u,un) are then treated as the total obs of u around
! (m,n).  
!
   allocate(Yb(nol,ne),Yb_m(nol),Yo(nol),Ro(nol,nol),rho(nol))
   call global2local(ub_g,vb_g,zb_g,ub_gm,vb_gm,zb_gm,uo_g,vo_g,zo_g,           &
                     uo_gm,vo_gm,zo_gm,uo,vo,zo,Xb,Xb_m,Yb,Yb_m,Yo,np,nol,      &
    	             olon,olat,m,n,nx,ny,ng,ne,no,nxl,rscale,rho,debug)
!
! define the local observational error covariance matrix R
!
   call observational_err_cov_mtx(Ro,nol,obs_err_u,obs_err_v,obs_err_z)
!
! calling the LETKF core now.
!
    call letkf_main(Xa,Ro,Xb,Xb_m,Yb,Yb_m,Yo,np,nol,ne,ifactor,rho,da_flag,debug)
    if (debug.ge.1) then
     print*,'Local Yb(:1) is'
     write(*,'(10E11.3)')(Yb(j,1),j=1,nol)
     print*,'Local Yo is'
     write(*,'(10E11.3)')(Yo(j),j=1,nol)
     print*,'Local Xb(:1) is'
     write(*,'(10E11.3)')(Xb(j,1),j=1,np)
     print*,'Local Xa(:1) is'
     write(*,'(10E11.3)')(Xa(j,1),j=1,np)
     read*
    endif
!
! project from local patch to global patch at the point $i$, which is corresponding
! to point (m,n) on the 2D grid
!
    call local2global(Xa,np,ne,ua_g,va_g,za_g,nx,ny,ng,m,n,nxl,i)
!
! deallocate the arrays
!
    deallocate(Yb,Yb_m,Yo,Ro,rho)
10  continue

  enddo grid_loop
!
! convert from 1D arrays to 2D arrays
!
  call vector2array(ua_g,va_g,za_g,ua,va,za,ng,ne,nx,ny) 
!
! output analysis x and updated Pa. Note that bmatrix.dat
! will be over-written by the updated Pa as this file
! will be discarded after performing KF.
!
  ofile='ana_000.dat'
  ie            = 0
9 continue
  ie            = ie + 1
  if (ie.le.9) then
   write(ofile(7:7),'(1I1)')ie
  elseif (ie.le.99) then
   write(ofile(6:7),'(1I2)')ie
  elseif (ie.le.999) then
   write(ofile(5:7),'(1I3)')ie
  else
   print*,'letkf.exe: Stop as too many ensemble members'
   stop
  endif
  if (debug.ge.1) print*,'letkf.exe:  Open output file is: ',ofile(1:20)
  open(12,file=ofile,status='unknown')
  write(12,'(6E13.6)')((ua(i,j,ie),i=1,nx),j=1,ny)
  write(12,'(6E13.6)')((va(i,j,ie),i=1,nx),j=1,ny)
  write(12,'(6E13.6)')((za(i,j,ie),i=1,nx),j=1,ny)
  close(12)
  if (ie.lt.ne) goto 9
!
! printout for viewing
!
  do ie          = 1,ne
   write(92,rec=irec)((ua(i,j,ie),i=1,nx),j=1,ny)
   irec          = irec + 1
  enddo
  do ie          = 1,ne 
   write(92,rec=irec)((va(i,j,ie),i=1,nx),j=1,ny)
   irec          = irec + 1
  enddo
  do ie          = 1,ne
   write(92,rec=irec)((za(i,j,ie),i=1,nx),j=1,ny)
   irec          = irec + 1
  enddo
  print*,'letkf.exe: LETKF finishes safely...!'
  end


  subroutine local_obs(m,n,nx,ny,nxl,olon,olat,no,nol,debug)
  implicit none
  integer m,n,no,nxl,nol,nx,ny
  real olon(no),olat(no)
  real xlon
  integer i,j,k,istart,iend,jstart,jend,debug
!
! defining the start and end position of the local patch around ig
!
  istart         = m - nxl/2
  iend           = m + nxl/2
  jstart         = n - nxl/2
  jend           = n + nxl/2
  if (debug.ge.1) print*,'debug: istart,iend,jstart,jend',istart,iend,jstart,jend
  if (istart.lt.1) then
   istart        = 1
   iend   	 = nxl
   nol           = 0 
   return
  endif
  if (iend.gt.nx)  then
   istart        = nx - nxl + 1
   iend          = nx
   nol           = 0
   return
  endif
  if (jstart.lt.1) then
   jstart        = 1
   jend          = nxl
   nol           = 0
   return
  endif
  if (jend.gt.ny)  then
   jstart        = ny - nxl + 1
   jend          = ny
   nol           = 0
   return
  endif
  if (debug.ge.1) print*,'debug fixed: istart,iend,jstart,jend',istart,iend,jstart,jend
!
! find the local obs within the local obs patch.
!
  nol            = 0
  do i           = 1,no   
   if (olon(i).ge.istart.and.olon(i).le.iend.and. &
       olat(i).ge.jstart.and.olat(i).le.jend) then
    nol          = nol + 1
   endif
  enddo  
  nol            = nol*3
  return
  end

  subroutine global2local(ub_g,vb_g,zb_g,ub_gm,vb_gm,zb_gm,uo_g,vo_g,zo_g,          &
                          uo_gm,vo_gm,zo_gm,uo,vo,zo,Xb,Xb_m,Yb,Yb_m,Yo,np,nol,     &
		          olon,olat,m,n,nx,ny,ng,ne,no,nxl,rscale,rho,debug)
  implicit none
  integer             :: nx,ny,no,np,nol,ng,ne,m,n,nxl
  real, intent(in)    :: ub_g(ng,ne),vb_g(ng,ne),zb_g(ng,ne)
  real, intent(in)    :: ub_gm(ng),vb_gm(ng),zb_gm(ng)
  real, intent(in)    :: uo_g(no,ne),vo_g(no,ne),zo_g(no,ne)
  real, intent(in)    :: uo_gm(no),vo_gm(no),zo_gm(no)
  real, intent(in)    :: uo(no),vo(no),zo(no)
  real, intent(in)    :: olat(no),olon(no),rscale
  real, intent(out)   :: Xb(np,ne),Xb_m(np),Yb(nol,ne),Yb_m(nol),Yo(nol),rho(nol)
  integer             :: i,j,k,istart,iend,jstart,jend,debug,nlocal,id
!
! defining the start and end position of the local patch around ig
!
  istart         = m - nxl/2
  iend           = m + nxl/2
  jstart         = n - nxl/2
  jend           = n + nxl/2
  if (debug.ge.1) print*,'debug: istart,iend,jstart,jend',istart,iend,jstart,jend
  if (istart.lt.1) then
   istart        = 1
   iend  	 = nxl
   print*,'global2local: bnd points are not allowed...stop'
   stop
  endif
  if (iend.gt.nx)  then
   istart        = nx - nxl + 1
   iend          = nx
   print*,'global2local: bnd points are not allowed...stop'
   stop
  endif
  if (jstart.lt.1) then
   jstart        = 1
   jend          = nxl
   print*,'global2local: bnd points are not allowed...stop'
   stop
  endif
  if (jend.gt.ny)  then
   jstart        = ny - nxl + 1
   jend          = ny
   print*,'global2local: bnd points are not allowed...stop'
   stop
  endif
  if (debug.ge.1) print*,'debug: fixed istart,iend,jstart,jend',istart,iend,jstart,jend	
!
! asign the global array to local patch
!
  k              = 0
  do i           = istart,iend
   do j          = jstart,jend
    k            = k + 1
    id           = (j-1)*nx+i
    Xb(k,1:ne)         = ub_g(id,1:ne) 
    Xb(np/3+k,1:ne)    = vb_g(id,1:ne)
    Xb(2*np/3+k,1:ne)  = zb_g(id,1:ne)  
    Xb_m(k)            = ub_gm(id)
    Xb_m(np/3+k)       = vb_gm(id)
    Xb_m(2*np/3+k)     = zb_gm(id)
   enddo
  enddo
!
! find and assign the global obs ensemble to the local obs patch.
! Note that the cyclinc boundary is applied only for +/- 5 points
! at each end of the domain.
!
  k              = 0
  do i           = 1,no   
   if (olon(i).ge.istart.and.olon(i).le.iend.and. &
       olat(i).ge.jstart.and.olat(i).le.jend) then
    k            = k + 1
    Yb(k,1:ne)         = uo_g(i,1:ne)
    Yb(nol/3+k,1:ne)   = vo_g(i,1:ne)
    Yb(2*nol/3+k,1:ne) = zo_g(i,1:ne)
    Yb_m(k)            = uo_gm(i)
    Yb_m(nol/3+k)      = vo_gm(i)
    Yb_m(2*nol/3+k)    = zo_gm(i)
    Yo(k)              = uo(i)
    Yo(nol/3+k)        = vo(i)
    Yo(2*nol/3+k)      = zo(i)
    rho(k)             = exp(-((olon(i)-m)**2. + (olat(i)-n)**2.)/rscale**2)
    rho(nol/3+k)       = rho(k)
    rho(2*nol/3+k)     = rho(k)
   endif
  enddo
  if (debug.ge.1) print*,'k,nol/3',k,nol/3
  if (k.ne.nol/3) then
   print*,'letkf.exe: The number of obs in local patch does not match the array shape'
   stop
  endif
  return
  end
  
  subroutine local2global(Xa,np,ne,ua_g,va_g,za_g,nx,ny,ng,m,n,nxl,ig)
  implicit none
  integer nx,ny,ne,np,ng,m,n,ig,nxl
  real Xa(np,ne),ua_g(ng,ne),va_g(ng,ne),za_g(ng,ne)
  integer istart,iend,jstart,jend,ic,debug
  debug          = 0
!
! defining the start and end position of the local patch around ig
!
  istart         = m - nxl/2
  iend           = m + nxl/2
  jstart         = n - nxl/2
  jend           = n + nxl/2
  if (debug.ge.1) print*,'debug: istart,iend,jstart,jend',istart,iend,jstart,jend
  if (istart.lt.1) then
   istart        = 1
   iend   		 = nxl
  endif
  if (iend.gt.nx)  then
   istart        = nx - nxl + 1
   iend          = nx
  endif
  if (jstart.lt.1) then
   jstart        = 1
   jend          = nxl
  endif
  if (jend.gt.ny)  then
   jstart        = ny - nxl + 1
   jend          = ny
  endif
  if (debug.ge.1) print*,'debug: fixed istart,iend,jstart,jend',istart,iend,jstart,jend
!
! note that np = nxl*nxl*3
!
  ic             = nxl*nxl/2 + 1
  ua_g(ig,1:ne)  = Xa(ic,1:ne)
  va_g(ig,1:ne)  = Xa(np/3+ic,1:ne)
  za_g(ig,1:ne)  = Xa(2*np/3+ic,1:ne)   
  return
  end

  subroutine letkf_main(Xa,R,Xb,Xb_m,Yb,Yb_m,Yo,nxl,nol,ne,ifactor,rho,da_flag,debug)
  use common
  use common_mtx
  implicit none
  integer nxl,nol,ne,debug
  real, intent(in)  :: R(nol,nol),Xb(nxl,ne),Xb_m(nxl),rho(nol)
  real, intent(in)  :: Yb(nol,ne),Yb_m(nol),Yo(nol),ifactor
  real, intent(out) :: Xa(nxl,ne)  
  integer           :: da_flag
  real              :: C(ne,nol),tem1(ne,ne),dy(nol)
  real              :: Wa(ne,ne),Pat(ne,ne),wa_m(ne),Rinv(nol,nol)
  integer           :: i1,j1,k1,i,j,k
!
! compute the R inverted first, and couple with the localization
! factor
!
  call compute_Rinv(R,Rinv,nol,da_flag)
  do i1         = 1,nol
   Rinv(i1,i1)  = Rinv(i1,i1)*rho(i1)
  enddo
  if (debug.ge.2) then 
   print*,'R matrix'
   do i1        = 1,nol
    write(0,'(20F6.3)')(R(i1,j1),j1=1,nol)
   enddo
  endif 
!
! compute matrix C (step 4 in H06)
!
  C             = matmul(transpose(Yb),Rinv)
  if (debug.ge.2) then 
   print*,'step 4: C matrix'
   do i1        = 1,ne
     write(0,'(20F6.3)')(C(i1,j1),j1=1,nol)
   enddo
  endif
!
! compute matrix Pa tilde (step 5 in H06)
!
  tem1          = matmul(C,Yb)
  do j          = 1,ne
   tem1(j,j)    = tem1(j,j) + (ne-1)/(1+ifactor)
  enddo
  call mtx_inv(ne,tem1,Pat)
  if (debug.ge.2) then
   print*,'step 5: Yb matrix'
   do i1        = 1,nol
    write(0,'(20F6.3)')(Yb(i1,j1),j1=1,20)
   enddo
   print*,'step 5: tem1=C*Yb matrix'
   do i1        = 1,20
    write(0,'(20F6.3)')(tem1(i1,j1),j1=1,20)
   enddo
  endif
!
! compute matrix Wa (step 6 in H06) 
!
  tem1          = (ne-1)*Pat
  call mtx_sqrt(ne,tem1,Wa)
  if (debug.ge.2) then
   print*,'step 6: Wa matrix'
   do i1        = 1,20
    write(0,'(20F6.3)')(Wa(i1,j1),j1=1,20)
   enddo
  endif
!
! compute wa_m and add it to Wa matrix (step 7 in H06)
!
  dy            = Yo - Yb_m
  if (debug.ge.2) print*,'step 7: ok dy'
  wa_m          = matmul(matmul(Pat,C),dy)
  if (debug.ge.2) print*,'step 7: ok wa_m'
  do j          = 1,ne
   Wa(:,j)      = Wa(:,j) + wa_m(:)
  enddo
  if (debug.ge.2) print*,'step 7: ok Wa'
!
! Computing the analysis ensemble members (step 8 in H06)
!
  Xa            = matmul(Xb,Wa)
  if (debug.ge.2) print*,'step 8: ok Xa'
  do j          = 1,ne
   Xa(:,j)      = Xa(:,j) + Xb_m(:)
  enddo
  return
  end

  subroutine observation_operator(nx,no,olon,H)
  implicit none
  integer nx,no,i,j,m,n,k,io,jo
  real H(no,nx),olon(no)
  H         = 0.
  do i      = 1,no
   j        = int(olon(i)+0.001)  
   H(i,j)   = 1.
  enddo
  return
  end

  subroutine localiztion_operator(nx,ny,no,nv,olon,olat,lopt)
  implicit none
  integer nx,ny,no,nv,i,j,m,n
  real lopt(nv,no),olat(no),olon(no),rscale,radi
  rscale    = 10.
  do i      = 1,nv
   do j     = 1,no
    m       = mod(i,nx)
    n       = i/nx + 1
    if (m.eq.0) then
     m      = nx
     n      = n - 1
    endif
    radi    = sqrt((olon(j)-m)**2. + (olat(j)-n)**2.)
    lopt(i,j) = exp(-(radi/rscale)**2)
   enddo
  enddo
  print*,'Checking Pat matrix'
  do i       = 1,20
   write(*,'(20F6.2)')(lopt(i,j),j=1,20)
  enddo
  return
  end

  subroutine background_err_cov_mtx(B,nx,bvar,rscale)
  implicit none
  integer i,j,m,n,nx
  real B(nx,nx),rscale,radi,bvar
  do i      = 1,nx
   do j     = 1,nx
    radi    = sqrt((j-i)**2. + (j-i)**2.)
    B(i,j)  = bvar*bvar*exp(-(radi/rscale)**2)
   enddo
  enddo
  return
  end


  subroutine compute_Htilde(H,no,nv,Xf,ne,Ht)
  implicit none
  integer ne,nv,no
  real H(no,nv),Xf(nv,ne),Ht(no,ne)
  integer i,j,k,m,n
  Ht       = matmul(H,Xf)
  return
  end

  subroutine observational_err_cov_mtx(R,no,uvar,vvar,zvar)
  implicit none
  integer no,i
  real R(no,no),uvar,vvar,zvar
  R        = 0
  do i     = 1,no
   if (i.le.no/3) then								  
    R(i,i) = uvar**2
   elseif (i.le.2*no/3.and.i.gt.no/3) then
    R(i,i) = vvar**2 
   else 
    R(i,i) = zvar**2 
   endif
  enddo
  return
  end


  subroutine compute_Rinv(R,Rinv,no,da_flag)
  implicit none
  integer no,i,da_flag
  real R(no,no),Rinv(no,no)
  Rinv     = 0.
  do i     = 1,no
   if (da_flag.eq.0) then
    Rinv(i,i) = 0.
   elseif (da_flag.eq.1) then
    if (i.le.no/3) Rinv(i,i) = 1/R(i,i)
   elseif (da_flag.eq.2) then
    if (i.gt.no/3.and.i.le.2*no/3) Rinv(i,i) = 1/R(i,i)
   elseif (da_flag.eq.3) then
    if (i.gt.2*no/3) Rinv(i,i) = 1/R(i,i)
   else
    Rinv(i,i) = 1/R(i,i)
   endif
  enddo
  return
  end
  
  subroutine obs_increment(H,no,nv,xfm,po,ne,obs_inc)
  implicit none
  integer no,ne,nv
  real H(no,nv),xfm(nv),po(no),obs_inc(nv)
  real tem
  integer i,j
  do i     = 1,no
   tem     = 0.
   do j    = 1,nv
    tem    = tem + H(i,j)*xfm(j)
   enddo
   obs_inc(i) = po(i) - tem
  enddo
  return
  end

  subroutine analysis_mean(K,lopt,nv,no,xfm,obs_inc,xam)
  implicit none
  integer no,nv
  real K(nv,no),xfm(nv),obs_inc(no),xam(nv),lopt(nv,no)
  integer i,j
  do i     = 1,nv
   xam(i)  = xfm(i)
   do j    = 1,no
    xam(i) = xam(i) + lopt(i,j)*K(i,j)*obs_inc(j)
   enddo
  enddo
  return
  end

  subroutine convert_vector_array(Xa,nx,ny,nv,pa)
  implicit none
  integer nv,nx,ny
  real Xa(nv),pa(nx,ny)
  integer i,j,k,m,n
  do i   = 1,nx
   do j  = 1,ny
    m    = (j-1)*nx + i
    pa(i,j)  = Xa(m)
   enddo
  enddo
  return
  end

  subroutine convert_vector_array1(Xa,nv,pa,nx,ny)
  implicit none
  integer nv,nx,ny
  real Xa(nv),pa(nx,ny)
  integer i,j,k,m,n
  do i       = 1,nx
   do j      = 1,ny
    m        = (j-1)*nx + i
    pa(i,j)  = Xa(m)
   enddo
  enddo
  return
  end

  subroutine convert_array_vector(ub,vb,zb,nx,ny,ne,ub_g,vb_g,zb_g,ng)
  implicit none
  integer nx,ny,ne,ng
  real ub(nx,ny,ne),vb(nx,ny,ne),zb(nx,ny,ne)
  real ub_g(ng,ne),vb_g(ng,ne),zb_g(ng,ne)
  integer i,j,k,m,n,nv
  do i           = 1,ng
   n             = i/nx + 1
   m             = mod(i,nx)
   if (m.eq.0) then
    m            = nx
    n            = n - 1
   endif
   ub_g(i,:)     = ub(m,n,:)
   vb_g(i,:)     = vb(m,n,:)
   zb_g(i,:)     = zb(m,n,:)
  enddo
  return
  end

  SUBROUTINE global_observation(uo_g,vo_g,zo_g,uo_gm,vo_gm,zo_gm,  &
                                ub,vb,zb,olat,olon,ne,nx,ny,ng,no)
  IMPLICIT NONE
  INTEGER nx,ny,ng,no,ne
  REAL ub(nx,ny,ne),vb(nx,ny,ne),zb(nx,ny,ne)
  REAL uo_g(no,ne),vo_g(no,ne),zo_g(no,ne)
  REAL uo_gm(no),vo_gm(no),zo_gm(no)
  REAL olat(no),olon(no)
  INTEGER i,j,k,m,n
!
! this is the first step in Hunt et al, for which H is simple = I because
! observations are given at the grid point. 
!
  DO k           = 1,ne
   DO i          = 1,no
    m            = nint(olon(i))
    n            = nint(olat(i))
    uo_g(i,k)    = ub(m,n,k)
    vo_g(i,k)    = vb(m,n,k)
    zo_g(i,k)    = zb(m,n,k)
   ENDDO
  ENDDO
!
! compute the $\ba Y$ according to H05 
!
  uo_gm(:)       = 0.
  vo_gm(:)       = 0.
  zo_gm(:)       = 0.
  do k           = 1,ne
   uo_gm(:)      = uo_gm(:) + uo_g(:,k)
   vo_gm(:)      = vo_gm(:) + vo_g(:,k)
   zo_gm(:)      = zo_gm(:) + zo_g(:,k)
  enddo 
  uo_gm(:)       = uo_gm(:)/ne
  vo_gm(:)       = vo_gm(:)/ne
  zo_gm(:)       = zo_gm(:)/ne
!
! substract the mean of the total global obs from the total global obs
! to create pertrubation global obs
!
  do k           = 1,ne
   uo_g(:,k)     = uo_g(:,k) - uo_gm(:)
   vo_g(:,k)     = vo_g(:,k) - vo_gm(:)
   zo_g(:,k)     = zo_g(:,k) - zo_gm(:)
  enddo
  RETURN
  END

  SUBROUTINE vector2array(ua_g,va_g,za_g,ua,va,za,ng,ne,nx,ny) 
  IMPLICIT NONE
  INTEGER nx,ny,ne,ng
  REAL ua_g(ng,ne),va_g(ng,ne),za_g(ng,ne)
  REAL ua(nx,ny,ne),va(nx,ny,ne),za(nx,ny,ne)
  INTEGER i,j,k,m,n
  DO i           = 1,nx
   DO j          = 1,ny
    m            = (j-1)*nx + i
	ua(i,j,:)    = ua_g(m,:)
    va(i,j,:)    = va_g(m,:)
    za(i,j,:)    = za_g(m,:)
   ENDDO
  ENDDO
  RETURN
  END


  SUBROUTINE input_namelist(debug,obs_err_u,obs_err_v,obs_err_z,model_flag, &
                            ini_flag,ifactor,rscale,nme,no,ne,nx,ny,nxl,da_flag)
  INCLUDE "../registry/swe.inc"
  RETURN
  END

