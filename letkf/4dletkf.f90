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
!    - 08 Sep 2011: revise for 4d-letkf and parallelize the code. Note
!                   that the code is not totally optimized in terms of
!                   memory yet. The focus is only on parallel computing.
!    - 09 Sep 2011: parallelization works well. Need to impletment 
!                   4D-LETKF next.
!    - 12 Sep 2011: implement 4d_letkf  
!    - 14 Sep 2011: modulize for better organization
!
! REFERENCE:
!    - Hunt et al., 2005: arxiv:physics/0511236
!    - Fertig et al. 2007: Tellus A, 96-100
!
! AUTHOR: 
!      Chanh Q. Kieu, Lecturer
!      Lab for Climate and Weather Research
!      Vietnam National University, Hanoi
!      email: chanhkq@vnu.edu.vn
!
! COPYRIGHT: (C) 2011
!
!=========================================================================
!
  program LETKF_4D_swe
  use common
  use common_mtx
  use common_mpi
  use module_interface
  use module_4dletkf
  use module_utils
  implicit none
  include 'mpif.h'
!
! MPI variables
!
  integer              :: ierr,nprocs,myrank,irank        
  integer              :: ista,iend,iprev,inext    
  integer              :: ureq1,vreq1,zreq1
  integer              :: istatus(MPI_STATUS_SIZE)
  integer, allocatable :: itype(:),ureq(:),vreq(:),zreq(:)
  integer              :: i,j,m,n,k,ie,it     
!
! initialize MPI communications
!
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
!
! reading namelist now
!
  call input_namelist(debug,obs_err_u,obs_err_v,obs_err_z, & 
                      model_flag,ini_flag,ifactor,rscale,  &
                      nme,no,ne,nx,ny,nxl,nt,da_flag)
!
! initialized model parameters
!
  call initilize_arrays(myrank)
!
! reading observation data
!
  call get_obs(myrank)
!
! reading forecast files
!
  call get_bgd(myrank)
!
! convert from 2d arrays to 1d column vectors and create the 
! background first for the analysis
!
  call convert_array_vector(ub,vb,zb,nx,ny,ne,nt,ub_g,vb_g,zb_g,ng)
  ua_g(:,:)      = ub_g(:,:,1)
  va_g(:,:)      = vb_g(:,:,1)
  za_g(:,:)      = zb_g(:,:,1)
! 
! compute the global pertrubation background observation ensemble 
! from the global background. (step 1 in H04)
! 
  call global_observation(uo_g,vo_g,zo_g,uo_gm,vo_gm,zo_gm,     &
                          ub,vb,zb,olat,olon,ne,nx,ny,ng,no,nt)
  deallocate(ub,vb,zb)
!
! compute the global mean and ensemble of pertubation background.
! Note that the global background will not be needed from now on,
! so the xb_g arrays will now store perturbation, not the total.
! This is step 2 in H04.
!
  call global_background_perturbation(ub_g,vb_g,zb_g,ub_gm,     &
                                      vb_gm,zb_gm,ng,ne,nt)
!
! indexing the MPI jobs
!
  call para_range(1,ng,nprocs,myrank,ista,iend)
  iprev         = myrank-1
  inext         = myrank+1
  if (myrank.eq.0)        iprev = MPI_PROC_NULL
  if (myrank.eq.nprocs-1) inext = MPI_PROC_NULL
  print*,'myrank = ',myrank,ista,iend
!
! Now loop over all grid point with the corresponding
! local patch.
!
  grid_loop: do i = ista,iend
   print*,'global grid loop @ i =',i,myrank
!
! compute first the grid location and var.
!
   n             = i/nx + 1
   m             = mod(i,nx)
   if (m.eq.0) then 
    n            = n - 1
    m            = nx
   endif             
   if (debug.ge.1) print*,'Working at the local patch:',i,m,n
!
! finding the number of local observation within each patch. 
! Note that the size of the local patch in 2D array will be 
! (m-nxl/2,m+nxl/2) x (n-nxl/2,n+nxl/2). This step is merely 
! for defining the local size 
!
   call local_obs(m,n,nx,ny,nxl,nt,olon,olat,no,nol,debug)
   if (debug.ge.1) print*,'Number of local obs is',nol,myrank
   if (nol.lt.1) goto 10
!
! now we will project from global to local patch at each grid 
! point (m,n). The cross correlation can be handled separately 
! to create a new set of, e.g., un from z. (u,un) are then 
! treated as the total obs of u around (m,n).  
!
   allocate(Yb(nol,ne,nt),Yb_m(nol,nt),Yo(nol,nt),Ro(nol,nol))
   allocate(Ro(nol,nol),rho(nol))
   call global2local(ub_g,vb_g,zb_g,ub_gm,vb_gm,zb_gm,uo_g,     &
                     vo_g,zo_g,uo_gm,vo_gm,zo_gm,uo,vo,zo,Xb,   &
                     Xb_m,Yb,Yb_m,Yo,np,nol,olon,olat,m,n,nx,ny,&
                     ng,ne,no,nxl,nt,rscale,rho,debug)
!
! define the local observational error covariance matrix R. Note 
! that obs are assumed to be the same at all slices, i.e.. $H_l 
! = I \forall l \in (1,...,L)$ (see Fertig et al. (2007)). Also, 
! all obs are of the same errors this experiment.
!
   call observational_err_cov_mtx(Ro,nol,obs_err_u,obs_err_v,obs_err_z)
!
! calling the LETKF core now.
!
    call letkf_main(Xa,Ro,Xb,Xb_m,Yb,Yb_m,Yo,np,nol,ne,nt,      &
                    ifactor,rho,da_flag,debug)
    call letkf_check(myrank)
!
! project from local patch to global patch at the point $i$, which 
! corresponds to point (m,n) on the 2D grid
!
    call local2global(Xa,np,ne,ua_g,va_g,za_g,nx,ny,ng,m,n,nxl,i)
!
! deallocate the arrays
!
    deallocate(Yb,Yb_m,Yo,Ro,rho)
10  continue
  enddo grid_loop
  allocate(itype(0:nprocs-1))
  allocate(ureq(0:nprocs-1))  
  allocate(vreq(0:nprocs-1))
  allocate(zreq(0:nprocs-1))
  do irank       = 0,nprocs - 1
   call para_range(1,ng,nprocs,irank,ista,iend)
   call para_type_block2(1,ng,1,ista,iend,1,ne,MPI_REAL,itype(irank))
  enddo
  if (myrank.eq.0) then
!
! gather the updates from all cores (non-contiguous memory)
!
   do irank      = 1,nprocs-1
    call MPI_IRECV(ua_g,1,itype(irank),irank,1,MPI_COMM_WORLD,ureq(irank),ierr)
    call MPI_IRECV(va_g,1,itype(irank),irank,1,MPI_COMM_WORLD,vreq(irank),ierr)
    call MPI_IRECV(za_g,1,itype(irank),irank,1,MPI_COMM_WORLD,zreq(irank),ierr)
   enddo
   do irank      = 1,nprocs-1
    call MPI_WAIT(ureq(irank),istatus,ierr)
    call MPI_WAIT(vreq(irank),istatus,ierr)
    call MPI_WAIT(zreq(irank),istatus,ierr)
   enddo
!
! convert from 1D arrays to 2D arrays
!
   call vector2array(ua_g,va_g,za_g,ua,va,za,ng,ne,nx,ny) 
!
! output analysis  
!
   call write_analysis(myrank)
  else
   call MPI_ISEND(ua_g,1,itype(myrank),0,1,MPI_COMM_WORLD,ureq1,ierr)
   call MPI_ISEND(va_g,1,itype(myrank),0,1,MPI_COMM_WORLD,vreq1,ierr)
   call MPI_ISEND(za_g,1,itype(myrank),0,1,MPI_COMM_WORLD,zreq1,ierr)
   call MPI_WAIT(ureq1,istatus,ierr)
   call MPI_WAIT(vreq1,istatus,ierr)
   call MPI_WAIT(zreq1,istatus,ierr)
  endif
  call MPI_FINALIZE(ierr)
  if (myrank.eq.0) print*,'4dletkf.exe: LETKF finishes safely...!'
  end
