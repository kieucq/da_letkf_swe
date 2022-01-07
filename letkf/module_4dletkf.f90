MODULE module_4dletkf
!
! NOTE: This module contains the core of the 4dletkf system   
!
! HISTORY:  
!    - 02 Feb 2010: Updated from the LETKF for L40 model.
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
!=================================================================
!
CONTAINS

  subroutine global_background_perturbation(ub_g,vb_g,zb_g,ub_gm, &
                                            vb_gm,zb_gm,ng,ne,nt)
  implicit none   
  integer :: ng,ne,nt,i,j,k
  real    :: ub_g(ng,ne,nt),vb_g(ng,ne,nt),zb_g(ng,ne,nt)
  real    :: ub_gm(ng,nt),vb_gm(ng,nt),zb_gm(ng,nt)
  ub_gm          = 0
  vb_gm          = 0
  zb_gm          = 0
  do i           = 1,ne
   ub_gm         = ub_gm + ub_g(:,i,:)
   vb_gm         = vb_gm + vb_g(:,i,:)
   zb_gm         = zb_gm + zb_g(:,i,:)
  enddo
  ub_gm          = ub_gm/ne 
  vb_gm          = vb_gm/ne
  zb_gm          = zb_gm/ne
  do i           = 1,ne
   ub_g(:,i,:)   = ub_g(:,i,:) - ub_gm(:,:)
   vb_g(:,i,:)   = vb_g(:,i,:) - vb_gm(:,:)
   zb_g(:,i,:)   = zb_g(:,i,:) - zb_gm(:,:)
  enddo
  return
  end subroutine global_background_perturbation


  subroutine local_obs(m,n,nx,ny,nxl,nt,olon,olat,no,nol,debug)
  implicit none
  integer m,n,no,nxl,nol,nx,ny,nt
  real olon(no,nt),olat(no,nt)
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
   if (olon(i,1).ge.istart.and.olon(i,1).le.iend.and. &
       olat(i,1).ge.jstart.and.olat(i,1).le.jend) then
    nol          = nol + 1
   endif
  enddo  
  nol            = nol*3
  return
  end subroutine local_obs

  subroutine global2local(ub_g,vb_g,zb_g,ub_gm,vb_gm,zb_gm,uo_g,vo_g,zo_g,          &
                          uo_gm,vo_gm,zo_gm,uo,vo,zo,Xb,Xb_m,Yb,Yb_m,Yo,np,nol,     &
		          olon,olat,m,n,nx,ny,ng,ne,no,nxl,nt,rscale,rho,debug)
  implicit none
  integer             :: nx,ny,no,np,nol,ng,ne,m,n,nxl,nt
  real, intent(in)    :: ub_g(ng,ne,nt),vb_g(ng,ne,nt),zb_g(ng,ne,nt)
  real, intent(in)    :: ub_gm(ng,nt),vb_gm(ng,nt),zb_gm(ng,nt)
  real, intent(in)    :: uo_g(no,ne,nt),vo_g(no,ne,nt),zo_g(no,ne,nt)
  real, intent(in)    :: uo_gm(no,nt),vo_gm(no,nt),zo_gm(no,nt)
  real, intent(in)    :: uo(no,nt),vo(no,nt),zo(no,nt)
  real, intent(in)    :: olat(no,nt),olon(no,nt),rscale
  real, intent(out)   :: Xb(np,ne),Xb_m(np),rho(nol)
  real, intent(out)   :: Yb(nol,ne,nt),Yb_m(nol,nt),Yo(nol,nt)
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
    Xb(k,1:ne)         = ub_g(id,1:ne,1) 
    Xb(np/3+k,1:ne)    = vb_g(id,1:ne,1)
    Xb(2*np/3+k,1:ne)  = zb_g(id,1:ne,1)  
    Xb_m(k)            = ub_gm(id,1)
    Xb_m(np/3+k)       = vb_gm(id,1)
    Xb_m(2*np/3+k)     = zb_gm(id,1)
   enddo
  enddo
!
! find and assign the global obs ensemble to the local obs patch.
! Note that the cyclinc boundary is applied only for +/- 5 points
! at each end of the domain.
!
  k              = 0
  do i           = 1,no   
   if (olon(i,1).ge.istart.and.olon(i,1).le.iend.and. &
       olat(i,1).ge.jstart.and.olat(i,1).le.jend) then
    k            = k + 1
    Yb(k,:,:)            = uo_g(i,:,:)
    Yb(nol/3+k,:,:)      = vo_g(i,:,:)
    Yb(2*nol/3+k,:,:)    = zo_g(i,:,:)
    Yb_m(k,:)            = uo_gm(i,:)
    Yb_m(nol/3+k,:)      = vo_gm(i,:)
    Yb_m(2*nol/3+k,:)    = zo_gm(i,:)
    Yo(k,:)              = uo(i,:)
    Yo(nol/3+k,:)        = vo(i,:)
    Yo(2*nol/3+k,:)      = zo(i,:)
    rho(k)               = exp(-((olon(i,1)-m)**2. + (olat(i,1)-n)**2.)/rscale**2)
    rho(nol/3+k)         = rho(k)
    rho(2*nol/3+k)       = rho(k)
   endif
  enddo
  if (debug.ge.1) print*,'k,nol/3',k,nol/3
  if (k.ne.nol/3) then
   print*,'4dletkf.exe: The number of obs in local patch does not match the array shape'
   stop
  endif
  return
  end subroutine global2local
  
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
  end subroutine local2global

  subroutine letkf_main(Xa,R,Xb,Xb_m,Yb,Yb_m,Yo,nxl,nol,ne,nt,ifactor,rho,da_flag,debug)
  use common
  use common_mtx
  implicit none
  integer nxl,nol,ne,debug,nt
  real, intent(in)  :: R(nol,nol),Xb(nxl,ne),Xb_m(nxl),rho(nol)
  real, intent(in)  :: Yb(nol,ne,nt),Yb_m(nol,nt),Yo(nol,nt),ifactor
  real, intent(out) :: Xa(nxl,ne)  
  integer           :: da_flag
  real              :: C(ne,nol,nt),tem1(ne,ne),dy(nol)
  real              :: temc(ne,nol),temy(nol,ne)
  real              :: Wa(ne,ne),Pat(ne,ne),wa_m(ne),Rinv(nol,nol)
  integer           :: i1,j1,k1,i,j,k,it
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
! compute matrcies $C_l \forall l \in (1,...,L)$ (advance from 
! step 4 in H06 to the coefficients of Eq. (6) in F07)
!
  do it         = 1,nt
   temy(:,:)    = Yb(:,:,it)
   temc         = matmul(transpose(temy),Rinv)
   C(:,:,it)    = temc(:,:)
  enddo
  if (debug.ge.2) then 
   print*,'step 4: C matrix'
   do i1        = 1,ne
    write(0,'(20F6.3)')(C(i1,j1,1),j1=1,nol)
   enddo
  endif
!
! compute matrix Pa tilde (step 5 in H06)
!
  tem1          = 0.
  do it         = 1,nt
   temy(:,:)    = Yb(:,:,it)
   temc(:,:)    = C(:,:,it)      
   tem1         = tem1 + matmul(temc,temy)
  enddo
  do j          = 1,ne
   tem1(j,j)    = tem1(j,j) + (ne-1)/(1+ifactor)
  enddo
  call mtx_inv(ne,tem1,Pat)
  if (debug.ge.2) then
   print*,'step 5: Yb matrix'
   do i1        = 1,nol
    write(0,'(20F6.3)')(Yb(i1,j1,1),j1=1,20)
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
  do it         = 1,nt
   dy(:)        = Yo(:,it) - Yb_m(:,it)
   temc(:,:)    = C(:,:,it)
   wa_m         = matmul(matmul(Pat,temc),dy)
  enddo
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
  end subroutine letkf_main

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
  end subroutine observation_operator

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
  end subroutine localiztion_operator

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
  end subroutine compute_Htilde

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
  end subroutine observational_err_cov_mtx


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
  end subroutine compute_Rinv
  
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
  end subroutine obs_increment

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
  end subroutine analysis_mean

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
  end subroutine convert_vector_array

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
  end subroutine convert_vector_array1

  subroutine convert_array_vector(ub,vb,zb,nx,ny,ne,nt,ub_g,vb_g,zb_g,ng)
  implicit none
  integer nx,ny,ne,ng,nt
  real ub(nx,ny,ne,nt),vb(nx,ny,ne,nt),zb(nx,ny,ne,nt)
  real ub_g(ng,ne,nt),vb_g(ng,ne,nt),zb_g(ng,ne,nt)
  integer i,j,k,m,n,nv
  do i           = 1,ng
   n             = i/nx + 1
   m             = mod(i,nx)
   if (m.eq.0) then
    m            = nx
    n            = n - 1
   endif
   ub_g(i,:,:)   = ub(m,n,:,:)
   vb_g(i,:,:)   = vb(m,n,:,:)
   zb_g(i,:,:)   = zb(m,n,:,:)
  enddo
  return
  end subroutine convert_array_vector

  SUBROUTINE global_observation(uo_g,vo_g,zo_g,uo_gm,vo_gm,zo_gm,  &
                                ub,vb,zb,olat,olon,ne,nx,ny,ng,no,nt)
  IMPLICIT NONE
  integer nx,ny,ng,no,ne,nt
  REAL ub(nx,ny,ne,nt),vb(nx,ny,ne,nt),zb(nx,ny,ne,nt)
  REAL uo_g(no,ne,nt),vo_g(no,ne,nt),zo_g(no,ne,nt)
  REAL uo_gm(no,nt),vo_gm(no,nt),zo_gm(no,nt)
  REAL olat(no,nt),olon(no,nt)
  integer i,j,k,m,n
!
! this is the first step in Hunt et al, for which H is simple = I because
! observations are given at the grid point. Also, it will be assumed that
! obs at all time slices are exactly the same. Otherwise, there will be
! separate H operator for each time slice.  
!
  DO i           = 1,no
   m             = nint(olon(i,1))
   n             = nint(olat(i,1))
   uo_g(i,:,:)   = ub(m,n,:,:)
   vo_g(i,:,:)   = vb(m,n,:,:)
   zo_g(i,:,:)   = zb(m,n,:,:)
  ENDDO
!
! compute the $\ba Y$ according to H05 
!
  uo_gm          = 0.
  vo_gm          = 0.
  zo_gm          = 0.
  do k           = 1,ne
   uo_gm(:,:)    = uo_gm(:,:) + uo_g(:,k,:)
   vo_gm(:,:)    = vo_gm(:,:) + vo_g(:,k,:)
   zo_gm(:,:)    = zo_gm(:,:) + zo_g(:,k,:)
  enddo 
  uo_gm          = uo_gm/ne
  vo_gm          = vo_gm/ne
  zo_gm          = zo_gm/ne
!
! substract the mean of the total global obs from the total global obs
! to create pertrubation global obs
!
  do k           = 1,ne
   uo_g(:,k,:)   = uo_g(:,k,:) - uo_gm(:,:)
   vo_g(:,k,:)   = vo_g(:,k,:) - vo_gm(:,:)
   zo_g(:,k,:)   = zo_g(:,k,:) - zo_gm(:,:)
  enddo
  RETURN
  END

  SUBROUTINE vector2array(ua_g,va_g,za_g,ua,va,za,ng,ne,nx,ny) 
  IMPLICIT NONE
  integer nx,ny,ne,ng
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
  END SUBROUTINE vector2array


  SUBROUTINE input_namelist(debug,obs_err_u,obs_err_v,obs_err_z,model_flag, &
                            ini_flag,ifactor,rscale,nme,no,ne,nx,ny,nxl,nt,da_flag)
  INCLUDE "../registry/swe.inc"
  RETURN
  END SUBROUTINE input_namelist
END MODULE module_4dletkf

