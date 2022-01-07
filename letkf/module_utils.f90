MODULE module_utils 
USE module_interface
!=======================================================================
!
! NOTE:    MPI interface for 4d-letkf communication
!
! HISTORY: 08 Sep 2011: created by CK
!         
!
!=======================================================================
CONTAINS
  SUBROUTINE initilize_arrays(myrank)
  IMPLICIT NONE
  INTEGER :: myrank
  mis           = -99999.
  irec          = 1
  if (myrank.eq.0) then
   print*,'4dletkf.exe: input namelist is'
   print*,'4dletkf.exe: debug         =',debug
   print*,'4dletkf.exe: obs_err_u     =',obs_err_u
   print*,'4dletkf.exe: obs_err_v     =',obs_err_v
   print*,'4dletkf.exe: obs_err_z     =',obs_err_z
   print*,'4dletkf.exe: model_flag    =',model_flag
   print*,'4dletkf.exe: ini_flag      =',ini_flag
   print*,'4dletkf.exe: ifactor       =',ifactor
   print*,'4dletkf.exe: rscale        =',rscale
   print*,'4dletkf.exe: da_flag       =',da_flag
   print*,'4dletkf.exe: nme           =',nme
   print*,'4dletkf.exe: no            =',no
   print*,'4dletkf.exe: ne            =',ne
   print*,'4dletkf.exe: nxl           =',nxl
   print*,'4dletkf.exe: nx            =',nx
   print*,'4dletkf.exe: ny            =',ny
   print*,'4dletkf.exe: nt            =',nt
   if (debug.ge.1) read*
  endif
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
  allocate(ub(nx,ny,ne,nt),vb(nx,ny,ne,nt),zb(nx,ny,ne,nt))
  allocate(ua(nx,ny,ne),va(nx,ny,ne),za(nx,ny,ne))
  allocate(ub_g(ng,ne,nt),vb_g(ng,ne,nt),zb_g(ng,ne,nt))
  allocate(ub_gm(ng,nt),vb_gm(ng,nt),zb_gm(ng,nt))
  allocate(ua_g(ng,ne),va_g(ng,ne),za_g(ng,ne))
  allocate(ua_gm(ng),va_gm(ng),za_gm(ng))
  allocate(Xb(np,ne),Xb_m(np),Xa(np,ne),Xa_m(np))
  allocate(uo(no,nt),vo(no,nt),zo(no,nt))
  allocate(olat(no,nt),olon(no,nt))
  allocate(uo_g(no,ne,nt),vo_g(no,ne,nt),zo_g(no,ne,nt))
  allocate(uo_gm(no,nt),vo_gm(no,nt),zo_gm(no,nt))
  END SUBROUTINE initilize_arrays

  SUBROUTINE get_obs(myrank)
  IMPLICIT NONE
  INTEGER :: i,j,k,it,myrank
  it            = 0
5 ofile='obs_000.dat'
  it            = it + 1
  if (it.le.9) then
   write(ofile(7:7),'(1I1)')it
  elseif (it.le.99) then
   write(ofile(6:7),'(1I2)')it
  elseif (it.le.999) then
   write(ofile(5:7),'(1I3)')it
  else
   print*,'4dletkf.exe: Stop as too many 4d time interval'
   stop
  endif
  if (myrank.eq.0) print*,'4dletkf.exe: obs file to be opened: ',ofile(1:12)
  open(90,file=ofile,status='old')
  i              = 1
2 continue
  read(90,*,end=3)olon(i,it),olat(i,it),uo(i,it),vo(i,it),zo(i,it)
  if (debug.ge.2) then
   if (i.eq.1) print*,'4dletkf.exe: Checking the obs'
   print*,i,olon(i,it),olat(i,it),uo(i,it),vo(i,it),zo(i,it)
  endif
  i             = i + 1
  if (i.le.no) goto 2
  close(90)
  goto 4
3 print*,'4dletkf.exe: There is not enough obs data as needed by no...stop'
  stop
4 continue
  close(90)
  if (it.lt.nt) goto 5
!
! quick check of obs data
!
  if (myrank.eq.0) then
   allocate(tem2d1(nx,ny))
   do it         = 1,nt
   tem2d1        = mis
   do i          = 1,no
    tem2d1(nint(olon(i,it)),nint(olat(i,it))) = uo(i,it)
   enddo
   write(92,rec=irec)((tem2d1(i,j),i=1,nx),j=1,ny)
   irec          = irec + 1
   enddo
   do it         = 1,nt
   tem2d1        = mis
   do i          = 1,no
    tem2d1(nint(olon(i,it)),nint(olat(i,it))) = vo(i,it)
   enddo
   write(92,rec=irec)((tem2d1(i,j),i=1,nx),j=1,ny)
   irec          = irec + 1
   enddo
   do it         = 1,nt
   tem2d1        = mis
   do i          = 1,no
    tem2d1(nint(olon(i,it)),nint(olat(i,it))) = zo(i,it)
   enddo
   write(92,rec=irec)((tem2d1(i,j),i=1,nx),j=1,ny)
   irec          = irec + 1
   enddo
   deallocate(tem2d1)
  endif
  END SUBROUTINE get_obs

  SUBROUTINE get_bgd(myrank)
  IMPLICIT NONE
  INTEGER :: i,j,k,it,ie,myrank  
  ifile='bgd_000_t000.dat'
  it            = 0
7 continue
  it            = it + 1
  if (it.le.9) then
   write(ifile(12:12),'(1I1)')it
  elseif (it.le.99) then
   write(ifile(11:12),'(1I2)')it
  elseif (it.le.999) then
   write(ifile(10:12),'(1I3)')it
  else
   print*,'4dletkf.exe: Stop as too many 4d time interval'
   stop
  endif
  write(ifile(5:5),'(1I1)')0
  write(ifile(6:6),'(1I1)')0
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
   print*,'4dletkf.exe: Stop as too many ensemble members'
   stop
  endif
  if (myrank.eq.0) print*,'4dletkf.exe: input forecast file: ',ifile(1:20)
  open(91,file=ifile)
  read(91,'(6E13.6)')((ub(i,j,ie,it),i=1,nx),j=1,ny)
  read(91,'(6E13.6)')((vb(i,j,ie,it),i=1,nx),j=1,ny)
  read(91,'(6E13.6)')((zb(i,j,ie,it),i=1,nx),j=1,ny)
  close(91)
  if (ie.lt.ne) goto 8
  if (it.lt.nt) goto 7
  if (myrank.eq.0) then
   do ie          = 1,ne
    write(92,rec=irec)((ub(i,j,ie,1),i=1,nx),j=1,ny)
    irec          = irec + 1
   enddo
   do ie          = 1,ne
    write(92,rec=irec)((vb(i,j,ie,1),i=1,nx),j=1,ny)
    irec          = irec + 1
   enddo
   do ie          = 1,ne
    write(92,rec=irec)((zb(i,j,ie,1),i=1,nx),j=1,ny)
    irec          = irec + 1
   enddo
  endif
  END SUBROUTINE get_bgd

  SUBROUTINE letkf_check(myrank)
  IMPLICIT NONE
  INTEGER :: i,j,k,myrank
  if (debug.ge.1.and.myrank.eq.0) then
   print*,'Local Yb(:1,1) is'
   write(*,'(10E11.3)')(Yb(j,1,1),j=1,nol)
   print*,'Local Yo is'
   write(*,'(10E11.3)')(Yo(j,1),j=1,nol)
   print*,'Local Xb(:1) is'
   write(*,'(10E11.3)')(Xb(j,1),j=1,np)
   print*,'Local Xa(:1) is'
   write(*,'(10E11.3)')(Xa(j,1),j=1,np)
   read*
  endif
  END SUBROUTINE letkf_check

  SUBROUTINE write_analysis(myrank)
  IMPLICIT NONE
  INTEGER :: i,j,k,ie,myrank
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
   print*,'4dletkf.exe: Stop as too many ensemble members'
   stop
  endif
  if (debug.ge.1) print*,'4dletkf.exe:  Open output file is: ',ofile(1:20)
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
  END SUBROUTINE write_analysis
END MODULE module_utils
