!
! Note:
!         This program is for performing an analysis of the output
!         from the ctl, truth, and assimilation runs.
!
! History: Created Feb 6, 2009
!
! Author: Chanh Q. Kieu
!
!===================================================================
  PROGRAM analysis
  IMPLICIT NONE
  INTEGER                               :: nx,ny,ne 
  INTEGER                               :: nt,iday,ihour,imin
  REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: ua,va,za
  REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: ut,vt,zt
  REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: ub,vb,zb
  REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: uo,vo,zo
  REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: uc,vc,zc
  REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: uspa,vspa,zspa
  REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: uspb,vspb,zspb
  REAL, ALLOCATABLE, DIMENSION(:)       :: rmsa,rmsb,rmso,rmsc,tu,tv,tz,spa,spb
  REAL, ALLOCATABLE, DIMENSION(:)       :: lono,lato
  REAL                                  :: time,tfcst,restart,tem,mis
  CHARACTER*50                          :: bfile,cfile,tfile,ofile,afile
  CHARACTER*50                          :: spdaf,spdbf
  INTEGER                               :: i,j,k,itime,irec,ntime,debug,no
  mis       = -99999
  call input_namelist(debug,nx,ny,ne,restart,tfcst,no)
  nt        = int(tfcst/restart) + 1
  print*,'ana.exe: restart   = ',restart
  print*,'ana.exe: nt        = ',nt
  print*,'ana.exe: nx        = ',nx
  print*,'ana.exe: ny        = ',ny
  print*,'ana.exe: ne        = ',ne
  print*,'ana.exe: no        = ',no
  allocate(ut(nx,ny,nt),vt(nx,ny,nt),zt(nx,ny,nt))
  allocate(ua(nx,ny,nt),va(nx,ny,nt),za(nx,ny,nt))
  allocate(ub(nx,ny,nt),vb(nx,ny,nt),zb(nx,ny,nt))
  allocate(uc(nx,ny,nt),vc(nx,ny,nt),zc(nx,ny,nt))
  allocate(uo(nx,ny,nt),vo(nx,ny,nt),zo(nx,ny,nt))
  allocate(uspa(nx,ny,nt),vspa(nx,ny,nt),zspa(nx,ny,nt))
  allocate(uspb(nx,ny,nt),vspb(nx,ny,nt),zspb(nx,ny,nt))
  allocate(rmsa(nt),rmsb(nt),rmso(nt),rmsc(nt),tu(no),tv(no),tz(no))
  allocate(lato(no),lono(no),spa(nt),spb(nt))
  OPEN(91,file='ana.dat',access='direct',form='unformatted',recl=nx*ny*4)
  OPEN(92,file='ana.txt')
  WRITE(92,'(A5,10A16)')'Cycle','RMSa','RMSb','RMSo','RMSf','Anlys spread','Bgnd spread'
  itime     = 1
  irec      = 1
  time      = 0
  uo        = mis
  vo        = mis
  zo        = mis
19 continue
!
! reading the truth
!
  iday      = ifix (time/86400.)
  ihour     = ifix ((time-iday*86400.)/3600.)
  imin      = ifix ((time-iday*86400.-ihour*3600.)/60.)
  if (debug.eq.1) print*,time,imin,ihour,iday
  tfile     = 'tru_00:00:00.dat'
  bfile     = 'bgd_00:00:00.dat'
  ofile     = 'obs_00:00:00.dat'
  afile     = 'ana_00:00:00.dat'
  cfile     = 'ctl_00:00:00.dat'
  spdaf     = 'spa_00:00:00.dat'
  spdbf     = 'spb_00:00:00.dat'
  if (iday.lt.10) then
   write(ofile(6:6),'(1I1)')iday
   write(afile(6:6),'(1I1)')iday
   write(tfile(6:6),'(1I1)')iday
   write(bfile(6:6),'(1I1)')iday
   write(cfile(6:6),'(1I1)')iday
   write(spdaf(6:6),'(1I1)')iday
   write(spdbf(6:6),'(1I1)')iday
  elseif (iday.lt.100) then
   write(ofile(5:6),'(1I2)')iday
   write(afile(5:6),'(1I2)')iday
   write(tfile(5:6),'(1I2)')iday
   write(bfile(5:6),'(1I2)')iday
   write(cfile(5:6),'(1I2)')iday
   write(spdaf(5:6),'(1I2)')iday
   write(spdbf(5:6),'(1I2)')iday
  else
   print*,'truth.exe: day string length is too long...stop'
   stop
  endif
  if (ihour.lt.10) then
   write(ofile(9:9),'(1I1)')ihour
   write(afile(9:9),'(1I1)')ihour
   write(tfile(9:9),'(1I1)')ihour
   write(bfile(9:9),'(1I1)')ihour
   write(cfile(9:9),'(1I1)')ihour
   write(spdaf(9:9),'(1I1)')ihour
   write(spdbf(9:9),'(1I1)')ihour
  elseif (iday.lt.100) then
   write(ofile(8:9),'(1I2)')ihour
   write(afile(8:9),'(1I2)')ihour
   write(tfile(8:9),'(1I2)')ihour
   write(bfile(8:9),'(1I2)')ihour
   write(cfile(8:9),'(1I2)')ihour
   write(spdaf(8:9),'(1I2)')ihour
   write(spdbf(8:9),'(1I2)')ihour
  else
   print*,'truth.exe: hour string length is too long...stop'
   stop
  endif
  if (imin.lt.10) then
   write(ofile(12:12),'(1I1)')imin
   write(afile(12:12),'(1I1)')imin
   write(tfile(12:12),'(1I1)')imin
   write(bfile(12:12),'(1I1)')imin
   write(cfile(12:12),'(1I1)')imin
   write(spdaf(12:12),'(1I1)')imin
   write(spdbf(12:12),'(1I1)')imin 
  elseif (iday.lt.100) then
   write(ofile(11:12),'(1I2)')imin
   write(afile(11:12),'(1I2)')imin
   write(tfile(11:12),'(1I2)')imin
   write(bfile(11:12),'(1I2)')imin
   write(cfile(11:12),'(1I2)')imin
   write(spdaf(11:12),'(1I2)')imin
   write(spdbf(11:12),'(1I2)')imin 
  else
   print*,'truth.exe: minute string length is too long...stop'
   stop
  endif 
  IF (debug.eq.1) PRINT*,'ana.exe: Open truth file is:  ',tfile(1:30)
  IF (debug.eq.1) PRINT*,'ana.exe: Open analysis file is: ',afile(1:30)
  IF (debug.eq.1) PRINT*,'ana.exe: Open background file is:  ',bfile(1:30)
  IF (debug.eq.1) PRINT*,'ana.exe: Open observation file is:  ',ofile(1:30)
  IF (debug.eq.1) PRINT*,'ana.exe: Open analysis spread is:  ',spdaf(1:30)
  IF (debug.eq.1) PRINT*,'ana.exe: Open background spread is:  ',spdbf(1:30)
  OPEN(71,file=tfile,status='old')
  OPEN(73,file=bfile,status='old')
  OPEN(74,file=ofile,status='old')
  OPEN(75,file=afile,status='old')
  OPEN(76,file=cfile,status='old')
  OPEN(77,file=spdaf,status='old')
  OPEN(78,file=spdbf,status='old')
  if (debug.eq.1) print*,'ana.exe: done opening data file'
  READ(71,'(6E13.6)')((ut(i,j,itime),i=1,nx),j=1,ny)
  READ(71,'(6E13.6)')((vt(i,j,itime),i=1,nx),j=1,ny)
  READ(71,'(6E13.6)')((zt(i,j,itime),i=1,nx),j=1,ny)
  READ(73,'(6E13.6)')((ub(i,j,itime),i=1,nx),j=1,ny)
  READ(73,'(6E13.6)')((vb(i,j,itime),i=1,nx),j=1,ny)
  READ(73,'(6E13.6)')((zb(i,j,itime),i=1,nx),j=1,ny)
  READ(75,'(6E13.6)')((ua(i,j,itime),i=1,nx),j=1,ny)
  READ(75,'(6E13.6)')((va(i,j,itime),i=1,nx),j=1,ny)
  READ(75,'(6E13.6)')((za(i,j,itime),i=1,nx),j=1,ny)
  READ(76,'(6E13.6)')((uc(i,j,itime),i=1,nx),j=1,ny)
  READ(76,'(6E13.6)')((vc(i,j,itime),i=1,nx),j=1,ny)
  READ(76,'(6E13.6)')((zc(i,j,itime),i=1,nx),j=1,ny)
  READ(77,'(6E13.6)')((uspa(i,j,itime),i=1,nx),j=1,ny)
  READ(77,'(6E13.6)')((vspa(i,j,itime),i=1,nx),j=1,ny)
  READ(77,'(6E13.6)')((zspa(i,j,itime),i=1,nx),j=1,ny)
  READ(78,'(6E13.6)')((uspb(i,j,itime),i=1,nx),j=1,ny)
  READ(78,'(6E13.6)')((vspb(i,j,itime),i=1,nx),j=1,ny)
  READ(78,'(6E13.6)')((zspb(i,j,itime),i=1,nx),j=1,ny)
  IF (no.eq.nx*ny) THEN
   DO k     = 1,no
    i       = mod(k,nx)
    j       = k/nx + 1
    if (i.eq.0) then
     i      = nx
     j      = j - 1
    endif
    READ(74,*)tem,tem,uo(i,j,itime),vo(i,j,itime),zo(i,j,itime)
   ENDDO
  ELSE
   DO k     = 1,no
    READ(74,*)lono(k),lato(k),tu(k),tv(k),tz(k)
    i       = nint(lono(k))
    j       = nint(lato(k))
    uo(i,j,itime) = tu(k)
    vo(i,j,itime) = tv(k)
    zo(i,j,itime) = tz(k)
   ENDDO
  ENDIF
  CLOSE(71)
  CLOSE(73)
  CLOSE(74)
  CLOSE(75)
  if (debug.eq.1) print*,'ana.exe: done reading data'

  WRITE(91,rec=irec)((ut(i,j,itime),i=1,nx),j=1,ny)
  irec     = irec + 1
  WRITE(91,rec=irec)((vt(i,j,itime),i=1,nx),j=1,ny)
  irec     = irec + 1
  WRITE(91,rec=irec)((zt(i,j,itime),i=1,nx),j=1,ny)
  irec     = irec + 1
  WRITE(91,rec=irec)((uo(i,j,itime),i=1,nx),j=1,ny)
  irec     = irec + 1
  WRITE(91,rec=irec)((vo(i,j,itime),i=1,nx),j=1,ny)
  irec     = irec + 1
  WRITE(91,rec=irec)((zo(i,j,itime),i=1,nx),j=1,ny)
  irec     = irec + 1
  WRITE(91,rec=irec)((ub(i,j,itime),i=1,nx),j=1,ny)
  irec     = irec + 1
  WRITE(91,rec=irec)((vb(i,j,itime),i=1,nx),j=1,ny)
  irec     = irec + 1
  WRITE(91,rec=irec)((zb(i,j,itime),i=1,nx),j=1,ny)
  irec     = irec + 1
  WRITE(91,rec=irec)((ua(i,j,itime),i=1,nx),j=1,ny)
  irec     = irec + 1
  WRITE(91,rec=irec)((va(i,j,itime),i=1,nx),j=1,ny)
  irec     = irec + 1
  WRITE(91,rec=irec)((za(i,j,itime),i=1,nx),j=1,ny)
  irec     = irec + 1
  WRITE(91,rec=irec)((uc(i,j,itime),i=1,nx),j=1,ny)
  irec     = irec + 1
  WRITE(91,rec=irec)((vc(i,j,itime),i=1,nx),j=1,ny)
  irec     = irec + 1
  WRITE(91,rec=irec)((zc(i,j,itime),i=1,nx),j=1,ny)
  irec     = irec + 1
  WRITE(91,rec=irec)((uspa(i,j,itime),i=1,nx),j=1,ny)
  irec     = irec + 1
  WRITE(91,rec=irec)((vspa(i,j,itime),i=1,nx),j=1,ny)
  irec     = irec + 1
  WRITE(91,rec=irec)((zspa(i,j,itime),i=1,nx),j=1,ny)
  irec     = irec + 1
  WRITE(91,rec=irec)((uspb(i,j,itime),i=1,nx),j=1,ny)
  irec     = irec + 1
  WRITE(91,rec=irec)((vspb(i,j,itime),i=1,nx),j=1,ny)
  irec     = irec + 1
  WRITE(91,rec=irec)((zspb(i,j,itime),i=1,nx),j=1,ny)
  irec     = irec + 1
  if (debug.eq.1) print*,'ana.exe: done writing out for checking data'
!
! Compute the stardard error devidation
!
  rmso(itime)  = 0.
  rmsb(itime)  = 0.
  rmsa(itime)  = 0.
  rmsc(itime)  = 0.
  spa(itime)   = 0.
  spb(itime)   = 0.
  DO i         = 1,nx
   DO j        = 1,ny
   rmsa(itime) = rmsa(itime) + 0.5*(ua(i,j,itime)-ut(i,j,itime))**2      &
                             + 0.5*(va(i,j,itime)-vt(i,j,itime))**2      &
                             + 9.81*abs(za(i,j,itime)-zt(i,j,itime)) 
   rmsb(itime) = rmsb(itime) + 0.5*(ub(i,j,itime)-ut(i,j,itime))**2      &
                             + 0.5*(vb(i,j,itime)-vt(i,j,itime))**2      &
                             + 9.81*abs(zb(i,j,itime)-zt(i,j,itime))
   rmsc(itime) = rmsc(itime) + 0.5*(uc(i,j,itime)-ut(i,j,itime))**2      &
                             + 0.5*(vc(i,j,itime)-vt(i,j,itime))**2      &
                             + 9.81*abs(zc(i,j,itime)-zt(i,j,itime))
   spa(itime)  = spa(itime)  + 0.5*uspa(i,j,itime)**2 + 0.5*vspa(i,j,itime)**2  &
                             + 9.81*zspa(i,j,itime)
   spb(itime)  = spb(itime)  + 0.5*uspb(i,j,itime)**2 + 0.5*vspb(i,j,itime)**2  &
                             + 9.81*zspb(i,j,itime)
   IF (no.eq.nx*ny)                                                      &
   rmso(itime) = rmso(itime) + 0.5*(uo(i,j,itime)-ut(i,j,itime))**2      &
                             + 0.5*(vo(i,j,itime)-vt(i,j,itime))**2      &
                             + 9.81*abs(zo(i,j,itime)-zt(i,j,itime))
   ENDDO
  ENDDO
  rmsa(itime)  = sqrt(rmsa(itime)/nx/ny)
  rmsb(itime)  = sqrt(rmsb(itime)/nx/ny)
  rmso(itime)  = sqrt(rmso(itime)/nx/ny)
  rmsc(itime)  = sqrt(rmsc(itime)/nx/ny)
  spa(itime)   = sqrt(spa(itime)/nx/ny)
  spb(itime)   = sqrt(spb(itime)/nx/ny)
  IF (debug.eq.1) PRINT*,'ana.exe: Checking error now',itime,rmsa(itime),rmsb(itime),rmso(itime),rmsc(itime)
  WRITE(92,'(I5,10F16.4)')itime,rmsa(itime),rmsb(itime),rmso(itime),rmsc(itime),spa(itime),spb(itime)
!
! advance and loop now
!
  time     = time + restart*3600
  itime    = itime + 1
  IF (itime.le.nt) GOTO 19
  PRINT*,'ana.exe: Program ends perfectly'
  END

  SUBROUTINE input_namelist(debug,nx,ny,ne,restart,tfcst,no)
  INCLUDE "../registry/swe.inc"
  RETURN
  END

