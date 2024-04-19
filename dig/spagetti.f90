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
  REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: ut,vt,zt
  REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: uf,vf,zf
  REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: uc,vc,zc
  REAL, ALLOCATABLE, DIMENSION(:)       :: lono,lato
  REAL, ALLOCATABLE, DIMENSION(:)       :: ztmin, zcmin
  REAL, ALLOCATABLE, DIMENSION(:,:)     :: zfmin
  REAL                                  :: time,tfcst,restart,tem,mis
  REAL                                  :: zt_min, zc_min, zf_min
  CHARACTER*50                          :: bfile,cfile,tfile,ofile,afile
  CHARACTER*50                          :: spdaf,spdbf
  INTEGER                               :: i,j,k,itime,irec,ntime,debug,no
  INTEGER                               :: imem
  mis       = -99999
  call input_namelist(debug,nx,ny,ne,restart,tfcst,no)
  nt        = int(tfcst/restart) + 1
  print*,'ana.exe: restart   = ',restart
  print*,'ana.exe: nt        = ',nt
  print*,'ana.exe: nx        = ',nx
  print*,'ana.exe: ny        = ',ny
  print*,'ana.exe: ne        = ',ne
  print*,'ana.exe: no        = ',no
  allocate(uf(nx,ny,ne),vf(nx,ny,ne),zf(nx,ny,ne))
  allocate(uc(nx,ny,nt),vc(nx,ny,nt),zc(nx,ny,nt))
  allocate(ut(nx,ny,nt),vt(nx,ny,nt),zt(nx,ny,nt))
  allocate(ztmin(nt),zcmin(nt),zfmin(nt,ne))
  OPEN(91,file='spa.dat',access='direct',form='unformatted',recl=nx*ny)
  OPEN(92,file='spa.txt')
  WRITE(92,'(A5,10A16)')'Cycle','RMSa','RMSb','RMSo','RMSf','Anlys spread','Bgnd spread'
  itime     = 1
  irec      = 1
  time      = 0
19 continue
!
! reading the truth/ctl
!
  iday      = ifix (time/86400.)
  ihour     = ifix ((time-iday*86400.)/3600.)
  imin      = ifix ((time-iday*86400.-ihour*3600.)/60.)
  if (debug.eq.1) print*,time,imin,ihour,iday
  tfile     = 'tru_00:00:00.dat'
  cfile     = 'ctl_00:00:00.dat'
  bfile     = './mem_000/fsc_00:00:00.dat'
  if (iday.lt.10) then
   write(tfile(6:6),'(1I1)')iday
   write(cfile(6:6),'(1I1)')iday
   write(bfile(16:16),'(1I1)')iday
  elseif (iday.lt.100) then
   write(tfile(5:6),'(1I2)')iday
   write(cfile(5:6),'(1I2)')iday
   write(bfile(15:6),'(1I2)')iday
  else
   print*,'truth.exe: day string length is too long...stop'
   stop
  endif
  if (ihour.lt.10) then
   write(tfile(9:9),'(1I1)')ihour
   write(cfile(9:9),'(1I1)')ihour
   write(bfile(19:19),'(1I1)')ihour
  elseif (iday.lt.100) then
   write(tfile(8:9),'(1I2)')ihour
   write(cfile(8:9),'(1I2)')ihour
   write(bfile(18:19),'(1I2)')ihour
  else
   print*,'truth.exe: hour string length is too long...stop'
   stop
  endif
  if (imin.lt.10) then
   write(tfile(12:12),'(1I1)')imin
   write(cfile(12:12),'(1I1)')imin
   write(bfile(22:22),'(1I1)')imin
  elseif (iday.lt.100) then
   write(tfile(11:12),'(1I2)')imin
   write(cfile(11:12),'(1I2)')imin
   write(bfile(21:22),'(1I2)')imin
  else
   print*,'truth.exe: minute string length is too long...stop'
   stop
  endif 
  !PRINT*,'ana.exe: Open truth file is:  ',tfile(1:30)
  !PRINT*,'ana.exe: Open CTL file is:    ',cfile(1:30)
  OPEN(71,file=tfile,status='old')
  OPEN(76,file=cfile,status='old')
  READ(71,'(6E13.6)')((ut(i,j,itime),i=1,nx),j=1,ny)
  READ(71,'(6E13.6)')((vt(i,j,itime),i=1,nx),j=1,ny)
  READ(71,'(6E13.6)')((zt(i,j,itime),i=1,nx),j=1,ny)
  READ(76,'(6E13.6)')((uc(i,j,itime),i=1,nx),j=1,ny)
  READ(76,'(6E13.6)')((vc(i,j,itime),i=1,nx),j=1,ny)
  READ(76,'(6E13.6)')((zc(i,j,itime),i=1,nx),j=1,ny)
  CLOSE(71)
  CLOSE(76)
!
! reading all member in
!
  do imem = 1,ne
     !PRINT*,'Open forecast file is:    ',bfile(1:50)
     if (imem.lt.10) then
        write(bfile(9:9),'(1I1)')imem
     else
        write(bfile(8:9),'(1I2)')imem
     endif
     OPEN(72,file=bfile,status='old')
     READ(72,'(6E13.6)')((uf(i,j,imem),i=1,nx),j=1,ny)
     READ(72,'(6E13.6)')((vf(i,j,imem),i=1,nx),j=1,ny)
     READ(72,'(6E13.6)')((zf(i,j,imem),i=1,nx),j=1,ny)
     call zmin(zf(:,:,imem),nx,ny,zf_min)
     zfmin(itime,imem) = zf_min
     CLOSE(72)
  enddo
!
! Compute the stardard error devidation
!
  call zmin(zt(:,:,itime),nx,ny,zt_min)
  call zmin(zc(:,:,itime),nx,ny,zc_min)
  ztmin(itime)   = zt_min
  zcmin(itime)   = zc_min
  WRITE(*,'(1I5,9F12.1)')itime,zc_min,zt_min,zfmin(itime,1),zfmin(itime,10),zfmin(itime,20),zfmin(itime,30)
!
! advance and loop now
!
  time     = time + restart*3600
  itime    = itime + 1
  IF (itime.le.nt) GOTO 19
  PRINT*,'Program ends perfectly'
  END

  SUBROUTINE input_namelist(debug,nx,ny,ne,restart,tfcst,no)
  INCLUDE "../../registry/swe.inc"
  RETURN
  END

  SUBROUTINE zmin(a,nx,ny,amin)
  INTEGER nx,ny
  REAL a(nx,ny),amin,t1
  amin = 9999999
  DO i = 1,nx
    DO j = 1,ny
     IF (amin.gt.a(i,j).and.j.le.11) THEN
       amin = a(i,j)
     ENDIF 
    ENDDO
  ENDDO 
  RETURN
  END 

