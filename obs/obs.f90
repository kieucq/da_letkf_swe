!
! This program is for creating an ensemble of  
! observation that are perturbed from the truth
!
! Author: Chanh Q. Kieu
!
! Date: June 1, 2009
!
!===================================================================
!
  PROGRAM obs
  IMPLICIT NONE
  INTEGER,PARAMETER :: nx=114,ny=17,nv=nx*ny 
  REAL              :: obs_err_u,obs_err_v,obs_err_z
  REAL              :: u(nx,ny),v(nx,ny),z(nx,ny)
  REAL(8)           :: rnd(nv)
  REAL              :: ur(nx,ny),vr(nx,ny),zr(nx,ny)
  CHARACTER*100     :: ofile,ifile,temc
  REAL              :: time,dt,tfcst,restart
  INTEGER           :: iday,ihour,imin,obs_flag,ne,itime,ntime
  INTEGER           :: id,jd,i,j,k,irec,debug,no,nt,out_flag,nprint
  REAL              :: t_window,dt_window
  INTEGER           :: icen,jcen
  OPEN(91,file='obs.dat',access='direct',form='unformatted',recl=nx*ny*4)
!
! reading namelist
!
  CALL input_namelist(debug,tfcst,restart,obs_err_u,obs_err_v, &
                      obs_err_z,no,ne,obs_flag,t_window,dt_window,icen,jcen)
  ntime     = int(tfcst*3600/dt_window) + 1
  nt        = nint(t_window/dt_window) + 1                                                         
  out_flag  = nint(restart*3600/dt_window)
  print*,'dt_window = ',dt_window
  print*,'obs.exe: nt             =',nt
  print*,'obs.exe: out_flag       =',out_flag
  PRINT*,'obs.exe: restart        =',restart
  PRINT*,'obs.exe: tfcst          =',tfcst
  PRINT*,'obs.exe: obs_error_u is: ',obs_err_u
  PRINT*,'obs.exe: obs_error_v is: ',obs_err_v
  PRINT*,'obs.exe: obs_error_z is: ',obs_err_z
  PRINT*,'obs.exe: number of ensemble is: ',ne
  PRINT*,'obs.exe: number of obs is: ',no
  PRINT*,'obs.exe: obs_flag: ',obs_flag
  IF (debug.eq.1) READ*
!
! start up output files
!
  ifile     = './tru_00:00:00.dat'
  ofile     = './obs_00:00:00.dat'
  irec      = 1
  itime     = 1
  time      = 0.
  nprint    = 0
19 continue
!
! reading the truth
!
  if (mod(itime,out_flag).eq.1) nprint = 1
  if (nprint.gt.0.and.nprint.le.nt) then
   print*,'time,itime,nprint:',time,itime,nprint
   iday      = ifix (time/86400.)
   ihour     = ifix ((time-iday*86400.)/3600.)
   imin      = ifix ((time-iday*86400.-ihour*3600.)/60.)
   if (debug.ge.1) print*,time,iday,ihour,imin
   ifile   = 'tru_00:00:00.dat'
   if (iday.lt.10) then
    write(ifile(6:6),'(1I1)')iday
   elseif (iday.lt.100) then
    write(ifile(5:6),'(1I2)')iday
   else
    print*,'obs.exe: day string length is too long...stop'
    stop
   endif
   if (ihour.lt.10) then
    write(ifile(9:9),'(1I1)')ihour
   elseif (iday.lt.100) then
    write(ifile(8:9),'(1I2)')ihour
   else
    print*,'obs.exe: hour string length is too long...stop'
    stop
   endif
   if (imin.lt.10) then
    write(ifile(12:12),'(1I1)')imin
   elseif (iday.lt.100) then
    write(ifile(11:12),'(1I2)')imin
   else
    print*,'obs.exe: minute string length is too long...stop'
    stop
   endif 
   if (debug.ge.1) print*,'obs.exe: input truth  is: ',ifile(1:20)
   OPEN(71,file=ifile,status='old')
   READ(71,'(6F13.6)')((u(i,j),i=1,nx),j=1,ny)
   READ(71,'(6F13.6)')((v(i,j),i=1,nx),j=1,ny)
   READ(71,'(6F13.6)')((z(i,j),i=1,nx),j=1,ny)
   CLOSE(71)
   WRITE(91,rec=irec)((u(i,j),i=1,nx),j=1,ny)
   irec     = irec + 1
   WRITE(91,rec=irec)((v(i,j),i=1,nx),j=1,ny)
   irec     = irec + 1
   WRITE(91,rec=irec)((z(i,j),i=1,nx),j=1,ny)
   irec     = irec + 1
!
! generate random perturbation now
!
   CALL com_randn(nx*ny,rnd)
   ur       = RESHAPE(rnd,(/nx,ny/))*obs_err_u
   CALL com_randn(nx*ny,rnd)
   vr       = RESHAPE(rnd,(/nx,ny/))*obs_err_v
   CALL com_randn(nx*ny,rnd)
   zr       = RESHAPE(rnd,(/nx,ny/))*obs_err_z
   ur       = u + ur
   vr       = v + vr
   zr       = z + zr
!
! Output the observation
!
   ofile   = 'obs_00:00:00.dat'
   if (iday.lt.10) then
    write(ofile(6:6),'(1I1)')iday
   elseif (iday.lt.100) then
    write(ofile(5:6),'(1I2)')iday
   else
    print*,'obs.exe: day string length is too long...stop'
    stop
   endif
   if (ihour.lt.10) then
    write(ofile(9:9),'(1I1)')ihour
   elseif (iday.lt.100) then
    write(ofile(8:9),'(1I2)')ihour
   else
    print*,'obs.exe: hour string length is too long...stop'
    stop
   endif
   if (imin.lt.10) then
    write(ofile(12:12),'(1I1)')imin
   elseif (iday.lt.100) then
    write(ofile(11:12),'(1I2)')imin
   else
    print*,'obs.exe: minute string length is too long...stop'
    stop
   endif
   IF (debug.eq.0) PRINT*,'obs.exe: output file is:  ',ofile(1:30)
   OPEN(72,file=ofile,status='unknown')
   IF (obs_flag.eq.0) THEN
     PRINT*,'Generating obs at all grid points obs_flag = 0',no
     IF (no.gt.nx*ny) THEN
       PRINT*,'Num of obs points > number of grid points for this opt...stop'
       STOP 
     ENDIF
     DO i     = 1,no
       jd     = i/nx + 1
       id     = mod(i,nx)
       if (id.eq.0) then 
        id    = nx
        jd    = jd - 1
       endif 
       WRITE(72,'(2I5,3E12.4)')id,jd,ur(id,jd),vr(id,jd),zr(id,jd)
     ENDDO
   ELSEIF (obs_flag.eq.1) THEN
     PRINT*,'Generating obs around TC center icen,jcen = ',icen,jcen  
     i = 1
     DO id = icen-nint(sqrt(no*1.0)/2),icen+nint(sqrt(no*1.0)/2)
       DO jd = jcen-nint(sqrt(no*1.0)/2),jcen+nint(sqrt(no*1.0)/2)
         IF (i.le.no.and.id.ge.1.and.id.le.nx.and.jd.ge.1.and.jd.le.ny) THEN
           WRITE(72,'(2I5,3E12.4)')id,jd,ur(id,jd),vr(id,jd),zr(id,jd)
           i = i + 1
         ENDIF
       ENDDO
     ENDDO
     PRINT*,'Actual num of obs points with obs_flag = 1 is: ',i-1,no  
   ELSE
     PRINT*,'obs.exe: not support for other obs_flag yet...stop'
     STOP
   ENDIF
!
! write out obs data
!
   IF (obs_flag.eq.1) THEN
     DO i = 1,nx
       DO j = 1,ny
         IF (i.lt.icen-nint(sqrt(no*1.0)/2).or.i.gt.icen+nint(sqrt(no*1.0)/2).or. &
             j.lt.jcen-nint(sqrt(no*1.0)/2).or.j.gt.jcen+nint(sqrt(no*1.0)/2)) THEN
            ur(i,j) = -99999
            vr(i,j) = -99999
            zr(i,j) = -99999
         ENDIF
       ENDDO
     ENDDO
   ENDIF
   WRITE(91,rec=irec)((ur(i,j),i=1,nx),j=1,ny)
   irec     = irec + 1
   WRITE(91,rec=irec)((vr(i,j),i=1,nx),j=1,ny)
   irec     = irec + 1
   WRITE(91,rec=irec)((zr(i,j),i=1,nx),j=1,ny)
   irec     = irec + 1
   CLOSE(72)
   nprint   = nprint + 1
  else
   nprint   = 0
  endif
  time      = time + dt_window
  itime     = itime + 1
  IF (itime.le.ntime) GOTO 19
  PRINT*,'obs.exe: program ends perfectly'
  END  


  SUBROUTINE input_namelist(debug,tfcst,restart,obs_err_u,obs_err_v,&
                            obs_err_z,no,ne,obs_flag,t_window,dt_window, &
                            icen,jcen)
  INCLUDE "../registry/swe.inc"
  RETURN
  END


  SUBROUTINE com_randn(ndim,var)
  USE mt19937
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ndim
  REAL(8),INTENT(OUT) :: var(1:ndim)
  REAL(8) :: rnd(2),pi
  INTEGER :: idate(8)
  INTEGER :: i,iseed
  LOGICAL,SAVE :: first=.true.
  pi = 4.*atan(1.)

  IF (first) THEN
    CALL DATE_AND_TIME(VALUES=idate)
    iseed = idate(8) + idate(7)*1000
    CALL init_genrand(iseed)
    first=.false.
  END IF

  IF( MOD(ndim,2)==0 ) THEN
    DO i=1,ndim/2
      rnd(1) = genrand_res53()
      rnd(2) = genrand_res53()
      var(i*2-1) = sqrt( -2.0d0 * log( rnd(1) ) ) * sin( 2.0d0*pi*rnd(2) )
      var(i*2) = sqrt( -2.0d0 * log( rnd(1) ) ) * cos( 2.0d0*pi*rnd(2) )
    END DO
  ELSE
    DO i=1,(ndim-1)/2
      rnd(1) = genrand_res53()
      rnd(2) = genrand_res53()
      var(i*2-1) = sqrt( -2.0d0 * log( rnd(1) ) ) * sin( 2.0d0*pi*rnd(2) )
      var(i*2) = sqrt( -2.0d0 * log( rnd(1) ) ) * cos( 2.0d0*pi*rnd(2) )
    END DO
    rnd(1) = genrand_res53()
    rnd(2) = genrand_res53()
    var(ndim) = sqrt( -2.0d0 * log( rnd(1) ) ) * sin( 2.0d0*pi*rnd(2) )
  END IF
  RETURN
END SUBROUTINE com_randn
