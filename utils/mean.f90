!
! This program is to compute the ensemble of outputs from the model
! with the perturbed forcings. This can be included into KF codes
! but putting it into another program will ease the same KF code
! for any further implemenetation
!
! Author: Chanh Q. Kieu
!
!==================================================================
    PROGRAM ensemble_mean
    IMPLICIT NONE
    INTEGER                             :: nx,ny,ne
    CHARACTER*100                       :: ofile,ifile
    REAL, ALLOCATABLE, DIMENSION(:,:)   :: um,vm,zm
    REAL, ALLOCATABLE, DIMENSION(:,:)   :: us,vs,zs
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ua,va,za
    REAL                                :: t1,t2,t3,t4
    INTEGER                             :: id,ie,i,j,k,debug
    ofile         = 'mean.dat'
    ifile         = 'mem_000.dat'
    CALL input_namelist(debug,nx,ny,ne)
    PRINT*,'emean.exe: nx dim                =',nx
    PRINT*,'emean.exe: ny dim                =',ny
    PRINT*,'emean.exe: number of ensemble ne =',ne
    ALLOCATE(ua(nx,ny,ne),va(nx,ny,ne),za(nx,ny,ne))
    ALLOCATE(um(nx,ny),vm(nx,ny),zm(nx,ny))
    ALLOCATE(us(nx,ny),vs(nx,ny),zs(nx,ny))
    DO id         = 1,ne
     IF (id.lt.10) THEN
      WRITE(ifile(7:7),'(1I1)')id
     ELSEIF (id.lt.100) THEN
      WRITE(ifile(6:7),'(1I2)')id
     ELSE
      WRITE(ifile(5:7),'(1I3)')id
     ENDIF
     IF (debug.eq.1) PRINT*,'emean.exe: Open file is:  ',ifile(1:30)
     OPEN(10,file=ifile,status='old')
     READ(10,'(6E13.6)')((ua(i,j,id),i=1,nx),j=1,ny)
     READ(10,'(6E13.6)')((va(i,j,id),i=1,nx),j=1,ny)
     READ(10,'(6E13.6)')((za(i,j,id),i=1,nx),j=1,ny)
     CLOSE(10)
    ENDDO
!
! compute state vector mean
!
    um             = 0.
    vm             = 0.
    zm             = 0.
    DO i           = 1,ne
     um(:,:)       = um(:,:) + ua(:,:,i)
     vm(:,:)       = vm(:,:) + va(:,:,i)
     zm(:,:)       = zm(:,:) + za(:,:,i)
    ENDDO
    um             = um/ne
    vm             = vm/ne
    zm             = zm/ne
    OPEN(11,file=ofile)
    WRITE(11,'(6E13.6)')((um(i,j),i=1,nx),j=1,ny)
    WRITE(11,'(6E13.6)')((vm(i,j),i=1,nx),j=1,ny)
    WRITE(11,'(6E13.6)')((zm(i,j),i=1,nx),j=1,ny)
    CLOSE(11)
!
! compute the ensemble spread
!
    DO i           = 1,ne
     us(:,:)       = us(:,:) + (ua(:,:,i)-um(:,:))**2
     vs(:,:)       = vs(:,:) + (va(:,:,i)-vm(:,:))**2
     zs(:,:)       = zs(:,:) + (za(:,:,i)-zm(:,:))**2
    ENDDO
    us             = sqrt(us/ne)
    vs             = sqrt(vs/ne)
    zs             = sqrt(zs/ne)
    OPEN(12,file='spread.dat')
    WRITE(12,'(6E13.6)')((us(i,j),i=1,nx),j=1,ny)
    WRITE(12,'(6E13.6)')((vs(i,j),i=1,nx),j=1,ny)
    WRITE(12,'(6E13.6)')((zs(i,j),i=1,nx),j=1,ny)
    CLOSE(12)
    END

    SUBROUTINE input_namelist(debug,nx,ny,ne)
    INCLUDE "../registry/swe.inc"
    RETURN
    END


