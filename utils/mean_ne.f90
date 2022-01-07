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
    INTEGER, PARAMETER                  :: n = 40
    CHARACTER*100                       :: ofile,ifile
    INTEGER                             :: ne
    REAL, ALLOCATABLE, DIMENSION(:)     :: xm
    REAL, ALLOCATABLE, DIMENSION(:,:)   :: x,Bm,Q
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: B
    REAL                                :: t1,t2,t3,t4
    INTEGER                             :: id,ie,i,j,k,debug
    ofile         = 'mean.dat'
    ifile         = 'mem_000.dat'
    CALL input_namelist(debug,ne)
    PRINT*,'mean_ne.exe: number of ensemble ne =',ne
    ALLOCATE(xm(n),x(n,ne),Q(n,n))
    DO id         = 1,ne
     IF (id.lt.10) THEN
      WRITE(ifile(7:7),'(1I1)')id
     ELSEIF (id.lt.100) THEN
      WRITE(ifile(6:7),'(1I2)')id
     ELSE
      WRITE(ifile(5:7),'(1I3)')id
     ENDIF
     IF (debug.eq.1) PRINT*,'mean_ne.exe: Open file is:  ',ifile(1:30)
     OPEN(10,file=ifile,status='old')
     READ(10,*)(x(i,id),i=1,n)
     CLOSE(10)
    ENDDO
!
! compute state vector mean
!
    xm             = 0.
    DO i           = 1,ne
     xm(:)         = xm(:) + x(:,i)
    ENDDO
    xm(:)          = xm(:)/ne
    OPEN(11,file=ofile)
    WRITE(11,*)(xm(i),i=1,n)
    CLOSE(11)
!
! compute covariance matrix for model errors
!
    DO i           = 1,n
     DO j          = 1,n
      Q(i,j)       = 0.
      DO k         = 1,ne
       Q(i,j)      = Q(i,j) + (x(i,k)-xm(i))*(x(j,k)-xm(j))
      ENDDO
      Q(i,j)       = Q(i,j)/ne
     ENDDO
    ENDDO
    END

    SUBROUTINE input_namelist(debug,ne)
    INCLUDE "../L40.inc"
    RETURN
    END


