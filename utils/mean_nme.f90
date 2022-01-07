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
    INTEGER                             :: nme
    REAL, ALLOCATABLE, DIMENSION(:)     :: xm
    REAL, ALLOCATABLE, DIMENSION(:,:)   :: x,Bm,Q
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: B
    REAL                                :: t1,t2,t3,t4
    INTEGER                             :: id,ie,i,j,k,debug
    ofile         = 'bgd.dat'
    ifile         = 'bgd_000.dat'
    CALL input_namelist(debug,nme)
    PRINT*,'emean.exe: number of ensemble nme =',nme
    ALLOCATE(xm(n),x(n,nme),B(n,n,nme),Bm(n,n),Q(n,n))
    DO id         = 1,nme
     IF (id.lt.10) THEN
      WRITE(ifile(7:7),'(1I1)')id
     ELSEIF (id.lt.100) THEN
      WRITE(ifile(6:7),'(1I2)')id
     ELSE
      WRITE(ifile(5:7),'(1I3)')id
     ENDIF
     IF (debug.eq.1) PRINT*,'emean.exe: Open file is:  ',ifile(1:30)
     OPEN(10,file=ifile,status='old')
     READ(10,*)(x(i,id),i=1,n)
     CLOSE(10)
    ENDDO
!
! compute state vector mean
!
    xm             = 0.
    DO i           = 1,nme
     xm(:)         = xm(:) + x(:,i)
    ENDDO
    xm(:)          = xm(:)/nme
    OPEN(11,file=ofile)
    WRITE(11,*)(xm(i),i=1,n)
    CLOSE(11)
!
! compute covariance matrix for model errors
!
    DO i           = 1,n
     DO j          = 1,n
      Q(i,j)       = 0.
      DO k         = 1,nme
       Q(i,j)      = Q(i,j) + (x(i,k)-xm(i))*(x(j,k)-xm(j))
      ENDDO
      Q(i,j)       = Q(i,j)/nme
     ENDDO
    ENDDO
    ofile          = 'mmatrix.dat'
    OPEN(11,file=ofile)
    DO i           = 1,n
     WRITE(11,*)(Q(i,j),j=1,n)
    ENDDO
    CLOSE(11)
!
! compute mean bakcground covariance matrix
!
    ofile         = 'tlmodel.dat'
    ifile         = 'tlmodel_000.dat'
    DO id         = 1,nme
     IF (id.lt.10) THEN
      WRITE(ifile(11:11),'(1I1)')id
     ELSEIF (id.lt.100) THEN
      WRITE(ifile(10:11),'(1I2)')id
     ELSE
      WRITE(ifile(9:11),'(1I3)')id
     ENDIF
     IF (debug.eq.1) PRINT*,'emean.exe: Open background file is:  ',ifile(1:30)
     OPEN(10,file=ifile,status='old')
     DO i          = 1,n
      READ(10,*)(B(i,j,id),j=1,n)
     ENDDO
     CLOSE(10)
    ENDDO
    Bm             = 0.
    DO i           = 1,nme
     Bm(:,:)       = Bm(:,:) + B(:,:,i)
    ENDDO
    Bm(:,:)        = Bm(:,:)/nme
    OPEN(11,file=ofile)
    DO i           = 1,n
     WRITE(11,*)(Bm(i,j),j=1,n)
    ENDDO
    CLOSE(11)
    END

    SUBROUTINE input_namelist(debug,nme)
    INCLUDE "../L40.inc"
    RETURN
    END


