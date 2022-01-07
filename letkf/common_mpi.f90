MODULE common_mpi
!=======================================================================
!
! NOTE:    MPI interface for indexing
!
! HISTORY: 08 Sep 2011: created by CK
!         
!
!=======================================================================
CONTAINS

      SUBROUTINE para_range(n1, n2, nprocs, irank, ista, iend)
      IMPLICIT none
      INTEGER      :: n1, n2, nprocs, irank, ista, iend
      INTEGER      :: iwork1,iwork2
      iwork1 = (n2 - n1 + 1) / nprocs
      iwork2 = MOD(n2 - n1 + 1, nprocs)
      ista = irank * iwork1 + n1 + MIN(irank, iwork2)
      iend = ista + iwork1 - 1
      IF (iwork2 > irank) iend = iend + 1
      END SUBROUTINE para_range

      SUBROUTINE para_type_block2a(imin, imax, ilen, jlen, ioldtype, inewtype)
      IMPLICIT none
      INCLUDE 'mpif.h'
      INTEGER       :: imin, imax, ilen, jlen, ioldtype, inewtype,ierr
      CALL MPI_TYPE_VECTOR(jlen, ilen, imax - imin + 1,ioldtype, inewtype, ierr)
      CALL MPI_TYPE_COMMIT(inewtype, ierr)
      END SUBROUTINE para_type_block2a

      SUBROUTINE para_type_block2(imin, imax, jmin, ista, iend, jsta, jend, ioldtype, inewtype)
      IMPLICIT none
      INCLUDE 'mpif.h'
      INTEGER       :: imin, imax, jmin, ista, iend, jsta, jend, ioldtype, inewtype
      INTEGER       :: iblock(2), idisp(2), itype(2), ierr,isize,ilen,jlen,itemp
      CALL MPI_TYPE_EXTENT(ioldtype, isize, ierr)
      ilen = iend - ista + 1
      jlen = jend - jsta + 1
      CALL MPI_TYPE_VECTOR(jlen, ilen, imax - imin + 1, ioldtype, itemp, ierr)
      iblock(1) = 1
      iblock(2) = 1
      idisp(1) = 0
      idisp(2) = ((imax-imin+1) * (jsta-jmin) + (ista-imin)) * isize
      itype(1) = MPI_LB
      itype(2) = itemp
      CALL MPI_TYPE_STRUCT(2, iblock, idisp, itype, inewtype, ierr)
      CALL MPI_TYPE_COMMIT(inewtype, ierr)
      END SUBROUTINE para_type_block2

END MODULE common_mpi
