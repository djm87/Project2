!============================================================================================!
! CUDA_Kernels                                                                               !
! Author      : Daniel Savage                                                                !
! Date        : 6/30/2013                                                                    !
! Description : CUDA kernels that build the exponential matrix and  multiple it by           !
!               by the supervector                                                           !
! Information : Matrix multiply: http://www.cise.ufl.edu/~sahni/papers/gpuMatrixMultiply.pdf !
!             : CUDA Fortran: http://www.pgroup.com/doc/pgicudaforug.pdf                     !
!============================================================================================!
MODULE CUDA_Kernels
  IMPLICIT NONE

CONTAINS
  
!! NOTES:
! Do ISOLER and IIPT actually do anything? 
!
!
!
!
attributes(global) SUBROUTINE CPSOLV()
	USE Variables




    RETURN
    
  END SUBROUTINE CPSOLV
  
END MODULE CUDA_Kernels
