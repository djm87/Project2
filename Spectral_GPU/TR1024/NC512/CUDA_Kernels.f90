!============================================================================================!
! CUDA_Kernels                                                                               !
! Author      : Daniel Savage                                                                !
! Date        :                                                                              !
! Description : CUDA kernels that run the texture operation across multiple streams          !
!               On GPU device.                                                               !
! Information : Matrix multiply: http://www.cise.ufl.edu/~sahni/papers/gpuMatrixMultiply.pdf !
!             : CUDA Fortran: http://www.pgroup.com/doc/pgicudaforug.pdf                     !
!============================================================================================!

MODULE CUDA_Kernels
  USE Variables
  IMPLICIT NONE
  REAL(RP), DIMENSION(3,3,NCRYS), DEVICE :: QMATT, RMATT
  	REAL(RP), DIMENSION(NCRYS,9), DEVICE :: G
CONTAINS
!=================================================================================
!================================================================================= 
  
 attributes(global) SUBROUTINE build_SC(SS1, SS2, SS3, SS4, P1, P2, P3, P4, SC, n, m)
	USE cudadevice
    IMPLICIT NONE
    
	!Global device memory used to initially store date from host before it is accessed by shared memory
    REAL(RP), DIMENSION(n, 2*m), DEVICE :: SC ! SC is of size(m,2*m) because all of the imaginary numbers are stored with an offset of m
    REAL(RP), DIMENSION(m),    DEVICE :: SS1, SS2, SS3, SS4
    REAL(RP), DIMENSION(n),    DEVICE :: P1, P2, P3, P4
   
    INTEGER, VALUE :: n, m
    INTEGER :: tx, ty, bx, by
    INTEGER :: i, j, k

	!The shared memory is declared 
    REAL(RP), DIMENSION(blkSzY), SHARED :: SS1_s, SS2_s, SS3_s, SS4_s
    REAL(RP), DIMENSION(blkSzY), SHARED :: SS12_s, SS22_s, SS32_s, SS42_s
    REAL(RP), DIMENSION(blkSzX), SHARED :: P1_s, P2_s, P3_s, P4_s
    REAL(RP), DIMENSION(blkSzX), SHARED :: P12_s, P22_s, P32_s, P42_s

    REAL(RP) :: mV_re, mV_im, mV 

    !tx, ty, bx, by are the thread indices and block sizes used in the kernel. The block sizes can be changed from the kernel
	!launch. This code is flexible and thread block sizes can be changed in variables by changing blkSzy and blkSzx.
	tx = threadidx%x !will be 1-16
    ty = threadidx%y !will be 1-16
    bx = blockdim%x  !will be 16
    by = blockdim%y  !will be 16
	
	!i and j go from 0 to n and 0 to m respectively
    i  = (blockIdx%x-1) * bx * 2 + tx
    j  = (blockIdx%y-1) * by * 2 + ty

    
    IF (  (i<=n).AND.(j<=m) ) THEN
       
       !The following section reads in 16 elements and an offset of another 16 elements. THis was done
	   !to optimize the amount of work being done by each thread. The simple version of this is only reading
	   !in P1_s and SS1_s multiply and storing them where i  = (blockIdx%x-1) * bx + tx and 
	   ! j  = (blockIdx%y-1) * by + ty.
	   P1_s(tx) = P1(i)
       P2_s(tx) = P2(i)
       P3_s(tx) = P3(i)
       P4_s(tx) = P4(i)

       P12_s(tx) = P1(i+bx)
       P22_s(tx) = P2(i+bx)
       P32_s(tx) = P3(i+bx)
       P42_s(tx) = P4(i+bx)

       SS1_s(ty) = SS1(j)
       SS2_s(ty) = SS2(j)
       SS3_s(ty) = SS3(j)
       SS4_s(ty) = SS4(j)

       SS12_s(ty) = SS1(j+bx)
       SS22_s(ty) = SS2(j+bx)
       SS32_s(ty) = SS3(j+bx)
       SS42_s(ty) = SS4(j+bx)
       
	   ! The next four blocks multiply elements of the supermat with corresponding euler angles. 
       mV  = SS1_s(ty)  * P1_s(tx)  + SS2_s(ty)  * P2_s(tx)  + SS3_s(ty)  * P3_s(tx)  + SS4_s(ty)  * P4_s(tx)
         
		  !Decomposes exponential into sin and cos
          mV_re =  COS(mV) 
          mV_im = -SIN(mV)       
       ! The SC matrix is updated in global device memory from shared in a chunk blkSzX x blkSzY
       SC(i,j)   = mV_re
       SC(i,m+j) = mV_im

       mV  = SS1_s(ty)  * P12_s(tx)  + SS2_s(ty)  * P22_s(tx)  + SS3_s(ty)  * P32_s(tx)  + SS4_s(ty)  * P42_s(tx)
         
          mV_re =  COS(mV)
          mV_im = -SIN(mV)       

       SC(i+bx,j)   = mV_re
       SC(i+bx,m+j) = mV_im

       mV  = SS12_s(ty)  * P1_s(tx)  + SS22_s(ty)  * P2_s(tx)  + SS32_s(ty)  * P3_s(tx)  + SS42_s(ty)  * P4_s(tx)
     
          mV_re =  COS(mV)
          mV_im = -SIN(mV)       

       SC(i,j+bx)   = mV_re
       SC(i,m+j+bx) = mV_im

       mV  = SS12_s(ty)  * P12_s(tx)  + SS22_s(ty)  * P22_s(tx)  + SS32_s(ty)  * P32_s(tx)  + SS42_s(ty)  * P42_s(tx)       
       
          mV_re =  COS(mV)
          mV_im = -SIN(mV)       
       
       SC(i+bx,j+bx)   = mV_re
       SC(i+bx,m+j+bx) = mV_im

    END IF

    RETURN
  END SUBROUTINE build_SC
!=================================================================================
!================================================================================= 
  attributes(global) SUBROUTINE mmul( A, B, C, N, M, L)!, offset )
   USE Variables
   IMPLICIT NONE
    
    REAL(RP), DEVICE :: A(N,M), B(M,L), C(N,L) !These large matrices are declared in device global memory (This can be accessed by host)
    INTEGER, VALUE :: N, M, L!, offset 
    INTEGER :: i, j, k, km
    INTEGER :: tx, ty, bx, by 
     REAL(RP) :: temp1, temp2
	!Sub arrays are located in shared memory which can only be accessed by thread blocks (fast memory access)
	REAL(RP), SHARED :: As1(16,17), As2(16,17), Bs(16,9) 
   
    !tx, ty, bx, by are the thread indices and block sizes used in the kernel. The block sizes can be changed from the kernel
	!launch; however, this code must be rewritten specifically for different block sizes!
    tx = threadidx%x !will be 1-16
    ty = threadidx%y !will be 1-9
    bx = blockdim%x  !will be 16
    by = blockdim%y  !will be 9
	
	!i goes from 0 to N in sections of 32. The reason for not declaring a block size of 32 is that the threads perform better 
    !when a certain level of work is given to them. See below to understand why i is like this. 	
    i = (blockidx%x-1) * 2 * bx + tx! + offset!The reason 
	
	!j goes from 0 to L over and over again and is used to read in elements of the Bs matrix and to fill the C matrix
    j = (blockidx%y-1) * by + ty 
    
    temp1 = 0._RP !Holds the sum of As1*Bs
	temp2 = 0._RP !Holds the sum of As2*Bs
	
	IF ( (i<=N).AND.(j<=L) ) THEN !This is the size of C
		DO k = 1, M, bx !steps from one 32 x 16 block to the next 32 x 16 block until the end of As is reached. This loop does one thread block! 
			IF (ty <= 8) THEN !Since As is reading only a 32 x 16 chunk of A at a time, we don't allow ty=9 to be used in this step
			As1(tx,ty) = A(i, ty + k - 1) !Fills As1 with a 16 x 8 chunk
			As1(tx,ty + 8) = A(i, 8 + ty + k - 1) !Fills As1 with a 16 x 8 chunk making a full As1 of 16 x 16
			
			As2(tx,ty) = A(i + bx, ty + k - 1) !Fills As2 with a 16 x 8 chunk offset from as1 by bx
			As2(tx,ty + 8) = A(i + bx, 8 + ty + k - 1) !Fills As2 with a 16 x 8 chunk making a full As2 of 16 x 16
			
			Bs(tx,ty) = B(k-1 + tx,j) !Fills Bs with a 16 x 8 chunk of B
			ELSE
			Bs(tx,ty) = B(k-1 + tx,j) !Finishes filling Bs with a 16 x 9 chunk of B and uses ty=9 to do this
			END IF
			CALL syncthreads() !This is used to make sure all parts of getting the submatrices are ready to be used in calculations
			
			DO km = 1, 16 !Performs matrix multiply
				temp1 = temp1 + As1(tx,km) * Bs(km,ty)
				temp2 = temp2 + As2(tx,km) * Bs(km,ty)
			ENDDO
			
			CALL syncthreads() 
		ENDDO
		
		C(i,j) = temp1 !Fills the C matrix with a completed value, also transfers data from shared memory to global device memory. 
		C(i + bx,j) = temp2 !Fills the C matrix with a completed value offset by bx. So a 16 x 9 thread block calculates a 32 x 9 chunk of C
		 
	END IF
    RETURN
    
  END SUBROUTINE mmul
!=================================================================================
!=================================================================================
    attributes(global) SUBROUTINE QMAT(Phi2, PHI, Phi1)
	    IMPLICIT NONE
		
		REAL(RP), DIMENSION(NCRYS), DEVICE ::PHI, Phi1,Phi2
		INTEGER :: i
		INTEGER :: tx, bx
		REAL(RP), SHARED :: to_rad
		!Total shared memory used per block is 6*256*8 = 12288 out of 49152
		REAL(RP), DIMENSION(256), SHARED :: SPhi1, CPhi1, SPHI, CPHI, SPhi2, CPhi2 !Try with different shared memory configuration
		tx = threadidx%x !will be 1-256
		bx = blockdim%x  !will be 256
		
		i = (blockidx%x-1) * bx + tx 

		to_rad = (PI / 180.d0)
		
		IF (i<=NCRYS) THEN 

			!change start	
			SPhi1(tx) = SIN(to_rad*(Phi1(i)*360/120/m_fine))
			CPhi1(tx) = COS(to_rad*(Phi1(i)*360/120/m_fine))
			SPHI(tx)  = SIN(to_rad*(PHI(i)*360/120/m_fine))
			CPHI(tx)  = COS(to_rad*(PHI(i)*360/120/m_fine))
			SPhi2(tx) = SIN(to_rad*(Phi2(i)*360/120/m_fine))
			CPhi2(tx) = COS(to_rad*(Phi2(i)*360/120/m_fine))
			!change end
			
			QMATT(1,1,i) = CPhi1(tx)*CPhi2(tx)-SPhi1(tx)*CPHI(tx)*SPhi2(tx)
			QMATT(1,2,i) = -CPhi1(tx)*SPhi2(tx)-SPhi1(tx)*CPHI(tx)*CPhi2(tx)
			QMATT(1,3,i) = SPhi1(tx)*SPHI(tx)
			QMATT(2,1,i) = SPhi1(tx)*CPhi2(tx)+CPhi1(tx)*CPHI(tx)*SPhi2(tx)
			QMATT(2,2,i) = -SPhi1(tx)*SPhi2(tx)+CPhi1(tx)*CPHI(tx)*CPhi2(tx)
			QMATT(2,3,i) = -CPhi1(tx)*SPHI(tx)
			QMATT(3,1,i) = SPHI(tx)*SPhi2(tx)
			QMATT(3,2,i) = SPHI(tx)*CPhi2(tx)
			QMATT(3,3,i) = CPHI(tx)    
			
						
		ENDIF
		RETURN
			 
    END SUBROUTINE QMAT
  
  !=================================================================================
  !=================================================================================
    attributes(global) SUBROUTINE G_flag0(Q_p_sam, PHI, Phi1, Phi2, &
        			QMAT_11, QMAT_12, QMAT_13, QMAT_21, QMAT_22, QMAT_23, &
					QMAT_31, QMAT_32, QMAT_33)
    IMPLICIT NONE
   
	REAL(RP), DIMENSION(9), DEVICE :: Q_p_sam
	REAL(RP), DIMENSION(NCRYS), DEVICE ::PHI, Phi1,Phi2
	REAL(RP), DIMENSION(NCRYS), DEVICE :: QMAT_11, QMAT_12, QMAT_13, QMAT_21, &
					QMAT_22, QMAT_23, QMAT_31, QMAT_32, QMAT_33
    INTEGER :: i
    INTEGER :: tx, bx
	!Total shared memory used per block is 9*256*8 + 9*8 + 256*9*8 +256*3*8= 43080 out of 49152
    REAL(RP), SHARED :: Gsub(256,9), sQ_p_sam(9), angles(256,3), to_deg, twoPI
	REAL(RP), DIMENSION(256), SHARED :: QMATs_11, QMATs_12, QMATs_13, QMATs_21, QMATs_22, &
										QMATs_23, QMATs_31, QMATs_32, QMATs_33
   
	tx = threadidx%x !will be 1-256
    bx = blockdim%x  !will be 256
    
    i = (blockidx%x-1) * bx + tx 

	sQ_p_sam = Q_p_sam
	!change start
	to_deg = 1._rp /(PI / 180._rp)
	!change end
	twoPI = 2._rp*PI
	 angles = 0._rp
	
	IF(i<=NCRYS) THEN
		QMATs_11(tx) = QMAT_11(i) 
		QMATs_12(tx) = QMAT_12(i)
		QMATs_13(tx) = QMAT_13(i)
		QMATs_21(tx) = QMAT_21(i)
		QMATs_22(tx) = QMAT_22(i)
		QMATs_23(tx) = QMAT_23(i)
		QMATs_31(tx) = QMAT_31(i)
		QMATs_32(tx) = QMAT_32(i)
		QMATs_33(tx) = QMAT_33(i)
	
		Gsub(tx,1) = sQ_p_sam(1)*QMATs_11(tx) + sQ_p_sam(2)*QMATs_21(tx) + sQ_p_sam(3)*QMATs_31(tx)
		Gsub(tx,2) = sQ_p_sam(4)*QMATs_11(tx) + sQ_p_sam(5)*QMATs_21(tx) + sQ_p_sam(6)*QMATs_31(tx)
		Gsub(tx,3) = sQ_p_sam(7)*QMATs_11(tx) + sQ_p_sam(8)*QMATs_21(tx) + sQ_p_sam(9)*QMATs_31(tx)
		Gsub(tx,4) = sQ_p_sam(1)*QMATs_12(tx) + sQ_p_sam(2)*QMATs_22(tx) + sQ_p_sam(3)*QMATs_32(tx)
		Gsub(tx,5) = sQ_p_sam(4)*QMATs_12(tx) + sQ_p_sam(5)*QMATs_22(tx) + sQ_p_sam(6)*QMATs_32(tx)
		Gsub(tx,6) = sQ_p_sam(7)*QMATs_12(tx) + sQ_p_sam(8)*QMATs_22(tx) + sQ_p_sam(9)*QMATs_32(tx)
		Gsub(tx,7) = sQ_p_sam(1)*QMATs_13(tx) + sQ_p_sam(2)*QMATs_23(tx) + sQ_p_sam(3)*QMATs_33(tx)
		Gsub(tx,8) = sQ_p_sam(4)*QMATs_13(tx) + sQ_p_sam(5)*QMATs_23(tx) + sQ_p_sam(6)*QMATs_33(tx)
		Gsub(tx,9) = sQ_p_sam(7)*QMATs_13(tx) + sQ_p_sam(8)*QMATs_23(tx) + sQ_p_sam(9)*QMATs_33(tx)		 

		IF (abs(Gsub(tx,9)) >0.9999) THEN
			angles(tx,1) = acos(Gsub(tx,1))
			angles(tx,2) = ACOS(Gsub(tx,9))
			angles(tx,3) = 0._rp

			IF (Gsub(tx,2) < 0) THEN
			   angles(tx,1) = twoPI - angles(tx,1)
			END IF
		ELSE
			angles(tx,1) = ATAN2(Gsub(tx,7),-1*Gsub(tx,8))
			angles(tx,2) = ACOS(Gsub(tx,9))
			angles(tx,3) = ATAN2(Gsub(tx,3),Gsub(tx,6))
		ENDIF

		IF (angles(tx,1) < 0) THEN
			angles(tx,1) = angles(tx,1)+ twoPI
		ENDIF
		 
		IF (angles(tx,2) < 0) THEN
			angles(tx,2) = angles(tx,2)+ twoPI
		ENDIF
		 
		IF (angles(tx,3) < 0) THEN
			angles(tx,3) = angles(tx,3)+ twoPI
		ENDIF
       !change, start,  transfering of euler angles in radians into r,s,t,q
		angles(tx,:) = NINT( to_deg * angles(tx,:) * 120 * m_fine / 360 )
	   !change, end   

		Phi1(i) = angles(tx,1)
		PHI(i)  = angles(tx,2)
		Phi2(i) = angles(tx,3)
		
	ENDIF
	RETURN	
  END SUBROUTINE G_flag0 
  
  !=================================================================================
  !=================================================================================
    attributes(global) SUBROUTINE RMAT(StressMat,W_app_pr, fact_DT, DT)
    IMPLICIT NONE
	REAL(RP), DEVICE :: StressMat(NCRYS,9), W_app_pr(3,3)
	REAL(RP), VALUE :: DT, fact_DT	
    INTEGER :: i
    INTEGER :: tx, bx
    REAL(RP), DIMENSION(256), SHARED :: ang, axis_1, axis_2, axis_3, w21_recs, w31_recs, w32_recs 
	REAL(RP), SHARED :: W_app_pr21, W_app_pr31, W_app_pr32 
    tx = threadidx%x !will be 1-256
    bx = blockdim%x  !will be 256
    
    i = (blockidx%x-1) * bx + tx 
  
	W_app_pr21 = W_app_pr(2,1)
	W_app_pr31 = W_app_pr(3,1)
	W_app_pr32 = W_app_pr(3,2)
	
	IF (i<=NCRYS) THEN
		
		w21_recs(tx) = StressMat(i,6)*fact_DT + W_app_pr21*DT
		w31_recs(tx) = StressMat(i,7)*fact_DT + W_app_pr31*DT
		w32_recs(tx) = StressMat(i,8)*fact_DT + W_app_pr32*DT
		
        
		
		ang(tx) = SQRT(w21_recs(tx)**2 + w31_recs(tx)**2 + w32_recs(tx)**2)
	
		IF (ang(tx) == 0) THEN
			axis_1(tx) = 1._rp 
			axis_2(tx) = 0._rp
			axis_3(tx) = 0._rp
		ELSE
			axis_1(tx) =  w32_recs(tx)/ang(tx)
			axis_2(tx) = -w31_recs(tx)/ang(tx)
			axis_3(tx) =  w21_recs(tx)/ang(tx)
		ENDIF

	
		RMATT(1,1,i) = (1._rp - axis_1(tx)**2)*COS(ang(tx)) + axis_1(tx)**2
		RMATT(1,2,i) = axis_1(tx)*axis_2(tx)*(1._rp - COS(ang(tx))) + axis_3(tx)*SIN(ang(tx))
		RMATT(1,3,i) = axis_1(tx)*axis_3(tx)*(1._rp - COS(ang(tx))) - axis_2(tx)*SIN(ang(tx))
		RMATT(2,1,i) = axis_1(tx)*axis_2(tx)*(1._rp - COS(ang(tx))) - axis_3(tx)*SIN(ang(tx))
		RMATT(2,2,i) = (1._rp - axis_2(tx)**2)*COS(ang(tx)) + axis_2(tx)**2
		RMATT(2,3,i) = axis_2(tx)*axis_3(tx)*(1._rp - COS(ang(tx))) + axis_1(tx)*SIN(ang(tx))
		RMATT(3,1,i) = axis_1(tx)*axis_3(tx)*(1._rp - COS(ang(tx))) + axis_2(tx)*SIN(ang(tx))
		RMATT(3,2,i) = axis_2(tx)*axis_3(tx)*(1._rp - COS(ang(tx))) - axis_1(tx)*SIN(ang(tx))
		RMATT(3,3,i) = (1._rp - axis_3(tx)**2)*COS(ang(tx)) + axis_3(tx)**2
		
    ENDIF
    RETURN    

  END SUBROUTINE RMAT
  
   
  !=================================================================================
  !=================================================================================
    attributes(global) SUBROUTINE G_flag(Q_p_sam,Phi2, PHI, Phi1, flag, n_steps)
    IMPLICIT NONE
	
	INTEGER, VALUE, INTENT(IN)  :: flag, n_steps
	REAL(RP), DIMENSION(9), DEVICE :: Q_p_sam

	REAL(RP), DIMENSION(NCRYS),   DEVICE :: PHI, Phi1,Phi2
    INTEGER :: i
    INTEGER :: tx, bx
	!Total shared memory used is 9*256*8 + 9*8 + 256*9*8 +256*3*8 + 9*8= 45128 out of 49152
    REAL(RP), SHARED :: sQ_p_sam(9), angles(256,3), GG(9), twoPI
    REAL(RP), DIMENSION(256), SHARED :: QMATs_11, QMATs_12, QMATs_13, QMATs_21, &
					QMATs_22, QMATs_23, QMATs_31, QMATs_32, QMATs_33
    REAL(RP), DIMENSION(256), SHARED :: Rs_11, Rs_12, Rs_13, Rs_21, Rs_22, Rs_23, &
					Rs_31, Rs_32, Rs_33
	REAL(RP) :: to_rad, to_deg
	tx = threadidx%x !will be 1-256
    bx = blockdim%x  !will be 256
  
    i = (blockidx%x-1) * bx + tx 
    
	sQ_p_sam = Q_p_sam
	to_rad = 1.e-03_rp * (PI / 180._rp)

    to_deg = 1._rp / to_rad
	twoPI = 2._rp*PI

	 
	IF (i<=NCRYS) THEN
	
		QMATs_11(tx) = QMATT(1,1,i) 
		QMATs_12(tx) = QMATT(1,2,i)
		QMATs_13(tx) = QMATT(1,3,i)
		QMATs_21(tx) = QMATT(2,1,i)
		QMATs_22(tx) = QMATT(2,2,i)
		QMATs_23(tx) = QMATT(2,3,i)
		QMATs_31(tx) = QMATT(3,1,i)
		QMATs_32(tx) = QMATT(3,2,i)
		QMATs_33(tx) = QMATT(3,3,i)
		
		Rs_11(tx) = RMATT(1,1,i)
		Rs_12(tx) = RMATT(1,2,i)
		Rs_13(tx) = RMATT(1,3,i)
		Rs_21(tx) = RMATT(2,1,i)
		Rs_22(tx) = RMATT(2,2,i)
		Rs_23(tx) = RMATT(2,3,i)
		Rs_31(tx) = RMATT(3,1,i) 
		Rs_32(tx) = RMATT(3,2,i)
		Rs_33(tx) = RMATT(3,3,i)
        
		G(i,1) = Rs_11(tx)*QMATs_11(tx) + Rs_21(tx)*QMATs_21(tx) + Rs_31(tx)*QMATs_31(tx)
		G(i,2) = Rs_12(tx)*QMATs_11(tx) + Rs_22(tx)*QMATs_21(tx) + Rs_32(tx)*QMATs_31(tx)
		G(i,3) = Rs_13(tx)*QMATs_11(tx) + Rs_23(tx)*QMATs_21(tx) + Rs_33(tx)*QMATs_31(tx)
		G(i,4) = Rs_11(tx)*QMATs_12(tx) + Rs_21(tx)*QMATs_22(tx) + Rs_31(tx)*QMATs_32(tx)
		G(i,5) = Rs_12(tx)*QMATs_12(tx) + Rs_22(tx)*QMATs_22(tx) + Rs_32(tx)*QMATs_32(tx)
		G(i,6) = Rs_13(tx)*QMATs_12(tx) + Rs_23(tx)*QMATs_22(tx) + Rs_33(tx)*QMATs_32(tx)
		G(i,7) = Rs_11(tx)*QMATs_13(tx) + Rs_21(tx)*QMATs_23(tx) + Rs_31(tx)*QMATs_33(tx)
		G(i,8) = Rs_12(tx)*QMATs_13(tx) + Rs_22(tx)*QMATs_23(tx) + Rs_32(tx)*QMATs_33(tx)
		G(i,9) = Rs_13(tx)*QMATs_13(tx) + Rs_23(tx)*QMATs_23(tx) + Rs_33(tx)*QMATs_33(tx)		 
		
		IF (flag == n_steps) THEN
		
				GG(:) = G(i,:)
				G(i,1) = sQ_p_sam(1)*GG(1) + sQ_p_sam(4)*GG(2) + sQ_p_sam(7)*GG(3)
				G(i,2) = sQ_p_sam(2)*GG(1) + sQ_p_sam(5)*GG(2) + sQ_p_sam(8)*GG(3)
				G(i,3) = sQ_p_sam(3)*GG(1) + sQ_p_sam(6)*GG(2) + sQ_p_sam(9)*GG(3)
				G(i,4) = sQ_p_sam(1)*GG(4) + sQ_p_sam(4)*GG(5) + sQ_p_sam(7)*GG(6) 
				G(i,5) = sQ_p_sam(2)*GG(4) + sQ_p_sam(5)*GG(5) + sQ_p_sam(8)*GG(6)
				G(i,6) = sQ_p_sam(3)*GG(4) + sQ_p_sam(6)*GG(5) + sQ_p_sam(9)*GG(6)
				G(i,7) = sQ_p_sam(1)*GG(7) + sQ_p_sam(4)*GG(8) + sQ_p_sam(7)*GG(9)
				G(i,8) = sQ_p_sam(2)*GG(7) + sQ_p_sam(5)*GG(8) + sQ_p_sam(8)*GG(9)
				G(i,9) = sQ_p_sam(3)*GG(7) + sQ_p_sam(6)*GG(8) + sQ_p_sam(9)*GG(9)

		ENDIF
		
		IF (abs(G(i,9)) >0.9999) THEN
			angles(tx,1) = ACOS(G(i,1))
			angles(tx,2) = ACOS(G(i,9))
			angles(tx,3) = 0._rp

			IF (G(i,2) < 0) THEN
			   angles(tx,1) = twoPI - angles(tx,1)
			END IF
		ELSE
			angles(tx,1) = ATAN2(G(i,7),-1*G(i,8))
			angles(tx,2) = ACOS(G(i,9))
			angles(tx,3) = ATAN2(G(i,3),G(i,6))
		ENDIF

		IF (angles(tx,1) < 0) THEN
			angles(tx,1) = angles(tx,1)+ twoPI
		ENDIF
		 
		IF (angles(tx,2) < 0) THEN
			angles(tx,2) = angles(tx,2)+ twoPI
		ENDIF
		 
		IF (angles(tx,3) < 0) THEN
			angles(tx,3) = angles(tx,3)+ twoPI
		ENDIF
		
		angles(tx,:) = NINT( to_deg * angles(tx,:) )
		
		Phi1(i) = angles(tx,1)
		PHI(i)  = angles(tx,2)
		Phi2(i) = angles(tx,3)
		
    ENDIF	
    RETURN
  END SUBROUTINE G_flag 

END MODULE CUDA_Kernels 
