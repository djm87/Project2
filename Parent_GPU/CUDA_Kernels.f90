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
  USE Variables
 
  
  REAL(RP), DIMENSION(NSLIP,NC), DEVICE :: SLPRA1,SLPRA2
  REAL(RP), DIMENSION(6,NC), DEVICE :: STRESS2BE
  !REAL(RP), DIMENSION(NC,6,6), DEVICE :: TJINV_1,TJINV_2
  !INTEGER, DEVICE ::LU_pivot(6*NC), LU_info(NC),INV_info(NC)
CONTAINS
  
!! NOTES:
! After program runs look into fast memory usage for main variable
!
!
!
!
attributes(global) SUBROUTINE CPSOLV(C6A,SSYSMAT,KUPD,STRS1PK,DG2, &
					STRV, ISTEP)!, SLPRA1, SLPRA2) 
	USE Variables
	IMPLICIT NONE
    !USE cublas_device
	REAL(RP), DEVICE :: C6A(6,NSLIP,*),SSYSMAT(3,3,NSLIP,*), &
				STRS1PK(6,*), DG2(12,NC), STRV(6,NC)!, SLPRA1(NSLIP,NC), &
				!SLPRA2(NSLIP,NC)
	INTEGER, VALUE :: KUPD,ISTEP
	INTEGER :: tx, bx, zz!, istat
	!type(cublasHandle) :: ch1
	!Local memory just over 83MB at 16384 crystals
	REAL(RP):: C6P6A(6,6,NSLIP), SLPRESA2(NSLIP),DRFCE2(NSLIP), &
				SLIPSIGN(NSLIP),ABSDGMAX,ABSDG(NSLIP),DELTPK(6), & 
				TPKOLD(6), EROR(6),DGDTPKA(NSLIP), ABSDELTPK(6), &
				ERORDELTPK,HB(NSLIP),SLPRESA2TMP(NSLIP),  &
				EORRSLPRES,DGMAXDGDO,EXPN,DMAX,DGDO, &
				DGMAX, TJINV(6,6) 	
	
	INTEGER :: IITR, ICORR,START_TIME,STOP_TIME
        INTEGER, SHARED :: counter

	INTEGER :: I,II,III,j,jj,k,kk
	
	!Shared memory 
	!REAL(RP), SHARED :: 
	
	tx = threadidx%x !will be 1-256
	bx = blockdim%x
	
	zz = (blockidx%x-1) * bx + tx

	IITR = 0 
	DGDO = GDO*DT 
	EXPN= 1.0D0/EXPM 
	DMAX = 0.05D0
		  
	C6P6A = 0.0D0
	
	
    START_TIME = 0
	STOP_TIME = 0
	
	IF(ISTEP == 1) THEN
		SLPRA1(:,zz) = SO
		SLPRA2(:,zz) = SO
		STRESS2BE(:,zz) = 0.0D0
	END IF

	!INITIALISE SLIP RESISTANCE AND GET CONSTANT C6P6A
	DO I = 1,NSLIP
		SLPRESA2(I)= SLPRA1(I,zz) 
		DO II = 1,6 
			C6P6A(II,1,I) = C6A(II,I,zz)*SSYSMAT(1,1,I,zz)            
			C6P6A(II,2,I) = C6A(II,I,zz)*SSYSMAT(2,2,I,zz)
			C6P6A(II,3,I) = C6A(II,I,zz)*SSYSMAT(3,3,I,zz)
			C6P6A(II,4,I) = C6A(II,I,zz)*(SSYSMAT(1,2,I,zz)+ &
				SSYSMAT(2,1,I,zz))
			C6P6A(II,5,I) = C6A(II,I,zz)*(SSYSMAT(2,3,I,zz)+ &
			        SSYSMAT(3,2,I,zz))
			C6P6A(II,6,I) = C6A(II,I,zz)*(SSYSMAT(1,3,I,zz)+ &
										SSYSMAT(3,1,I,zz))
		ENDDO
	ENDDO   
	 
	!START OF SOLUTION PROCEDURE
	111  CONTINUE 

	DRFCE2 = 0.0D0 
	 
	DO JJ = 1,NSLIP
		DRFCE2(JJ) = DRFCE2(JJ)+STRESS2BE(1,zz)*SSYSMAT(1,1,JJ,zz)
		DRFCE2(JJ) = DRFCE2(JJ)+STRESS2BE(2,zz)*SSYSMAT(2,2,JJ,zz)
		DRFCE2(JJ) = DRFCE2(JJ)+STRESS2BE(3,zz)*SSYSMAT(3,3,JJ,zz)
		DRFCE2(JJ) = DRFCE2(JJ)+STRESS2BE(4,zz)*(SSYSMAT(1,2,JJ,zz)+ &
										 SSYSMAT(2,1,JJ,zz))
		DRFCE2(JJ) = DRFCE2(JJ)+STRESS2BE(5,zz)*(SSYSMAT(2,3,JJ,zz)+ &
										 SSYSMAT(3,2,JJ,zz))
		DRFCE2(JJ) = DRFCE2(JJ)+STRESS2BE(6,zz)*(SSYSMAT(1,3,JJ,zz)+ &
										 SSYSMAT(3,1,JJ,zz))

		SLIPSIGN(JJ) = SIGN(1.0D0,DRFCE2(JJ))
	ENDDO
    
	ABSDGMAX = 0.0D0
	DO I = 1,NSLIP          
		ABSDG(I) = 0.0D0
		DG2(I,zz) = 0.0D0
		IF(DRFCE2(I).NE.0.0D0)THEN
			ABSDG(I) = DGDO*ABS(DRFCE2(I)/SLPRESA2(I))**EXPN
			DG2(I,zz) = ABSDG(I)*SLIPSIGN(I)
			ABSDGMAX = MAX(ABSDGMAX,ABSDG(I))
		ENDIF
	ENDDO
    
	IF (ABSDGMAX.GT.DMAX.AND.IITR.NE.0.AND.ICORR.LE.30)THEN
		DO J=1,6
			DELTPK(J) = DELTPK(J)*0.25
			STRESS2BE(J,zz) = TPKOLD(J)-DELTPK(J)
		ENDDO
		ICORR = ICORR + 1
		GO TO 111
	ENDIF
	
	ICORR = 0      

	DO J=1,6 
		TPKOLD(J) = STRESS2BE(J,zz)
	ENDDO
    
	DO J=1,6
		EROR(J) =  STRESS2BE(J,zz) - STRV(J,zz)
	ENDDO
    
	DO I=1,NSLIP
		DO II=1,6
			EROR(II) = EROR(II)+DG2(I,zz)*C6A(II,I,zz)
		ENDDO
	ENDDO

	DO K=1,6
		DO KK=1,6
			TJINV(KK,K) = 0.0D0
		ENDDO
	ENDDO
	TJINV(1,1) = 1.0D0
	TJINV(2,2) = 1.0D0
	TJINV(3,3) = 1.0D0
	TJINV(4,4) = 1.0D0
	TJINV(5,5) = 1.0D0
	TJINV(6,6) = 1.0D0    

	DO I = 1,NSLIP
		IF (DG2(I,zz).NE.0.0D0) THEN
		DGDTPKA(I) = DG2(I,zz)*EXPN/DRFCE2(I)
			DO K = 1,6
				DO KK = 1,6
					TJINV(KK,K) = TJINV(KK,K) +C6P6A(KK,K,I)* DGDTPKA(I)
				ENDDO
			ENDDO
		ENDIF
	ENDDO

	!GET STRESS CORRECTION AND CORRECT PK STRESS
	!IF(zz==1) START_TIME = clock()
	call CPMATINV(TJINV,6,6)    
       ! IF(zz==1) STOP_TIME = clock()
        !IF(zz==1) print *,"time INV=",(STOP_TIME-START_TIME)/705.500
        !IF(zz==1) counter  = counter + 1
       ! IF(zz==1) print *,"Count=",counter	
	!Call syncthreads()

	!IF(zz == (blockidx%x)*bx) THEN  
	!	istat = cublasCreate_v2(ch1)
    !
	!	istat = cublasSgetrfBatched_v2(ch1,6,TJINV_1(1+(blockidx%x-1)*bx:256+(blockidx%x-1)*bx,1,1),6, &
	!			LU_pivot(1+(blockidx%x-1)*bx),LU_info(1+(blockidx%x-1)*bx),bx)
    !
	!	istat = cublasDgetriBatched_v2(ch1,6,TJINV_1(1+(blockidx%x-1)*bx:256+(blockidx%x-1)*bx,1,1),6, &
	!			LU_pivot(1+(blockidx%x-1)*bx),TJINV_2(1+(blockidx%x-1)*bx:256+(blockidx%x-1)*bx,1,1),6, &
	!			INV_info(1+(blockidx%x-1)*bx),bx)
	!	istat = cublasDestroy_v2(ch1)
	!END IF
	
	DO II = 1,6
		DELTPK(II)=TJINV(II,1)*EROR(1)+TJINV(II,2)*EROR(2)+ &
			TJINV(II,3)*EROR(3)+TJINV(II,4)*EROR(4)+ &
			TJINV(II,5)*EROR(5)+TJINV(II,6)*EROR(6)
	ENDDO

	DO JJ = 1,6
		ABSDELTPK(JJ) = ABS(DELTPK(JJ))
		STRESS2BE(JJ,zz)= STRESS2BE(JJ,zz) - DELTPK(JJ)
	ENDDO
	
    
	ERORDELTPK = 0.0D0
	DO J = 1,6
		ERORDELTPK = MAX(ERORDELTPK,ABSDELTPK(J))
	ENDDO

	IITR = IITR + 1
	IF(IITR.EQ.201)RETURN
	IF(ERORDELTPK.GT.SO/10000000.0)GO TO 111
	IF(IITR.LE.3)GO TO 111
	
    !IF(zz ==30 .and. ISTEP == 1) print *,"STRESS2BE(6)11 = ", STRESS2BE(6)
	!IF(zz ==30 .and. ISTEP == 1) print *,"tests"
	
	DO I = 1,NSLIP
		DRFCE2(I) = 0.0D0
	ENDDO

	DO JJ = 1,NSLIP
		DRFCE2(JJ) = DRFCE2(JJ)+STRESS2BE(1,zz)*SSYSMAT(1,1,JJ,zz)
		DRFCE2(JJ) = DRFCE2(JJ)+STRESS2BE(2,zz)*SSYSMAT(2,2,JJ,zz)
		DRFCE2(JJ) = DRFCE2(JJ)+STRESS2BE(3,zz)*SSYSMAT(3,3,JJ,zz)
		DRFCE2(JJ) = DRFCE2(JJ)+STRESS2BE(4,zz)*(SSYSMAT(1,2,JJ,zz)+ &
					 SSYSMAT(2,1,JJ,zz))
		DRFCE2(JJ) = DRFCE2(JJ)+STRESS2BE(5,zz)*(SSYSMAT(2,3,JJ,zz)+ &
					 SSYSMAT(3,2,JJ,zz))
		DRFCE2(JJ) = DRFCE2(JJ)+STRESS2BE(6,zz)*(SSYSMAT(1,3,JJ,zz)+ &
					 SSYSMAT(3,1,JJ,zz))
	ENDDO

	DO I = 1,NSLIP         
	ABSDG(I) = 0.0D0
	DG2(I,zz) = 0.0D0
		IF(DRFCE2(I).NE.0.0D0)THEN
			ABSDG(I) = DGDO*ABS(DRFCE2(I)/SLPRESA2(I))**EXPN
			DG2(I,zz) =SIGN(ABSDG(I),DRFCE2(I))
		ENDIF
	ENDDO
		  

	DO I = 1,NSLIP
		IF (SLPRESA2(I).GT.SS)THEN
			HB(I) = 0.0D0
		ELSE
			HB(I) = HO*(1.0D0 - SLPRESA2(I)/SS)**EXPA
		ENDIF
		SLPRESA2TMP(I) = SLPRA1(I,zz)
	ENDDO

	DO II = 1,NSLIP
		DO I = 1,NSLIP
			SLPRESA2TMP(I) = SLPRESA2TMP(I)+HB(II)*ABSDG(II)
		ENDDO
	ENDDO

	EORRSLPRES = 0.0D0
	DO I = 1,NSLIP
		EORRSLPRES = DMAX1(DABS(SLPRESA2TMP(I)-SLPRESA2(I)),EORRSLPRES)
		SLPRESA2(I) = SLPRESA2TMP(I)
	ENDDO

	IF(EORRSLPRES.GT.SO/10000000.0)GO TO 111
		  
	!END OF SLIP RESISTANCE LOOP

	STRS1PK(1,zz)=STRESS2BE(1,zz)
	STRS1PK(2,zz)=STRESS2BE(2,zz)
	STRS1PK(3,zz)=STRESS2BE(3,zz)
	STRS1PK(4,zz)=STRESS2BE(4,zz)
	STRS1PK(5,zz)=STRESS2BE(5,zz)
	STRS1PK(6,zz)=STRESS2BE(6,zz)

	!IF(zz ==30 .and. ISTEP == 1) print *,"STRS1PK11 = ", STRS1PK(6,zz)
	!IF(zz ==31 .and. ISTEP == 1) print *,"STRS1PK12 = ", STRS1PK(6,zz)
	!IF(zz ==38.and. ISTEP == 1) print *,"STRS1PK13 = ", STRS1PK(6,zz)
	DGMAX = 0.0D0

	IF(KUPD.EQ.1)THEN
		DO I = 1,NSLIP
			SLPRA2(I,zz)=SLPRESA2(I)
			SLPRA1(I,zz)=SLPRESA2(I)
			DGMAX = MAX(DGMAX,ABSDG(I))
		ENDDO
	ENDIF	

  RETURN  
  END SUBROUTINE CPSOLV
  
  attributes(DEVICE) SUBROUTINE CPMATINV(A,N,ID)
    
      IMPLICIT REAL(RP) (A-H,O-Z)
      DIMENSION IPIVOT(1000),IND(1000,2),PIVOT(1000)                        
      DIMENSION A(ID,1),B(500,1)                                     
            
!     SEARCH FOR PIVOT ELEMENT                                          
      M=0
      DO 30 J=1,N                                                       
   30 IPIVOT (J)=0                                                      
      DO 470 I=1,N                                                      
      AMAX=0.0E0                                                        
      DO 150 J=1,N                                                      
      IF(IPIVOT(J).EQ.1) GO TO 150                                      
      DO 140 K=1,N                                                      
      IF(IPIVOT(K).GT.1) GO TO 590                                      
      IF(IPIVOT(K).EQ.1) GO TO 140                                      
      IF(ABS(AMAX).GE.ABS(A(J,K))) GO TO 140                          
      IROW=J                                                            
      ICOLUM=K                                                          
      AMAX=A(J,K)                                                       
  140 CONTINUE                                                          
  150 CONTINUE                                                          
      IPIVOT(ICOLUM)=IPIVOT(ICOLUM)+1                                   
!     INTERCHANGE ROWS TO PUT PIVOT ELEMENT ON DIAGONAL                 
      IF(IROW.EQ.ICOLUM) GO TO 280                                      
!      DETERM=-DETERM                                                    
      DO 220 L=1,N                                                      
      SWAP=A(IROW,L)                                                    
      A(IROW,L)=A(ICOLUM,L)                                             
  220 A(ICOLUM,L)=SWAP                                                  
      M=0
      IF(M.LE.0) GO TO 280                                              
      DO 270 L=1,M                                                      
      SWAP=B(IROW,L)                                                    
      B(IROW,L)=B(ICOLUM,L)                                             
  270 B(ICOLUM,L)=SWAP                                                  
  280 IND(I,1)=IROW                                                   
      IND(I,2)=ICOLUM                                                 
      PIVOT(I)=A(ICOLUM,ICOLUM)                                         
!      DETERM=DETERM*PIVOT(I)                                            
!     DIVIDE PIVOT ROW BY PIVOT ELEMENT                                 
      A(ICOLUM,ICOLUM)=1.0E0                                            
      DO 340 L=1,N                                                      
  340 A(ICOLUM,L)=A(ICOLUM,L)/PIVOT(I)                                  
      IF(M.LE.0) GO TO 380                                              
      DO 370 L=1,M                                                      
  370 B(ICOLUM,L)=B(ICOLUM,L)/PIVOT(I)                                  
!     REDUCE NON-PIVOT ROWS                                             
  380 DO 470 L1=1,N                                                     
      IF(L1.EQ.ICOLUM) GO TO 470                                        
      T=A(L1,ICOLUM)                                                    
      A(L1,ICOLUM)=0.0E0                                                
      DO 430 L=1,N                                                      
  430 A(L1,L)=A(L1,L)-A(ICOLUM,L)*T                                     
      IF(M.LE.0) GO TO 470                                              
      DO 460 L=1,M                                                      
  460 B(L1,L)=B(L1,L)-B(ICOLUM,L)*T                                     
  470 CONTINUE                                                          
!     INTERCHANGE COLUMNS                                               
      DO 580 I=1,N                                                      
      L=N+1-I                                                           
      IF(IND(L,1).EQ.IND(L,2)) GO TO 580                            
      JROW=IND(L,1)                                                   
      JCOLUM=IND(L,2)                                                 
      DO 570 K=1,N                                                      
      SWAP=A(K,JROW)                                                    
      A(K,JROW)=A(K,JCOLUM)                                             
      A(K,JCOLUM)=SWAP                                                  
  570 CONTINUE                                                          
  580 CONTINUE                                                          
  590 RETURN                                                            
    END SUBROUTINE CPMATINV    
  
END MODULE CUDA_Kernels
