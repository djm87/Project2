
    PROGRAM POLYPL
	  
	  USE CUDA_Kernels
	  USE Variables
	  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
      

      REAL*8 PHI1(NC), PHI(NC), PHI2(NC)
 
      CHARACTER FILENAME*200
      CHARACTER FILENAME2*200

     ! call cpu_time(toti)
	
      FILENAME2 = 'texture.txt'
      FILENAME = 'stress.txt'
	  
      OPEN(5, FILE='euler.inp')
    


      DO IC = 1,NC
          READ(5,*)PHI1(IC) ,PHI(IC),PHI2(IC)
      ENDDO
      CLOSE(5) 

      

      CALL POLYPLMDL(FILENAME,FILENAME2,IPOINT,NSTEP,IDEF, &
                       DT,TEPS,DEFGRDXX,DEFGRDXY,DEFGRDXZ, &
                       DEFGRDYX,DEFGRDYY,DEFGRDYZ,DEFGRDZX, &
                       DEFGRDZY,DEFGRDZZ,PHI1,PHI,PHI2,NC)



      
      END PROGRAM POLYPL


      SUBROUTINE POLYPLMDL(FILENAME,FILENAME2,IPOINT,NSTEP,IDEF, &
                       DT,TEPS,DEFGRDXX,DEFGRDXY,DEFGRDXZ, &
                       DEFGRDYX,DEFGRDYY,DEFGRDYZ,DEFGRDZX, &
                       DEFGRDZY,DEFGRDZZ,PHI1,PHI,PHI2,NC)
      USE Variables
	  USE cudafor
	  USE CUDA_Kernels
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
      
      DIMENSION STRS2FEM(6),ESV(6,NC),B6A(6,NSLIP,NC),&
               TEMP3(3,3,NSLIP,NC),TEMP2(3,3),TEMP(3,3),FFFF(3,3,NC)
   
      DIMENSION PHI1(*), PHI(*), PHI2(*)
     
      REAL*8,DIMENSION(:,:,:),ALLOCATABLE:: DEFGRDN,DEFGRDC,QC2S
      REAL*8,DIMENSION(:,:,:),ALLOCATABLE:: DEFPI1,DEFPI2
      REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE, PINNED:: SSYSMAT
	  REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE, DEVICE ::SSYSMAT_d
      REAL*8,DIMENSION(:,:,:),ALLOCATABLE::ELSTIF      
      REAL*8,DIMENSION(:,:),ALLOCATABLE, PINNED ::STRS1PK
	  REAL*8,DIMENSION(:,:),ALLOCATABLE, DEVICE ::STRS1PK_d
	  REAL*8,DIMENSION(:,:),ALLOCATABLE, PINNED ::DG2
	  REAL*8,DIMENSION(:,:),ALLOCATABLE, DEVICE ::DG2_d
	  REAL*8,DIMENSION(:,:,:),ALLOCATABLE, PINNED ::C6A
	  REAL*8,DIMENSION(:,:,:),ALLOCATABLE, DEVICE ::C6A_d 
	  REAL*8,DIMENSION(:,:),ALLOCATABLE, DEVICE ::STRV_d
	  REAL*8,DIMENSION(:,:),ALLOCATABLE, PINNED ::STRV
	  !REAL*8,DIMENSION(:,:),ALLOCATABLE, DEVICE ::SLPRA1_d, SLPRA2_d
	  !REAL*8,DIMENSION(:,:),ALLOCATABLE, PINNED ::SLPRA1,SLPRA2
      REAL*8 SCAUCHYAVG(6)       
	  REAL :: ti, tf, cpu_t(3), toti, totf, cuda_t(6)
	  INTEGER(KIND=CUDA_STREAM_KIND) :: stream(2) ! This initializes CUDA stream. 
	  INTEGER :: iostat
	  INTEGER :: imode
	  
	  CHARACTER FILENAME*200
      CHARACTER FILENAME2*200
	  
	  TYPE(dim3) :: grid_CPSOLV, block_CPSOLV
      TYPE(cudaEvent) :: beg, fin
      TYPE (cudaDeviceProp) :: prop
	  
	  iostat = cudaSetDevice(0)
	  iostat = cudaGetDevice(j)
	  iostat = cudaGetDeviceProperties(prop, j)
	  write(*,"(' Device: ', a,/)") trim(prop%name)
	  iostat = cudaRuntimeGetVersion(j)
	  write(*,"(' CUDA Ver: ', I4,/)") j

	  cpu_t  = 0.0D0
	  cuda_t = 0.0D0

	  iostat = cudaEventCreate(beg)
	  iostat = cudaEventCreate(fin)
	  
      OPEN(7, FILE=FILENAME2) 
      OPEN(9, FILE=FILENAME)                 

      ALLOCATE (SSYSMAT(3,3,12,NC),SSYSMAT_d(3,3,12,NC)) 
      ALLOCATE (ELSTIF(6,6,NC))
      ALLOCATE (DEFPI1(3,3,NC))  
      ALLOCATE (DEFPI2(3,3,NC))  
      ALLOCATE (DG2(NSLIP,NC),DG2_d(NSLIP,NC))
      ALLOCATE (DEFGRDN(3,3,NC))    
      ALLOCATE (DEFGRDC(3,3,NC))       
      ALLOCATE (QC2S(3,3,NC))        
      ALLOCATE (STRS1PK(6,NC),STRS1PK_d(6,NC),STRV(6,NC),STRV_d(6,NC))
      ALLOCATE (C6A(6,NSLIP,NC),C6A_d(6,NSLIP,NC))
      !ALLOCATE (SLPRA1(NSLIP,NC),SLPRA1_d(NSLIP,NC))
	  !ALLOCATE (SLPRA2(NSLIP,NC),SLPRA2_d(NSLIP,NC))
      PI = 4.0D0*DATAN(1.0D0)
      DO IC = 1,NC
!        PHI1(IC) = PHI1(IC)*PI/180.0D0 
!        PHI(IC) = PHI(IC)*PI/180.0D0 
!        PHI2(IC) = PHI2(IC)*PI/180.0D0 
      ENDDO

     
      !Does all the updating happening in CPSETUP?
      CALL CPSETUP(SSYSMAT,ELSTIF,DEFPI1,DEFPI2, &
                        DEFGRDN,DEFGRDC,NC, &
                        PHI1,PHI,PHI2,QC2S) 
      
      TEXOUT = DT*NSTEP
      TAU = 0.0D0
      
      call cpu_time(toti)


      do ISTEP = 1,NSTEP   

      TAU = TAU + DT 


      EDOT = 0.001D0
      
      if (idef.eq.1) then
      
        QFAC = 0.0D0
    
        DEFGRDXX = DEXP(EDOT*TAU)
        DEFGRDXY = 0.0D0
        DEFGRDXZ = 0.0D0
        DEFGRDYX = 0.0D0
        DEFGRDYY = DEXP(-QFAC*EDOT*TAU)
        DEFGRDYZ = 0.0D0
        DEFGRDZX = 0.0D0
        DEFGRDZY = 0.0D0
        DEFGRDZZ = DEXP(-(1-QFAC)*EDOT*TAU)  

      elseif(idef.eq.2)then
        DEFGRDXX = 0.0D0
        DEFGRDXY = EDOT*TAU
        DEFGRDXY = 0.0D0
        DEFGRDXZ = 0.0D0
        DEFGRDYX = 0.0D0
        DEFGRDYY = 0.0D0
        DEFGRDYZ = 0.0D0
        DEFGRDZX = 0.0D0
        DEFGRDZY = 0.0D0
      endif 

      EQVST = EDOT*DSQRT((2*(1+QFAC*QFAC+(1-QFAC)*(1-QFAC)))/3)
      STRAINLF = DABS(DLOG(DEFGRDXX))

      DEF2(1,1) = DEFGRDXX
      DEF2(2,1) = DEFGRDYX
      DEF2(3,1) = DEFGRDZX            
      DEF2(1,2) = DEFGRDXY
      DEF2(2,2) = DEFGRDYY
      DEF2(3,2) = DEFGRDZY           
      DEF2(1,3) = DEFGRDXZ
      DEF2(2,3) = DEFGRDYZ
      DEF2(3,3) = DEFGRDZZ  

      SCAUCHYAVG(1) = 0.0D0
      SCAUCHYAVG(2) = 0.0D0
      SCAUCHYAVG(3) = 0.0D0
      SCAUCHYAVG(4) = 0.0D0
      SCAUCHYAVG(5) = 0.0D0
      SCAUCHYAVG(6) = 0.0D0
	  
	call cpu_time(ti)  
    DO IC = 1, NC
	!IF(IC ==1 .and. ISTEP == 2) print *,"DEFPI1 = ", DEFPI1(2,1,IC)
		CALL MATMATMLT(DEF2(1,1), DEF2(1,1), TEMP(1,1), 2) 

		CALL MATMATMLT(TEMP(1,1), DEFPI1(1,1,IC), TEMP2(1,1), 1) 

		CALL MATMATMLT(DEFPI1(1,1,IC), TEMP2(1,1), FFFF(1,1,IC), 2)
		
		ESV(1,IC)=FFFF(1,1,IC)-1.0D0
		ESV(2,IC)=FFFF(2,2,IC)-1.0D0
		ESV(3,IC)=FFFF(3,3,IC)-1.0D0
		ESV(4,IC)=FFFF(1,2,IC)*2.0D0
		ESV(5,IC)=FFFF(2,3,IC)*2.0D0
		ESV(6,IC)=FFFF(1,3,IC)*2.0D0
		
		STRV(1,IC) = 0.5D0*(ELSTIF(1,1,IC)*ESV(1,IC) + ELSTIF(1,2,IC)*ESV(2,IC) + &
					  ELSTIF(1,3,IC)*ESV(3,IC) + ELSTIF(1,4,IC)*ESV(4,IC) + &
					  ELSTIF(1,5,IC)*ESV(5,IC) + ELSTIF(1,6,IC)*ESV(6,IC))
		STRV(2,IC) = 0.5D0*(ELSTIF(2,1,IC)*ESV(1,IC) + ELSTIF(2,2,IC)*ESV(2,IC) + &
					  ELSTIF(2,3,IC)*ESV(3,IC) + ELSTIF(2,4,IC)*ESV(4,IC) + &
					  ELSTIF(2,5,IC)*ESV(5,IC) + ELSTIF(2,6,IC)*ESV(6,IC))
		STRV(3,IC) = 0.5D0*(ELSTIF(3,1,IC)*ESV(1,IC) + ELSTIF(3,2,IC)*ESV(2,IC) + &
					  ELSTIF(3,3,IC)*ESV(3,IC) + ELSTIF(3,4,IC)*ESV(4,IC) + &
					  ELSTIF(3,5,IC)*ESV(5,IC) + ELSTIF(3,6,IC)*ESV(6,IC))
		STRV(4,IC) = 0.5D0*(ELSTIF(4,1,IC)*ESV(1,IC) + ELSTIF(4,2,IC)*ESV(2,IC) + &
					  ELSTIF(4,3,IC)*ESV(3,IC) + ELSTIF(4,4,IC)*ESV(4,IC) + &
					  ELSTIF(4,5,IC)*ESV(5,IC) + ELSTIF(4,6,IC)*ESV(6,IC))
		STRV(5,IC) = 0.5D0*(ELSTIF(5,1,IC)*ESV(1,IC) + ELSTIF(5,2,IC)*ESV(2,IC) + &
					  ELSTIF(5,3,IC)*ESV(3,IC) + ELSTIF(5,4,IC)*ESV(4,IC) + &
					  ELSTIF(5,5,IC)*ESV(5,IC) + ELSTIF(5,6,IC)*ESV(6,IC))
		STRV(6,IC) = 0.5D0*(ELSTIF(6,1,IC)*ESV(1,IC) + ELSTIF(6,2,IC)*ESV(2,IC) + &
					  ELSTIF(6,3,IC)*ESV(3,IC) + ELSTIF(6,4,IC)*ESV(4,IC) + &
					  ELSTIF(6,5,IC)*ESV(5,IC) + ELSTIF(6,6,IC)*ESV(6,IC)) 
					  
		DO I = 1,NSLIP
			
			CALL MATMATMLT(FFFF(1,1,IC), SSYSMAT(1,1,I,IC), TEMP3(1,1,I,IC), 1)             
			
			B6A(1,I,IC) = TEMP3(1,1,I,IC)+TEMP3(1,1,I,IC)
			B6A(2,I,IC) = TEMP3(2,2,I,IC)+TEMP3(2,2,I,IC)
			B6A(3,I,IC) = TEMP3(3,3,I,IC)+TEMP3(3,3,I,IC)
			B6A(4,I,IC) = (TEMP3(1,2,I,IC)+TEMP3(2,1,I,IC))*2.0D0
			B6A(5,I,IC) = (TEMP3(2,3,I,IC)+TEMP3(3,2,I,IC))*2.0D0       
			B6A(6,I,IC) = (TEMP3(1,3,I,IC)+TEMP3(3,1,I,IC))*2.0D0             

			
			C6A(1,I,IC)=0.5D0*(ELSTIF(1,1,IC)*B6A(1,I,IC)+ELSTIF(1,2,IC)*B6A(2,I,IC)+ &
						ELSTIF(1,3,IC)*B6A(3,I,IC)+ELSTIF(1,4,IC)*B6A(4,I,IC) + &
						ELSTIF(1,5,IC)*B6A(5,I,IC)+ELSTIF(1,6,IC)*B6A(6,I,IC))
			C6A(2,I,IC)=0.5D0*(ELSTIF(2,1,IC)*B6A(1,I,IC)+ELSTIF(2,2,IC)*B6A(2,I,IC)+ &
						ELSTIF(2,3,IC)*B6A(3,I,IC)+ELSTIF(2,4,IC)*B6A(4,I,IC) + &
						ELSTIF(2,5,IC)*B6A(5,I,IC)+ELSTIF(2,6,IC)*B6A(6,I,IC))
			C6A(3,I,IC)=0.5D0*(ELSTIF(3,1,IC)*B6A(1,I,IC)+ELSTIF(3,2,IC)*B6A(2,I,IC)+ &
						ELSTIF(3,3,IC)*B6A(3,I,IC)+ELSTIF(3,4,IC)*B6A(4,I,IC) + &
						ELSTIF(3,5,IC)*B6A(5,I,IC)+ELSTIF(3,6,IC)*B6A(6,I,IC))
			C6A(4,I,IC)=0.5D0*(ELSTIF(4,1,IC)*B6A(1,I,IC)+ELSTIF(4,2,IC)*B6A(2,I,IC)+ &
						ELSTIF(4,3,IC)*B6A(3,I,IC)+ELSTIF(4,4,IC)*B6A(4,I,IC) + &
						ELSTIF(4,5,IC)*B6A(5,I,IC)+ELSTIF(4,6,IC)*B6A(6,I,IC))
			C6A(5,I,IC)=0.5D0*(ELSTIF(5,1,IC)*B6A(1,I,IC)+ELSTIF(5,2,IC)*B6A(2,I,IC)+ &
						ELSTIF(5,3,IC)*B6A(3,I,IC)+ELSTIF(5,4,IC)*B6A(4,I,IC) + &
						ELSTIF(5,5,IC)*B6A(5,I,IC)+ELSTIF(5,6,IC)*B6A(6,I,IC))
			C6A(6,I,IC)=0.5D0*(ELSTIF(6,1,IC)*B6A(1,I,IC)+ELSTIF(6,2,IC)*B6A(2,I,IC)+ &
						ELSTIF(6,3,IC)*B6A(3,I,IC)+ELSTIF(6,4,IC)*B6A(4,I,IC) + &
						ELSTIF(6,5,IC)*B6A(5,I,IC)+ELSTIF(6,6,IC)*B6A(6,I,IC))
			
		ENDDO
		!IF(IC ==1 .and. ISTEP == 2) print *,"C6A = ", C6A(6,6,IC)
		!IF(IC ==1 .and. ISTEP == 2) print *,"SSYSMAT = ", SSYSMAT(2,3,2,IC)
    ENDDO
    call cpu_time(tf)
	cpu_t(1) = cpu_t(1) + tf-ti
!===================================================================================
!=================================================================================== 

	block_CPSOLV = dim3(256, 1,1)
	grid_CPSOLV = dim3(NC/256, 1,1) 
	iostat = cudaEventRecord(beg, 0)
	SSYSMAT_d = SSYSMAT 
	C6A_d = C6A
    STRV_d = STRV

	iostat = cudaEventRecord(fin, 0)
    iostat = cudaThreadSynchronize()
    iostat = cudaEventElapsedTime(tf, beg, fin)	
	cuda_t(1) = cuda_t(1) + tf
	
	iostat = cudaEventRecord(beg, 0)	
	
		CALL CPSOLV<<<grid_CPSOLV, block_CPSOLV>>>(C6A_d,SSYSMAT_d,1, &
							STRS1PK_d,DG2_d, STRV_d, ISTEP)
		ir = cudaGetLastError()
		IF(ir) print *, cudaGetErrorString( ir ) 	
	iostat = cudaEventRecord(fin, 0)
    iostat = cudaThreadSynchronize()
    iostat = cudaEventElapsedTime(tf, beg, fin)	
	cuda_t(2) = cuda_t(2) + tf
	
	iostat = cudaEventRecord(beg, 0)

	STRS1PK = STRS1PK_d 
	DG2= DG2_d
	
    iostat = cudaEventRecord(fin, 0)
    iostat = cudaThreadSynchronize()
    iostat = cudaEventElapsedTime(tf, beg, fin)
    cuda_t(3) = cuda_t(3) + tf
!===================================================================================
!=================================================================================== 
    	
	DO IC = 1, NC
   		! IF(IC ==NC .and. ISTEP == 1) print *,"STRS1PK_2 = ", STRS1PK(1,IC)
		! IF(IC ==NC .and. ISTEP == 2) print *,"STRS1PK_2 = ", STRS1PK(1,IC)
		! IF(IC ==NC .and. ISTEP == 3) print *,"STRS1PK_2 = ", STRS1PK(1,IC)
  !  IF(IC ==1 .and. ISTEP == 1) print *,"DG2 = ", DG2(6,IC)   
		call cpu_time(ti)
		CALL DEFPIUPD(IC,IIPT,IC,SSYSMAT,DEFPI1,DEFPI2, &
		                 STRS2FEM,STRS1PK,DG2)
						 
		!IF(IC ==NC .and. ISTEP == 1) print *,"Stress2BE1 = ", STRESS2BE(1,IC)
		!IF(IC ==1 .and. ISTEP == 1) print *,"STRS2FEM_1 = ", STRS2FEM(3)
		!IF(IC ==1 .and. ISTEP == 2) print *,"STRS2FEM_2 = ", STRS2FEM(3)
		!IF(IC ==1 .and. ISTEP == 3) print *,"STRS2FEM_3 = ", STRS2FEM(3)
		call cpu_time(tf)  
		cpu_t(2) = cpu_t(2) + tf-ti
	!	IF(IC ==1 .and. ISTEP == 1) print *,"DEFPI1 = ", DEFPI1(2,1,IC)
	    !IF(IC ==200 .and. ISTEP == 1) print *,"STRS2FEM = ", STRS2FEM(3)
		!SCAUCHYAVG is zeroed each ISTEP
		SCAUCHYAVG(1) = SCAUCHYAVG(1)+STRS2FEM(1)/dble(NC)
		!IF(IC ==200 .and. ISTEP == 2) print *,"STRS2FEM = ", STRS2FEM(3)
		SCAUCHYAVG(2) = SCAUCHYAVG(2)+STRS2FEM(2)/dble(NC)
		SCAUCHYAVG(3) = SCAUCHYAVG(3)+STRS2FEM(3)/dble(NC)
		SCAUCHYAVG(4) = SCAUCHYAVG(4)+STRS2FEM(4)/dble(NC)
		SCAUCHYAVG(5) = SCAUCHYAVG(5)+STRS2FEM(5)/dble(NC)
		SCAUCHYAVG(6) = SCAUCHYAVG(6)+STRS2FEM(6)/dble(NC)
		
    ENDDO ! IC
      !IF(ISTEP == 1) print *,"STRS2FEM = ", STRS2FEM(3)
      if(TAU.EQ.TEXOUT)then
	  call cpu_time(ti)
      CALL EVOLTEXTURE(DEFPI1,DEFPI2,DEFGRDN,DEFGRDC,NC,QC2S, &
                      PHI1,PHI,PHI2)
	  call cpu_time(tf)
	  cpu_t(3) = cpu_t(3) + tf-ti
      WRITE(7,*)TAU
      DO IC = 1,NC
          WRITE(7,559)PHI1(IC)*180.0D0/PI ,PHI(IC)*180.0D0/PI, &
                     PHI2(IC)*180.0D0/PI,1.0D0/dble(NC)
      ENDDO
      endif

559   format(4f15.5) 
561   format(7f15.5) 


      write(9,561) STRAINLF, & !DABS(EQVST*TAU), 
                  SCAUCHYAVG(1),SCAUCHYAVG(2),SCAUCHYAVG(3), &
                  SCAUCHYAVG(4),SCAUCHYAVG(5),SCAUCHYAVG(6)
      IF(ISTEP == NSTEP) THEN 
		print *,"",STRAINLF, & !DABS(EQVST*TAU), 
                  SCAUCHYAVG(1),SCAUCHYAVG(2),SCAUCHYAVG(3), &
                  SCAUCHYAVG(4),SCAUCHYAVG(5),SCAUCHYAVG(6)
	  END IF
      ENDDO ! ISTEP

		call cpu_time(totf)
		cpu_t3=totf-toti
		
  WRITE(*,'(A70)')               'Time (seconds) and Percentage Time                                '
  WRITE(*,'(A70)')               '======================================================================='
 !WRITE(*,'(A39,4X,2(1PE13.4))') 'Host to Device Transfer           ', cuda_t(1)*1.e-3_rp, cuda_t(1)*1.e-3_rp/(totf-toti)*100
  WRITE(*,'(A39,4X,2(1PE13.4))') 'CPSOLV      	                    ', cuda_t(2)*1.e-3_rp, cuda_t(2)*1.e-3_rp/(totf-toti)*100
  WRITE(*,'(A39,4X,2(1PE13.4))') 'Device to Host Transfer           ', cuda_t(3)*1.e-3_rp, cuda_t(3)*1.e-3_rp/(totf-toti)*100
  WRITE(*,'(A39,4X,2(1PE13.4))') 'C Star                            ', cpu_t(1), cpu_t(1)/(totf-toti)*100
  WRITE(*,'(A39,4X,2(1PE13.4))') 'DEFPIUPD                          ', cpu_t(2), cpu_t(2)/(totf-toti)*100
  WRITE(*,'(A39,4X,2(1PE13.4))') 'Texture                           ', cpu_t(3), cpu_t(3)/(totf-toti)*100
  WRITE(*,'(A70)')               '======================================================================='

  WRITE(*,'(A39,4X,1PE13.4)')    'Total                                 ', cpu_t3
		

      CLOSE(7)        
      CLOSE(9)
      iostat = cudaEventDestroy(beg)
      iostat = cudaEventDestroy(fin)
	  

      DEALLOCATE (SSYSMAT(3,3,12,NC),SSYSMAT_d(3,3,12,NC)) 
      DEALLOCATE (ELSTIF(6,6,NC))
      DEALLOCATE (DEFPI1(3,3,NC))  
      DEALLOCATE (DEFPI2(3,3,NC))  
      DEALLOCATE (DG2(NSLIP,NC),DG2_d(NSLIP,NC))
      DEALLOCATE (DEFGRDN(3,3,NC))    
      DEALLOCATE (DEFGRDC(3,3,NC))       
      DEALLOCATE (QC2S(3,3,NC))        
      DEALLOCATE (STRS1PK(6,NC),STRS1PK_d(6,NC),STRV(6,NC),STRV_d(6,NC))
      DEALLOCATE (C6A(6,NSLIP,NC),C6A_d(6,NSLIP,NC))

      RETURN
      END


!---+----+----+----+----+----+----+----+----+----+----+----+----+----+--    
      SUBROUTINE CPSETUP(SSYSMAT,ELSTIF,DEFPI1,DEFPI2, &
                        DEFGRDN,DEFGRDC,NUMEL, &
                        PHI1,PHI,PHI2,QC2S)
	  USE Variables
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N) 

!      COMMON /FIRSTTIME/ IFST 
      COMMON /ELBRICK/ NIPT
      DIMENSION TEMP(3,3),SSYSMAT(3,3,12,NUMEL),ELSTIF(6,6,NUMEL)
      DIMENSION DEFPI1(3,3,NUMEL),DEFPI2(3,3,NUMEL)
      DIMENSION DEFGRDN(3,3,NUMEL), &
               DEFGRDC(3,3,NUMEL), &
               PHI1(NUMEL),PHI(NUMEL),PHI2(NUMEL), &
               QC2S(3,3,NUMEL)
     
      CHARACTER*80 IUSRVL
      COMMON /IUSR/ IUSRVL(10)
      
      REAL*8,DIMENSION(:,:,:),ALLOCATABLE:: QTC2S,SLPMATC
      REAL*8,DIMENSION(:,:),ALLOCATABLE:: SNOR, SDIR
      
      C11 = 168400.0D0
      C12 = 121400.0D0
      C44 = 75400.0D0
     
      ALLOCATE (QTC2S(3,3,NUMEL))      

      QC2S(1,1,:)=DCOS(PHI1)*DCOS(PHI2)-DSIN(PHI1)*DCOS(PHI)*DSIN(PHI2)
      QC2S(1,2,:)=-DCOS(PHI1)*DSIN(PHI2)-DSIN(PHI1)*DCOS(PHI)*DCOS(PHI2)
      QC2S(1,3,:)=DSIN(PHI1)*DSIN(PHI)
      QC2S(2,1,:)=DSIN(PHI1)*DCOS(PHI2)+DCOS(PHI1)*DCOS(PHI)*DSIN(PHI2)
      QC2S(2,2,:)=-DSIN(PHI1)*DSIN(PHI2)+DCOS(PHI1)*DCOS(PHI)*DCOS(PHI2)
      QC2S(2,3,:)=-DCOS(PHI1)*DSIN(PHI)
      QC2S(3,1,:)=DSIN(PHI)*DSIN(PHI2)
      QC2S(3,2,:)=DSIN(PHI)*DCOS(PHI2)
      QC2S(3,3,:)=DCOS(PHI)
      
      QTC2S(1,1,:)=QC2S(1,1,:)
      QTC2S(1,2,:)=QC2S(2,1,:)
      QTC2S(1,3,:)=QC2S(3,1,:)
      QTC2S(2,1,:)=QC2S(1,2,:)
      QTC2S(2,2,:)=QC2S(2,2,:)
      QTC2S(2,3,:)=QC2S(3,2,:)
      QTC2S(3,1,:)=QC2S(1,3,:)
      QTC2S(3,2,:)=QC2S(2,3,:)
      QTC2S(3,3,:)=QC2S(3,3,:)
      
              
      ELSTIF(1,1,:) = C12+2.0D0*C44+(C11-C12-2.0D0*C44)* &
                     (QC2S(1,1,:)*QC2S(1,1,:)*QC2S(1,1,:)*QC2S(1,1,:)+ &
                      QC2S(1,2,:)*QC2S(1,2,:)*QC2S(1,2,:)*QC2S(1,2,:)+ &
                      QC2S(1,3,:)*QC2S(1,3,:)*QC2S(1,3,:)*QC2S(1,3,:))
     
      ELSTIF(2,2,:) = C12+2.0D0*C44+(C11-C12-2.0D0*C44)* &
                     (QC2S(2,1,:)*QC2S(2,1,:)*QC2S(2,1,:)*QC2S(2,1,:)+ &
                      QC2S(2,2,:)*QC2S(2,2,:)*QC2S(2,2,:)*QC2S(2,2,:)+ &
                      QC2S(2,3,:)*QC2S(2,3,:)*QC2S(2,3,:)*QC2S(2,3,:))

      ELSTIF(3,3,:) = C12+2.0D0*C44+(C11-C12-2.0D0*C44)* &
                     (QC2S(3,1,:)*QC2S(3,1,:)*QC2S(3,1,:)*QC2S(3,1,:)+ &
                      QC2S(3,2,:)*QC2S(3,2,:)*QC2S(3,2,:)*QC2S(3,2,:)+ &
                      QC2S(3,3,:)*QC2S(3,3,:)*QC2S(3,3,:)*QC2S(3,3,:))

      ELSTIF(1,2,:) = C12+(C11-C12-2.0D0*C44)* &
                    (QC2S(1,1,:)*QC2S(1,1,:)*QC2S(2,1,:)*QC2S(2,1,:)+ &
                     QC2S(1,2,:)*QC2S(1,2,:)*QC2S(2,2,:)*QC2S(2,2,:)+ &
                     QC2S(1,3,:)*QC2S(1,3,:)*QC2S(2,3,:)*QC2S(2,3,:))

      ELSTIF(1,3,:) = C12+(C11-C12-2.0D0*C44)* &
                     (QC2S(1,1,:)*QC2S(1,1,:)*QC2S(3,1,:)*QC2S(3,1,:)+ &
                      QC2S(1,2,:)*QC2S(1,2,:)*QC2S(3,2,:)*QC2S(3,2,:)+ &
                      QC2S(1,3,:)*QC2S(1,3,:)*QC2S(3,3,:)*QC2S(3,3,:))

      ELSTIF(2,3,:) = C12+(C11-C12-2.0D0*C44)* &
                    (QC2S(2,1,:)*QC2S(2,1,:)*QC2S(3,1,:)*QC2S(3,1,:)+ &
                     QC2S(2,2,:)*QC2S(2,2,:)*QC2S(3,2,:)*QC2S(3,2,:)+ &
                     QC2S(2,3,:)*QC2S(2,3,:)*QC2S(3,3,:)*QC2S(3,3,:))
     
      ELSTIF(2,1,:)=ELSTIF(1,2,:)
      ELSTIF(3,1,:)=ELSTIF(1,3,:)
      ELSTIF(3,2,:)=ELSTIF(2,3,:)
          
      ELSTIF(4,4,:) = C44+(C11-C12-2.0D0*C44)* &
                     (QC2S(1,1,:)*QC2S(2,1,:)*QC2S(1,1,:)*QC2S(2,1,:)+ &
                      QC2S(1,2,:)*QC2S(2,2,:)*QC2S(1,2,:)*QC2S(2,2,:)+ &
                      QC2S(1,3,:)*QC2S(2,3,:)*QC2S(1,3,:)*QC2S(2,3,:))
     
      ELSTIF(6,6,:) = C44+(C11-C12-2.0D0*C44)* &
                     (QC2S(1,1,:)*QC2S(3,1,:)*QC2S(1,1,:)*QC2S(3,1,:)+ &
                      QC2S(1,2,:)*QC2S(3,2,:)*QC2S(1,2,:)*QC2S(3,2,:)+ &
                      QC2S(1,3,:)*QC2S(3,3,:)*QC2S(1,3,:)*QC2S(3,3,:))

      ELSTIF(5,5,:) = C44+(C11-C12-2.0D0*C44)* &
                     (QC2S(2,1,:)*QC2S(3,1,:)*QC2S(2,1,:)*QC2S(3,1,:)+ &
                      QC2S(2,2,:)*QC2S(3,2,:)*QC2S(2,2,:)*QC2S(3,2,:)+ &
                      QC2S(2,3,:)*QC2S(3,3,:)*QC2S(2,3,:)*QC2S(3,3,:))
    
      ELSTIF(4,6,:) = (C11-C12-2*C44)* &
                     (QC2S(1,1,:)*QC2S(2,1,:)*QC2S(1,1,:)*QC2S(3,1,:)+ &
                      QC2S(1,2,:)*QC2S(2,2,:)*QC2S(1,2,:)*QC2S(3,2,:)+ &
                      QC2S(1,3,:)*QC2S(2,3,:)*QC2S(1,3,:)*QC2S(3,3,:))     
     
      ELSTIF(4,5,:) = (C11-C12-2*C44)* &
                     (QC2S(1,1,:)*QC2S(2,1,:)*QC2S(2,1,:)*QC2S(3,1,:)+ &
                      QC2S(1,2,:)*QC2S(2,2,:)*QC2S(2,2,:)*QC2S(3,2,:)+ &
                      QC2S(1,3,:)*QC2S(2,3,:)*QC2S(2,3,:)*QC2S(3,3,:))     

      ELSTIF(5,6,:) = (C11-C12-2*C44)* &
                     (QC2S(2,1,:)*QC2S(3,1,:)*QC2S(1,1,:)*QC2S(3,1,:)+ &
                      QC2S(2,2,:)*QC2S(3,2,:)*QC2S(1,2,:)*QC2S(3,2,:)+ &
                      QC2S(2,3,:)*QC2S(3,3,:)*QC2S(1,3,:)*QC2S(3,3,:))


      ELSTIF(5,4,:)=ELSTIF(4,5,:)
      ELSTIF(6,4,:)=ELSTIF(4,6,:)
      ELSTIF(6,5,:)=ELSTIF(5,6,:)
        
      ELSTIF(1,4,:) = (C11-C12-2*C44)* &
                     (QC2S(1,1,:)*QC2S(1,1,:)*QC2S(1,1,:)*QC2S(2,1,:)+ &
                      QC2S(1,2,:)*QC2S(1,2,:)*QC2S(1,2,:)*QC2S(2,2,:)+ &
                      QC2S(1,3,:)*QC2S(1,3,:)*QC2S(1,3,:)*QC2S(2,3,:))


      ELSTIF(1,6,:) = (C11-C12-2*C44)* &
                     (QC2S(1,1,:)*QC2S(1,1,:)*QC2S(1,1,:)*QC2S(3,1,:)+ &
                      QC2S(1,2,:)*QC2S(1,2,:)*QC2S(1,2,:)*QC2S(3,2,:)+ &
                      QC2S(1,3,:)*QC2S(1,3,:)*QC2S(1,3,:)*QC2S(3,3,:))


      ELSTIF(1,5,:) = (C11-C12-2*C44)* &
                     (QC2S(1,1,:)*QC2S(1,1,:)*QC2S(2,1,:)*QC2S(3,1,:)+ &
                      QC2S(1,2,:)*QC2S(1,2,:)*QC2S(2,2,:)*QC2S(3,2,:)+ &
                      QC2S(1,3,:)*QC2S(1,3,:)*QC2S(2,3,:)*QC2S(3,3,:))


      ELSTIF(2,4,:) = (C11-C12-2*C44)* &
                     (QC2S(2,1,:)*QC2S(2,1,:)*QC2S(1,1,:)*QC2S(2,1,:)+ &
                      QC2S(2,2,:)*QC2S(2,2,:)*QC2S(1,2,:)*QC2S(2,2,:)+ &
                      QC2S(2,3,:)*QC2S(2,3,:)*QC2S(1,3,:)*QC2S(2,3,:))


      ELSTIF(2,6,:) = (C11-C12-2*C44)* &
                     (QC2S(2,1,:)*QC2S(2,1,:)*QC2S(1,1,:)*QC2S(3,1,:)+ &
                      QC2S(2,2,:)*QC2S(2,2,:)*QC2S(1,2,:)*QC2S(3,2,:)+ &
                      QC2S(2,3,:)*QC2S(2,3,:)*QC2S(1,3,:)*QC2S(3,3,:))


      ELSTIF(2,5,:) = (C11-C12-2*C44)* &
                     (QC2S(2,1,:)*QC2S(2,1,:)*QC2S(2,1,:)*QC2S(3,1,:)+ &
                      QC2S(2,2,:)*QC2S(2,2,:)*QC2S(2,2,:)*QC2S(3,2,:)+ &
                      QC2S(2,3,:)*QC2S(2,3,:)*QC2S(2,3,:)*QC2S(3,3,:))


      ELSTIF(3,4,:) = (C11-C12-2*C44)* &
                     (QC2S(3,1,:)*QC2S(3,1,:)*QC2S(1,1,:)*QC2S(2,1,:)+ &
                      QC2S(3,2,:)*QC2S(3,2,:)*QC2S(1,2,:)*QC2S(2,2,:)+ &
                      QC2S(3,3,:)*QC2S(3,3,:)*QC2S(1,3,:)*QC2S(2,3,:))



      ELSTIF(3,6,:) = (C11-C12-2*C44)* &
                     (QC2S(3,1,:)*QC2S(3,1,:)*QC2S(1,1,:)*QC2S(3,1,:)+ &
                      QC2S(3,2,:)*QC2S(3,2,:)*QC2S(1,2,:)*QC2S(3,2,:)+ &
                      QC2S(3,3,:)*QC2S(3,3,:)*QC2S(1,3,:)*QC2S(3,3,:))



      ELSTIF(3,5,:) = (C11-C12-2*C44)*&
                     (QC2S(3,1,:)*QC2S(3,1,:)*QC2S(2,1,:)*QC2S(3,1,:)+ &
                      QC2S(3,2,:)*QC2S(3,2,:)*QC2S(2,2,:)*QC2S(3,2,:)+ &
                      QC2S(3,3,:)*QC2S(3,3,:)*QC2S(2,3,:)*QC2S(3,3,:))
     
     
      ELSTIF(4,1,:)=ELSTIF(1,4,:)
      ELSTIF(4,2,:)=ELSTIF(2,4,:)
      ELSTIF(4,3,:)=ELSTIF(3,4,:)
      ELSTIF(5,1,:)=ELSTIF(1,5,:)
      ELSTIF(5,2,:)=ELSTIF(2,5,:)
      ELSTIF(5,3,:)=ELSTIF(3,5,:)
      ELSTIF(6,1,:)=ELSTIF(1,6,:)
      ELSTIF(6,2,:)=ELSTIF(2,6,:)
      ELSTIF(6,3,:)=ELSTIF(3,6,:)

      ALLOCATE (SLPMATC(3,3,NSLIP))
      ALLOCATE (SNOR(NSLIP,3))      
      ALLOCATE (SDIR(NSLIP,3))      

      SNOR(1,1)=1.0D0  
      SNOR(1,2)=1.0D0
      SNOR(1,3)=1.0D0
      SDIR(1,1)=1.0D0
      SDIR(1,2)=-1.0D0
      SDIR(1,3)=0.0D0
      SNOR(2,1)=1.0D0
      SNOR(2,2)=1.0D0
      SNOR(2,3)=1.0D0
      SDIR(2,1)=1.0D0  
      SDIR(2,2)= 0.0D0  
      SDIR(2,3)=-1.0D0
      SNOR(3,1)=1.0D0  
      SNOR(3,2)= 1.0D0  
      SNOR(3,3)=1.0D0
      SDIR(3,1)=0.0D0  
      SDIR(3,2)= 1.0D0  
      SDIR(3,3)=-1.0D0
      SNOR(4,1)=-1.0D0  
      SNOR(4,2)= 1.0D0  
      SNOR(4,3)=1.0D0
      SDIR(4,1)=1.0D0   
      SDIR(4,2)= 0.0D0   
      SDIR(4,3)=1.0D0
      SNOR(5,1)=-1.0D0  
      SNOR(5,2)= 1.0D0  
      SNOR(5,3)=1.0D0
      SDIR(5,1)=1.0D0  
      SDIR(5,2)= 1.0D0   
      SDIR(5,3)=0.0D0
      SNOR(6,1)=-1.0D0  
      SNOR(6,2)= 1.0D0  
      SNOR(6,3)=1.0D0
      SDIR(6,1)=0.0D0  
      SDIR(6,2)= 1.0D0   
      SDIR(6,3)=-1.0D0        
      SNOR(7,1)=1.0D0  
      SNOR(7,2)= -1.0D0  
      SNOR(7,3)=1.0D0
      SDIR(7,1)=1.0D0  
      SDIR(7,2)= 0.0D0   
      SDIR(7,3)=-1.0D0 
      SNOR(8,1)=1.0D0  
      SNOR(8,2)= -1.0D0  
      SNOR(8,3)=1.0D0
      SDIR(8,1)=0.0D0  
      SDIR(8,2)= 1.0D0   
      SDIR(8,3)=1.0D0
      SNOR(9,1)=1.0D0  
      SNOR(9,2)= -1.0D0  
      SNOR(9,3)=1.0D0
      SDIR(9,1)=1.0D0  
      SDIR(9,2)= 1.0D0   
      SDIR(9,3)=0.0D0        
      SNOR(10,1)=1.0D0  
      SNOR(10,2)= 1.0D0  
      SNOR(10,3)=-1.0D0
      SDIR(10,1)=1.0D0  
      SDIR(10,2)= -1.0D0   
      SDIR(10,3)=0.0D0        
      SNOR(11,1)=1.0D0  
      SNOR(11,2)= 1.0D0  
      SNOR(11,3)=-1.0D0
      SDIR(11,1)=1.0D0  
      SDIR(11,2)= 0.0D0   
      SDIR(11,3)=1.0D0        
      SNOR(12,1)=1.0D0  
      SNOR(12,2)= 1.0D0  
      SNOR(12,3)=-1.0D0
      SDIR(12,1)=0.0D0  
      SDIR(12,2)= 1.0D0   
      SDIR(12,3)=1.0D0
     

      DO I = 1,NSLIP
           SDIRNORM = DSQRT(SDIR(I,1)**2+SDIR(I,2)**2+SDIR(I,3)**2)
           SNORNORM = DSQRT(SNOR(I,1)**2+SNOR(I,2)**2+SNOR(I,3)**2)
           
           SDIR(I,1) = SDIR(I,1)/SDIRNORM
           SDIR(I,2) = SDIR(I,2)/SDIRNORM
           SDIR(I,3) = SDIR(I,3)/SDIRNORM             
           SNOR(I,1) = SNOR(I,1)/SNORNORM
           SNOR(I,2) = SNOR(I,2)/SNORNORM
           SNOR(I,3) = SNOR(I,3)/SNORNORM

           SLPMATC(1,1,I) = SDIR(I,1)*SNOR(I,1)
           SLPMATC(1,2,I) = SDIR(I,1)*SNOR(I,2)
           SLPMATC(1,3,I) = SDIR(I,1)*SNOR(I,3)
           SLPMATC(2,1,I) = SDIR(I,2)*SNOR(I,1)
           SLPMATC(2,2,I) = SDIR(I,2)*SNOR(I,2)
           SLPMATC(2,3,I) = SDIR(I,2)*SNOR(I,3)
           SLPMATC(3,1,I) = SDIR(I,3)*SNOR(I,1)
           SLPMATC(3,2,I) = SDIR(I,3)*SNOR(I,2)
           SLPMATC(3,3,I) = SDIR(I,3)*SNOR(I,3)
      ENDDO
      
      DO J =1,NUMEL
        DO I = 1,NSLIP
        
      CALL MATMATMLT(SLPMATC(1,1,I), QTC2S(1,1,J), TEMP(1,1), 1)     

      CALL MATMATMLT(QC2S(1,1,J), TEMP(1,1), SSYSMAT(1,1,I,J), 1)
     
        ENDDO
      ENDDO
      
      DO I = 1,NUMEL
    
           DEFPI1(1,1,I) = 1.0D0
           DEFPI1(2,2,I) = 1.0D0
           DEFPI1(3,3,I) = 1.0D0
           DEFPI1(1,2,I) = 0.0D0
           DEFPI1(2,1,I) = 0.0D0
           DEFPI1(1,3,I) = 0.0D0                      
           DEFPI1(3,1,I) = 0.0D0                                 
           DEFPI1(2,3,I) = 0.0D0                      
           DEFPI1(3,2,I) = 0.0D0       
           
           DEFPI2(1,1,I) = 1.0D0
           DEFPI2(2,2,I) = 1.0D0
           DEFPI2(3,3,I) = 1.0D0
           DEFPI2(1,2,I) = 0.0D0
           DEFPI2(2,1,I) = 0.0D0
           DEFPI2(1,3,I) = 0.0D0                      
           DEFPI2(3,1,I) = 0.0D0                                 
           DEFPI2(2,3,I) = 0.0D0                      
           DEFPI2(3,2,I) = 0.0D0  
           
           DEFGRDN(1,1,I) = 1.0D0
           DEFGRDN(2,2,I) = 1.0D0
           DEFGRDN(3,3,I) = 1.0D0
           DEFGRDN(1,2,I) = 0.0D0
           DEFGRDN(2,1,I) = 0.0D0
           DEFGRDN(1,3,I) = 0.0D0                      
           DEFGRDN(3,1,I) = 0.0D0                                 
           DEFGRDN(2,3,I) = 0.0D0                      
           DEFGRDN(3,2,I) = 0.0D0  
           
           DEFGRDC(1,1,I) = 1.0D0
           DEFGRDC(2,2,I) = 1.0D0
           DEFGRDC(3,3,I) = 1.0D0
           DEFGRDC(1,2,I) = 0.0D0
           DEFGRDC(2,1,I) = 0.0D0
           DEFGRDC(1,3,I) = 0.0D0                      
           DEFGRDC(3,1,I) = 0.0D0                                 
           DEFGRDC(2,3,I) = 0.0D0                      
           DEFGRDC(3,2,I) = 0.0D0 
      
                                      

          

      ENDDO 

      DEALLOCATE(QTC2S)
      DEALLOCATE(SNOR)
      DEALLOCATE(SDIR)
      DEALLOCATE(SLPMATC)      
	
      RETURN
      END
      
      
!---+----+----+----+----+----+----+----+----+----+----+----+----+----+--      
      SUBROUTINE MATMATMLT(AMAT, BMAT, CMAT, IFLG)
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
      DIMENSION AMAT(3,3), BMAT(3,3), CMAT(3,3)
      
      IF(IFLG.EQ.1)THEN
      CMAT(1,1) = AMAT(1,1)*BMAT(1,1) + AMAT(1,2)*BMAT(2,1) + &
                     AMAT(1,3)*BMAT(3,1)
      CMAT(1,2) = AMAT(1,1)*BMAT(1,2) + AMAT(1,2)*BMAT(2,2) + &
                     AMAT(1,3)*BMAT(3,2)
      CMAT(1,3) = AMAT(1,1)*BMAT(1,3) + AMAT(1,2)*BMAT(2,3) + &
                     AMAT(1,3)*BMAT(3,3)
      CMAT(2,1) = AMAT(2,1)*BMAT(1,1) + AMAT(2,2)*BMAT(2,1) + &
                     AMAT(2,3)*BMAT(3,1)
      CMAT(2,2) = AMAT(2,1)*BMAT(1,2) + AMAT(2,2)*BMAT(2,2) + &
                     AMAT(2,3)*BMAT(3,2)
      CMAT(2,3) = AMAT(2,1)*BMAT(1,3) + AMAT(2,2)*BMAT(2,3) + &
                     AMAT(2,3)*BMAT(3,3)
      CMAT(3,1) = AMAT(3,1)*BMAT(1,1) + AMAT(3,2)*BMAT(2,1) + &
                     AMAT(3,3)*BMAT(3,1)
      CMAT(3,2) = AMAT(3,1)*BMAT(1,2) + AMAT(3,2)*BMAT(2,2) + &
                     AMAT(3,3)*BMAT(3,2)
      CMAT(3,3) = AMAT(3,1)*BMAT(1,3) + AMAT(3,2)*BMAT(2,3) + &
                     AMAT(3,3)*BMAT(3,3)
      ENDIF
      
      
      IF(IFLG.EQ.2)THEN
      CMAT(1,1) = AMAT(1,1)*BMAT(1,1) + AMAT(2,1)*BMAT(2,1) + &
                     AMAT(3,1)*BMAT(3,1)
      CMAT(1,2) = AMAT(1,1)*BMAT(1,2) + AMAT(2,1)*BMAT(2,2) + &
                     AMAT(3,1)*BMAT(3,2)
      CMAT(1,3) = AMAT(1,1)*BMAT(1,3) + AMAT(2,1)*BMAT(2,3) + &
                     AMAT(3,1)*BMAT(3,3)
      CMAT(2,1) = AMAT(1,2)*BMAT(1,1) + AMAT(2,2)*BMAT(2,1) + &
                     AMAT(3,2)*BMAT(3,1)
      CMAT(2,2) = AMAT(1,2)*BMAT(1,2) + AMAT(2,2)*BMAT(2,2) + &
                     AMAT(3,2)*BMAT(3,2)
      CMAT(2,3) = AMAT(1,2)*BMAT(1,3) + AMAT(2,2)*BMAT(2,3) + &
                     AMAT(3,2)*BMAT(3,3)
      CMAT(3,1) = AMAT(1,3)*BMAT(1,1) + AMAT(2,3)*BMAT(2,1) + &
                     AMAT(3,3)*BMAT(3,1)
      CMAT(3,2) = AMAT(1,3)*BMAT(1,2) + AMAT(2,3)*BMAT(2,2) + &
                     AMAT(3,3)*BMAT(3,2)
      CMAT(3,3) = AMAT(1,3)*BMAT(1,3) + AMAT(2,3)*BMAT(2,3) + &
                     AMAT(3,3)*BMAT(3,3)
      ENDIF
      
      IF(IFLG.EQ.3)THEN
      CMAT(1,1) = AMAT(1,1)*BMAT(1,1) + AMAT(1,2)*BMAT(1,2) + &
                     AMAT(1,3)*BMAT(1,3)
      CMAT(1,2) = AMAT(1,1)*BMAT(2,1) + AMAT(1,2)*BMAT(2,2) + &
                     AMAT(1,3)*BMAT(2,3)
      CMAT(1,3) = AMAT(1,1)*BMAT(3,1) + AMAT(1,2)*BMAT(3,2) + &
                     AMAT(1,3)*BMAT(3,3)
      CMAT(2,1) = AMAT(2,1)*BMAT(1,1) + AMAT(2,2)*BMAT(1,2) + &
                     AMAT(2,3)*BMAT(1,3)
      CMAT(2,2) = AMAT(2,1)*BMAT(2,1) + AMAT(2,2)*BMAT(2,2) + &
                     AMAT(2,3)*BMAT(2,3)
      CMAT(2,3) = AMAT(2,1)*BMAT(3,1) + AMAT(2,2)*BMAT(3,2) + &
                     AMAT(2,3)*BMAT(3,3)
      CMAT(3,1) = AMAT(3,1)*BMAT(1,1) + AMAT(3,2)*BMAT(1,2) + &
                     AMAT(3,3)*BMAT(1,3)
      CMAT(3,2) = AMAT(3,1)*BMAT(2,1) + AMAT(3,2)*BMAT(2,2) + &
                     AMAT(3,3)*BMAT(2,3)
      CMAT(3,3) = AMAT(3,1)*BMAT(3,1) + AMAT(3,2)*BMAT(3,2) + &
                     AMAT(3,3)*BMAT(3,3)
      ENDIF
      
      RETURN
      END      
      
!---+----+----+----+----+----+----+----+----+----+----+----+----+----+--
!
      SUBROUTINE CPMATINV(A,N,ID)
      IMPLICIT REAL*8 (A-H,O-Z)
!      INCLUDE 'IMPL.INC'      
!                                                                       
!     ******************************************************************
!     THIS SUBROUTINE COMPUTES THE INVERSE AND DETERMINANT OF MATRIX A *
!     OF ORDER N,BY THE GAUSS-JORDAN METHOD, A-INVERSE REPLACES A ,AND *
!     THE DETERMINANT OF A IS PLACED IN DETERM. IF M=1 THE VECTOR B    *
!     CONTAINS THE CONSTANT VECTOR WHEN MATINV IS CALLED, AND THIS IS  *
!     REPLACED WITH THE SOLUTION VECTOR IF M=0,NO SIMULTANEOUS         *
!     EQUATION SOLUTION IS CALLED FOR, AND B IS NOT PERTINENT. N IS NOT*
!     TO EXCEED 100.                                                   *
!      A--IS THE MATRIX OF COEFFICIENTS OR THE MATRIX TO BE INVERTED.  *
!     A CONTANS A-INVERSE AFTER EXECUTION.                            *
!      N-- IS THE ORDER OF THE SQUARE MATRIX, A.                       *
!      B--IS THE MATRIX CONTANING COLUMN VECTORS OF CONSTANTS (EACH    *
!         COLUMN VECTOR IS ASSOCIATED WITH A IN THE FOLLOWING          *
!         MANNER--AX=B.).                                              *
!      M--IS THE NUMBER OF CRITERION VECTORS (COLUMN VECTORS OF        *
!         SIMULTANEOUS SOLUTIONS) TO BE CALCULATED AND STORED IN B.    *
!      M=0 RESULTS IN THE COMPUTATION OF ONLY THE INVERSE AND          *
!          DETERMINANT OF A.                                           *
!      DETERM--CONTANS THE VALUE OF THE DETERMINANT AFTER EXECUTION.   *
!     ******************************************************************
     
      DIMENSION IPIVOT(1000),IND(1000,2),PIVOT(1000)                        
      DIMENSION A(ID,1),B(500,1)                                     

!     INITIALIZATION                                                    
!      DETERM=1.0E0                                                      
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
      END    
      
      
!---+----+----+----+----+----+----+----+----+----+----+----+----+----+--
      SUBROUTINE DEFPIUPD(IC,IIPT,IELM,SSYSMAT,DEFPI1,DEFPI2, &
                      STRS2FEM,STRS1PK,DG2)
	  USE Variables
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)

      COMMON /CAUCHYST/ STRESK(3,3)
      COMMON /JACTMP/TANG1(3,3) 
      COMMON /ELBRICK/ NIPT
      
      DIMENSION SSYSMAT(3,3,NSLIP,*),DEFPI1(3,3,*), &
               DEFPI2(3,3,*),STRS1PK(6,*),DG2(NSLIP,*)	
      DIMENSION DEFPI2TEMP(3,3),STRESPK(3,3), &
               TEMP(3,3),STRS2FEM(6), DETDEFS2

      DVELP(1,1) = 0.0D0
      DVELP(1,2) = 0.0D0
      DVELP(1,3) = 0.0D0
      DVELP(2,1) = 0.0D0
      DVELP(2,2) = 0.0D0
      DVELP(2,3) = 0.0D0
      DVELP(3,1) = 0.0D0
      DVELP(3,2) = 0.0D0
      DVELP(3,3) = 0.0D0
      DO I = 1,NSLIP
        DVELP(1,1) = DVELP(1,1)+DG2(I,IC)*SSYSMAT(1,1,I,IC)
        DVELP(1,2) = DVELP(1,2)+DG2(I,IC)*SSYSMAT(1,2,I,IC)
        DVELP(1,3) = DVELP(1,3)+DG2(I,IC)*SSYSMAT(1,3,I,IC)

        DVELP(2,1) = DVELP(2,1)+DG2(I,IC)*SSYSMAT(2,1,I,IC)
        DVELP(2,2) = DVELP(2,2)+DG2(I,IC)*SSYSMAT(2,2,I,IC)
        DVELP(2,3) = DVELP(2,3)+DG2(I,IC)*SSYSMAT(2,3,I,IC)
        DVELP(3,1) = DVELP(3,1)+DG2(I,IC)*SSYSMAT(3,1,I,IC)
        DVELP(3,2) = DVELP(3,2)+DG2(I,IC)*SSYSMAT(3,2,I,IC)
        DVELP(3,3) = DVELP(3,3)+DG2(I,IC)*SSYSMAT(3,3,I,IC)
      ENDDO        
      TEMP(1,1) = 1.0D0-DVELP(1,1)
      TEMP(2,2) = 1.0D0-DVELP(2,2)
      TEMP(3,3) = 1.0D0-DVELP(3,3)
      TEMP(1,2) = -DVELP(1,2)
      TEMP(2,1) = -DVELP(2,1)
      TEMP(3,1) = -DVELP(3,1)
      TEMP(1,3) = -DVELP(1,3)
      TEMP(3,2) = -DVELP(3,2)
      TEMP(2,3) = -DVELP(2,3)      
      
      CALL MATMATMLT(DEFPI1(1,1,IELM),TEMP(1,1),DEFPI2TEMP(1,1),1)

!
      DET = DEFPI2TEMP(1,1)*(DEFPI2TEMP(2,2)*DEFPI2TEMP(3,3)- &
                       DEFPI2TEMP(2,3)*DEFPI2TEMP(3,2))- &
           DEFPI2TEMP(1,2)*(DEFPI2TEMP(2,1)*DEFPI2TEMP(3,3)- &
                       DEFPI2TEMP(3,1)*DEFPI2TEMP(2,3))+ &
           DEFPI2TEMP(1,3)*(DEFPI2TEMP(2,1)*DEFPI2TEMP(3,2)- &
                       DEFPI2TEMP(3,1)*DEFPI2TEMP(2,2))
!
        IF(DABS(DET-1.0D0).GT.1.0D-3)THEN
          WRITE(*,*)'DET OF DEFPI2TEMP NOT GOOD = ',DET,'IC = ',IC
        ENDIF

      IF (DET.GT.0.0D0) THEN
        DEFPI2TEMP(1,1) = DEFPI2TEMP(1,1)/(DET**(1.0D0/3.0D0))
        DEFPI2TEMP(1,2) = DEFPI2TEMP(1,2)/(DET**(1.0D0/3.0D0))
        DEFPI2TEMP(1,3) = DEFPI2TEMP(1,3)/(DET**(1.0D0/3.0D0))
        DEFPI2TEMP(2,1) = DEFPI2TEMP(2,1)/(DET**(1.0D0/3.0D0))
        DEFPI2TEMP(2,2) = DEFPI2TEMP(2,2)/(DET**(1.0D0/3.0D0))
        DEFPI2TEMP(2,3) = DEFPI2TEMP(2,3)/(DET**(1.0D0/3.0D0))
        DEFPI2TEMP(3,1) = DEFPI2TEMP(3,1)/(DET**(1.0D0/3.0D0))
        DEFPI2TEMP(3,2) = DEFPI2TEMP(3,2)/(DET**(1.0D0/3.0D0))
        DEFPI2TEMP(3,3) = DEFPI2TEMP(3,3)/(DET**(1.0D0/3.0D0))
      ENDIF

      CALL MATMATMLT(DEF2(1,1), DEFPI2TEMP(1,1), DEFS2(1,1), 1)


      STRESPK(1,1)=STRS1PK(1,IC) 
      STRESPK(1,2)=STRS1PK(4,IC) 
      STRESPK(2,3)=STRS1PK(5,IC)
      STRESPK(2,1)=STRS1PK(4,IC) 
      STRESPK(2,2)=STRS1PK(2,IC) 
      STRESPK(1,3)=STRS1PK(6,IC)
      STRESPK(3,2)=STRS1PK(5,IC) 
      STRESPK(3,1)=STRS1PK(6,IC) 
      STRESPK(3,3)=STRS1PK(3,IC)

      CALL MATMATMLT(STRESPK(1,1), DEFS2(1,1), TANG1(1,1), 3)


      CALL MATMATMLT(DEFS2(1,1), TANG1(1,1), STRESK(1,1) , 1)
      

      DETDEFS2=DEFS2(1,1)*(DEFS2(2,2)*DEFS2(3,3)-DEFS2(2,3)*DEFS2(3,2))- &
              DEFS2(1,2)*(DEFS2(2,1)*DEFS2(3,3)-DEFS2(3,1)*DEFS2(2,3))+ &
              DEFS2(1,3)*(DEFS2(2,1)*DEFS2(3,2)-DEFS2(3,1)*DEFS2(2,2))

            STRESK(1,1) = STRESK(1,1)/DETDEFS2
            STRESK(1,2) = STRESK(1,2)/DETDEFS2
            STRESK(1,3) = STRESK(1,3)/DETDEFS2
            STRESK(2,1) = STRESK(2,1)/DETDEFS2
            STRESK(2,2) = STRESK(2,2)/DETDEFS2
            STRESK(2,3) = STRESK(2,3)/DETDEFS2
            STRESK(3,1) = STRESK(3,1)/DETDEFS2
            STRESK(3,2) = STRESK(3,2)/DETDEFS2
            STRESK(3,3) = STRESK(3,3)/DETDEFS2
            
      DETDEFS2 = 1.0D0/DETDEFS2
            
      DO  I = 1,3
        DO  II = 1,3            
          DEFPI2(I,II,IELM) = DEFPI2TEMP(I,II)
          DEFPI1(I,II,IELM) = DEFPI2TEMP(I,II)
        ENDDO
      ENDDO
 
      STRS2FEM(1)=STRESK(1,1) 
      STRS2FEM(2)=STRESK(2,2) 
      STRS2FEM(3)=STRESK(3,3)
      STRS2FEM(4)=STRESK(1,2) 
      STRS2FEM(6)=STRESK(1,3) 
      STRS2FEM(5)=STRESK(2,3) 
      
      RETURN
      END   
            
                  
!---+----+----+----+----+----+----+----+----+----+----+----+----+----+--      
      SUBROUTINE EVOLTEXTURE(DEFPI1,DEFPI2,DEFGRDN,DEFGRDC,NUMEL,QC2S, &
                            PHI1OUT,PHIOUT,PHI2OUT)
	  USE Variables
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
       
      COMMON /FIRSTTIME/ IFST 
      COMMON /ELBRICK/ NIPT     
 
      COMMON/VECT5/ F1,F2,F3,F4,F5,F6,F7,F8,F9, & 
     	 WXX,WYY,WZZ,DETI,DUM5(45)
      COMMON/VECT5MK/UN2C(3,3),RN2C(3,3)
      DIMENSION DEFPI1(3,3,NUMEL),DEFPI2(3,3,NUMEL)
      DIMENSION DEFGRDN(3,3,*),DEFGRDC(3,3,*),TEMP(3,3), &
               TXMT(3,3),ANGLES(3),QC2S(3,3,NUMEL), &
               PHI1OUT(NUMEL),PHIOUT(NUMEL),PHI2OUT(NUMEL)   
       
      PI = 4.0D0*DATAN(1.0D0)      
!      OPEN (29,FILE = 'ODFOUT.TXT',ACCESS= 'APPEND',STATUS = 'UNKNOWN')           
      DO J = 1,NUMEL

          CALL MATMATMLT(DEF2(1,1), DEFPI2(1,1,J), TEMP(1,1), 1)
          
      	    F1 = TEMP(1,1)
      	    F2 = TEMP(2,2)
      	    F3 = TEMP(3,3)      	     
      	    F4 = TEMP(1,2)
      	    F5 = TEMP(1,3)
      	    F6 = TEMP(2,1)      	    
      	    F7 = TEMP(2,3)
      	    F8 = TEMP(3,1)
      	    F9 = TEMP(3,2)      	   
          
            CALL GETVRT(Q11,Q12,Q13,Q21,Q22,Q23,Q31,Q32,Q33)
            
            CALL MATMATMLT(RN2C(1,1), QC2S(1,1,J), TXMT(1,1), 1)
            
            IF(ABS(TXMT(3,3)) == 1.0D0)THEN
              ANGLES(1) = DACOS(TXMT(1,1))
              ANGLES(2) = DACOS(TXMT(3,3))
              ANGLES(3) = 0.0

              IF(TXMT(2,1)<0.0)THEN
                ANGLES(1) = 2*PI - ANGLES(1)
              ENDIF
            ELSE
              ANGLES(1) = DATAN2(TXMT(1,3),-1.0D0*TXMT(2,3))
              ANGLES(2) = DACOS(TXMT(3,3))
              ANGLES(3) = DATAN2(TXMT(3,1),TXMT(3,2))
            ENDIF

!           MAKE SURE ANGLES ARE POSITIVE
                IF(ANGLES(1) <0.0D0 )THEN
                    ANGLES(1) = ANGLES(1)+ 2.0D0*PI
                ENDIF
                IF(ANGLES(2) <0.0D0 )THEN
                    ANGLES(2) = ANGLES(2)+ 2.0D0*PI
                ENDIF
                IF(ANGLES(3) <0.0D0 )THEN
                    ANGLES(3) = ANGLES(3)+ 2.0D0*PI
                ENDIF
                
!         WRITE(29,559)180.0D0/PI*ANGLES(1),180.0D0/PI*ANGLES(2),
!     +          180.0D0/PI*ANGLES(3)
	PHI1OUT(J) = ANGLES(1)
	PHIOUT(J) = ANGLES(2)
        PHI2OUT(J) = ANGLES(3)
      
      ENDDO
!      CLOSE(29)
559   FORMAT(3F15.9) 
      
      RETURN
      END
                        
                        
                        
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE GETVRT(Q11,Q12,Q13,Q21,Q22,Q23,Q31,Q32,Q33)
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
      COMMON/VECT5/ F11,F22,F33,F12,F13,F21, & 
       F23,F31,F32,C11,C22,C33, &
       C12,C13,C23,C1212,C2323,C1313, &
       ALI1,ALI2,ALI3,UI1,UI2,UI3, &
       S11,S22,S33,S12,S13,S23, &
       UI11,UI22,UI33,UI12,UI13,UI23, &
       SI11,SI22,SI33,SI12,SI13,SI23, &
       AMDA1,AMDA2,AMDA3,R,Q,SDETM, &
       XMOD,UI1S,SCL1,SCL2,SCL0,CT3,ST3, &
       C2313,C1223,C1213
     
      COMMON/VECT5MK/UN2C(3,3),RN2C(3,3)

      DATA O3/.333333333333333333333333D0/

      XXX = DABS(F12) + DABS(F21) + DABS(F13) + DABS(F31) + DABS(F23) + &
           DABS(F32) 
      IF(XXX.LT.1.D-12) THEN
        Q11 = 1.D0
        Q22 = 1.D0
        Q33 = 1.D0
        Q12 = 0.D0
        Q21 = 0.D0
        Q13 = 0.D0
        Q31 = 0.D0
        Q23 = 0.D0
        Q32 = 0.D0
        GOTO 200
      ENDIF

      ROOT3=DSQRT(3.D0)
      C11=F11*F11+F21*F21+F31*F31
      C12=F11*F12+F21*F22+F31*F32
      C13=F11*F13+F21*F23+F31*F33
      C23=F12*F13+F22*F23+F32*F33
      C22=F12*F12+F22*F22+F32*F32
   10 C33=F13*F13+F23*F23+F33*F33

      C1212=C12*C12
      C1313=C13*C13
      C2323=C23*C23
      C2313=C23*C13
      C1223=C12*C23
      C1213=C12*C13
      S11=C22*C33-C2323
      UI1=O3*(C11+C22+C33)  
!      UI2=S11+C11*C22+C33*C11-C1212-C1313
!      UI3=C11*S11+C12*(C2313-C12*C33)
!     1      +C13*(C1223-C22*C13)
      
      UI2=C11*C11+C22*C22+C33*C33+C12*C12+C12*C12+C13*C13+C13*C13+ &
            C23*C23+C23*C23
      UI2=0.5D0*((C11+C22+C33)*(C11+C22+C33)-UI2)
      UI3=C11*(C22*C33-C23*C23)-C12*(C12*C33-C23*C13)+ &
     	     C13*(C12*C23-C22*C13)
     
                
   20 UI1S=UI1*UI1

      Q    =DSQRT(-DMIN1(O3*UI2-UI1S,0.D0))
      R    =0.5D0*(UI3-UI1*UI2)+UI1*UI1S

      XMOD =Q*Q*Q
   30 SCL1 =.5D0+DSIGN(.5D0,XMOD-1.D-30)

      SCL2 =.5D0+DSIGN(.5D0,XMOD-DABS(R))
      SCL0 =DMIN1(SCL1,SCL2)
      SCL1 =1.0D0-SCL0

      XXX = XMOD+SCL1
      IF(XXX.EQ.0.) XXX = 1.D-20
      XXX=R/XXX
      IF(DABS(XXX).GT.1.D0) XXX = DABS(XXX)/XXX
      SDETM=DACOS(XXX)*O3
!      SDETM=ACOS(R/(XMOD+SCL1))*O3
      Q  =SCL0*Q
      CT3=Q*DCOS(SDETM)
   40 ST3=Q*ROOT3*DSIN(SDETM)

      SDETM=SCL1*DSQRT(DMAX1(0.0D0,R))
      AMDA1=2.000*(CT3+SDETM)+UI1
      AMDA1=DMAX1(0.D0,AMDA1)
      AMDA1=DSQRT(AMDA1)
      AMDA2=-CT3+ST3-SDETM+UI1
      AMDA3=-CT3-ST3-SDETM+UI1
      AMDA2=DMAX1(0.D0,AMDA2)      
      AMDA3=DMAX1(0.D0,AMDA3)      
      AMDA2=DSQRT(AMDA2)
   50 AMDA3=DSQRT(AMDA3)
!      AMDA1=SQRT(2.000*(CT3+SDETM)+UI1)
!      AMDA2=SQRT(-CT3+ST3-SDETM+UI1)
!   50 AMDA3=SQRT(-CT3-ST3-SDETM+UI1)

      SDETM=AMDA1*AMDA2
      ALI1=AMDA1+AMDA2+AMDA3
      ALI2= SDETM+AMDA2*AMDA3+AMDA3*AMDA1
   60 ALI3= SDETM*AMDA3/ALI1

      S11= C11+ALI3
      S22= C22+ALI3
      S33= C33+ALI3
      S12= C2313-C12*S33
      S13= C1223-S22*C13
      S23=-C2323+S22*S33
   70 DD = S11*S23+C12*S12+C13*S13
      SDETM=ALI1*DD
      IF(SDETM.EQ.0.) THEN
         SDETM = 1.D-20
      ENDIF
      SDETM=1./SDETM

      C11=C11+ALI2
      C22=C22+ALI2
      C33=C33+ALI2
      SI11=SDETM*S23
      SI12=SDETM*S12
      SI13=SDETM*S13
      SI22=SDETM*( S11*S33-C1313)
      SI23=SDETM*(-S11*C23+C1213)
   90 SI33=SDETM*( S11*S22-C1212)

      S12=C12*SI12
      S13=C13*SI13
      S23=C23*SI23
      UI11=C11*SI11+S12+S13
      UI22=S12+C22*SI22+S23
      UI33=S13+S23+C33*SI33
      UI12=C11*SI12+C12*SI22+C13*SI23
      UI13=C11*SI13+C12*SI23+C13*SI33
  100 UI23=C12*SI13+C22*SI23+C23*SI33

      Q11=F11*UI11+F12*UI12+F13*UI13
      Q12=F11*UI12+F12*UI22+F13*UI23
      Q13=F11*UI13+F12*UI23+F13*UI33
      Q21=F21*UI11+F22*UI12+F23*UI13
      Q22=F21*UI12+F22*UI22+F23*UI23
      Q23=F21*UI13+F22*UI23+F23*UI33
      Q31=F31*UI11+F32*UI12+F33*UI13
      Q32=F31*UI12+F32*UI22+F33*UI23
  120 Q33=F31*UI13+F32*UI23+F33*UI33
      
      UNC11=Q11*F11+Q21*F21+Q31*F31
      UNC22=Q12*F12+Q22*F22+Q32*F32
      UNC33=Q13*F13+Q23*F23+Q33*F33
      UNC12=Q11*F12+Q21*F22+Q31*F32
      UNC23=Q12*F13+Q22*F23+Q32*F33
  130 UNC13=Q11*F13+Q21*F23+Q31*F33
  
      UN2C(1,1) = UNC11
      UN2C(2,2) = UNC22
      UN2C(3,3) = UNC33
      UN2C(1,2) = UNC12
      UN2C(2,1) = UNC12
      UN2C(1,3) = UNC13
      UN2C(3,1) = UNC13
      UN2C(2,3) = UNC23
      UN2C(3,2) = UNC23
      
 200  RN2C(1,1) = Q11
      RN2C(2,2) = Q22
      RN2C(3,3) = Q33
      RN2C(1,2) = Q12
      RN2C(2,1) = Q21
      RN2C(1,3) = Q13
      RN2C(3,1) = Q31
      RN2C(2,3) = Q23
      RN2C(3,2) = Q32
      
      RETURN
      END
                              
