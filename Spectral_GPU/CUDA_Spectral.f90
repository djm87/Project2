PROGRAM SCP_FCC
  USE cudafor
  USE cublas
  USE CUDA_Kernels
  USE Variables
  !IMPLICIT NONE
  
  !!!! DEVICE VARS !!!!
  REAL(RP), DIMENSION(:,:), ALLOCATABLE         :: FoMatRE
  REAL(RP), DIMENSION(:,:), ALLOCATABLE         :: FoMatIM
  REAL(RP), DIMENSION(:,:), ALLOCATABLE, DEVICE :: FoMat_d !This allocates FoMat_d in the device global memory (as seen by this document)
 
  REAL(RP), DIMENSION(:,:), ALLOCATABLE, DEVICE :: StressMat_d !This allocates StressMat_d in the device global memory (as seen by this document)
  REAL(RP), DIMENSION(:,:), ALLOCATABLE, PINNED :: StressMat !Pin memory is allocate and doesn't move around in the memory (always the same address). Good practice when moving data between host and device
  
  REAL(RP), DIMENSION(:,:), ALLOCATABLE         :: SC
  REAL(RP), DIMENSION(:,:), ALLOCATABLE, DEVICE :: SC_d

  REAL(RP), DIMENSION(:), ALLOCATABLE, DEVICE :: s_dot_d
  REAL(RP), DIMENSION(:,:), ALLOCATABLE, DEVICE ::  W_app_pr_d
  REAL(RP), DIMENSION(:,:), ALLOCATABLE, PINNED ::  W_app_pr
  
  REAL(RP), DIMENSION(:), ALLOCATABLE :: sum_gamma_dot

  
  REAL(RP), DIMENSION(:), ALLOCATABLE, PINNED :: Phi1, PHI, Phi2, termsdim4
  REAL(RP), DIMENSION(:), ALLOCATABLE, DEVICE :: Phi1_d, PHI_d, Phi2_d, termsdim4_d 

  REAL(RP), DIMENSION(:,:), ALLOCATABLE, PINNED :: Super_set_out 
  REAL(RP), DIMENSION(:), ALLOCATABLE :: Super_set_120_1, Super_set_120_2, Super_set_120_3, Super_set_120_4
  REAL(RP), DIMENSION(:), ALLOCATABLE, DEVICE :: Super_set_1_d, Super_set_2_d, Super_set_3_d, Super_set_4_d
  

  REAL(RP), DIMENSION(:), ALLOCATABLE, DEVICE :: Q_p_sam_d
  
  !REAL(RP), DIMENSION(:,:,:), ALLOCATABLE, DEVICE :: QMATT, RMATT

			 
  !!!!!!!!!!!!!!!!!!!!!!!

  INTEGER, PARAMETER::  NMAX = 3
  INTEGER :: ind_max, ind_min, ind_mid, n_steps, klm, nrot !change, added klm, nrot,sizeFm
  INTEGER :: i, j, k, ij, jk
  INTEGER :: imode
  INTEGER :: ICRYS, JSET
  
  REAL(RP), DIMENSION(3,3) :: L, D_app_sa, W_app_sa, V_eig,&
                              STRESS_TENSOR_PR, STRESS_TENSOR_SA

  REAL(RP), DIMENSION(:),   ALLOCATABLE :: S11, S22, S12, S13, S23
  REAL(RP), DIMENSION(:,:), ALLOCATABLE :: STRESS_VEC
  REAL(RP) :: eps_dot, DT, fact_DT, theta, Total_Strain
  REAL     :: ti, tf, cpu_t(3), toti, totf
  REAL     :: cuda_t(7)

  TYPE(dim3) :: grid_sz, block_sz, blocks_texture, threads_texture
  TYPE(dim3) :: grid_sz_1D, block_sz_1D
  TYPE(cudaEvent) :: beg, fin 
  INTEGER(kind=cuda_stream_kind) :: streams(16)
  INTEGER :: iostat, offset
  TYPE (cudaDeviceProp) :: prop
  iostat = cudaSetDevice(0)
  iostat = cudaGetDevice(j)
  iostat = cudaGetDeviceProperties(prop, j)
  write(*,"(' Device: ', a,/)") trim(prop%name)
  iostat = cudaRuntimeGetVersion(j)
  write(*,"(' CUDA Ver: ', I4,/)") j
  
  
  
  cpu_t  = 0
  cuda_t = 0
  
  iostat = cudaEventCreate(beg)
  iostat = cudaEventCreate(fin)

    DO i =1,16
	  iostat = cudaStreamCreate(streams(i))
	ENDDO

  
  ALLOCATE ( SC(NCRYS,SUP_IND), SC_d(NCRYS,2*SUP_IND))
  ALLOCATE ( FoMatRE(SUP_IND,9),   FoMatIM(SUP_IND,9), &
             FoMat_d(2*SUP_IND,9) )
  ALLOCATE ( StressMat(NCRYS,6), StressMat_d(NCRYS,9) )
  ALLOCATE ( Phi1(NCRYS), PHI(NCRYS), Phi2(NCRYS), termsdim4(NCRYS), &
             Phi1_d(NCRYS), PHI_d(NCRYS), Phi2_d(NCRYS), termsdim4_d(NCRYS) )
  ALLOCATE ( s_dot_d(NCRYS), W_app_pr_d(3,3), W_app_pr(3,3) )
  ALLOCATE ( sum_gamma_dot(NCRYS) )
  ALLOCATE ( Super_set_out(SUP_IND,4), Super_set_120_1(SUP_IND), Super_set_120_2(SUP_IND), Super_set_120_3(SUP_IND), Super_set_120_4(SUP_IND), &
             Super_set_1_d(SUP_IND), Super_set_2_d(SUP_IND), Super_set_3_d(SUP_IND), Super_set_4_d(SUP_IND) )
  !ALLOCATE ( QMATT(3,3,NCRYS), RMATT(3,3,NCRYS))

  ALLOCATE ( Q_p_sam_d(9))			 
  
  block_sz = dim3(blkSzX,blkSzY,1) ! Used in the build_SC CUDA kernel and values of dim3 are declared in Variables
  grid_sz  = dim3(grdSzX,grdSzY,1)
  blocks_texture = dim3(NCRYS/256 , 1 ,1)
  threads_texture = dim3(256,1,1)
  
  CALL CPU_TIME(ti) !Starts recording Preamble

  ones = 1._rp

  Total_time = Total_time*SQRT(2.0_rp)
  sizeF= (/ 120._rp, 120._rp, 120._rp, 120._rp /)

  
   L(1,1) =0.001d0
   L(1,2)=.0
   L(1,3)=.0
   L(2,1)=.0
   L(2,2)=.0
   L(2,3)=.0
   L(3,1)=.0
   L(3,2)=.0
   L(3,3)= -0.001d0
  
  OPEN (UNIT = 10, FILE = 'Fo_ns_11.txt', status='OLD')
  OPEN (UNIT = 20, FILE = 'Fo_ns_22.txt', status='OLD')
  OPEN (UNIT = 30, FILE = 'Fo_ss_12.txt', status='OLD')
  OPEN (UNIT = 40, FILE = 'Fo_ss_13.txt', status='OLD')
  OPEN (UNIT = 50, FILE = 'Fo_ss_23.txt', status='OLD')
  OPEN (UNIT = 60, FILE = 'Fo_ws_12.txt', status='OLD')
  OPEN (UNIT = 70, FILE = 'Fo_ws_13.txt', status='OLD')
  OPEN (UNIT = 80, FILE = 'Fo_ws_23.txt', status='OLD')
  OPEN (UNIT = 90, FILE = 'Fo_GD.txt', status='OLD')
  !change, Super_set_120
  OPEN (UNIT = 99, FILE = 'Super_set_120.txt', status='OLD')
  !change, Super_set_120

  DO i = 1,SUP_IND
     READ(10,*) Fo_ns_11(i)
     READ(20,*) Fo_ns_22(i)
     READ(30,*) Fo_ss_12(i)
     READ(40,*) Fo_ss_13(i)
     READ(50,*) Fo_ss_23(i)
     READ(60,*) Fo_ws_12(i)
     READ(70,*) Fo_ws_13(i)
     READ(80,*) Fo_ws_23(i)
     READ(90,*) Fo_GD(i)
     !change, Super_set_120_1 instead of Super_set_1 ...
     READ(99,*) Super_set_120_1(i), Super_set_120_2(i), Super_set_120_3(i), Super_set_120_4(i)
     !change, end 
  ENDDO

  CLOSE(UNIT=10)
  CLOSE(UNIT=20)
  CLOSE(UNIT=30)
  CLOSE(UNIT=40)
  CLOSE(UNIT=50)
  CLOSE(UNIT=60)
  CLOSE(UNIT=70)
  CLOSE(UNIT=80)
  CLOSE(UNIT=90)
  CLOSE(UNIT=99)

  !OPEN (UNIT = 33, FILE = 'euler.txt', status='OLD')
  !DO i = 1, NCRYS
  !   READ(33, *) Phi1(i), PHI(i), Phi2(i)
  !ENDDO
  !CLOSE(33)
    Phi1=0.0
    PHI=0.0
    Phi2=0.0
     !change, start, Super_set transformation according to m_fine
    Super_set_out(:,1)=Super_set_120_1
    Super_set_out(:,2)=Super_set_120_2
    Super_set_out(:,3)=Super_set_120_3
    Super_set_out(:,4)=Super_set_120_4
    do ij = 1,1024
        do jk = 1,4
            if (Super_set_out(ij,jk).gt.61) then
                Super_set_out(ij,jk) = Super_set_out(ij,jk) + m_fine*120 - 120
            end if
        end do
    end do


   !change, end
  
  D_app_sa = 0.5_rp * (L + TRANSPOSE(L))
  W_app_sa = 0.5_rp * (L - TRANSPOSE(L))

   !change, start (eigenvalues for new L)
   call Jacobi(D_app_sa, 3, D_eig, V_eig, nrot)

   D_app_sa = 0.5d0 * (L + transpose(L)) !comment: because Jacobi subroutine changes D_app_sa


   V_eig(1,3) =   V_eig(2,1)*V_eig(3,2) - V_eig(3,1)*V_eig(2,2)
   V_eig(2,3) = - V_eig(1,1)*V_eig(3,2) + V_eig(3,1)*V_eig(1,2)
   V_eig(3,3) =   V_eig(1,1)*V_eig(2,2) - V_eig(2,1)*V_eig(1,2)
print *,"V_eig",V_eig
print *,"D_eig",D_eig
  !change, end
  
  eps_dot = SQRT( SUM(D_eig*D_eig) ) 
  
  !
  Total_Strain = eps_dot*Total_Time
  DT = Strain_Increment/eps_dot
  n_steps =  NINT(Total_Time/DT)

  ALLOCATE(S11(1:n_steps))
  ALLOCATE(S22(1:n_steps))
  ALLOCATE(S12(1:n_steps))
  ALLOCATE(S13(1:n_steps))
  ALLOCATE(S23(1:n_steps))
  ALLOCATE(STRESS_VEC(1:n_steps,1:6))

  D_unit = D_eig/eps_dot
  
  DT = Total_Time/n_steps
  
  theta = DATAN2( -(-2.d0*D_unit(1) - D_unit(3)), -dsqrt(3.0D0)*D_unit(3) )
  IF ( theta < 0 ) THEN
     theta = theta + 2*PI
  ENDIF

  Q_p_sam = RESHAPE((V_eig), (/9/))
print *,"Q_p_sam",Q_p_sam
  W_app_pr = MATMUL(TRANSPOSE(V_eig), MATMUL(W_app_sa,V_eig))

  fact_1 = 2._rp / PRODUCT(sizeF)
  fact_2 = fact_1 * 1.e-02_rp / REAL(NCRYS,KIND=RP) * SIGN((eps_dot/DGo)**mo, eps_dot/DGo)
  fact_DT = fact_1 * eps_dot/DGo * DT
  fact_dot = fact_DT * ho / Ss**exp_a
  
	
  !to_rad = 1.e-03_rp * (PI / 180._rp)
  !to_deg = 1._rp / to_rad

  ! 
  !Phi2(:) = NINT( to_deg * Phi2(:) )
  ! PHI(:) = NINT( to_deg *  PHI(:) )
  !Phi1(:) = NINT( to_deg * Phi1(:) )

  CALL TEXTURE(Phi1, PHI, Phi2, 0, 0, 0, Q_p_sam, 0, n_steps)
   print *,"PHIS",Phi1(100), PHI(100),Phi2(100) !Good
  sizeF = sizeF * REAL(m_fine,KIND=RP)

  fact_s = 2._rp*PI / sizeF(1)

  !change, start
  term_s = NINT((180*theta/PI)*120*m_fine/360)
  termsdim4 = spread(term_s, 1, NCRYS)
  theta = term_s*360/120/m_fine*PI/180
  !change, end
  
  s_dot = So * ones

  s_dot_d = s_dot
  
  Super_set_out(:,1) = fact_s * ( Super_set_out(:,1)  - 1._rp )
  Super_set_out(:,2) = fact_s * ( Super_set_out(:,2)  - 1._rp )
  Super_set_out(:,3) = fact_s * ( Super_set_out(:,3)  - 1._rp )
  Super_set_out(:,4) = fact_s * ( Super_set_out(:,4)  - 1._rp )
  print *,"Super_set_1", Super_set_out(115,3)
  print *,"fact_s",fact_s
  CALL CPU_TIME(tf) ! End Preamble
  
  cpu_t(3) = tf-ti

  PRINT*, 'Building FoMat'
  FoMatRE = 0._rp
  FoMatIM = 0._rp
  FoMatRE(:,1) = Fo_ns_11
  FoMatRE(:,2) = Fo_ns_22
  FoMatRE(:,3) = Fo_ss_12
  FoMatRE(:,4) = Fo_ss_13
  FoMatRE(:,5) = Fo_ss_23
  FoMatRE(:,6) = Fo_ws_12
  FoMatRE(:,7) = Fo_ws_13
  FoMatRE(:,8) = Fo_ws_23
  FoMatRE(:,9) = Fo_GD

  FoMatIM(:,1) = AIMAG(Fo_ns_11)
  FoMatIM(:,2) = AIMAG(Fo_ns_22) 
  FoMatIM(:,3) = AIMAG(Fo_ss_12) 
  FoMatIM(:,4) = AIMAG(Fo_ss_13) 
  FoMatIM(:,5) = AIMAG(Fo_ss_23)
  FoMatIM(:,6) = AIMAG(Fo_ws_12) 
  FoMatIM(:,7) = AIMAG(Fo_ws_13) 
  FoMatIM(:,8) = AIMAG(Fo_ws_23)
  FoMatIM(:,9) = AIMAG(Fo_GD)
  
  FoMat_d(1:SUP_IND,1:9)           = FoMatRE
  FoMat_d(SUP_IND+1:2*SUP_IND,1:9) = FoMatIM
   
  CALL CPU_TIME(toti)
  
  iostat = cudaEventRecord(beg, 0)

	 !Initially toss Phi*_d onto device where it stays until the end
     Phi2_d = Phi2
     PHI_d  = PHI
     Phi1_d = Phi1
     W_app_pr_d = W_app_pr
     Q_p_sam_d = Q_p_sam
     Super_set_1_d = Super_set_out(:,1)
     Super_set_2_d = Super_set_out(:,2)
     Super_set_3_d = Super_set_out(:,3)
     Super_set_4_d = Super_set_out(:,4)
     termsdim4_d = termsdim4

     iostat = cudaEventRecord(fin, 0)
     iostat = cudaThreadSynchronize()
     iostat = cudaEventElapsedTime(tf, beg, fin)
     cuda_t(1) = cuda_t(1) + tf

    
    DO j = 1,n_steps
    
   iostat = cudaEventRecord(beg, 0)
	 !Launches a CUDA Kernel from CUDA_Kernels
     CALL build_SC<<<grid_sz,block_sz>>>(Super_set_1_d, Super_set_2_d, Super_set_3_d, Super_set_4_d, &
                                          Phi2_d, PHI_d, Phi1_d, termsdim4_d, SC_d, NCRYS, SUP_IND)
     

     iostat = cudaEventRecord(fin, 0)
     iostat = cudaThreadSynchronize()
     iostat = cudaEventElapsedTime(tf, beg, fin)
     cuda_t(2) = cuda_t(2) + tf

     iostat = cudaEventRecord(beg, 0)
	 
     !Launches a CUDA Kernel from CUDA_Kernels
	 DO i = 1,1
	 offset = (i-1)* NCRYS/1
     CALL mmul<<<dim3(NCRYS/16,1,1),dim3(16,9,1),0,streams(i)>>>( SC_d, FoMat_d(:,1:9), StressMat_d(:,1:9), NCRYS, 2*SUP_IND, 9)!, offset)
	 
		 iostat = cudaMemcpy2DAsync(StressMat(offset + 1,1), NCRYS , StressMat_d(offset + 1,1), NCRYS, NCRYS,5, cudaMemcpyDeviceToHost, streams(i)) 
		 iostat = cudaMemcpyAsync(StressMat(offset + 1,9) , StressMat_d(offset + 1,9) , NCRYS, streams(i))
	 
	 ENDDO
    
     iostat = cudaEventRecord(fin, 0)
     iostat = cudaThreadSynchronize() ! Makes sure all of mmul is complete and returned before moving on.
     iostat = cudaEventElapsedTime(tf, beg, fin)
     cuda_t(3) = cuda_t(3) + tf
     
     

	
	!WRITE(*,'(3F20.5)'), StressMat(1,9), StressMat(2,9), &
    !   StressMat(3,9) 
	
     sum_gamma_dot = StressMat(:,9)
     if(j.eq.1)print *,"sum_gamma_dot", sum_gamma_dot(100)
     call cpu_time(ti)

     s_dot = s_dot + fact_dot * ( Ss - s_dot )**exp_a * sum_gamma_dot

     call cpu_time(tf)
     cuda_t(5) = cuda_t(5) + tf-ti
     call cpu_time(ti)
     
     
     S11(j) = DOT_PRODUCT(StressMat(:,1),s_dot)
     S22(j) = DOT_PRODUCT(StressMat(:,2),s_dot)
     S12(j) = DOT_PRODUCT(StressMat(:,3),s_dot)
     S13(j) = DOT_PRODUCT(StressMat(:,4),s_dot)
     S23(j) = DOT_PRODUCT(StressMat(:,5),s_dot)

     call cpu_time(tf)
     cuda_t(6) = cuda_t(6) + tf-ti
     
     CALL CPU_TIME(ti)

     STRESS_TENSOR_PR = RESHAPE((/S11(j), S12(j), S13(j), &
                                  S12(j), S22(j), S23(j), &
                                  S13(j), S23(j), - ( S11(j) + S22(j) )/), (/3,3/))
     STRESS_TENSOR_SA = fact_2 * MATMUL(V_eig,MATMUL(STRESS_TENSOR_PR,TRANSPOSE(V_eig)))

     STRESS_VEC(j,1) = STRESS_TENSOR_SA(1,1)
     STRESS_VEC(j,2) = STRESS_TENSOR_SA(2,2)
     STRESS_VEC(j,3) = STRESS_TENSOR_SA(3,3)
     STRESS_VEC(j,4) = STRESS_TENSOR_SA(1,2)
     STRESS_VEC(j,5) = STRESS_TENSOR_SA(1,3) 
     STRESS_VEC(j,6) = STRESS_TENSOR_SA(2,3)
    
     CALL CPU_TIME(tf)

     cpu_t(2) = cpu_t(2) + tf-ti

	 
	 iostat = cudaEventRecord(beg, 0)
	 CALL QMAT<<<blocks_texture,threads_texture>>>(Phi2_d, PHI_d, Phi1_d)
		
    
	 CALL RMAT<<<blocks_texture,threads_texture>>>(StressMat_d,W_app_pr_d, fact_DT, DT)
	
				
	 CALL G_flag<<<blocks_texture,threads_texture>>>(Q_p_sam_d, Phi2_d, PHI_d, Phi1_d, j, n_steps)
	 
     iostat = cudaEventRecord(fin, 0) 
	 iostat = cudaThreadSynchronize()
	 iostat = cudaEventElapsedTime(tf, beg, fin)
	 cuda_t(7) = cuda_t(7) + tf 
	
	
		!WRITE(*,'(3F20.5)'), Phi1(1)/1000._rp, PHI(1)/1000._rp, &
       !   Phi2(1)/1000._rp 
  ENDDO
   
  
     iostat = cudaEventRecord(beg, 0)
	 
     Phi2 = Phi2_d
     PHI  = PHI_d
     Phi1 = Phi1_d

     iostat = cudaEventRecord(fin, 0)
     iostat = cudaThreadSynchronize()
     iostat = cudaEventElapsedTime(tf, beg, fin)
     cuda_t(1) = cuda_t(1) + tf
  
  
  
  CALL CPU_TIME(totf)
  
  OPEN (unit=55,file='stressSCP.dat')
  OPEN (unit=56,file='texture.dat')
  WRITE(55,'(7F11.5)') 0._rp,0._rp,0._rp,0._rp,0._rp,0._rp,0._rp
  DO I = 1,n_steps
     WRITE(*,'(7F11.5)')  REAL(i,KIND=RP)*1.0_rp/REAL(n_steps,KIND=RP),STRESS_VEC(i,:)
     WRITE(55,'(7F11.5)') REAL(i,KIND=RP)*1.0_rp/REAL(n_steps,KIND=RP),STRESS_VEC(i,:)
  ENDDO
  CLOSE (55)
  WRITE(56,*) 'TEXTURE AT STRAIN =    1.0000'
  WRITE(56,*) ' 0.947   1.221   0.864 '
  WRITE(56,*) '  90.00   90.00   90.00  '
  WRITE(56,*) 'B   8192'
  DO I = 1,NCRYS
     WRITE(56,'(4F20.5)'), Phi1(I)/1000._rp, PHI(I)/1000._rp, &
          Phi2(I)/1000._rp,1.0_rp/REAL(NCRYS,KIND=RP)
  END DO
  CLOSE (56)

  WRITE(*,'(A70)')               'Time (milliseconds) and Percentage Time                                '
  WRITE(*,'(A70)')               '======================================================================='
  WRITE(*,'(A39,4X,1PE13.4)')    'Preamble                              ', cpu_t(3)*1.e3_rp
  WRITE(*,'(A39,4X,2(1PE13.4))') 'Phi Transferred to Device             ', cuda_t(1), cuda_t(1)/(totf-toti)*1.e-1_rp
  WRITE(*,'(A39,4X,2(1PE13.4))') 'Build_SC Matrix Kernel                ', cuda_t(2), cuda_t(2)/(totf-toti)*1.e-1_rp
  WRITE(*,'(A39,4X,2(1PE13.4))') 'Mmul and StressMat Transfer           ', cuda_t(3), cuda_t(3)/(totf-toti)*1.e-1_rp
  WRITE(*,'(A39,4X,2(1PE13.4))') 'Texture Subroutine                    ', cuda_t(7), cuda_t(7)/(totf-toti)*1.e-1_rp
  WRITE(*,'(A39,4X,2(1PE13.4))') 'Sdot                                  ', cuda_t(5)*1.e3_rp, cuda_t(5)/(totf-toti)*1.e2_rp
  WRITE(*,'(A39,4X,2(1PE13.4))') 'Dot Product of Stress Matrix and Sdot ', cuda_t(6)*1.e3_rp, cuda_t(6)/(totf-toti)*1.e2_rp
  WRITE(*,'(A39,4X,2(1PE13.4))') 'Stress Tensor to Stress Vector        ', cpu_t(2)*1.e3_rp, cpu_t(2)/(totf-toti)*1.e2_rp
  WRITE(*,'(A70)')               '======================================================================='
  WRITE(*,'(A39,4X,1PE13.4)')    'Steps                                 ', (totf-toti)*1.e3_rp

  WRITE(*,'(A39,4X,1PE13.4)')    'Total                                 ', SUM(cpu_t)*1.e3_rp + SUM(cuda_t) 

  OPEN(unit=14, file='Timing.out')
  WRITE(14,'(A5,4X,1PE16.8)') 'Preamble ', cpu_t(3)*1.e3_rp 
  WRITE(14,'(A5,4X,2(1PE16.8))') 'Phi trans', cuda_t(1), cuda_t(1)/(totf-toti)*1.e-1_rp
  WRITE(14,'(A5,4X,2(1PE16.8))') 'SC Mat   ', cuda_t(2), cuda_t(2)/(totf-toti)*1.e-1_rp
  WRITE(14,'(A5,4X,2(1PE16.8))') 'MatMul   ', cuda_t(3), cuda_t(3)/(totf-toti)*1.e-1_rp
  WRITE(14,'(A5,4X,2(1PE16.8))') 'Str trans', cuda_t(4), cuda_t(4)/(totf-toti)*1.e-1_rp
  WRITE(14,'(A5,4X,2(1PE16.8))') 'Texture  ', cuda_t(7), cuda_t(7)/(totf-toti)*1.e-1_rp
  WRITE(14,'(A5,4X,2(1PE16.8))') 'Sdot     ', cuda_t(5), cuda_t(5)/(totf-toti)*1.e-1_rp
  WRITE(14,'(A5,4X,2(1PE16.8))') 'MatVec   ', cuda_t(6)*1.e3_rp, cuda_t(6)/(totf-toti)*1.e2_rp
  WRITE(14,'(A5,4X,2(1PE16.8))') 'Epilogue ', cpu_t(2)*1.e3_rp, cpu_t(2)/(totf-toti)*1.e2_rp
  WRITE(14,'(A5,4X,1PE16.8)') 'Steps    ', (totf-toti)*1.e3_rp
  WRITE(14,'(A5,4X,1PE16.8)') 'Total    ', SUM(cpu_t)*1.e3_rp + SUM(cuda_t) 
  CLOSE(14)
  

  iostat = cudaEventDestroy(beg)
  iostat = cudaEventDestroy(fin)
  DO i=1,16
     iostat = cudaStreamDestroy(streams(i))
  ENDDO
  

  DEALLOCATE( S11, S22, S12, S13, S23, STRESS_VEC )
  DEALLOCATE( SC, SC_d )
  DEALLOCATE( FoMatRE, FoMatIM, FoMat_d )
  DEALLOCATE( Phi1, PHI, Phi2, termsdim4, Phi1_d, PHI_d, Phi2_d, termsdim4_d )
  DEALLOCATE( StressMat, StressMat_d, s_dot_d )
  DEALLOCATE( sum_gamma_dot, W_app_pr_d, W_app_pr )
  DEALLOCATE( Super_set_out,Super_set_120_1, Super_set_120_2, Super_set_120_3, Super_set_120_4, &
              Super_set_1_d, Super_set_2_d, Super_set_3_d, Super_set_4_d )

  DEALLOCATE ( Q_p_sam_d)
  

   
END PROGRAM SCP_FCC


SUBROUTINE TEXTURE(Phi1, PHI, Phi2, w21_rec, w31_rec, w32_rec, Q_p_sam, flag, n_steps)

  USE Variables
  IMPLICIT NONE
  
  REAL(RP), DIMENSION(NCRYS) :: SPhi1, CPhi1, SPHI, CPHI, SPhi2, CPhi2
  REAL(RP), DIMENSION(NCRYS) :: w21_rec, w31_rec, w32_rec
  REAL(RP), DIMENSION(NCRYS,9) :: G
  REAL(RP) :: GG(9)
  REAL(RP) :: angles(NCRYS,3), Q_p_sam(9)!, to_rad,to_deg
  INTEGER, INTENT(IN) :: flag, n_steps
  REAL(RP), DIMENSION(NCRYS) :: QMAT_11, QMAT_12, QMAT_13, QMAT_21, QMAT_22, QMAT_23, QMAT_31, QMAT_32, QMAT_33
  REAL(RP), DIMENSION(NCRYS) :: R_11, R_12, R_13, R_21, R_22, R_23, R_31, R_32, R_33
  REAL(RP), DIMENSION(NCRYS) :: ang, axis_1, axis_2, axis_3
  REAL(RP), DIMENSION(NCRYS) :: Phi1, PHI, Phi2
  INTEGER :: i
  REAL(RP) :: twoPI

  twoPI = 2._rp*PI
  !Change start 
  to_rad = (PI / 180.d0)

  to_deg = 1.d0 / to_rad
  !change end
  angles = 0._rp

  ! The function receives the Euler angles in degrees and converts them to radians
  ! COORDINATE TRANSFORMATION MATRIX CRYSTAL TO SAMPLE

  SPhi1 = SIN(Phi1)
  CPhi1 = COS(Phi1)
  SPHI  = SIN(PHI)
  CPHI  = COS(PHI)
  SPhi2 = SIN(Phi2)
  CPhi2 = COS(Phi2)

  QMAT_11 = CPhi1*CPhi2-SPhi1*CPHI*SPhi2
  QMAT_12 = -CPhi1*SPhi2-SPhi1*CPHI*CPhi2
  QMAT_13 = SPhi1*SPHI 
  QMAT_21 = SPhi1*CPhi2+CPhi1*CPHI*SPhi2
  QMAT_22 = -SPhi1*SPhi2+CPhi1*CPHI*CPhi2
  QMAT_23 = -CPhi1*SPHI
  QMAT_31 = SPHI*SPhi2
  QMAT_32 = SPHI*CPhi2
  QMAT_33 = CPHI

  IF (flag == 0) THEN
     DO i = 1,NCRYS
        G(i,1) = Q_p_sam(1)*QMAT_11(i) + Q_p_sam(2)*QMAT_21(i) + Q_p_sam(3)*QMAT_31(i)
        G(i,2) = Q_p_sam(4)*QMAT_11(i) + Q_p_sam(5)*QMAT_21(i) + Q_p_sam(6)*QMAT_31(i)
        G(i,3) = Q_p_sam(7)*QMAT_11(i) + Q_p_sam(8)*QMAT_21(i) + Q_p_sam(9)*QMAT_31(i)
        G(i,4) = Q_p_sam(1)*QMAT_12(i) + Q_p_sam(2)*QMAT_22(i) + Q_p_sam(3)*QMAT_32(i)
        G(i,5) = Q_p_sam(4)*QMAT_12(i) + Q_p_sam(5)*QMAT_22(i) + Q_p_sam(6)*QMAT_32(i)
        G(i,6) = Q_p_sam(7)*QMAT_12(i) + Q_p_sam(8)*QMAT_22(i) + Q_p_sam(9)*QMAT_32(i)
        G(i,7) = Q_p_sam(1)*QMAT_13(i) + Q_p_sam(2)*QMAT_23(i) + Q_p_sam(3)*QMAT_33(i)
        G(i,8) = Q_p_sam(4)*QMAT_13(i) + Q_p_sam(5)*QMAT_23(i) + Q_p_sam(6)*QMAT_33(i)
        G(i,9) = Q_p_sam(7)*QMAT_13(i) + Q_p_sam(8)*QMAT_23(i) + Q_p_sam(9)*QMAT_33(i)
     END DO

     WHERE(ABS(G(:,9)) == 1)
        angles(:,1) = ACOS(G(:,1))
        angles(:,2) = ACOS(G(:,9))
        angles(:,3) = 0._rp

        WHERE(G(:,2) < 0)
           angles(:,1) = twoPI - angles(:,1)
        endwhere
     ELSEWHERE
        angles(:,1) = ATAN2(G(:,7),-1*G(:,8))
        angles(:,2) = ACOS(G(:,9))
        angles(:,3) = ATAN2(G(:,3),G(:,6))
     endwhere

     WHERE(angles(:,1) < 0)
        angles(:,1) = angles(:,1)+ twoPI
     endwhere
     WHERE(angles(:,2) < 0)
        angles(:,2) = angles(:,2)+ twoPI
     endwhere
     WHERE(angles(:,3) < 0)
        angles(:,3) = angles(:,3)+ twoPI
     endwhere

     !change, start,  transfering of euler angles in radians into r,s,t
     angles = NINT( to_deg*angles*120*m_fine/360 )
     !change, end

     Phi1(:) = angles(:,1)
     PHI(:)  = angles(:,2)
     Phi2(:) = angles(:,3)

     RETURN
  ENDIF
  
  !
  ang = SQRT(w21_rec**2 + w31_rec**2 + w32_rec**2)
  WHERE(ang == 0)
     axis_1 = 1._rp 
     axis_2 = 0._rp
     axis_3 = 0._rp
  ELSEWHERE
     axis_1 =  w32_rec/ang
     axis_2 = -w31_rec/ang
     axis_3 =  w21_rec/ang
  endwhere

  R_11 = (1._rp - axis_1**2)*COS(ang) + axis_1**2
  R_12 = axis_1*axis_2*(1._rp - COS(ang)) + axis_3*SIN(ang)
  R_13 = axis_1*axis_3*(1._rp - COS(ang)) - axis_2*SIN(ang)
  R_21 = axis_1*axis_2*(1._rp - COS(ang)) - axis_3*SIN(ang)
  R_22 = (1._rp - axis_2**2)*COS(ang) + axis_2**2
  R_23 = axis_2*axis_3*(1._rp - COS(ang)) + axis_1*SIN(ang)
  R_31 = axis_1*axis_3*(1._rp - COS(ang)) + axis_2*SIN(ang)
  R_32 = axis_2*axis_3*(1._rp - COS(ang)) - axis_1*SIN(ang)
  R_33 = (1._rp - axis_3**2)*COS(ang) + axis_3**2

  ! ind(1,:) is the first column of "g"; Vect_Rot is transpose so that is why they have the same indices. This is realy rotation

!!$  !omp parallel shared(G, R_11, R_12, R_13, R_31, R_32, R_33, R_21, R_22, R_23)&
!!$  !omp shared(QMAT_11, QMAT_12, QMAT_13, QMAT_31, QMAT_32, QMAT_33, QMAT_21, QMAT_22, QMAT_23) &
!!$  !omp private(i)
!!$  !omp do
!!$  DO i = 1,NCRYS
!!$     G(i,1) = R_11(i)*QMAT_11(i) + R_21(i)*QMAT_21(i) + R_31(i)*QMAT_31(i) 
!!$     G(i,2) = R_12(i)*QMAT_11(i) + R_22(i)*QMAT_21(i) + R_32(i)*QMAT_31(i)
!!$     G(i,3) = R_13(i)*QMAT_11(i) + R_23(i)*QMAT_21(i) + R_33(i)*QMAT_31(i)
!!$     G(i,4) = R_11(i)*QMAT_12(i) + R_21(i)*QMAT_22(i) + R_31(i)*QMAT_32(i)
!!$     G(i,5) = R_12(i)*QMAT_12(i) + R_22(i)*QMAT_22(i) + R_32(i)*QMAT_32(i)
!!$     G(i,6) = R_13(i)*QMAT_12(i) + R_23(i)*QMAT_22(i) + R_33(i)*QMAT_32(i)
!!$     G(i,7) = R_11(i)*QMAT_13(i) + R_21(i)*QMAT_23(i) + R_31(i)*QMAT_33(i)
!!$     G(i,8) = R_12(i)*QMAT_13(i) + R_22(i)*QMAT_23(i) + R_32(i)*QMAT_33(i)
!!$     G(i,9) = R_13(i)*QMAT_13(i) + R_23(i)*QMAT_23(i) + R_33(i)*QMAT_33(i)
!!$  END DO
!!$  !omp end do
!!$  !omp end parallel
!!$  
  G(:,1) = R_11*QMAT_11 + R_21*QMAT_21 + R_31*QMAT_31 
  G(:,2) = R_12*QMAT_11 + R_22*QMAT_21 + R_32*QMAT_31
  G(:,3) = R_13*QMAT_11 + R_23*QMAT_21 + R_33*QMAT_31
  G(:,4) = R_11*QMAT_12 + R_21*QMAT_22 + R_31*QMAT_32
  G(:,5) = R_12*QMAT_12 + R_22*QMAT_22 + R_32*QMAT_32
  G(:,6) = R_13*QMAT_12 + R_23*QMAT_22 + R_33*QMAT_32
  G(:,7) = R_11*QMAT_13 + R_21*QMAT_23 + R_31*QMAT_33
  G(:,8) = R_12*QMAT_13 + R_22*QMAT_23 + R_32*QMAT_33
  G(:,9) = R_13*QMAT_13 + R_23*QMAT_23 + R_33*QMAT_33

  ! In the last time increment transform the texture in the sample frame
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  IF (flag == n_steps) THEN

     DO i = 1,NCRYS
        GG(:) = G(i,:)
        G(i,1) = Q_p_sam(1)*GG(1) + Q_p_sam(4)*GG(2) + Q_p_sam(7)*GG(3)
        G(i,2) = Q_p_sam(2)*GG(1) + Q_p_sam(5)*GG(2) + Q_p_sam(8)*GG(3)
        G(i,3) = Q_p_sam(3)*GG(1) + Q_p_sam(6)*GG(2) + Q_p_sam(9)*GG(3)
        G(i,4) = Q_p_sam(1)*GG(4) + Q_p_sam(4)*GG(5) + Q_p_sam(7)*GG(6) 
        G(i,5) = Q_p_sam(2)*GG(4) + Q_p_sam(5)*GG(5) + Q_p_sam(8)*GG(6)
        G(i,6) = Q_p_sam(3)*GG(4) + Q_p_sam(6)*GG(5) + Q_p_sam(9)*GG(6)
        G(i,7) = Q_p_sam(1)*GG(7) + Q_p_sam(4)*GG(8) + Q_p_sam(7)*GG(9)
        G(i,8) = Q_p_sam(2)*GG(7) + Q_p_sam(5)*GG(8) + Q_p_sam(8)*GG(9)
        G(i,9) = Q_p_sam(3)*GG(7) + Q_p_sam(6)*GG(8) + Q_p_sam(9)*GG(9)
     END DO

  ENDIF
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !function angles = GRmatmult(newpos)

  !angles = zeros(size(newpos,1),3);
  !angles(:,1) = datan2(newpos(:,7),-1*newpos(:,8))
  !angles(:,2) = dacos(newpos(:,9))
  !angles(:,3) = datan2(newpos(:,3),newpos(:,6))

  WHERE(ABS(G(:,9)) == 1)
     angles(:,1) = ACOS(G(:,1))
     angles(:,2) = ACOS(G(:,9))
     angles(:,3) = 0._rp

     WHERE(G(:,2) < 0)
        angles(:,1) = twoPI - angles(:,1)
     endwhere
  ELSEWHERE
     angles(:,1) = ATAN2(G(:,7),-G(:,8))
     angles(:,2) = ACOS(G(:,9))
     angles(:,3) = ATAN2(G(:,3),G(:,6))
  endwhere
  
  !% Make sure angles are positive
  WHERE(angles(:,1) < 0)
     angles(:,1) = angles(:,1) + twoPI
  endwhere
  WHERE(angles(:,2) < 0)
     angles(:,2) = angles(:,2) + twoPI
  endwhere
  WHERE(angles(:,3) < 0)
     angles(:,3) = angles(:,3) + twoPI
  endwhere
  
    !change, start,  transfering of euler angles in radians into r,s,t,q
    angles = NINT( to_deg*angles*120*m_fine/360 )
    !change, end

  Phi1(:) = angles(:,1)
  PHI(:)  = angles(:,2)
  Phi2(:) = angles(:,3)

  RETURN
END SUBROUTINE TEXTURE

!*************************************************************
!* This subroutine computes all eigenvalues and eigenvectors *
!* of a real symmetric square matrix A(N,N). On output, ele- *
!* ments of A above the diagonal are destroyed. D(N) returns *
!* the eigenvalues of matrix A. V(N,N) contains, on output,  *
!* the eigenvectors of A by columns. THe normalization to    *
!* unity is made by main program before printing results.    *
!* NROT returns the number of Jacobi matrix rotations which  *
!* were required.                                            *
!* --------------------------------------------------------- *
!* Ref.:"NUMERICAL RECIPES, Cambridge University Press, 1986,*
!*       chap. 11, pages 346-348" [BIBLI 08].                *
!*************************************************************
Subroutine Jacobi(A,N,D,V,NROT)
integer N,NROT
real*8  A(1:N,1:N),D(1:N),V(1:N,1:N)
real*8, pointer :: B(:), Z(:)
real*8  c,g,h,s,sm,t,tau,theta,tresh

allocate(B(1:100),stat=ialloc)
allocate(Z(1:100),stat=ialloc)

  do ip=1, N    !initialize V to identity matrix
    do iq=1, N
      V(ip,iq)=0.d0 
    end do
      V(ip,ip)=1.d0
  end do  
  do ip=1, N
    B(ip)=A(ip,ip)
    D(ip)=B(ip)
    Z(ip)=0.d0    
  end do
  NROT=0
  do i=1, 50
    sm=0.d0
    do ip=1, N-1     !sum off-diagonal elements
      do iq=ip+1, N
        sm=sm+DABS(A(ip,iq))
      end do
    end do
    if(sm==0.d0) return  !normal return
    if(i.lt.4) then
      tresh=0.2d0*sm**2
    else
      tresh=0.d0
    end if
    do ip=1, N-1
      do iq=ip+1, N
        g=100.d0*DABS(A(ip,iq))
! after 4 sweeps, skip the rotation if the off-diagonal element is small
        if((i.gt.4).and.(DABS(D(ip))+g.eq.DABS(D(ip))) &
		.and.(DABS(D(iq))+g.eq.DABS(D(iq)))) then
		  A(ip,iq)=0.d0
        else if(DABS(A(ip,iq)).gt.tresh) then
	  h=D(iq)-D(ip)
	  if(DABS(h)+g.eq.DABS(h)) then
	    t=A(ip,iq)/h
          else
	    theta=0.5d0*h/A(ip,iq)  
            t=1.d0/(DABS(theta)+DSQRT(1.d0+theta**2))
	    if(theta.lt.0.d0) t=-t
          end if
	  c=1.d0/DSQRT(1.d0+t**2)
	  s=t*c
          tau=s/(1.d0+c)
	  h=t*A(ip,iq)
	  Z(ip)=Z(ip)-h
	  Z(iq)=Z(iq)+h
	  D(ip)=D(ip)-h
	  D(iq)=D(iq)+h
	  A(ip,iq)=0.d0
	  do j=1, ip-1
	    g=A(j,ip)
	    h=A(j,iq)
	    A(j,ip)=g-s*(h+g*tau)
	    A(j,iq)=h+s*(g-h*tau)
          end do
	  do j=ip+1, iq-1
	    g=A(ip,j)
	    h=A(j,iq)
	    A(ip,j)=g-s*(h+g*tau)
	    A(j,iq)=h+s*(g-h*tau)
          end do		      
	  do j=iq+1, N
	    g=A(ip,j)
	    h=A(iq,j)
	    A(ip,j)=g-s*(h+g*tau)
	    A(iq,j)=h+s*(g-h*tau)
          end do		  
	  do j=1, N
	    g=V(j,ip)
	    h=V(j,iq)
	    V(j,ip)=g-s*(h+g*tau)
	    V(j,iq)=h+s*(g-h*tau)
          end do		  
          NROT=NROT+1
        end if !if ((i.gt.4)...
      end do !main iq loop
    end do !main ip loop
    do ip=1, N
      B(ip)=B(ip)+Z(ip)
      D(ip)=B(ip)
      Z(ip)=0.d0
    end do
  end do !main i loop
  pause ' 50 iterations !'
  return
END    

! end of file ujacobi.f90
