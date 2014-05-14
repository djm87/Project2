module SCP_VARIABLES
  implicit none

  integer, parameter :: SUP_IND = 1024, NCRYS = 1024, m_fine = 3000 !change, added m_fine
  double precision, parameter:: PI = ACOS(-1.0) 						 
  double complex, dimension(SUP_IND) :: Fo_ns_11, Fo_ns_22, Fo_ss_12, Fo_ss_13, Fo_ss_23, &
                                        Fo_ws_12, Fo_ws_13, Fo_ws_23, Fo_GD
  double complex, dimension(SUP_IND) :: Super_mat
  double precision, dimension(SUP_IND) :: Super_set_1, Super_set_2, Super_set_3, Super_set_4
  double precision, dimension(SUP_IND) :: Super_set_120_1, Super_set_120_2, Super_set_120_3, Super_set_120_4!change, added this row
  double precision, dimension(SUP_IND,4) :: Super_set_out !change, added this row
  double complex :: ii
  DOUBLE PRECISION :: euler(NCRYS,3), Q_p_sam(9)
  double precision :: Total_Time = 1000 ! User enters
  double precision, parameter:: So = 16.0
  double precision, parameter:: ho = 180.0
  double precision, parameter:: Ss = 148.0
  double precision, parameter:: exp_a = 2.25

  double precision, parameter:: Strain_Increment = 0.01
  double precision, parameter:: mo = 0.01 ! The rate sensitivity parameter
  double precision, parameter:: DGo = 0.001 ! Reff. slip rate
  double precision, parameter:: s_ref = 100.0 ! Reff. slip resistance
end module SCP_VARIABLES

program SCP_FCC
  use SCP_VARIABLES
  implicit none
real :: t, t_1, t_2, t_3
  integer, PARAMETER::  NMAX = 3
  double precision, dimension(NCRYS) ::  Phi1, PHI, Phi2, termsdim4
  INTEGER          ind_max, ind_min, ind_mid, n_steps, i, j, IMODE, III, klm, nrot !change, added klm, nrot
  DOUBLE PRECISION D_eig(3), D_unit(3), eps_dot, DT, theta, Total_Strain, sizeF(4)
  INTRINSIC        ABS, MAX

  double precision, dimension(NCRYS) :: ones, s_dot, sum_gamma_dot, &
                                           stress_11, stress_22, stress_12, stress_13, stress_23, &
                                           w_12, w_13, w_23
  double precision, dimension(SUP_IND) :: SC

  integer :: ICRYS, JSET, k, p,ij, jk,mk,lk !change, added  k, p, ij, jk,mk,lk	
  double precision, dimension(3,3) :: L, D_app_sa, W_app_sa, V_eig, W_app_pr, &
                                      STRESS_TENSOR_PR, STRESS_TENSOR_SA

  double precision, dimension(:), ALLOCATABLE :: S11, S22, S12, S13, S23, strain
  double precision, dimension(:,:), ALLOCATABLE :: STRESS_VEC

  double precision fact_1, fact_2, fact_DT, fact_dot, to_rad,to_deg, fact_s, term_s


  ii = (0.d0,1.d0)

  ones = 1.d0

  Total_time = Total_time*dsqrt(2.0D0)
  sizeF= (/120.0, 120.0, 120.0, 120.0/)
  
   L(1,1) =0.0005d0
   L(1,2)=.0
   L(1,3)=.0
   L(2,1)=.0
   L(2,2)=0.0005d0
   L(2,3)=.0
   L(3,1)=.0
   L(3,2)=.0
   L(3,3)= -0.001d0

  OPEN (UNIT = 1, FILE = 'Fo_ns_11.txt', status='OLD')
  OPEN (UNIT = 2, FILE = 'Fo_ns_22.txt', status='OLD')
  OPEN (UNIT = 3, FILE = 'Fo_ss_12.txt', status='OLD')
  OPEN (UNIT = 4, FILE = 'Fo_ss_13.txt', status='OLD')
  OPEN (UNIT = 5, FILE = 'Fo_ss_23.txt', status='OLD')
  OPEN (UNIT = 66, FILE = 'Fo_ws_12.txt', status='OLD')
  OPEN (UNIT = 7, FILE = 'Fo_ws_13.txt', status='OLD')
  OPEN (UNIT = 8, FILE = 'Fo_ws_23.txt', status='OLD')
  OPEN (UNIT = 9, FILE = 'Fo_GD.txt', status='OLD')
  !change, Super_set_120
  OPEN (UNIT = 10, FILE = 'Super_set_120.txt', status='OLD')
  !change, Super_set_120
 
  DO i = 1,SUP_IND
     READ(1,*) Fo_ns_11(i)
     READ(2,*) Fo_ns_22(i)
     READ(3,*) Fo_ss_12(i)
     READ(4,*) Fo_ss_13(i)
     READ(5,*) Fo_ss_23(i)
     READ(66,*) Fo_ws_12(i)
     READ(7,*) Fo_ws_13(i)
     READ(8,*) Fo_ws_23(i)
     READ(9,*) Fo_GD(i)
     !change, Super_set_120_1 instead of Super_set_1 ...
     READ(10,*) Super_set_120_1(i), Super_set_120_2(i), Super_set_120_3(i), Super_set_120_4(i)
     !change, end 
  ENDDO

     Fo_ns_11 = dconjg( Fo_ns_11 )
     Fo_ns_22 = dconjg( Fo_ns_22 )
     Fo_ss_12 = dconjg( Fo_ss_12 )
     Fo_ss_13 = dconjg( Fo_ss_13 )
     Fo_ss_23 = dconjg( Fo_ss_23 )
     Fo_ws_12 = dconjg( Fo_ws_12 )
     Fo_ws_13 = dconjg( Fo_ws_13 )
     Fo_ws_23 = dconjg( Fo_ws_23 )
     Fo_GD = dconjg( Fo_GD )

  close(UNIT=1)
  close(UNIT=2)
  close(UNIT=3)
  close(UNIT=4)
  close(UNIT=5)
  close(UNIT=66)
  close(UNIT=7)
  close(UNIT=8)
  close(UNIT=9)
  close(UNIT=10)
  

  OPEN (UNIT = 33, FILE = 'euler.txt', status='OLD')
  DO i = 1, NCRYS
     READ(33, *) Phi1(i), PHI(i), Phi2(i)
  ENDDO
  CLOSE(33)

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
  
  D_app_sa = 0.5d0 * (L + transpose(L))
  W_app_sa = 0.5d0 * (L - transpose(L))

!change, start (eigenvalues for new L)
call Jacobi(D_app_sa, 3, D_eig, V_eig, nrot)

D_app_sa = 0.5d0 * (L + transpose(L)) !comment: because Jacobi subroutine changes D_app_sa


V_eig(1,3) =   V_eig(2,1)*V_eig(3,2) - V_eig(3,1)*V_eig(2,2)
V_eig(2,3) = - V_eig(1,1)*V_eig(3,2) + V_eig(3,1)*V_eig(1,2)
V_eig(3,3) =   V_eig(1,1)*V_eig(2,2) - V_eig(2,1)*V_eig(1,2)

!change, end


  eps_dot = dsqrt(sum(D_eig*D_eig)) 
!
  Total_Strain = eps_dot*Total_Time
  DT = Strain_Increment/eps_dot
  n_steps = NINT(Total_Time/DT)

  ALLOCATE(S11(1:n_steps))
  ALLOCATE(S22(1:n_steps))
  ALLOCATE(S12(1:n_steps))
  ALLOCATE(S13(1:n_steps))
  ALLOCATE(S23(1:n_steps))
  ALLOCATE(STRESS_VEC(1:n_steps,1:6))
  ALLOCATE(strain(n_steps))

  D_unit = D_eig/eps_dot

  DT = Total_Time/n_steps

  !change, start, added minuses
  theta = datan2( -(-2.d0*D_unit(1) - D_unit(3)), -dsqrt(3.0D0)*D_unit(3) )
  !change, end
  if (theta<0) then
     !change, start, added 2 pi 
     theta = theta + 2*PI
     !change, end
  endif

  Q_p_sam = reshape((V_eig), (/9/))
print *,"Q_p_sam",Q_p_sam
  W_app_pr = matmul(transpose(V_eig), matmul(W_app_sa,V_eig))

  fact_1 = 2.d0 / product(sizeF)
  fact_2 = fact_1 * 1.d-02 / DBLE(NCRYS) * sign((eps_dot/DGo)**mo, eps_dot/DGo)
  fact_DT = fact_1 * eps_dot/DGo * DT
  fact_dot = fact_DT * ho / Ss**exp_a

 
  CALL TEXTURE(Phi1, PHI, Phi2, 0, 0, 0, 0, n_steps) ! 
 print *,"PHIS",Phi1(100), PHI(100),Phi2(100)
  sizeF = sizeF * DBLE(m_fine)

  fact_s = 2.d0*PI / sizeF(1)
   Super_set_1=fact_s * (Super_set_out(:,1)  - 1.d0 ) 
   Super_set_2=fact_s * (Super_set_out(:,2)  - 1.d0 )  
   Super_set_3=fact_s * (Super_set_out(:,3)  - 1.d0 ) 
   Super_set_4=fact_s * (Super_set_out(:,4)  - 1.d0 ) 
   print *,"Super_set_1", Super_set_3(115)
   print *,"fact_s",fact_s
  !change, start
  term_s = NINT((180*theta/PI)*120*m_fine/360)
  termsdim4 = spread(term_s, 1, NCRYS)
  theta = term_s*360/120/m_fine*PI/180
  !change, end
 
  s_dot = So * ones

 call cpu_time(t)
t_1 = t

  do j = 1,n_steps

     DO i = 1,NCRYS
         
        !change, start
        SC = Super_set_1 * Phi2(i) &
             + Super_set_2 * PHI(i) &
             + Super_set_3 * Phi1(i) &
             + Super_set_4 * termsdim4(i) 
        !change, end
        
        Super_mat = dcmplx( dcos( SC ), dsin( SC ) )

        sum_gamma_dot(i) = dot_product( Fo_GD,Super_mat )
		if ((j.eq.1).and.(i.eq.100)) print *,"sum_gamma_dot", sum_gamma_dot(100)
        stress_11(i) = dot_product( Fo_ns_11,Super_mat ) 
        stress_22(i) = dot_product( Fo_ns_22,Super_mat )
        stress_12(i) = dot_product( Fo_ss_12,Super_mat )
        stress_13(i) = dot_product( Fo_ss_13,Super_mat )
        stress_23(i) = dot_product( Fo_ss_23,Super_mat )

             w_12(i) = dot_product( Fo_ws_12,Super_mat )
             w_13(i) = dot_product( Fo_ws_13,Super_mat )
             w_23(i) = dot_product( Fo_ws_23,Super_mat )
     end do


     call TEXTURE (Phi1, PHI, Phi2, w_12*fact_DT + W_app_pr(2,1)*DT, &
                                    w_13*fact_DT + W_app_pr(3,1)*DT, &
                                    w_23*fact_DT + W_app_pr(3,2)*DT, j, n_steps)

     s_dot = s_dot + fact_dot * ( Ss - s_dot )**exp_a * sum_gamma_dot

     S11(j) = dot_product(stress_11,s_dot)
     S22(j) = dot_product(stress_22,s_dot)
     S12(j) = dot_product(stress_12,s_dot)
     S13(j) = dot_product(stress_13,s_dot)
     S23(j) = dot_product(stress_23,s_dot)

     STRESS_TENSOR_PR = RESHAPE((/S11(j), S12(j), S13(j), &
                                  S12(j), S22(j), S23(j), &
                                  S13(j), S23(j), - ( S11(j) + S22(j) )/), (/3,3/))
     STRESS_TENSOR_SA = fact_2 * matmul(V_eig,matmul(STRESS_TENSOR_PR,transpose(V_eig)))

     STRESS_VEC(j,1) = STRESS_TENSOR_SA(1,1)
     STRESS_VEC(j,2) = STRESS_TENSOR_SA(2,2)
     STRESS_VEC(j,3) = STRESS_TENSOR_SA(3,3)
     STRESS_VEC(j,4) = STRESS_TENSOR_SA(1,2)
     STRESS_VEC(j,5) = STRESS_TENSOR_SA(1,3)
     STRESS_VEC(j,6) = STRESS_TENSOR_SA(2,3)

  enddo



		call cpu_time(t)
t_2 = t-t_1
print*, t_2
!
  open (unit=55,file='stressSCP.dat')
  open (unit=56,file='texture.dat')
     write(55,'(7F11.5)') 0.,0.,0.,0.,0.,0.,0.
  DO I = 1,n_steps
     !write(*,'(7F11.5)') dble(i)*1.0d0/dble(n_steps),STRESS_VEC(i,:)
     write(55,'(7F11.5)') dble(i)*1.0d0/dble(n_steps),STRESS_VEC(i,:)
  enddo
  close (55)
 write(56,*) 'TEXTURE AT STRAIN =    1.0000'
 write(56,*) ' 0.947   1.221   0.864 '
 write(56,*) '  90.00   90.00   90.00  '
 write(56,*) 'B   1024'
        DO III = 1,NCRYS
 	  write(56,'(4F20.5)'), Phi1(III)*360/120/m_fine, PHI(III)*360/120/m_fine, Phi2(III)*360/120/m_fine,1.0d0/dble(NCRYS)
	END DO
print *,"(Phi1, PHI, Phi2) = ", Phi1(1:16), PHI(1:16), Phi2(1:16)
 close (56)

  DEALLOCATE(S11,S22,S12,S13,S23, STRESS_VEC, strain)

end program SCP_FCC


SUBROUTINE TEXTURE(Phi1, PHI, Phi2, w21_rec, w31_rec, w32_rec, flag, n_steps)

  use SCP_VARIABLES
  implicit none
  double precision, dimension(NCRYS) :: SPhi1, CPhi1, SPHI, CPHI, SPhi2, CPhi2
  double precision, dimension(NCRYS) :: w21_rec, w31_rec, w32_rec
  double precision, dimension(NCRYS,9) :: G
  double precision GG(9)
  double precision angles(NCRYS,3)
  integer, intent(in) :: flag, n_steps
  double precision, dimension(NCRYS) :: QMAT_11, QMAT_12, QMAT_13, QMAT_21, QMAT_22, QMAT_23, QMAT_31, QMAT_32, QMAT_33
  double precision, dimension(NCRYS) :: R_11, R_12, R_13, R_21, R_22, R_23, R_31, R_32, R_33
  double precision, dimension(NCRYS) :: ang, axis_1, axis_2, axis_3
  integer :: i
  double precision, dimension(NCRYS) :: Phi1, PHI, Phi2

  double precision twoPI, to_rad,to_deg

  to_rad = (PI / 180.d0)
  to_deg = 1.d0 / to_rad
    !change, start, transfering of r,s,t into euler angles in radians
    if (flag.ne.0) then
        Phi1 = Phi1*360/120/m_fine;
        Phi1 = Phi1*to_rad;
        PHI = PHI*360/120/m_fine;
        PHI = PHI*to_rad;
        Phi2 = Phi2*360/120/m_fine;
        Phi2 = Phi2*to_rad;
    end if
    !change, end
  twoPI = 2.d0*PI
  angles = 0.d0
  


  SPhi1 = DSIN(Phi1)
  CPhi1 = DCOS(Phi1)
  SPHI  = DSIN(PHI)
  CPHI  = DCOS(PHI)
  SPhi2 = DSIN(Phi2)
  CPhi2 = DCOS(Phi2)

  QMAT_11 = CPhi1*CPhi2-SPhi1*CPHI*SPhi2
  QMAT_12 = -CPhi1*SPhi2-SPhi1*CPHI*CPhi2
  QMAT_13 = SPhi1*SPHI 
  QMAT_21 = SPhi1*CPhi2+CPhi1*CPHI*SPhi2
  QMAT_22 = -SPhi1*SPhi2+CPhi1*CPHI*CPhi2
  QMAT_23 = -CPhi1*SPHI
  QMAT_31 = SPHI*SPhi2
  QMAT_32 = SPHI*CPhi2
  QMAT_33 = CPHI

  if (flag == 0) then
     do i = 1,NCRYS
        G(i,1) = Q_p_sam(1)*QMAT_11(i) + Q_p_sam(2)*QMAT_21(i) + Q_p_sam(3)*QMAT_31(i)
        G(i,2) = Q_p_sam(4)*QMAT_11(i) + Q_p_sam(5)*QMAT_21(i) + Q_p_sam(6)*QMAT_31(i)
        G(i,3) = Q_p_sam(7)*QMAT_11(i) + Q_p_sam(8)*QMAT_21(i) + Q_p_sam(9)*QMAT_31(i)
        G(i,4) = Q_p_sam(1)*QMAT_12(i) + Q_p_sam(2)*QMAT_22(i) + Q_p_sam(3)*QMAT_32(i)
        G(i,5) = Q_p_sam(4)*QMAT_12(i) + Q_p_sam(5)*QMAT_22(i) + Q_p_sam(6)*QMAT_32(i)
        G(i,6) = Q_p_sam(7)*QMAT_12(i) + Q_p_sam(8)*QMAT_22(i) + Q_p_sam(9)*QMAT_32(i)
        G(i,7) = Q_p_sam(1)*QMAT_13(i) + Q_p_sam(2)*QMAT_23(i) + Q_p_sam(3)*QMAT_33(i)
        G(i,8) = Q_p_sam(4)*QMAT_13(i) + Q_p_sam(5)*QMAT_23(i) + Q_p_sam(6)*QMAT_33(i)
        G(i,9) = Q_p_sam(7)*QMAT_13(i) + Q_p_sam(8)*QMAT_23(i) + Q_p_sam(9)*QMAT_33(i)
     end do

     where(abs(G(:,9)) == 1)
        angles(:,1) = dacos(G(:,1))
        angles(:,2) = dacos(G(:,9))
        angles(:,3) = 0.

        where(G(:,2) < 0)
           angles(:,1) = twoPI - angles(:,1)
        endwhere
     elsewhere
        angles(:,1) = datan2(G(:,7),-1*G(:,8))
        angles(:,2) = dacos(G(:,9))
        angles(:,3) = datan2(G(:,3),G(:,6))
     endwhere

     where(angles(:,1) < 0)
        angles(:,1) = angles(:,1)+ twoPI
     endwhere
     where(angles(:,2) < 0)
        angles(:,2) = angles(:,2)+ twoPI
     endwhere
     where(angles(:,3) < 0)
        angles(:,3) = angles(:,3)+ twoPI
     endwhere
     !change, start,  transfering of euler angles in radians into r,s,t
     angles = NINT( to_deg*angles*120*m_fine/360 )
     !change, end

     Phi1(:) = angles(:,1)
     PHI(:)  = angles(:,2)
     Phi2(:) = angles(:,3)

     return
  endif

  !
  ang = dsqrt(w21_rec**2 + w31_rec**2 + w32_rec**2)
  where(ang == 0)
     axis_1 = 1.d0 
     axis_2 = 0.
     axis_3 = 0.
  elsewhere
     axis_1 =  w32_rec/ang
     axis_2 = -w31_rec/ang
     axis_3 =  w21_rec/ang
  endwhere

  R_11 = (1.d0 - axis_1**2)*dcos(ang) + axis_1**2
  R_12 = axis_1*axis_2*(1.d0 - dcos(ang)) + axis_3*dsin(ang)
  R_13 = axis_1*axis_3*(1.d0 - dcos(ang)) - axis_2*dsin(ang)
  R_21 = axis_1*axis_2*(1.d0 - dcos(ang)) - axis_3*dsin(ang)
  R_22 = (1.d0 - axis_2**2)*dcos(ang) + axis_2**2
  R_23 = axis_2*axis_3*(1.d0 - dcos(ang)) + axis_1*dsin(ang)
  R_31 = axis_1*axis_3*(1.d0 - dcos(ang)) + axis_2*dsin(ang)
  R_32 = axis_2*axis_3*(1.d0 - dcos(ang)) - axis_1*dsin(ang)
  R_33 = (1.d0 - axis_3**2)*dcos(ang) + axis_3**2

  ! ind(1,:) is the first column of "g"; Vect_Rot is transpose so that is why they have the same indices. This is realy rotation
  do i = 1,NCRYS
     G(i,1) = R_11(i)*QMAT_11(i) + R_21(i)*QMAT_21(i) + R_31(i)*QMAT_31(i) 
     G(i,2) = R_12(i)*QMAT_11(i) + R_22(i)*QMAT_21(i) + R_32(i)*QMAT_31(i)
     G(i,3) = R_13(i)*QMAT_11(i) + R_23(i)*QMAT_21(i) + R_33(i)*QMAT_31(i)
     G(i,4) = R_11(i)*QMAT_12(i) + R_21(i)*QMAT_22(i) + R_31(i)*QMAT_32(i)
     G(i,5) = R_12(i)*QMAT_12(i) + R_22(i)*QMAT_22(i) + R_32(i)*QMAT_32(i)
     G(i,6) = R_13(i)*QMAT_12(i) + R_23(i)*QMAT_22(i) + R_33(i)*QMAT_32(i)
     G(i,7) = R_11(i)*QMAT_13(i) + R_21(i)*QMAT_23(i) + R_31(i)*QMAT_33(i)
     G(i,8) = R_12(i)*QMAT_13(i) + R_22(i)*QMAT_23(i) + R_32(i)*QMAT_33(i)
     G(i,9) = R_13(i)*QMAT_13(i) + R_23(i)*QMAT_23(i) + R_33(i)*QMAT_33(i)
  end do

  ! In the last time increment transform the texture in the sample frame
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if (flag == n_steps) then

     do i = 1,NCRYS
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
     end do

  endif
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !function angles = GRmatmult(newpos)

  !angles = zeros(size(newpos,1),3);
  !angles(:,1) = datan2(newpos(:,7),-1*newpos(:,8))
  !angles(:,2) = dacos(newpos(:,9))
  !angles(:,3) = datan2(newpos(:,3),newpos(:,6))

  where(abs(G(:,9)) == 1)
     angles(:,1) = dacos(G(:,1))
     angles(:,2) = dacos(G(:,9))
     angles(:,3) = 0.0

     where(G(:,2) < 0)
        angles(:,1) = twoPI - angles(:,1)
     endwhere
  elsewhere
     angles(:,1) = datan2(G(:,7),-G(:,8))
     angles(:,2) = dacos(G(:,9))
     angles(:,3) = datan2(G(:,3),G(:,6))
  endwhere

  !% Make sure angles are positive
  where(angles(:,1) < 0)
     angles(:,1) = angles(:,1) + twoPI
  endwhere
  where(angles(:,2) < 0)
     angles(:,2) = angles(:,2) + twoPI
  endwhere
  where(angles(:,3) < 0)
     angles(:,3) = angles(:,3) + twoPI
  endwhere


    !change, start,  transfering of euler angles in radians into r,s,t,q
    angles = NINT( to_deg*angles*120*m_fine/360 )
    !change, end
    IF(flag==n_steps) THEN
	!print *,"angles1 = ",angles(5,1),"angles2 = ",angles(5,2),"angles3 = ",angles(5,3)
	!print *,"angles1 = ",angles(6,1),"angles2 = ",angles(6,2),"angles3 = ",angles(6,3)
	!print *,"angles1 = ",angles(7,1),"angles2 = ",angles(7,2),"angles3 = ",angles(7,3)
	!pause
    ENDIF
  Phi1(:) = angles(:,1)
  PHI(:)  = angles(:,2)
  Phi2(:) = angles(:,3)

 ! print *,"Phi1(6) = ",Phi1(6), "Phi1(8)= ", Phi1(8), "Itter= ", flag
 ! pause
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
