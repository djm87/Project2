MODULE Variables

  IMPLICIT NONE
  
 INTEGER, PARAMETER, PUBLIC :: RP = KIND(0.D0)
  
  INTEGER, PARAMETER :: blkSz  = 16
  INTEGER, PARAMETER :: blkSzX = 16
  INTEGER, PARAMETER :: blkSzY = 16
  INTEGER :: grdSz, grdSzX, grdSzY
  
  INTEGER, PARAMETER :: SUP_IND = 1024, NCRYS = 8192, m_fine = 3000

  PARAMETER(grdSz  = NCRYS/blkSz)		
  PARAMETER(grdSzX = NCRYS/blkSzX)
  PARAMETER(grdSzY = SUP_IND/blkSzY)
  
  COMPLEX(RP), DIMENSION(SUP_IND) :: Fo_ns_11, Fo_ns_22, Fo_ss_12, Fo_ss_13, Fo_ss_23, &
                                     Fo_ws_12, Fo_ws_13, Fo_ws_23, Fo_GD
   
  REAL(RP), DIMENSION(3) :: D_eig, D_unit
  REAL(RP), DIMENSION(4) :: sizeF

  REAL(RP), DIMENSION(NCRYS) :: ones, s_dot
  REAL(RP), DIMENSION(NCRYS) :: stress_11, stress_22, stress_12, stress_13, stress_23
                        
  REAL(RP) :: fact_1, fact_2, fact_dot, to_rad,to_deg, fact_s, term_s  

  REAL(RP) :: euler(NCRYS,3), Q_p_sam(9) 
  REAL(RP) :: Total_Time = 1000._RP ! User enters originally 1000
  REAL(RP), PARAMETER:: So = 16.0_RP
  REAL(RP), PARAMETER:: ho = 180.0_RP
  REAL(RP), PARAMETER:: Ss = 148.0_RP
  REAL(RP), PARAMETER:: exp_a = 2.25_RP

  REAL(RP), PARAMETER:: Strain_Increment = 0.01_RP
  REAL(RP), PARAMETER:: mo = 0.01_RP ! The rate sensitivity parameter
  REAL(RP), PARAMETER:: DGo = 0.001_RP ! Reff. slip rate
  REAL(RP), PARAMETER:: s_ref = 100.0_RP ! Reff. slip resistance


  REAL(RP) :: PI, PI2
  PARAMETER( PI = ACOS(-1._RP) )
  PARAMETER( PI2 = PI * PI )


END MODULE Variables
