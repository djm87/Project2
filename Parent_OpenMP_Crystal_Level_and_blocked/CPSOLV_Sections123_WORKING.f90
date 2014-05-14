MODULE CPSOLV
    CONTAINS
    SUBROUTINE CPSOLV(ISOLER,DTIME,C6A,SLPRA1,SLPRA2, &
                    SSYSMAT,KUPD,STRS1PK, ISTEP)
    USE Variables
  !  USE advisor_annotate

    IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)

    DIMENSION SLPRESA2(NSLIP,NC),DGDTPKA(NSLIP), &
            EROR(6,NC),TJINV(6,6,NC),DELTPK(6,NC), &
            ABSDELTPK(6),ABSDG(NSLIP),SLIPSIGN(NSLIP), &
            SLPRESA2TMP(NSLIP),HB(NSLIP),TPKOLD(6,NC),EYE(6,6)
    DIMENSION SLPRA1(NSLIP,NC),SLPRA2(NSLIP,NC), &
            C6A(6,NSLIP,NC),SSYSMAT(3,3,NSLIP,NC),STRS1PK(6,*),INDX(6), &
            IFLAG(NC),IITR(NC),C6P6A(6,6,NSLIP,NC),ICORR(NC) 

    IFLAG = 0
    IITR = 0 
    
    EYE = 0.0D0
    DO I = 1,6
		EYE(I,I) = 1.0D0
    ENDDO
    DGDO = GDO*DTIME 
    EXPN= 1.0D0/EXPM 
    DMAX = 0.05D0

 	!$OMP PARALLEL DEFAULT(PRIVATE) &
    !$omp shared(SLPRESA2,SLPRA1,C6P6A,C6A,SSYSMAT)  
    !$omp do 
    DO IC = 1,NC !IC SECTION 1

          !INITIALISE SLIP RESISTANCE AND GET CONSTANT C6P6A
		  DO I = 1,NSLIP
			  SLPRESA2(I,IC)= SLPRA1(I,IC) 
			  DO II = 1,6 
				C6P6A(II,1,I,IC) = C6A(II,I,IC)*SSYSMAT(1,1,I,IC)            
				C6P6A(II,2,I,IC) = C6A(II,I,IC)*SSYSMAT(2,2,I,IC)
				C6P6A(II,3,I,IC) = C6A(II,I,IC)*SSYSMAT(3,3,I,IC)
				C6P6A(II,4,I,IC) = C6A(II,I,IC)*(SSYSMAT(1,2,I,IC)+ &
											SSYSMAT(2,1,I,IC))
				C6P6A(II,5,I,IC) = C6A(II,I,IC)*(SSYSMAT(2,3,I,IC)+ &
										   SSYSMAT(3,2,I,IC))
				C6P6A(II,6,I,IC) = C6A(II,I,IC)*(SSYSMAT(1,3,I,IC)+ &
											SSYSMAT(3,1,I,IC))
			  ENDDO
		  ENDDO   
    ENDDO !IC SECTION 1
 !$omp end do
 !$omp end parallel
  			
				!print *,"C6P6A =",C6P6A(:,1,1,1)
				!print *,"SSYSMAT =",SSYSMAT(:,1,1,1)
				!print *,"C6A =",C6A(:,1,1)
				!print *,"DG2 =",DG2(:,IC)
				!print *,"IITR =",IITR(IC)		
			!pause


   !!START OF SOLUTION PROCEDURE

    111 CONTINUE
	!IF (ISTEP .eq. 1) THEN
	!print *,"Stress2be(:,2)", STRESS2BE(:,2)
	!pause
	!endif
    IF (ALL(IFLAG .eq. 1)) THEN 

    ELSE
 	!$OMP PARALLEL DEFAULT(PRIVATE) &
    !$omp shared(IFLAG,SSYSMAT,STRESS2BE,DG2,DGDO,EXPN,DMAX,SLPRESA2,IITR,DELTPK,TPKOLD,EROR,EYE,STRV,C6A,C6P6A,TJINV,ICORR)  
    !$omp do 
        DO IC = 1, NC  !IC SECTION 2

            IF (IFLAG(IC) .eq. 1) THEN
            ELSE
                222 CONTINUE

			    DRFCE2 = 0.0D0 
			    DRFCE2(:) = DRFCE2(:)+STRESS2BE(1,IC)*SSYSMAT(1,1,:,IC)
			    DRFCE2(:) = DRFCE2(:)+STRESS2BE(2,IC)*SSYSMAT(2,2,:,IC)
			    DRFCE2(:) = DRFCE2(:)+STRESS2BE(3,IC)*SSYSMAT(3,3,:,IC)
			    DRFCE2(:) = DRFCE2(:)+STRESS2BE(4,IC)*(SSYSMAT(1,2,:,IC)+SSYSMAT(2,1,:,IC))
			    DRFCE2(:) = DRFCE2(:)+STRESS2BE(5,IC)*(SSYSMAT(2,3,:,IC)+SSYSMAT(3,2,:,IC))
			    DRFCE2(:) = DRFCE2(:)+STRESS2BE(6,IC)*(SSYSMAT(1,3,:,IC)+SSYSMAT(3,1,:,IC))
                
			    SLIPSIGN(:) = SIGN(1.0D0,DRFCE2(:))
		
		        ABSDGMAX = 0.0D0
		        DO I = 1,NSLIP          
			        ABSDG(I) = 0.0D0
			        DG2(I,IC) = 0.0D0
			        IF(DRFCE2(I).NE.0.0D0) THEN
			            ABSDG(I) = DGDO*DABS(DRFCE2(I)/SLPRESA2(I,IC))**EXPN
			            DG2(I,IC) = ABSDG(I)*SLIPSIGN(I)
			            ABSDGMAX = DMAX1(ABSDGMAX,ABSDG(I))
			        ENDIF
                ENDDO
                	!IF(IC .eq. 2) THEN			
				!print *,"ICORR =",ICORR(IC)
				!print *,"ABSDGMAX =",ABSDGMAX
				!print *,"DMAX =",DMAX
				!print *,"DG2 =",DG2(:,IC)
				!print *,"IITR =",IITR(IC)		
			!pause
		!ENDIF
                IF (ABSDGMAX.GT.DMAX.AND.IITR(IC).NE.0.AND.ICORR(IC).LE.30)THEN
                    DELTPK(:,IC) = DELTPK(:,IC)*0.25
                    STRESS2BE(:,IC) = TPKOLD(:,IC)-DELTPK(:,IC)
                    ICORR(IC) = ICORR(IC) + 1
		!IF(IC .eq. 2) THEN			
			!print *,"ICORR =",ICORR(IC)		
			!pause
		!ENDIF
                    GO TO 222
                ENDIF
		ICORR(IC)=0
		              

		        TPKOLD(:,IC) = STRESS2BE(:,IC)

                EROR(:,IC) =  STRESS2BE(:,IC) - STRV(:,IC)

		        DO I=1,NSLIP
				    EROR(:,IC) = EROR(:,IC)+DG2(I,IC)*C6A(:,I,IC)
		        ENDDO
			  
                TJINV(:,:,IC) = EYE
                

		        DO I = 1,NSLIP
			        IF (DG2(I,IC).NE.0.0D0) THEN
				        DGDTPKA(I) = DG2(I,IC)*EXPN/DRFCE2(I)
                        
				        DO K = 1,6
					        TJINV(:,K,IC) = TJINV(:,K,IC) +C6P6A(:,K,I,IC)* DGDTPKA(I)
                       		        ENDDO
                        
			        ENDIF
               		 ENDDO


            ENDIF

        ENDDO !IC SECTION 2
!$omp end parallel
	!IF (ISTEP .eq. 1) THEN
	!print *,"TJINV(:,1,2)", TJINV(:,1,2)
	!print *,"ERROR(:,2)", EROR(:,2)
	!print *,"Stress2be(:,2)", STRESS2BE(:,2)
	!pause
	!endif

 	!$OMP PARALLEL DEFAULT(PRIVATE) &
    !$omp shared(TJINV,EROR,DELTPK)  
    !$omp do 
        DO IC = 1,NC  !IC SECTION 3
            IF (IFLAG(IC) .eq. 1) THEN
            !ELSEIF(1.eq.2) THEN
            !    CALL CPMATINV(TJINV(:,:,IC),6,6)
            !
            !    DELTPK(:,IC)=TJINV(:,1,IC)*EROR(1,IC)+TJINV(:,2,IC)*EROR(2,IC)+ &
            !            TJINV(:,3,IC)*EROR(3,IC)+TJINV(:,4,IC)*EROR(4,IC)+ &
            !            TJINV(:,5,IC)*EROR(5,IC)+TJINV(:,6,IC)*EROR(6,IC)
            ELSE
                   
                call LUDCMP(TJINV(:,:,IC),6,INDX,ID,ICODE)
                
                IF(ICODE .eq. 1) THEN
                    print *,"Singular IC = ",IC
                    RETURN
                END IF
                
                call LUBKSB(TJINV(:,:,IC),6,INDX,EROR(:,IC))
	            DELTPK(:,IC)=EROR(:,IC)
		
             ENDIF
         ENDDO !IC SECTION 3
!$omp end parallel
	!IF (ISTEP .eq. 1) THEN
	!	print *,"EROR(:,2)", EROR(:,2)
	!pause
	!endif
        DO IC = 1, NC  !IC SECTION 4
            IF (IFLAG(IC) .eq. 1) THEN
            ELSE

			    ABSDELTPK = DABS(DELTPK(:,IC))
			    STRESS2BE(:,IC)= STRESS2BE(:,IC) - DELTPK(:,IC)
                
		        ERORDELTPK = 0.0D0
		        DO J = 1,6
			        ERORDELTPK = DMAX1(ERORDELTPK,ABSDELTPK(J))
		        ENDDO
        
		        IITR(IC) = IITR(IC) + 1
		!IF(IC .eq.1) THEN
		!print *,"IITR=",IITR(IC)
		!pause
		!ENDIF
		IF(IITR(IC).EQ.200)THEN 
                    ISOLER=1
                ENDIF
                IF(IITR(IC).EQ.201)THEN
                    RETURN
                ENDIF
                IF((ERORDELTPK.GT.SO/1.0E7) .or.(IITR(IC).LE.3) ) THEN
                ELSE
		        

			        DRFCE2 = 0.0D0
		            DRFCE2(:) = DRFCE2(:)+STRESS2BE(1,IC)*SSYSMAT(1,1,:,IC)
		            DRFCE2(:) = DRFCE2(:)+STRESS2BE(2,IC)*SSYSMAT(2,2,:,IC)
		            DRFCE2(:) = DRFCE2(:)+STRESS2BE(3,IC)*SSYSMAT(3,3,:,IC)
		            DRFCE2(:) = DRFCE2(:)+STRESS2BE(4,IC)*(SSYSMAT(1,2,:,IC)+SSYSMAT(2,1,:,IC))
		            DRFCE2(:) = DRFCE2(:)+STRESS2BE(5,IC)*(SSYSMAT(2,3,:,IC)+SSYSMAT(3,2,:,IC))
		            DRFCE2(:) = DRFCE2(:)+STRESS2BE(6,IC)*(SSYSMAT(1,3,:,IC)+SSYSMAT(3,1,:,IC))


		            DO I = 1,NSLIP         
			            ABSDG(I) = 0.0D0
			            DG2(I,IC) = 0.0D0
			            IF(DRFCE2(I).NE.0.0D0)THEN
			                ABSDG(I) = DGDO*DABS(DRFCE2(I)/SLPRESA2(I,IC))**EXPN
			                DG2(I,IC) = SIGN(ABSDG(I),DRFCE2(I))
			            ENDIF
		            ENDDO
		  

		            DO I = 1,NSLIP
			            IF (SLPRESA2(I,IC).GT.SS)THEN
				            HB(I) = 0.0D0
			            ELSE
				            HB(I) = HO*(1.0D0 - SLPRESA2(I,IC)/SS)**EXPA
			            ENDIF
			            SLPRESA2TMP(I) = SLPRA1(I,IC)
		            ENDDO

		            DO II = 1,NSLIP
			            SLPRESA2TMP(:) = SLPRESA2TMP(:)+HB(II)*ABSDG(II)
		            ENDDO

		            EORRSLPRES = 0.0D0
		            DO I = 1,NSLIP
		                EORRSLPRES = DMAX1(DABS(SLPRESA2TMP(I)-SLPRESA2(I,IC)),EORRSLPRES)
		                SLPRESA2(I,IC) = SLPRESA2TMP(I)
                    ENDDO
          
                
		            IF(EORRSLPRES.GT.SO/1.0E7)THEN 
                    ELSE
	                    !END OF SLIP RESISTANCE LOOP
	                    STRS1PK(:,IC)=STRESS2BE(:,IC)

	                    IF(KUPD.EQ.1)THEN
		                    DO I = 1,NSLIP
			                    SLPRA2(I,IC)=SLPRESA2(I,IC)
			                    SLPRA1(I,IC)=SLPRESA2(I,IC)
		                    ENDDO
                        ENDIF
                        IFLAG(IC)=1
                    ENDIF
                ENDIF
            ENDIF
        ENDDO !IC SECTION 4

    GO TO 111
    ENDIF

	   RETURN 

 
    
    END SUBROUTINE CPSOLV  
END MODULE CPSOLV
