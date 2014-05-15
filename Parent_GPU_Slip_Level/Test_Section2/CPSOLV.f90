MODULE CPSOLV
    !Notes: total flow cost is between 0.5 and 2 seconds depending on number of crystals
    CONTAINS
    SUBROUTINE CPSOLV(ISOLER,DTIME,C6A,SLPRA1,SLPRA2, &
                    SSYSMAT,KUPD,STRS1PK, ISTEP,cpu_t5,cpu_t6,cpu_t7)
    USE Variables
  !  USE advisor_annotate

    IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)

    DIMENSION SLPRESA2(NSLIP,NC),DGDTPKA(NSLIP), &
            EROR(6,NC),TJINV(6,6,NC),DELTPK(6,NC), &
            ABSDELTPK(6),ABSDG(NSLIP),SLIPSIGN(NSLIP), &
            SLPRESA2TMP(NSLIP),HB(NSLIP),TPKOLD(6,NC),EYE(6,6)
    DIMENSION SLPRA1(NSLIP,NC),SLPRA2(NSLIP,NC), &
            C6A(6,NSLIP,NC),SSYSMAT(3,3,NSLIP,NC),STRS1PK(6,NC),INDX(6), &
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
    
	IF(ISTEP==2) THEN
	OPEN(102, FILE='/home/administrator/Documents/DAN/itter_spectral_comp/Parent_GPU_Slip_Level/Test_Section2/C6A.txt')
    OPEN(103, FILE='/home/administrator/Documents/DAN/itter_spectral_comp/Parent_GPU_Slip_Level/Test_Section2/SSYSMAT.txt')
	OPEN(104, FILE='/home/administrator/Documents/DAN/itter_spectral_comp/Parent_GPU_Slip_Level/Test_Section2/C6P6A.txt')
	OPEN(105, FILE='/home/administrator/Documents/DAN/itter_spectral_comp/Parent_GPU_Slip_Level/Test_Section2/STRESS2BE.txt')
	OPEN(106, FILE='/home/administrator/Documents/DAN/itter_spectral_comp/Parent_GPU_Slip_Level/Test_Section2/SLPRESA2.txt')
	OPEN(107, FILE='/home/administrator/Documents/DAN/itter_spectral_comp/Parent_GPU_Slip_Level/Test_Section2/STRV.txt')
	OPEN(108, FILE='/home/administrator/Documents/DAN/itter_spectral_comp/Parent_GPU_Slip_Level/Test_Section2/TJINV.txt')
	
	ENDIF

    !$OMP PARALLEL DEFAULT(PRIVATE) &
    !$omp shared(SO,SS,EXPA,KUPD,NSLIP,NC,IFLAG,SSYSMAT,STRESS2BE,DG2,DGDO,EXPN) &   
    !$omp shared(DMAX,SLPRESA2,SLPRA1,SLPRA2,IITR,STRS1PK,DELTPK,TPKOLD,EROR,EYE,STRV) &
    !$omp shared(C6A,C6P6A,TJINV,ICORR,cpu_t5,cpu_t6,cpu_t7,IMFLAG)
	!$OMP SINGLE
	call cpu_time(ti)
	!$OMP END SINGLE  
    !$omp do 
    DO IC = 1,NC !IC SECTION 1

          !INITIALISE SLIP RESISTANCE AND GET CONSTANT C6P6A
		DO I = 1,NSLIP
			  SLPRESA2(I,IC)= SLPRA1(I,IC)
			  		IF(ISTEP==2)THEN
     					WRITE(106,'(F15.7)') SLPRESA2(I,IC)
						WRITE(103,'(F15.7)') SSYSMAT(1,1,I,IC)
						WRITE(103,'(F15.7)') SSYSMAT(2,2,I,IC)
						WRITE(103,'(F15.7)') SSYSMAT(3,3,I,IC)
						WRITE(103,'(F15.7)') SSYSMAT(1,2,I,IC)
						WRITE(103,'(F15.7)') SSYSMAT(2,3,I,IC)
						WRITE(103,'(F15.7)') SSYSMAT(1,3,I,IC)
						WRITE(103,'(F15.7)') SSYSMAT(2,1,I,IC)
						WRITE(103,'(F15.7)') SSYSMAT(3,2,I,IC)
						WRITE(103,'(F15.7)') SSYSMAT(3,1,I,IC)
					ENDIF
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
				IF(ISTEP==2)THEN
					WRITE(102,'(F15.7)') C6A(II,I,IC)

					WRITE(104,'(F15.7)') C6P6A(II,1,I,IC)
					WRITE(104,'(F15.7)') C6P6A(II,2,I,IC)
					WRITE(104,'(F15.7)') C6P6A(II,3,I,IC)
					WRITE(104,'(F15.7)') C6P6A(II,4,I,IC)
					WRITE(104,'(F15.7)') C6P6A(II,5,I,IC)
					WRITE(104,'(F15.7)') C6P6A(II,6,I,IC)
				ENDIF
			  ENDDO
		ENDDO   
		IF(ISTEP==2)THEN
			DO I = 1,6
				WRITE(105,'(F15.7)') STRESS2BE(I,IC)
				WRITE(107,'(F15.7)') STRV(I,IC)
			ENDDO
		ENDIF
    ENDDO !IC SECTION 1
	!$OMP END DO! NOWAIT
	!$OMP SINGLE
	call cpu_time(tf)
        cpu_t6=cpu_t6+tf-ti
	!$OMP END SINGLE
   !!START OF SOLUTION PROCEDURE
   IF(ISTEP==2)THEN
	CLOSE(102)
    CLOSE(103)
	CLOSE(104)
	CLOSE(105)
	CLOSE(106)
	CLOSE(107)

	ENDIF
	
    111 CONTINUE
    IF (any(IFLAG .eq. 0,1)) THEN 
	!ELSE
	!!$OMP SINGLE
	!call cpu_time(ti)
	!!$OMP END SINGLE 
    !$OMP DO 
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
                IF(ISTEP ==2)THEN
					print*, "DRFCE2", DRFCE2(1),DRFCE2(2),DRFCE2(3)
					pause
				ENDIF
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

                IF (ABSDGMAX.GT.DMAX.AND.IITR(IC).NE.0.AND.ICORR(IC).LE.30)THEN !Make full kernel case
                    DELTPK(:,IC) = DELTPK(:,IC)*0.25
                    STRESS2BE(:,IC) = TPKOLD(:,IC)-DELTPK(:,IC)
                    ICORR(IC) = ICORR(IC) + 1

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
				IF(ISTEP==2 .and. Itter==0)THEN
					DO I = 1,6
						DO II = 1,6
							WRITE(108,'(F15.7)') TJINV(I,II,IC)
						ENDDO
					ENDDO
				ENDIF
            ENDIF
        ENDDO !IC SECTION 2
		CLOSE(108)
	!$OMP END DO NOWAIT
	!!$OMP SINGLE
	!call cpu_time(tf)
        !cpu_t7=cpu_t7+tf-ti
	!!$OMP END SINGLE
	!!$OMP SINGLE
	!call cpu_time(ti)
	!!$OMP END SINGLE
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
	!$OMP END DO NOWAIT

	!!$OMP SINGLE
	!call cpu_time(tf)
       ! cpu_t5=cpu_t5+tf-ti
	!!$OMP END SINGLE

    !$OMP DO 
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

		IF(IITR(IC).EQ.200)THEN 
                    ISOLER=1
                ENDIF
                IF(IITR(IC).EQ.201)THEN
                    RETURN
                ENDIF
                IF(((ERORDELTPK.GT.SO/1.0E7) .or.(IITR(IC).LE.3))) THEN
	
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
       !$OMP END DO NOWAIT
    GO TO 111
    ELSE
    ENDIF
    !$omp end parallel
    RETURN 

 
    
    END SUBROUTINE CPSOLV  
END MODULE CPSOLV