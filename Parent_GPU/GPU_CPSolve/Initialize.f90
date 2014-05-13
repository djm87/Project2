Module Initialize
CONTAINS
!Launching a three dimensional grid of 216 theads
attributes(global) SUBROUTINE Initialize(C6P6A,C6A,SSYSMAT,NC,NSLIP)

 
    IMPLICIT NONE
	INTEGER, VALUE ::	NC, NSLIP
    REAL(8), DEVICE :: C6P6A(NSLIP,6,NC),C6A(NSLIP,6,NC),SSYSMAT(NSLIP,9,NC)
    INTEGER :: tx,ty,tz,by,bz,IC
	REAL(8), SHARED :: S_SSYS(12,9),S_C6A(12,6)
	
	tx = threadidx%x !will be 1-12?
	ty = threadidx%y !will be 1-6?
	tz = threadidx%z !will be 1-4?
	by = blockdim%z  !will be 4

	IC = (blockIdx%z-1)*bz+tz
	
	!Read in SSYSMAT to shared memory 

		S_SSYS(tx,1) = SSYSMAT(tx,1,IC)
		S_SSYS(tx,7) = SSYSMAT(tx,2,IC)
		S_SSYS(tx,9) = SSYSMAT(tx,3,IC)
		S_SSYS(tx,4) = SSYSMAT(tx,1,IC)
		S_SSYS(tx,2) = SSYSMAT(tx,2,IC)	
		S_SSYS(tx,8) = SSYSMAT(tx,3,IC)
		S_SSYS(tx,6) = SSYSMAT(tx,1,IC)
		S_SSYS(tx,5) = SSYSMAT(tx,2,IC)
		S_SSYS(tx,3) = SSYSMAT(tx,3,IC)

		S_C6A(tx,ty) = C6A(tx,ty,IC)
		!S_C6A(tx,1) = C6A(tx,1,IC)
		!S_C6A(tx,2) = C6A(tx,2,IC)
		!S_C6A(tx,3) = C6A(tx,3,IC)
		!S_C6A(tx,4) = C6A(tx,4,IC)
		!S_C6A(tx,5) = C6A(tx,5,IC)
		!S_C6A(tx,6) = C6A(tx,6,IC)
		
		C6P6A(tx,1,IC) = S_C6A(tx,1)*S_SSYS(tx,1)	
		C6P6A(tx,2,IC) = S_C6A(tx,2)*S_SSYS(tx,2)	
		C6P6A(tx,3,IC) = S_C6A(tx,3)*S_SSYS(tx,3)	
		C6P6A(tx,4,IC) = S_C6A(tx,4)*S_SSYS(tx,4)+S_SSYS(tx,7) 
		C6P6A(tx,5,IC) = S_C6A(tx,5)*S_SSYS(tx,5)+S_SSYS(tx,8)
		C6P6A(tx,6,IC) = S_C6A(tx,6)*S_SSYS(tx,6)+S_SSYS(tx,9)

	RETURN

END SUBROUTINE Initialize
END Module Initialize

