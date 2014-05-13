MODULE Variables

	IMPLICIT NONE

	INTEGER, PARAMETER, PUBLIC :: RP = KIND(0.D0)

	INTEGER, PARAMETER :: NSLIP = 12, NC = 2048, NSTEP = 1000

	INTEGER, PARAMETER :: idef = 1
	REAL(RP), PARAMETER :: GDO = 0.001D0, EXPM = 0.01D0, SO = 16.0D0, HO = 180.0D0, &
				SS = 148.0D0, EXPA = 2.250D0, DT = 1.D0
	
	
	!These are public variables used to replace the common blocks
	REAL(RP), PUBLIC :: STRESS2BE(6,NC),STRV(6,NC) !CPSTRESS
	REAL(RP), PUBLIC :: DEF1(3,3),DEF2(3,3),DEFS1(3,3),DEFS2(3,3) !CPDEF
	REAL(RP), PUBLIC :: DRFCE2(12),DG2(12,NC),DVELP(3,3) !RSS2
	!REAL(RP), PUBLIC :: 
	!REAL(RP), PUBLIC :: 
	!REAL(RP), PUBLIC :: 

END MODULE Variables
