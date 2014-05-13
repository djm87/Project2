 !compile with pgfortran SystemSolveTestProgram.f90
  program main

!====================================================================

! Program is test bed for system solver tests

!====================================================================

implicit none

integer, parameter :: n=6, bSz = 20000

double precision :: a(n,n),b(n),xsol(n),bhold(n),bans(n),Ahold(n,n),As(n,n),d,gflops

double precision :: timet,time, ti, tf

integer i,j,k, indx(n), code

  OPEN (UNIT = 10, FILE = 'A.txt', status='OLD')
  OPEN (UNIT = 20, FILE = 'b.txt', status='OLD')

  DO i = 1,n
	  READ(10,*) (A(i,j),j=1,n)
  ENDDO

  READ(20,*) (b(j),j=1,N)
 
  CLOSE(UNIT=10)
  CLOSE(UNIT=20)

data (bans(i), i=1,6) / -0.079450891595162, -0.022586335306939,1.347813177920734,-0.932308530681985,0.760671633368628,-0.258499798443604 /
!data (b(i), i=1,n) / 0.8143,    0.2435,    0.9293,    0.3500,    0.1966/
!data (bans(i), i=1,n) / -0.127805211200580,   0.066015737256944,   1.070218063207432,  -0.765355198529404,   0.710091617668916 /
! print a header and the original matrix
print *,"========================================================"
  write (*,200)
 Ahold = 0.D0
 bhold = 0.D0
  do i=1,n
     do j =1,n

        Ahold(i,j)=a(i,j)
        bhold(j)=b(j)

     end do
     write (*,201) (ahold(i,j),j=1,n)
  end do
  
 !write (*,204)
 !write (*,201) (bhold(j),j=1,n)

call cpu_time(ti)

DO i=1,bSz

  !call inverse(a,c,n) 

  a = Ahold !Used to refresh a and b since they are destroyed each iteration
  b= bhold

call LUDCMP(a,n,indx,d,code)                    !This is my solver 
!call ludcmp_mk(a,n,n,indx,d,code) 
 ! do j = 1,n

 !      write (*,201)  (a(j,k),k=1,n)

 ! end do
if(code.NE.0) THEN
	print *,'Singular or poorly conditioned matrix'
	return
end if
! call CPMATINV(a,n,n)
!DO j = 1,6
!	xsol(j)=a(j,1)*b(1)+a(j,2)*b(2)+ &
!		a(j,3)*b(3)+a(j,4)*b(4)+ &
!		a(j,5)*b(5)+a(j,6)*b(6)
!ENDDO
		  

call LUBKSB(a,n,indx,b)
!call lubksb_mk(a,n,n,indx,b)				     !This is my solver 

ENDDO

!II = MatMul(c,Ahold)

call cpu_time(tf)

  timet = tf-ti

  time = (tf-ti)/bSz

gflops = (0.33334 * n * n * n + n * n) * bSz/timet/1d9



! print the inverse matrix C = A^{-1} 
write (*,203)
write (*,201) (b(j),j=1,n)
write (*,205)
write (*,201) (bans(j),j=1,n)
write (*,202)


  do i = 1,n

       write (*,201)  (a(i,j),j=1,n)

  end do

print *,"Time individual    =",time

print *,"     Time total    =",timet
write (*,901) gflops
        
200 format (' Computing Inverse matrix ',/,/, &

            ' Matrix A')
201 format (6f12.10)
 901 format('    Performance     ='f8.3,' GFlops/s')
202 format (/,' LU')
203 format (/,' Solution Vector X')
204 format (/,' Vector b')
205 format (/,' Matlab Solution Vector X')
end
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
!  ***************************************************************
!  * Given an N x N matrix A, this routine replaces it by the LU *
!  * decomposition of a rowwise permutation of itself. A and N   *
!  * are input. INDX is an output vector which records the row   *
!  * permutation effected by the partial pivoting; D is output   *
!  * as -1 or 1, depending on whether the number of row inter-   *
!  * changes was even or odd, respectively. This routine is used *
!  * in combination with LUBKSB to solve linear equations or to  *
!  * invert a matrix. Return code is 1, if matrix is singular.   *
!  ***************************************************************
 Subroutine LUDCMP(A,N,INDX,D,CODE)
 PARAMETER(NMAX=100,TINY=1.5D-16)
 REAL*8  AMAX,DUM, SUM, A(N,N),VV(NMAX),D
 INTEGER CODE, INDX(N)

 D=1.D0
 CODE=0

 DO I=1,N
   AMAX=0.d0
   DO J=1,N
     IF (DABS(A(I,J)).GT.AMAX) AMAX=DABS(A(I,J))
   END DO ! j loop
   IF(AMAX.LT.TINY) THEN
     CODE = 1
     RETURN
   END IF
   VV(I) = 1.d0 / AMAX
 END DO ! i loop

 DO J=1,N
   DO I=1,J-1
     SUM = A(I,J)
     DO K=1,I-1
       SUM = SUM - A(I,K)*A(K,J) 
     END DO ! k loop
     A(I,J) = SUM
   END DO ! i loop
   AMAX = 0.d0
   DO I=J,N
     SUM = A(I,J)
     DO K=1,J-1
       SUM = SUM - A(I,K)*A(K,J) 
     END DO ! k loop
     A(I,J) = SUM
     DUM = VV(I)*DABS(SUM)
     IF(DUM.GE.AMAX) THEN
       IMAX = I
       AMAX = DUM
     END IF
   END DO ! i loop  
   
   IF(J.NE.IMAX) THEN
     DO K=1,N
       DUM = A(IMAX,K)
       A(IMAX,K) = A(J,K)
       A(J,K) = DUM
     END DO ! k loop
     D = -D
     VV(IMAX) = VV(J)
   END IF

   INDX(J) = IMAX
   IF(DABS(A(J,J)) < TINY) A(J,J) = TINY

   IF(J.NE.N) THEN
     DUM = 1.d0 / A(J,J)
     DO I=J+1,N
       A(I,J) = A(I,J)*DUM
     END DO ! i loop
   END IF 
 END DO ! j loop

		
 RETURN
 END subroutine LUDCMP


!  ******************************************************************
!  * Solves the set of N linear equations A . X = B.  Here A is     *
!  * input, not as the matrix A but rather as its LU decomposition, *
!  * determined by the routine LUDCMP. INDX is input as the permuta-*
!  * tion vector returned by LUDCMP. B is input as the right-hand   *
!  * side vector B, and returns with the solution vector X. A, N and*
!  * INDX are not modified by this routine and can be used for suc- *
!  * cessive calls with different right-hand sides. This routine is *
!  * also efficient for plain matrix inversion.                     *
!  ******************************************************************
 Subroutine LUBKSB(A,N,INDX,B)
 REAL*8  SUM, A(N,N),B(N)
 INTEGER INDX(N)

 II = 0	


 DO I=1,N
   LL = INDX(I)
   SUM = B(LL)
   B(LL) = B(I)
   IF(II.NE.0) THEN
     DO J=II,I-1
       SUM = SUM - A(I,J)*B(J)
     END DO ! j loop
   ELSE IF(SUM.NE.0.d0) THEN
     II = I
   END IF
   B(I) = SUM
 END DO ! i loop

 DO I=N,1,-1
   SUM = B(I)
!print *,"SUM, I, N ",SUM, I, N
   IF(I < N) THEN
     DO J=I+1,N
       SUM = SUM - A(I,J)*B(J)
!print *,"SUM = ",SUM
     END DO ! j loop
   END IF

 		
		 
   B(I) = SUM / A(I,I)
 END DO ! i loop

 END subroutine LUBKSB

      SUBROUTINE ludcmp_mk(a,n,np,indx,d,isingular)
      INTEGER n,np,indx(n),NMAX
      REAL d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER i,imax,j,k,isingular
      REAL aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
!
!        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
!
        if(aamax.eq.0.) then
        isingular=1
        return
        endif
!
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.

        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
!
        if(a(j,j).eq.0.) a(j,j)=TINY
!
        if(a(j,j).eq.0.) then
        isingular=1
        return
        endif
!
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
!
      isingular=0
!
      return
      END

      SUBROUTINE lubksb_mk(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
!  ***************************************************************
!  * Given an N x N matrix A, this routine replaces it by the LU *
!  * decomposition of a rowwise permutation of itself. A and N   *
!  * are input. INDX is an output vector which records the row   *
!  * permutation effected by the partial pivoting; D is output   *
!  * as -1 or 1, depending on whether the number of row inter-   *
!  * changes was even or odd, respectively. This routine is used *
!  * in combination with LUBKSB to solve linear equations or to  *
!  * invert a matrix. Return code is 1, if matrix is singular.   *
!  ***************************************************************
 Subroutine LUDCMP2(A,N,INDX,D,CODE)
 PARAMETER(NMAX=100,TINY=1.5D-16)
 REAL*8  AMAX,DUM, SUM, A(N,N),VV(NMAX),D
 INTEGER CODE, INDX(N)

 D=1.D0
 CODE=0

 DO I=1,N
   AMAX=0.d0
   DO J=1,N
     IF (DABS(A(I,J)).GT.AMAX) AMAX=DABS(A(I,J))
   END DO ! j loop
   IF(AMAX.LT.TINY) THEN
     CODE = 1
     RETURN
   END IF
   VV(I) = 1.d0 / AMAX
 END DO ! i loop

 DO J=1,N
   DO I=1,J-1
     SUM = A(I,J)
     DO K=1,I-1
       SUM = SUM - A(I,K)*A(K,J) 
     END DO ! k loop
     A(I,J) = SUM
   END DO ! i loop
   AMAX = 0.d0
   DO I=J,N
     SUM = A(I,J)
     DO K=1,J-1
       SUM = SUM - A(I,K)*A(K,J) 
     END DO ! k loop
     A(I,J) = SUM
     DUM = VV(I)*DABS(SUM)
     IF(DUM.GE.AMAX) THEN
       IMAX = I
       AMAX = DUM
     END IF
   END DO ! i loop  
   
   IF(J.NE.IMAX) THEN
     DO K=1,N
       DUM = A(IMAX,K)
       A(IMAX,K) = A(J,K)
       A(J,K) = DUM
     END DO ! k loop
     D = -D
     VV(IMAX) = VV(J)
   END IF

   INDX(J) = IMAX
   IF(DABS(A(J,J)) < TINY) A(J,J) = TINY

   IF(J.NE.N) THEN
     DUM = 1.d0 / A(J,J)
     DO I=J+1,N
       A(I,J) = A(I,J)*DUM
     END DO ! i loop
   END IF 
 END DO ! j loop

 RETURN
 END subroutine LUDCMP2


!  ******************************************************************
!  * Solves the set of N linear equations A . X = B.  Here A is     *
!  * input, not as the matrix A but rather as its LU decomposition, *
!  * determined by the routine LUDCMP. INDX is input as the permuta-*
!  * tion vector returned by LUDCMP. B is input as the right-hand   *
!  * side vector B, and returns with the solution vector X. A, N and*
!  * INDX are not modified by this routine and can be used for suc- *
!  * cessive calls with different right-hand sides. This routine is *
!  * also efficient for plain matrix inversion.                     *
!  ******************************************************************
 Subroutine LUBKSB2(A,N,INDX,B)
 REAL*8  SUM, A(N,N),B(N)
 INTEGER INDX(N)

 II = 0

 DO I=1,N
   LL = INDX(I)
   SUM = B(LL)
   B(LL) = B(I)
   IF(II.NE.0) THEN
     DO J=II,I-1
       SUM = SUM - A(I,J)*B(J)
     END DO ! j loop
   ELSE IF(SUM.NE.0.d0) THEN
     II = I
   END IF
   B(I) = SUM
 END DO ! i loop

 DO I=N,1,-1
   SUM = B(I)
   IF(I < N) THEN
     DO J=I+1,N
       SUM = SUM - A(I,J)*B(J)
     END DO ! j loop
   END IF
   B(I) = SUM / A(I,I)
 END DO ! i loop

 RETURN
 END subroutine LUBKSB2
