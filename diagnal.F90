module diagnal

	private
	
	public::diagnal_exact_real,diagm,diagnal_exact_complex,diagnal_tridiag_real
	
	private::rebakc,reducc,tql2,trbakc,tredc,ranff
	
	contains
!--------------------------public-------------------------------

	subroutine diagnal_exact_real(n,H,e,X)
		integer::n,i,j
		real(8),pointer::H(:,:),e(:),zi(:,:),zr(:,:)
		real(8),pointer,optional::X(:,:)
		allocate(zi(n,n))
		allocate(zr(n,n))
		call diagm(n,n,H,H,e,Zr,Zi)
		if(present(X)) then
			do i=1,n
				do j=1,n
					X(j,i)=zr(j,i)
				end do
			end do
		end if
		deallocate(zi)
		deallocate(zr)
	end subroutine
	
	
	subroutine diagnal_exact_complex(n,H,e,X)
		implicit none
		integer::n,i,j
		real(8),pointer::H(:,:),e(:),zi(:,:),zr(:,:)
		complex(8),pointer,optional::X(:,:)
		allocate(zi(n,n))
		allocate(zr(n,n))
		call diagm(n,n,H,H,e,Zr,Zi)
		if(present(X)) then
			do i=1,n
				do j=1,n
					X(j,i)=zr(j,i)+(0.0,1.0)*zi(j,i)
				end do
			end do
		end if
		deallocate(zi)
		deallocate(zr)
	end subroutine
	
	subroutine diagnal_tridiag_real(n,d,e,Z)
		implicit none
		integer,intent(in)::n
		real(8),dimension(n),intent(inout)::d,e
		real(8),dimension(n,n),intent(inout),optional::Z
		real(8),pointer,dimension(:,:)::X
		integer::ierr
		if(present(Z)) then
			call tql2(n,n,d,e,Z,ierr)
		else
			allocate(X(n,n))
			call tql2(n,n,d,e,X,ierr)
			deallocate(X)
		end if
		if(ierr/=0) stop 'tql2'
	end subroutine
	
	
!--------------------------private------------------------------
	
      SUBROUTINE DIAGM(NM,N,H,S,E,ZR,ZI)
!**********************************************************************
!
!     SOLVES THE COMPLEX HERMITIAN EIGENVALUE PROBLEM (H-E*S)Z=0
!     FOR MATRICES H AND S OF ORDER N.  ONLY THE LOWER TRIANGLE
!     OF THE HERMITIAN MATRICES NEED BE SUPPLIED.  IF A(I,J) IS A
!     REAL MATRIX AND B IS A HERMITIAN MATRIX, THEN B IS STORED
!     IN A IN THE FOLLOWING WAY (I.GE.J):
!       A(I,J) = REAL ( B(I,J) )
!       A(J,I) = IMAG ( B(I,J) )
!     SINCE THE DIAGONAL ELEMENTS OF B ARE REAL, THERE IS NO NEED
!     TO STORE THE ZERO IMAGINARY PARTS.
!
!                  M. WEINERT   JULY 1983
!
!     NM         FIRST DIMENSION OF ARRAYS
!     N          ORDER OF MATRICES
!     H          HAMILTONIAN MATRIX (OVERWRITTEN ON OUTPUT)
!     S          OVERLAP MATRIX (OVERWRITTEN ON OUTPUT)
!     E          EIGENVALUES
!     ZR,ZI      EIGENVECTORS (REAL, IMAGINARY PARTS)
!     E1,E2,TAU  WORK ARRAYS
!     TIME       TIME REQUIRED FOR DIAGONALIZATION
!
!     MODIFIED FOR USE IN BEST  FALL 1993 - SPRING 1994
!**********************************************************************
!
      IMPLICIT NONE
      INTEGER  I,IERR,J,N,NM
      REAL*8   T1,T2,TIME

!--->    HAMILTONIAN AND OVERLAP MATRICES
      REAL*8 S(NM,N),H(NM,N)
!--->   EIGENVALUES AND EIGENVECTORS
      REAL*8 E(N),ZR(NM,N),ZI(NM,N)
!--->   WORK ARRAYS
      REAL*8 E1(N),E2(N),TAU(2,N)
!
!--->    REDUCE THE GENERAL PROBLEM TO THE STANDARD PROBLEM
!      CALL REDUCC(NM,N,H,S)
!--->    REDUCE THE STANDARD PROBLEM TO REAL TRIDIAGONAL FORM
      CALL TREDC(NM,N,H,E,E1,E2,TAU)
!--->    FIND EIGENVALUES AND EIGENVECTORS OF REAL TRIADIAGONAL MATRIX
      DO I=1,N
         DO J=1,N
            ZR(J,I)=0.0D0
         ENDDO
         ZR(I,I)=1.0D0
      ENDDO
      CALL TQL2(NM,N,E,E1,ZR,IERR)
      IF(IERR.NE.0) STOP 'TQL2'
!--->    BACK-TRANSFORM THE EIGENVECTORS TO THE STANDARD PROBLEM
      CALL TRBAKC(NM,N,H,TAU,N,ZR,ZI)
!--->    BACK-TRANSFORM THE EIGENVECTORS TO THE ORIGINAL PROBLEM
!      CALL REBAKC(NM,N,N,S,ZR,ZI)
      RETURN
      END SUBROUTINE
      
      
      SUBROUTINE REBAKC(NM,N,M,B,ZR,ZI)
!**********************************************************************
!     COMPLEX VERSION OF THE ALGOL PROCEDURE REBAKA, LINEAR
!     ALGEBRA, VOL. II, 1971 BY WILKINSON AND REINSCH.
!     FORMS THE EIGENVECTORS OF THE GENERALIZED HERMITIAN
!     EIGENSYSTEM BY BACK-TRANSFORMING THOSE OF THE DERIVED
!     STANDARD MATRIX DETERMINED BY REDUCC.
!     INPUT:
!      NM     ROW DIMENSION OF THE 2-D ARRAYS
!      N      ORDER OF THE MATRIX SYSTEM
!      M      NUMBER OF EIGENVECTORS TO BACK-TRANSFORM
!      B      CONTAINS THE CHOLESKY DECOMPOSITION OBTAINED
!             IN REDUCC IN COMPACT STORAGE MODE.
!      ZR,ZI  CONTAIN THE REAL AND IMAGINARY PARTS OF THE
!             EIGENVECTORS TO BE BACK-TRANSFORMED IN THE
!             FIRST M COLUMNS.
!     OUTPUT:
!      ZR,ZI  CONTAIN THE BACK-TRANSFORMED EIGENVECTORS
!                 M. WEINERT   JULY 1983
!**********************************************************************
      IMPLICIT NONE
      INTEGER  I,J,K,M,N,NM
      REAL*8   XI,XR
      REAL*8   B(NM,N),ZR(NM,N),ZI(NM,N)

      DO J=1,M
         DO I=N,1,-1
            XR=ZR(I,J)
            XI=ZI(I,J)
            DO K=I+1,N
               XR=XR - B(K,I)*ZR(K,J) - B(I,K)*ZI(K,J)
               XI=XI - B(K,I)*ZI(K,J) + B(I,K)*ZR(K,J)
            ENDDO
            ZR(I,J)=XR/B(I,I)
            ZI(I,J)=XI/B(I,I)
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE
      
      
      SUBROUTINE REDUCC(NM,N,A,B)
!**********************************************************************
!     COMPLEX VERSION OF THE ALGOL PROCEDURE REDUC1, LINEAR
!     ALGEBRA, VOL. II, WILKINSON AND REINSCH, 1971.
!
!     REDUCTION OF THE GENERAL HERMITIAN EIGENVALUE PROBLEM
!     A*X=LAMDA*B*X TO THE EQUIVALENT PROBLEM P*Z=LAMDA*Z
!     USING CHOLESKY DECOMPOSITION.
!
!     THE PROCEDURE WILL FAIL IF B, PERHAPS DUE TO ROUNDING
!     ERRORS, IS NOT POSITIVE DEFINITE.
!
!     INPUT:
!      NM     ROW DIMENSION OF 2-D ARRAYS A AND B AS DECLARED IN
!             CALL ROUTINE
!      N      ORDER OF THE EIGENSYSTEM
!      A,B    THE LOWER TRIANGLE OF THE HERMITIAN MATRICES STORED
!             IN COMPACT MODE AS:  (I.GE.J)
!                A(I,J) = REAL ( (H(I,J) )
!                A(J,I) = IMAG ( (H(I,J) )
!
!     OUTPUT:
!      A      CONTAINS THE LOWER TRIANGLE OF THE REDUCED PROBLEM
!             STORED IN COMPACT MODE
!      B      CONTAINS THE LOWER TRIANGULAR CHOLESKY DECOMPOSITION
!             STORED IN COMPACT MODE
!
!                     M. WEINERT     JULY 1983
!
!     NOTE THAT A IS IN THE FORM REQUIRED IN TREDC AND B IS IN THE
!     FORM REQUIRED IN REBAKC.
!**********************************************************************
      IMPLICIT NONE

      INTEGER  I,J,K,N,NM
      REAL*8   A(NM,N),B(NM,N)
      REAL*8   XI,XR,Y,ZERO

      PARAMETER(ZERO=0.0D0)

!--->    FORM L IN LOWER TRIANGLE OF B
      DO J=1,N
         XR=B(J,J)
         DO K=1,J-1
            XR=XR - B(J,K)*B(J,K) - B(K,J)*B(K,J)
         ENDDO
         IF(XR.LE.ZERO) STOP 'REDUCC'
         Y=SQRT(XR)
         B(J,J)=Y
         DO I=J+1,N
            XR=B(I,J)
            XI=B(J,I)
            DO K=1,J-1
               XR=XR - B(I,K)*B(J,K) - B(K,I)*B(K,J)
               XI=XI - B(K,I)*B(J,K) + B(I,K)*B(K,J)
            ENDDO
            B(J,I)=XI/Y
            B(I,J)=XR/Y
         ENDDO
      ENDDO
!--->    FORM HERMITIAN CONJUGATE OF INV(L)*A
      DO J=1,N
         Y=B(J,J)
         XR=A(J,J)
         DO K=1,J-1
            XR=XR - B(J,K)*A(J,K) - B(K,J)*A(K,J)
         ENDDO
         A(J,J)=XR/Y
         DO I=J+1,N
            XR=A(I,J)
            XI=A(J,I)
            DO K=1,J-1
               XR=XR - B(J,K)*A(I,K) - B(K,J)*A(K,I)
               XI=XI + B(K,J)*A(I,K) - B(J,K)*A(K,I)
            ENDDO
            A(J,I)=XI/Y
            A(I,J)=XR/Y
         ENDDO
      ENDDO
!--->    PREMULTIPLY BY INV(L)
      DO I=1,N
         Y=B(I,I)
         DO J=1,I-1
            XR=A(I,J) - B(I,J)*A(J,J)
            XI=A(J,I) - B(J,I)*A(J,J)
            DO K=J+1,I-1
               XR=XR - B(I,K)*A(K,J) + B(K,I)*A(J,K)
               XI=XI - B(K,I)*A(K,J) - B(I,K)*A(J,K)
            ENDDO
            DO K=1,J-1
               XR=XR - B(I,K)*A(J,K) - B(K,I)*A(K,J)
               XI=XI - B(K,I)*A(J,K) + B(I,K)*A(K,J)
            ENDDO
            A(J,I)=XI/Y
            A(I,J)=XR/Y
         ENDDO
         XR=A(I,I)
         DO K=1,I-1
            XR=XR - B(I,K)*A(I,K) - B(K,I)*A(K,I)
         ENDDO
         A(I,I)=XR/Y
      ENDDO

      RETURN
      END SUBROUTINE
      
      
      
      

      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
!*********************************************************************
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,
!     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND
!     WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
!
!     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
!     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.
!     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
!     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS
!     FULL MATRIX TO TRIDIAGONAL FORM.
!
!     ON INPUT
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT.
!        N IS THE ORDER OF THE MATRIX.
!        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
!        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
!          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
!        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
!          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
!          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
!          THE IDENTITY MATRIX.
!
!      ON OUTPUT
!        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
!          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
!          UNORDERED FOR INDICES 1,2,...,IERR-1.
!        E HAS BEEN DESTROYED.
!        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
!          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
!          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
!          EIGENVALUES.
!        IERR IS SET TO
!          ZERO       FOR NORMAL RETURN,
!          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
!                     DETERMINED AFTER 30 ITERATIONS.
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
      REAL*8 D(N),E(N),Z(NM,N)
      REAL*8 B,C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2
!     .......... FIRST EXECUTABLE STATEMENT  TQL2 .........
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
      DO 100 I = 2, N
  100 E(I-1) = E(I)
      F = 0.0D0
      B = 0.0D0
      E(N) = 0.0D0
      DO 240 L = 1, N
         J = 0
         H = ABS(D(L)) + ABS(E(L))
         IF (B .LT. H) B = H
!     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
         DO 110 M = L, N
            IF (B + ABS(E(M)) .EQ. B) GO TO 120
!     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
!     .......... FORM SHIFT ..........
         L1 = L + 1
         L2 = L1 + 1
         G = D(L)
         P = (D(L1) - G) / (2.0D0 * E(L))
         R = SQRT( 1.0D0 + P**2 )
         D(L) = E(L) / (P + SIGN(R,P))
         D(L1) = E(L) * (P + SIGN(R,P))
         DL1 = D(L1)
         H = G - D(L)
         IF (L2 .GT. N) GO TO 145
         DO 140 I = L2, N
  140    D(I) = D(I) - H
  145    F = F + H
!     .......... QL TRANSFORMATION ..........
         P = D(M)
         C = 1.0D0
         C2 = C
         EL1 = E(L1)
         S = 0.0D0
         MML = M - L
!     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            C3 = C2
            C2 = C
            S2 = S
            I = M - II
            G = C * E(I)
            H = C * P
            IF (ABS(P) .LT. ABS(E(I))) GO TO 150
            C = E(I) / P
            R = SQRT(C*C+1.0D0)
            E(I+1) = S * P * R
            S = C / R
            C = 1.0D0 / R
            GO TO 160
  150       C = P / E(I)
            R = SQRT(C*C+1.0D0)
            E(I+1) = S * E(I) * R
            S = 1.0D0 / R
            C = C * S
  160       P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
!     .......... FORM VECTOR ..........
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
  200    CONTINUE
         P = -S * S2 * C3 * EL1 * E(L) / DL1
         E(L) = S * P
         D(L) = C * P
         IF (B + ABS(E(L)) .GT. B) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
!     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
  300 CONTINUE
      GO TO 1001
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END SUBROUTINE
      
      
      
      SUBROUTINE TRBAKC(NM,N,A,TAU,M,ZR,ZI)
!*********************************************************************
!     COMPLEX VERSION OF THE ALGOL PROCEDURE TRBAK3 IN
!     LINEAR ALGEBRA, VOL. II, 1971 BY WILKIINSON AND REINSCH.
!     FORMS THE EIGENVECTORS OF A COMPLEX HERMITIAN
!     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
!     REAL SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY TREDC.
!     INPUT:
!      NM     ROW DIMENSION OF 2-D ARRAYS
!      N      ORDER OF THE MATRIX SYSTEM
!      A      CONTAINS INFORMATION ABOUT THE UNITARY TRANSFORMATIONS
!             USED IN THE REDUCTION BY TREDC
!      TAU    CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATIONS
!      M      NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED
!      ZR     CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
!             IN ITS FIRST M COLUMNS
!     OUTPUT:
!      ZR, ZI CONTAIN THE REAL AND IMAGINARY PARTS, RESPECTIVELY,
!             OF THE TRANSFORMED EIGENVECTORS IN THE FIRST M COLUMNS
!
!     THE LAST COMPONENT OF EACH RETURNED VECTOR IS REAL AND THE
!     VECTOR EUCLIDEAN NORMS ARE PRESERVED.
!                M. WEINERT
!*********************************************************************
      IMPLICIT NONE
      INTEGER  I,J,K,M,N,NM
      REAL*8   A(NM,N),TAU(2,N),ZR(NM,M),ZI(NM,M)
      REAL*8   H,S,SI
!
      IF (M .EQ. 0) RETURN
!--->   TRANSFORM EIGENVECTORS OF THE REAL SYMMETRIC TRIDIAGONAL MATRIX
!--->   TO THOSE OF THE HERMITIAN TRIDIAGONAL MATRIX
      DO K=1,N
         DO J=1,M
            ZI(K,J)= -ZR(K,J)*TAU(2,K)
            ZR(K,J)=  ZR(K,J)*TAU(1,K)
         ENDDO
      ENDDO
!--->   APPLY THE HOUSEHOLDER TRANSFORMATIONS
      DO I=2,N
         H=A(I,I)
         IF(H.NE.0.0D0) THEN
            DO J=1,M
               S =0.0D0
               SI=0.0D0
               DO K=1,I-1
                  S =S +A(I,K)*ZR(K,J)
                  SI=SI+A(I,K)*ZI(K,J)
               ENDDO
               DO K=1,I-1
                  S =S -A(K,I)*ZI(K,J)
                  SI=SI+A(K,I)*ZR(K,J)
               ENDDO
               S =(S /H)/H
               SI=(SI/H)/H
               DO K=1,I-1
                  ZR(K,J)=ZR(K,J)-S *A(I,K)-SI*A(K,I)
                  ZI(K,J)=ZI(K,J)-SI*A(I,K)+S *A(K,I)
               ENDDO
            ENDDO
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE
      
      
      SUBROUTINE TREDC(NM,N,A,D,E,E2,TAU)
!*********************************************************************
!
!     COMPLEX VERSION OF THE ALGOL PROCEDURE TRED3, LINEAR
!     ALGEBRA, VOL. II, WILKINSON AND REINSCH, 1971.
!     REDUCES A COMPLEX HERMITIAN MATRIX, STORED AS A SINGLE
!     SQUARE ARRAY, TO A REAL SYMMETRIC TRIDIAGONAL MATRIX
!     USING UNITARY SIMILARITY TRANSFORMATIONS.
!
!     ON INPUT
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!           ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!           DIMENSION STATEMENT.
!        N  IS THE ORDER OF THE MATRIX.
!        A  CONTAINS THE LOWER TRIANGLE OF THE COMPLEX HERMITIAN INPUT
!           MATRIX.  THE REAL PARTS OF THE MATRIX ELEMENTS ARE STORED
!           IN THE FULL LOWER TRIANGLE OF A, AND THE IMAGINARY PARTS
!           ARE STORED IN THE TRANSPOSED POSITIONS OF THE STRICT UPPER
!           TRIANGLE OF A.  NO STORAGE IS REQUIRED FOR THE ZERO
!           IMAGINARY PARTS OF THE DIAGONAL ELEMENTS.
!
!     ON OUTPUT
!        A   CONTAINS INFORMATION ABOUT THE UNITARY TRANSFORMATIONS
!            USED IN THE REDUCTION.
!        D   CONTAINS THE DIAGONAL ELEMENTS OF THE THE TRIDIAGONAL
!            MATRIX.
!        E   CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
!            MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.
!        E2  CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
!            E2 AND E MUST BE SEPARATE LOCATIONS.
!        TAU CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATIONS.
!
!                    M. WEINERT  JUNE 1990
!*********************************************************************
      IMPLICIT NONE
      INTEGER  I,J,K,L,N,II,NM,JM1,JP1
      REAL*8   A(NM,N),D(N),E(N),E2(N),TAU(2,N)
      REAL*8   F,G,H,FI,GI,HH,SI,SCALE
!
      TAU(1,N) = 1.0D0
      TAU(2,N) = 0.0D0
!
      DO 300 I=N,2,-1
!--->    USE D AND E2 HAS TEMPORARY STORAGE
        DO K=1,I-1
           D(K) =A(I,K)
           E2(K)=A(K,I)
        ENDDO
!--->    SCALE ROWS
        SCALE=0.0D0
        DO K=1,I-1
           SCALE = SCALE + ABS(D(K)) + ABS(E2(K))
        ENDDO
!--->    IF SCALE IS TOO SMALL TO GUARANTEE ORTHOGONALITY,
!--->    TRANSFORMATION IS SKIPPED
        H=0.0D0
        IF (SCALE .EQ. 0.0D0) THEN
           TAU(1,I-1)=1.0D0
           TAU(2,I-1)=0.0D0
           E(I) =0.0D0
           E2(I)=0.0D0
           GO TO 200
        ENDIF
        DO K=1,I-1
           D(K) = D(K)/SCALE
           E2(K)=E2(K)/SCALE
        ENDDO
        H=0.0D0
        DO K=1,I-1
           H=H + D(K)*D(K) + E2(K)*E2(K)
        ENDDO
        E2(I)=SCALE*SCALE*H
        G=SQRT(H)
        E(I)=SCALE*G
        F=SQRT( D(I-1)**2 + E2(I-1)**2 )
!--->     FORM NEXT DIAGONAL ELEMENT
        IF (F.EQ.0.0D0) THEN
           TAU(1,I-1) = -TAU(1,I)
           SI=TAU(2,I)
           D(I-1)=G
           A(I,I-1)=SCALE*D(I-1)
        ELSE
           TAU(1,I-1)=(E2(I-1) * TAU(2,I) -  D(I-1) * TAU(1,I))/F
           SI        =( D(I-1) * TAU(2,I) + E2(I-1) * TAU(1,I))/F
           H = H + F * G
           G = 1.0D0 + G / F
           D(I-1) = G *  D(I-1)
           E2(I-1)= G * E2(I-1)
           A(I,I-1) = SCALE* D(I-1)
           A(I-1,I) = SCALE*E2(I-1)
        ENDIF
!--->     FORM ELEMENT OF A*U
        F = 0.0D0
        DO J=1,I-1
           G =0.0D0
           GI=0.0D0
           DO K=1,J-1
              G = G  + A(J,K) *  D(K)
              GI= GI - A(J,K) * E2(K)
           ENDDO
           DO K=1,J-1
              G = G  + A(K,J) * E2(K)
              GI= GI + A(K,J) *  D(K)
           ENDDO
           G = G  + A(J,J) *  D(J)
           GI= GI - A(J,J) * E2(J)
           DO K = J+1,I-1
              G = G  + A(K,J) *  D(K)
              GI= GI - A(K,J) * E2(K)
           ENDDO
           DO K = J+1,I-1
              G = G  - A(J,K) * E2(K)
              GI= GI - A(J,K) *  D(K)
           ENDDO
!--->    FORM ELEMENT OF P
           E(J) = G / H
           TAU(2,J) = GI / H
           F = F + E(J) * D(J) - TAU(2,J) * E2(J)
        ENDDO

        HH = F / (H + H)
!--->    FORM REDUCED A
        DO J=1,I-1
           F = D(J)
           G = E(J) - HH * F
           E(J) = G
           FI = -E2(J)
           GI = TAU(2,J) - HH * FI
           TAU(2,J) = -GI
           A(J,J) = A(J,J) - 2.0D0 * (F * G + FI * GI)
           DO K=1,J-1
              A(J,K) = A(J,K) - F * E(K) - G * D(K) + FI * TAU(2,K) + GI * E2(K)
           ENDDO
           DO K=1,J-1
              A(K,J) = A(K,J) - F * TAU(2,K) - G * E2(K) - FI * E(K) - GI * D(K)
           ENDDO
        ENDDO

        TAU(2,I-1) = -SI
  200   D(I) = A(I,I)
        A(I,I) = SCALE * SQRT(H)
  300 CONTINUE
      E(1)=0.0D0
      E2(1)=0.0D0
      D(1)=A(1,1)
      A(1,1)=0.0D0
      RETURN
      END SUBROUTINE
      
      
      FUNCTION RANFF(IDUM)
!
!     PARK AND MILLER 'MINIMAL' RANDOM NUMBER GENERATOR WITH 
!     BAYS-DURHAM SHUFFLE AND SAFEGUARDS. RETURNS A UNIFORM DEVIATE
!     IN (0.0,1.0). INITIALIZE WITH IDUM A NEGATIVE INTEGER, BUT
!     DO NOT CHANGE WITH CALLS. BASED ON THE ROUTINE RAN1 IN 
!     PRESS, ET AL. (2ND EDITION), P. 271.
!
      IMPLICIT NONE
      INTEGER  IA,IM,IQ,IR,NTAB,NDIV,J,K,IDUM
      REAL*8   RANFF,AM,EPS,RNMX

      PARAMETER(IA=16807,IM=2147483647,AM=1.D0/IM,IQ=127773,IR=2836)
      PARAMETER(NTAB=32,NDIV=1+(IM-1)/NTAB)
      PARAMETER(EPS=2.D0**(-46),RNMX=1.D0-EPS)
!
      INTEGER IV(NTAB),IY
      SAVE IV,IY
      DATA IY/0/
!--->    INITIALIZE AND MAKE SURE IDUM.NE.0
      IF(IDUM.LE.0.OR.IY.EQ.0) THEN
         IDUM=MAX(-IDUM,1)
         DO J=NTAB+8,1,-1
            K=IDUM/IQ
            IDUM=IA*(IDUM-K*IQ)-IR*K
            IF(IDUM.LT.0) IDUM=IDUM+IM
            IF(J.LE.NTAB) IV(J)=IDUM
         ENDDO
         IY=IV(1)
      ENDIF
!--->    GET RANDOM NUMBER AND REFILL SHUFFLE TABLE
      K=IDUM/IQ
      IDUM=IA*(IDUM-K*IQ)-IR*K
      IF(IDUM.LT.0) IDUM=IDUM+IM
      J=1+IY/NDIV
      IY=IV(J)
      IV(J)=IDUM
      RANFF=MIN(AM*IY,RNMX)
      RETURN
      END FUNCTION

end module diagnal
	
	
