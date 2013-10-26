module underwood_mod

	integer,private,parameter::INMAX=10000

contains

! Date: Tue, 21 Jul 92 15:49:01 -0700
! From: Cheryl Carey <carey@sccm.Stanford.EDU>
!*****************************************************************
!*****************************************************************
!           AN ITERATIVE BLOCK LANCZOS METHOD FOR THE
!              SOLUTION OF LARGE SPARSE SYMMETRIC
!                       EIGENPROBLEMS
!                            by
!                   Richard Ray Underwood
!
!
!     if further details concerning the content of the code are
!     required, you are advised to refer to R. R. Underwood's
!     PhD thesis - "An Iterative Block Lanczos Method for the
!     Solution of Large Sparse Symmetric Eigenproblems", 1975.
!*****************************************************************
!*****************************************************************
!
!
!
!
	SUBROUTINE MINVAL(N,Q,PINIT,R,MMAX,EPS,OP,M,D,X,IECODE)
!=================================================================
!     this subroutine is the main subroutine implementing the
!     iterative block lanczos method for computing the
!     eigenvalues and eigenvectors of symmetric matrices.
!=================================================================
		IMPLICIT NONE
		INTEGER N , Q , PINIT , R , MMAX , M
		REAL(8) EPS
		EXTERNAL OP
		REAL(8) D(Q) , X(N,Q)
		INTEGER IECODE
		REAL(8) E(Q) , C(Q,Q)
		REAL(8) U(INMAX) , V(INMAX) , ERRC
		INTEGER P , S , PS , K , ITER , IMM , NCONV
!
!------
!     description of parameters:
!     N:     integer variable. The order of the symmetric
!            matrix A whose eigenvalues and eigenvectors are
!            being computed. The value of N should be less than
!            or equal to 1296 and greater than or equal to  2.
!
!     Q:     integer variable. The number of vectors of length
!            N contained in the array X. The value of Q should
!            be less than or equal to 25, at least one greater
!            than the value of R and less than or equal to N.
!
!     PINIT: integer variable. The initial block size to be used
!            in the block lanczos method. If PINIT is negative,
!            then -PINIT is used for the block size and columns
!            M+L,...,M+(-PINIT) of the array X are assumed to be
!            initialized to the initial matrix used to start the
!            block lanczos method. If the subroutine terminates
!            with a value of M less than R, then PINIT is assigned
!            a value -P where P is the final block size chosen.
!            In this circumstance, columns M+1,...M+P will contain
!            the most recent set of eigenvector approximations
!            which can be used to restart the subroutine if desired.
!
!     R:     integer variable. The number of eigenvalues and
!            eigenvectors being computed. That is, MINVAL attempts
!            to compute accurate approximations to the R least
!            eigenvalues and eigenvectors of the matrix A. The value
!            of R should be greater than zero and less than Q.
!
!     MMAX:  integer variable. The maximum number of matrix-vector
!            products A*X (where X is a vector) that are allowed
!            during one call of this subroutine to complete its task
!            of computing R eigenvalues and eigenvectors. Unless the
!            problem indicates otherwise, MMAX should be given a very
!            large value.
!
!     EPS:   REAL(8) variable. Initially, EPS should contain a value
!            indicating the relative precision to which MINVAL will
!            attempt to compute the eigenvalues and eigenvectors of A.
!            For eigenvalues less in modulus than 1, EPS will be an
!            absolute tolerance. Because of the way this method works,
!            it may happen that the later eigenvalues cannot be
!            computed to the same relative precision as those less in
!            value.
!
!     OP:    subroutine name. The actual argument corresponding to OP
!            should be the name of a subroutine used to define the
!            matrix A. This subroutine should have three arguments
!            N, U, and V, say, where N is an integer variable giving
!            the order of A, and U and V are two one-dimensional
!            arrays of length N. If W denotes the vector of order N
!            stored in U, then the statement
!                    CALL OP(N,U,V)
!            should result in the vector A*W being computed and stored
!            in V. The contents of U can be modified by this call.
!
!     M:     integer variable. M gives the number of eigenvalues and
!            eigenvectors already computed. Thus, initially, M should
!            be zero. If M is greater than zero, then columns one
!            through M of the array X are assumed to contain the
!            computed approximations to the M least eigenvalues and
!            eigenvectors of A. On exit, M contains a value equal to
!            the total number of eigenvalues and eigenvectors
!            computed including any already computed when MINVAL was
!            entered. Thus, on exit, the first M elements of D and the
!            first M columns of X will contain approximations to the
!            M least eigenvalues of A.
!
!     D:     REAL(8) array. D contains the computed eigenvalues. D
!            should be a one-dimensional array with at least Q
!            elements.
!
!     X:     REAL(8) array. X contains the computed eigenvectors. X
!            should be an array containing at least N*Q elements. X
!            is used not only to store the eigenvectors computed by
!            MINVAL, but also as working storage for the block lanczos
!            method. On exit, the first N*M elements of X contain the
!            eigenvector approximations - the first vector in the
!            first N elements, the second in the second N elements,
!            etc...
!
!
!     IECODE:integer variable. The value of IECODE indicates whether
!            MINVAL terminated successfully, and if not, the reason
!            why.
!
!               IECODE=0 : successful termination.
!               IECODE=1 : the value of N is less than 2.
!               IECODE=2 : the value of N exceeds MAX.
!               IECODE=3 : the value of R is less than 1.
!               IECODE=4 : the value of Q is less than or equal to R.
!               IECODE=5 : the value of Q is greater than 25.
!               IECODE=6 : the value of Q exceeds N.
!               IECODE=7 : the value of MMAX was exceeded before R
!                          eigenvalues and eigenvectors were
!                          computed.
!
!      Note that the subroutine has been designed to allow initial
!      approximations to the eigenvectors corresponding to the least
!      eigenvalues to be utilised if they are known (by storing them
!      in X and assigning PINIT minus the value of their number).
!      Furthermore, it has also been designed to allow restarting if
!      it stops with IECODE=7. Thus, the user of this program can
!      restart it after examining any partial results without loss of
!      previous work.
!------
!
!------
!     check that the initial values of the subroutine parameters
!     are in range.
!------
		IF ( N<2 ) THEN
			IECODE = 1
			RETURN
		ELSEIF ( N>INMAX ) THEN
			IECODE = 2
			RETURN
		ELSEIF ( R<1 ) THEN
			IECODE = 3
			RETURN
		ELSEIF ( Q<=R ) THEN
			IECODE = 4
			RETURN
!		ELSEIF ( Q>25 ) THEN
!			IECODE = 5
!			RETURN
		ELSEIF ( Q>N ) THEN
			IECODE = 6
			RETURN
		ELSE
!
!------
!     choose initial values for the block size P, the number of
!     steps that the block lanczos method is carried out, and
!     choose an initial N-by-P orthonormal matrix X1 used to start
!     the block lanczos method.
!------
			P = PINIT
			IF ( P<0 ) P = -P
			S = (Q-M)/P
			IF ( S<=2 ) THEN
				S = 2
				P = (Q-M)/2
			ENDIF
			IF ( PINIT>=0 ) THEN
				DO K = M+1 , M+P
					CALL RANDOM(N,Q,K,X)
				ENDDO
			ENDIF
			IF ( M<=0 ) THEN
				CALL ORTHG(N,Q,M,P,C,X)
!
!------
!     rotate the initial N-by-P matrix X1 so that
!        X1'*A*X1=diag(D1,D2,...,DP)
!     where DI is stored in D(I), I=1,...,P.
!------
				CALL SECTN(N,Q,M,P,OP,X,C,D,U,V)
				ERRC = 0.D0
			ENDIF
			ITER = 0
			IMM = 0
!
!------
!     the main body of the subroutine starts here. IMM
!     counts the number of matrix-vector products computed,
!     which is the number of times the subroutine named by
!     OP is called. ERRC measures the accumulated error in
!     the eigenvalues and eigenvectors.
!------
!
20		IF ( M>=R ) THEN
!
!------
!     this is the end of the main body of the subroutine. Now set
!     the value of IECODE and EXIT.
!------
				IECODE = 0
				RETURN
			ELSEIF ( IMM>MMAX ) THEN
				IECODE = 7
				PINIT = -P
			ELSE
				ITER = ITER + 1
				PS = P*S
!
!------
!     BKLANC carries out the block lanczos method and stores
!     the resulting block tridiagonal matrix MS in C and the
!     N-by-PS orthonormal matrix XS in X. The initial N-by-P
!     orthonormal matrix is assumed to be stored in columns
!     M+1 through M+PS of X. The residuals for these vectors
!     and the eigenvalue approximations in D are computed and
!     stored in E.
!------
				CALL BKLANC(N,Q,M,P,S,OP,D,C,X,E,U,V)
!
!------
!     EIGEN solves the eigenproblem for MS, storing the eigenvalues
!     in elements M+1 through M+PS of D and the eigenvectors in the
!     first P*S rows and columns of C (overwriting MS, possibly).
!------
				CALL EIGEN(Q,M,P,PS,C,D)
!
!------
!     CNVTST determines how many of the eigenvalues and eigenvectors
!     have converged using the error estimates stored in E. The number
!     that have converged is stored in NCONV. If NCONV=0, then none
!     have converged.
!------
				CALL CNVTST(N,Q,M,P,ERRC,EPS,D,E,NCONV)
!
!------
!     PCH chooses new values for P and S, the block size and the
!     number of steps for the block lanczos subprogram, respectively.
!------
				CALL PCH(N,Q,M,R,NCONV,P,S)
!
!------
!     ROTATE computes the eigenvectors of the restricted matrix
!     using XS stored in X and the eigenvectors of MS stored in C.
!     These vectors serve both as eigenvector approximations and
!     to form the matrix used to start the block lanczos method in
!     the next iteration.
!------
				CALL ROTATE(N,Q,M,PS,NCONV+P,C,X)
!
				M = M + NCONV
				IMM = IMM + P*S
				!WRITE (*,"(' =>ITER,IMM,P,PS =',4I5)") ITER , IMM , P , PS
				GOTO 20
			ENDIF
		ENDIF
!
	END SUBROUTINE
!
!
!
!
	SUBROUTINE BKLANC(N,Q,M,P,S,OP,D,C,X,E,U,V)
!====================================================================
!     this subroutine implements the block lanczos method
!     with reorthogonalization. BKLANC computes a block
!     tridiagonal matrix MS which it stores in rows and
!     columns M+1 through M+P*S of the array C, and an
!     orthonormal matrix XS which it stores in columns M+1
!     through M+P*S of the N-by-Q array X. MS is a symmetric
!     matrix with P-by-P symmetric matrices M(1),...,M(S) on
!     its diagonal and P-by-P upper triangular matrices
!     R(2),...,R(S) along its lower diagonal. Since MS is 
!     symmetric and banded, only its lower triangle (P+1
!     diagonals) is stored in C. XS is composed of S N-by-P
!     orthonormal matrices X(1),...,X(S) where X(1) is given
!     and should be stored in columns M+1 through M+P of X.
!     Furthermore, X(1) is assumed to satisfy X(1)*A*X(1) =
!     diag(D(M+1),D(M+2),...,D(M+P)), and if M>0, then X(1) is
!     assumed to be orthogonal to the vectors stored in columns
!     1 through M of X. OP is the name of an external subroutine
!     used to define the matrix A. During the first step, the
!     subroutine ERR is called and the quantities EJ are computed
!     where EJ=||A*X1J-D(M+J)*X1J||, X1J is the J-th column of x(1),
!     and ||*|| denotes the Euclidean norm. EJ is stored in E(M+J),
!     J=1,2,...,P. U and V are auxilliary vectors used by OP.
!====================================================================
		IMPLICIT NONE
		INTEGER N,Q,M,P,S
		EXTERNAL OP
		REAL(8) D(Q),C(Q,Q),X(:,:)
		REAL(8) E(Q),U(N),V(N),T
		INTEGER MP1,MPPS,L,LL,LU,K,I,J,IT,KP1,IL
!
		MP1 = M+1
		MPPS = M+P*S
		DO L = 1 , S
			LL = M + (L-1)*P + 1
			LU = M + L*P
			DO K = LL , LU
				DO I = 1 , N
					U(I) = X(I,K)
				ENDDO
				CALL OP(N,U,V)
				IF ( L>1 ) THEN
					DO I = K , LU
						T = 0
						DO J = 1 , N
							T = T + V(J)*X(J,I)
						ENDDO
						C(I,K) = T
					ENDDO
					IT = K - P
					DO I = 1 , N
						T = 0
						DO J = IT , K
							T = T + X(I,J)*C(K,J)
						ENDDO
						IF ( K==LU ) THEN
							V(I) = V(I) - T
						ELSE
						KP1 = K + 1
						DO J = KP1 , LU
							T = T + X(I,J)*C(J,K)
						ENDDO
						ENDIF
					ENDDO
				ELSE
					DO I = K , LU
						C(I,K) = 0.D0
					ENDDO
					C(K,K) = D(K)
					DO I = 1 , N
						V(I) = V(I) - D(K)*X(I,K)
					ENDDO
				ENDIF
				IF ( L==S ) CYCLE
				DO I = 1 , N
					 X(I,K+P) = V(I)
				ENDDO
			ENDDO
			IF ( L==1 ) CALL ERR(N,Q,M,P,X,E)
			IF ( L==S ) CYCLE
			CALL ORTHG(N,Q,LU,P,C,X)
			IL = LU + 1
			IT = LU
			DO J = 1 , P
			IT = IT + 1
			DO I = IL , IT
				C(I,IT-P) = C(I,IT)
			ENDDO
			ENDDO
		ENDDO
!
	END SUBROUTINE
!
!
!
!
	SUBROUTINE PCH(N,Q,M,R,NCONV,P,S)
!====================================================================
!     based on the values of N, Q, M, R and NCONV, PCH chooses new
!     values for P and S, the block size and number of steps for the
!     block lanczos method. The strategy used here is to choose P to
!     be the smaller of the two following values: 
!            1) the previous block size
!     and,   2) the number of values left to be computed. S is chosen
!     as large as possible subject to the constraints imposed by the
!     limits of storage. In any event, S is greater than or equal to
!     2. N is the order of the problem and Q is the number of vectors
!     available for storing eigenvectors and applying the block
!     lanczos method. M is the number of eigenvalues and eigenvectors
!     that have already been computed and R is the required number.
!     Finally, NCONV is the number of eigenvalues and eigenvectors
!     that have converged in the current iteration.
!====================================================================
		IMPLICIT NONE
		INTEGER N , Q , M , R , NCONV , P , S
		INTEGER PT , ST , MT
		!
		MT = M + NCONV
		PT = R - MT
		IF ( PT>P ) PT = P
		IF ( PT>0 ) THEN
			ST = (Q-MT)/PT
		ELSE
		  P = 0
		  RETURN
		ENDIF
		!
		IF ( ST>2 ) THEN
		  P = PT
		  S = ST
		ELSE
			ST = 2
			PT = (Q-MT)/2
			P = PT
			S = ST
		ENDIF
!
	END SUBROUTINE
!
!
!
!
	SUBROUTINE ERR(N,Q,M,P,X,E)
!================================================
!     ERR COMPUTES THE EUCLIDEAN LENGTHS OF THE
!     VECTORS STORED IN THE COLUMNS M+P+1 THROUGH
!     M+P+P OF THE N-BY-Q ARRAY X AND STORES THEM
!     IN ELEMENTS M+1 THROUGH M+P OF E.
!================================================
		IMPLICIT NONE
!
		INTEGER N , Q , M , P
		REAL(8) X(N,Q) , E(Q) , T , DSQRT
		INTEGER MP1 , MPP , K , I
!
		MP1 = M + P + 1
		MPP = M + P + P
		DO K = MP1 , MPP
			T = 0.D0
			DO I = 1 , N
				T = T + X(I,K)**2
			ENDDO
			E(K-P) = DSQRT(T)
		ENDDO
!
	END SUBROUTINE
!
!
!
!
	SUBROUTINE CNVTST(N,Q,M,P,ERRC,EPS,D,E,NCONV)
!===========================================================
!     CNVTST determines which of the P eigenvalues stored
!     in elements M+1 through M+P of D have converged. ERRC
!     is a measure of the accumulated error in the M
!     previously computed eigenvalues and eigenvectors. ERRC
!     is updated if more approximations have converged. The 
!     norms of the residual vectors are stored in elements
!     M+1 through M+P of E. EPS is the precision to which we
!     are computing the approximations. Finally, NCONV is the
!     number that have converged. If NCONV=0, then none have
!     converged.
!============================================================
		IMPLICIT NONE
		INTEGER N , Q , M , P
		REAL(8) EPS
		REAL(8) D(Q) , E(Q)
		REAL(8) ERRC , T , DSQRT
		INTEGER NCONV , K , I
		REAL(8) , PARAMETER :: CHEPS = 2.22D-16
!
		K = 0
		DO I = 1 , P
			T = DABS(D(M+I))
			IF ( T<1.D0 ) T = 1.D0
			IF ( E(M+I)>T*(EPS+10.D0*N*CHEPS)+ERRC ) EXIT
			K = I
		ENDDO
		NCONV = K
		IF ( K==0 ) RETURN
		T = 0.D0
		DO I = 1 , K
			T = T + E(M+I)**2
		ENDDO
		ERRC = DSQRT(ERRC**2+T)
!
	END SUBROUTINE
!
!
!
!
	SUBROUTINE EIGEN(Q,M,P,PS,C,D)
!=======================================================
!     EIGEN solves the eigenproblem for the symmetric
!     matrix MS of order PS stored in rows and columns
!     M+1 through M+PS of C. The eigenvalues of MS are
!     stored in elements M+1 through M+PS of D and the
!     eigenvactors are stored in rows and columns 1 
!     through PS of C possibly overwriting MS. EIGEN
!     simply re-stores MS in a manner acceptable to the
!     subroutines TRED2 and TQL2. These two routines are
!     available through eispack.
!=======================================================
		IMPLICIT NONE
		INTEGER M , P , Q , PS , PP1
		REAL(8) C(Q,Q) , D(Q)
		REAL(8) DD(Q) , V(Q)
		INTEGER I , LIM , LM1 , J , IERR
!
		DO I = 1 , PS
			LIM = I - P
			IF ( I<=P ) LIM = 1
			IF ( LIM>1 ) THEN
				LM1 = LIM - 1
				DO J = 1 , LM1
					C(I,J) = 0.D0
				ENDDO
			ENDIF
			DO J = LIM , I
				C(I,J) = C(I+M,J+M)
			ENDDO
		ENDDO
!
!------
!     CALL EISPACK ROUTINES HERE
!------
		CALL TRED2(Q,PS,C,DD,V,C)
		CALL TQL2(Q,PS,DD,V,C,IERR)
!
		!WRITE (*,"(' => ORDER =',I4,/,('    EIGENVALUES =',10D10.4))") PS , (DD(I),I=1,PS)
		!PP1 = P + 1
		!DO J = 1 , PP1
			!WRITE (*,"(' => J =',I4,/,('    EIGENVECTORS =',10D10.4))") J , (C(I,J),I=1,PS)
		!ENDDO
		DO I = 1 , PS
			D(M+I) = DD(I)
		ENDDO
!
	END SUBROUTINE
!
!
!
!
	SUBROUTINE SECTN(N,Q,M,P,OP,X,C,D,U,V)
!==============================================================
!     SECTN transforms the N-by-P orthonormal matrix X1,
!     say, stored in columns M+1 through M+P of the N-by-Q
!     array X so that X1'*A*X1 = diag(D1,D2,...,DP), where
!     ' denotes transpose and A is a symmetric matrix of 
!     order N defined by the subroutine OP. The values D1,...
!     ,DP are stored in elements M+1 through M+P of D. SECTN
!     forms the matrix X1'*A*X1 = CP, storing CP in the array
!     C. The values D1,D2,...,DP and the eigenvectors QP of CP
!     are computed by EIGEN and stored in D and C respectively.
!     ROTATE then carries out the transformation X1<=X1*QP.
!==============================================================
		IMPLICIT NONE
		INTEGER N , Q , M , P
		EXTERNAL OP
		REAL(8) X(N,Q) , C(Q,Q) , D(Q) , U(N) , V(N) , T
		INTEGER ICOL1 , ICOL2 , I , J , K
!
		ICOL1 = M
		DO J = 1 , P
			ICOL1 = ICOL1 + 1
			DO I = 1 , N
				U(I) = X(I,ICOL1)
			ENDDO
			CALL OP(N,U,V)
			ICOL2 = M
			DO I = 1 , J
				ICOL2 = ICOL2 + 1
				T = 0.D0
				DO K = 1 , N
					T = T + V(K)*X(K,ICOL2)
				ENDDO
				C(ICOL1,ICOL2) = T
			ENDDO
		ENDDO
		CALL EIGEN(Q,M,P,P,C,D)
		CALL ROTATE(N,Q,M,P,P,C,X)
!
	END SUBROUTINE
!
!
!
!
	SUBROUTINE ROTATE(N,Q,M,PS,L,C,X)
!======================================================
!     ROTATE computes the first L columns of the matrix
!     XS*QS where XS is an N-by-PS orthonormal matrix 
!     stored in columns M+1 through M+PS of the N-by-Q
!     array X and QS is a PS-by-PS orthonormal matrix
!     stored in rows and columns 1 through Ps of the
!     array C. The result is stored in columns M+1
!     through M+L of X overwriting part of XS.
!======================================================
		IMPLICIT NONE
		INTEGER N , Q , M , PS , L
		REAL(8) C(Q,Q) , X(N,Q)
		REAL(8) V(Q) , T
		INTEGER I , J , K
!
		DO I = 1 , N
			DO K = 1 , L
				T = 0.D0
				DO J = 1 , PS
					T = T + X(I,M+J)*C(J,K)
				ENDDO
				V(K) = T
			ENDDO
			DO K = 1 , L
				X(I,M+K) = V(K)
			ENDDO
		ENDDO
!
	END SUBROUTINE
!
!
!
!
	SUBROUTINE ORTHG(N,Q,F,P,B,X)
!=========================================================
!     ORTHG reorthogonalizes the N-by-P matrix Z stored
!     in columns F+1 through F+P of the N-by-Q array X
!     with respect to the vectors stored in the first F
!     columns of X and then decomposes the resulting
!     matrix into the product of an N-by-P orthonormal
!     matrix XORTH, say, and a P-by-P upper triangular
!     matrix R. XORTH is stored over Z and the upper
!     triangle of R is stored in rows and columns F+1
!     through F+P of the Q-by-Q array B. A stable variant
!     of the Gram-Schmidt orthogonalization method is
!     utilised.
!=========================================================
		IMPLICIT NONE
		INTEGER N , Q , F , P
		REAL(8) B(Q,Q) , X(N,Q) , T , S , DSQRT
		INTEGER FP1 , FPP , K , KM1 , I , J
		LOGICAL ORIG
!
		IF ( P==0 ) RETURN
		FP1 = F + 1
		FPP = F + P
		DO K = FP1 , FPP
			ORIG = .TRUE.
			KM1 = K - 1
50		T = 0.D0
			IF ( KM1>=1 ) THEN
				DO I = 1 , KM1
					S = 0.D0
					DO J = 1 , N
						S = S + X(J,I)*X(J,K)
					ENDDO
					IF ( ORIG .AND. I>F ) B(I,K) = S
					T = T + S*S
					DO J = 1 , N
						X(J,K) = X(J,K) - S*X(J,I)
					ENDDO
				ENDDO
			ENDIF
			S = 0.D0
			DO J = 1 , N
				S = S + X(J,K)*X(J,K)
			ENDDO
			T = T + S
			IF ( S>T/100.D0 ) THEN
				S = DSQRT(S)
				B(K,K) = S
				IF ( S/=0 ) S = 1.D0/S
				DO J = 1 , N
					X(J,K) = S*X(J,K)
				ENDDO
			ELSE
				ORIG = .FALSE.
				GOTO 50
			ENDIF
		ENDDO
!
!
	END SUBROUTINE
!
!
!
!
	SUBROUTINE RANDOM(N,Q,L,X)
	use rand_mod
!======================================================
!     RANDOM computes and stores a sequence of N
!     pseudo-random numbers in the L-th column of the
!     N-by-Q array X. RANDOM generates two sequences of
!     pseudo-random numbers, filling an array with one
!     sequence and using the second to access the array 
!     in a random fashion.
!======================================================
		IMPLICIT NONE
		INTEGER N , Q , L , FT , X1 , F1 , F2 , X0
		INTEGER I , K
		REAL(8) X(N,Q)
		REAL(8) T(100)
		INTEGER , PARAMETER :: A = 6821 , C = 5327
		F1 = 71416
		F2 = 27183
		X0 = 5328
!
		!DO I = 1 , 100
		!	X1 = A*X0 + C
		!	IF ( X1>=10000 ) X1 = X1 - 10000
		!	T(I) = X1/9999.D0 - 0.5D0
		!	X0 = X1
		!ENDDO
		DO I = 1 , N
		!	FT = F1 + F2
		!	IF ( FT>=1000000 ) FT = FT - 1000000
		!	F1 = F2
		!	F2 = FT
		!	K = FT/1.D6*100 + 1
		!	X(I,L) = T(K)
			X(I,L) = random_n(0.d0,1.d2)
		!	X1 = A*X0 + C
		!	IF ( X1>=10000 ) X1 = X1 - 10000
		!	T(K) = X1/9999D0 - 0.5D0
		!	X0 = X1
		ENDDO
!
!
	END SUBROUTINE
!
!
!
!
      SUBROUTINE TRED2(NM,N,A,D,E,Z)
!
      IMPLICIT NONE
      INTEGER I,J,K,L,N,II,NM,JP1
      DOUBLE PRECISION A(NM,N),D(N),E(N),Z(NM,N)
      DOUBLE PRECISION F,G,H,HH,SCALE
!     ======================================================
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED2,
!     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!
!     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX TO A
!     SYMMETRIC TRIDIAGONAL MATRIX USING AND ACCUMULATING
!     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
!
!     ON INPUT
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT.
!
!        N IS THE ORDER OF THE MATRIX.
!
!        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE
!          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
!
!     ON OUTPUT
!
!        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.
!
!        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
!          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.
!
!        Z CONTAINS THE ORTHOGONAL TRANSFORMATION MATRIX
!          PRODUCED IN THE REDUCTION.
!
!        A AND Z MAY COINCIDE.  IF DISTINCT, A IS UNALTERED.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
!     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
!
!     THIS VERSION DATED AUGUST 1983.
!
!     ======================================================
!
      DO 100 I = 1, N
!
         DO 80 J = I, N
   80    Z(J,I) = A(J,I)
!
         D(I) = A(N,I)
  100 CONTINUE
!
      IF (N == 1) GO TO 510
!     .......... FOR I=N STEP -1 UNTIL 2 DO -- ..........
      DO 300 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0.0D0
         SCALE = 0.0D0
         IF (L < 2) GO TO 130
!     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + DABS(D(K))
!
         IF (SCALE /= 0.0D0) GO TO 140
  130    E(I) = D(L)
!
         DO 135 J = 1, L
            D(J) = Z(L,J)
            Z(I,J) = 0.0D0
            Z(J,I) = 0.0D0
  135    CONTINUE
!
         GO TO 290
!
  140    DO 150 K = 1, L
            D(K) = D(K) / SCALE
            H = H + D(K) * D(K)
  150    CONTINUE
!
         F = D(L)
         G = -DSIGN(DSQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         D(L) = F - G
!     .......... FORM A*U ..........
         DO 170 J = 1, L
  170    E(J) = 0.0D0
!
         DO 240 J = 1, L
            F = D(J)
            Z(J,I) = F
            G = E(J) + Z(J,J) * F
            JP1 = J + 1
            IF (L < JP1) GO TO 220
!
            DO 200 K = JP1, L
               G = G + Z(K,J) * D(K)
               E(K) = E(K) + Z(K,J) * F
  200       CONTINUE
!
  220       E(J) = G
  240    CONTINUE
!     .......... FORM P ..........
         F = 0.0D0
!
         DO 245 J = 1, L
            E(J) = E(J) / H
            F = F + E(J) * D(J)
  245    CONTINUE
!
         HH = F / (H + H)
!     .......... FORM Q ..........
         DO 250 J = 1, L
  250    E(J) = E(J) - HH * D(J)
!     .......... FORM REDUCED A ..........
         DO 280 J = 1, L
            F = D(J)
            G = E(J)
!
            DO 260 K = J, L
  260       Z(K,J) = Z(K,J) - F * E(K) - G * D(K)
!
            D(J) = Z(L,J)
            Z(I,J) = 0.0D0
  280    CONTINUE
!
  290    D(I) = H
  300 CONTINUE
!     .......... ACCUMULATION OF TRANSFORMATION MATRICES ..........
      DO 500 I = 2, N
         L = I - 1
         Z(N,L) = Z(L,L)
         Z(L,L) = 1.0D0
         H = D(I)
         IF (H == 0.0D0) GO TO 380
!
         DO 330 K = 1, L
  330    D(K) = Z(K,I) / H
!
         DO 360 J = 1, L
            G = 0.0D0
!
            DO 340 K = 1, L
  340       G = G + Z(K,I) * Z(K,J)
!
            DO 360 K = 1, L
               Z(K,J) = Z(K,J) - G * D(K)
  360    CONTINUE
!
  380    DO 400 K = 1, L
  400    Z(K,I) = 0.0D0
!
  500 CONTINUE
!
  510 DO 520 I = 1, N
         D(I) = Z(N,I)
         Z(N,I) = 0.0D0
  520 CONTINUE
!
      Z(N,N) = 1.0D0
      E(1) = 0.0D0
      END SUBROUTINE
!
!
!
!
      DOUBLE PRECISION FUNCTION PYTHAG (A, B)
!   ======================================================
!   BEGIN PROLOGUE  PYTHAG
!   SUBSIDIARY
!   PURPOSE  Compute the complex square root of a complex number without
!            destructive overflow or underflow.
!   LIBRARY   SLATEC
!   TYPE      SINGLE PRECISION (PYTHAG-S)
!   AUTHOR  (UNKNOWN)
!   DESCRIPTION
!
!     Finds sqrt(A**2+B**2) without overflow or destructive underflow
!
!   SEE ALSO  EISDOC
!   ROUTINES CALLED  (NONE)
!   REVISION HISTORY  (YYMMDD)
!   811101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!   END PROLOGUE  PYTHAG
!     ======================================================
      IMPLICIT NONE
      DOUBLE PRECISION A,B
      DOUBLE PRECISION P,Q,R,S,T
!   FIRST EXECUTABLE STATEMENT  PYTHAG
      P = MAX(ABS(A),ABS(B))
      Q = MIN(ABS(A),ABS(B))
      IF (Q == 0.0E0) GO TO 20
   10 CONTINUE
         R = (Q/P)**2
         T = 4.0E0 + R
         IF (T == 4.0E0) GO TO 20
         S = R/T
         P = P + 2.0E0*P*S
         Q = Q*S
      GO TO 10
   20 PYTHAG = P
      END FUNCTION
!
!
!
!
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
!
      IMPLICIT NONE
      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
      DOUBLE PRECISION D(N),E(N),Z(NM,N)
      DOUBLE PRECISION C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2,TST1,TST2
!     ======================================================
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
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT.
!
!        N IS THE ORDER OF THE MATRIX.
!
!        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
!
!        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
!          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
!
!        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
!          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
!          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
!          THE IDENTITY MATRIX.
!
!      ON OUTPUT
!
!        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
!          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
!          UNORDERED FOR INDICES 1,2,...,IERR-1.
!
!        E HAS BEEN DESTROYED.
!
!        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
!          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
!          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
!          EIGENVALUES.
!
!        IERR IS SET TO
!          ZERO       FOR NORMAL RETURN,
!          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
!                     DETERMINED AFTER 30 ITERATIONS.
!
!     CALLS PYTHAG FOR  DSQRT(A*A + B*B) .
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
!     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
!
!     THIS VERSION DATED AUGUST 1983.
!
!     ======================================================
!
      IERR = 0
      IF (N == 1) RETURN
!
      DO 100 I = 2, N
  100 E(I-1) = E(I)
!
      F = 0.0D0
      TST1 = 0.0D0
      E(N) = 0.0D0
!
      DO 240 L = 1, N
         J = 0
         H = DABS(D(L)) + DABS(E(L))
         IF (TST1 < H) TST1 = H
!     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
         DO 110 M = L, N
            TST2 = TST1 + DABS(E(M))
            IF (TST2 == TST1) GO TO 120
!     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
!
  120    IF (M == L) GO TO 220
  130    IF (J == 30) GO TO 1000
         J = J + 1
!     .......... FORM SHIFT ..........
         L1 = L + 1
         L2 = L1 + 1
         G = D(L)
         P = (D(L1) - G) / (2.0D0 * E(L))
         R = PYTHAG(P,1.0D0)
         D(L) = E(L) / (P + DSIGN(R,P))
         D(L1) = E(L) * (P + DSIGN(R,P))
         DL1 = D(L1)
         H = G - D(L)
         IF (L2 > N) GO TO 145
!
         DO 140 I = L2, N
  140    D(I) = D(I) - H
!
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
            R = PYTHAG(P,E(I))
            E(I+1) = S * R
            S = E(I) / R
            C = P / R
            P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
!     .......... FORM VECTOR ..........
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
!
  200    CONTINUE
!
         P = -S * S2 * C3 * EL1 * E(L) / DL1
         E(L) = S * P
         D(L) = C * P
         TST2 = TST1 + DABS(E(L))
         IF (TST2 > TST1) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
!     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
!
         DO 260 J = II, N
            IF (D(J) >= P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
!
         IF (K == I) GO TO 300
         D(K) = D(I)
         D(I) = P
!
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
!
  300 CONTINUE
!
      RETURN
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
      END SUBROUTINE
!
!
!
!
!
!
!
!


end module underwood_mod
