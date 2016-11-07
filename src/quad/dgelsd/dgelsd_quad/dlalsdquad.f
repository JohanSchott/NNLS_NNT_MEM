*> \brief \b dlalsdquad uses the singular value decomposition of A to solve the least squares problem.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download dlalsdquad + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlalsd.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlalsd.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlalsd.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE dlalsdquad( UPLO, SMLSIZ, N, NRHS, D, E, B, LDB, RCOND,
*                          RANK, WORK, IWORK, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDB, N, NRHS, RANK, SMLSIZ
*       real*16   RCOND
*       ..
*       .. Array Arguments ..
*       INTEGER            IWORK( * )
*       real*16   B( LDB, * ), D( * ), E( * ), WORK( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> dlalsdquad uses the singular value decomposition of A to solve the least
*> squares problem of finding X to minimize the Euclidean norm of each
*> column of A*X-B, where A is N-by-N upper bidiagonal, and X and B
*> are N-by-NRHS. The solution X overwrites B.
*>
*> The singular values of A smaller than RCOND times the largest
*> singular value are treated as zero in solving the least squares
*> problem; in this case a minimum norm solution is returned.
*> The actual singular values are returned in D in ascending order.
*>
*> This code makes very mild assumptions about floating point
*> arithmetic. It will work on machines with a guard digit in
*> add/subtract, or on those binary machines without guard digits
*> which subtract like the Cray XMP, Cray YMP, Cray C 90, or Cray 2.
*> It could conceivably fail on hexadecimal or decimal machines
*> without guard digits, but we know of none.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>         = 'U': D and E define an upper bidiagonal matrix.
*>         = 'L': D and E define a  lower bidiagonal matrix.
*> \endverbatim
*>
*> \param[in] SMLSIZ
*> \verbatim
*>          SMLSIZ is INTEGER
*>         The maximum size of the subproblems at the bottom of the
*>         computation tree.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>         The dimension of the  bidiagonal matrix.  N >= 0.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>         The number of columns of B. NRHS must be at least 1.
*> \endverbatim
*>
*> \param[in,out] D
*> \verbatim
*>          D is real*16 array, dimension (N)
*>         On entry D contains the main diagonal of the bidiagonal
*>         matrix. On exit, if INFO = 0, D contains its singular values.
*> \endverbatim
*>
*> \param[in,out] E
*> \verbatim
*>          E is real*16 array, dimension (N-1)
*>         Contains the super-diagonal entries of the bidiagonal matrix.
*>         On exit, E has been destroyed.
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is real*16 array, dimension (LDB,NRHS)
*>         On input, B contains the right hand sides of the least
*>         squares problem. On output, B contains the solution X.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>         The leading dimension of B in the calling subprogram.
*>         LDB must be at least max(1,N).
*> \endverbatim
*>
*> \param[in] RCOND
*> \verbatim
*>          RCOND is real*16
*>         The singular values of A less than or equal to RCOND times
*>         the largest singular value are treated as zero in solving
*>         the least squares problem. If RCOND is negative,
*>         machine precision is used instead.
*>         For example, if diag(S)*X=B were the least squares problem,
*>         where diag(S) is a diagonal matrix of singular values, the
*>         solution would be X(i) = B(i) / S(i) if S(i) is greater than
*>         RCOND*max(S), and X(i) = 0 if S(i) is less than or equal to
*>         RCOND*max(S).
*> \endverbatim
*>
*> \param[out] RANK
*> \verbatim
*>          RANK is INTEGER
*>         The number of singular values of A greater than RCOND times
*>         the largest singular value.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is real*16 array, dimension at least
*>         (9*N + 2*N*SMLSIZ + 8*N*NLVL + N*NRHS + (SMLSIZ+1)**2),
*>         where NLVL = max(0, INT(log_2 (N/(SMLSIZ+1))) + 1).
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension at least
*>         (3*N*NLVL + 11*N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>         = 0:  successful exit.
*>         < 0:  if INFO = -i, the i-th argument had an illegal value.
*>         > 0:  The algorithm failed to compute a singular value while
*>               working on the submatrix lying in rows and columns
*>               INFO/(N+1) through MOD(INFO,N+1).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup doubleOTHERcomputational
*
*> \par Contributors:
*  ==================
*>
*>     Ming Gu and Ren-Cang Li, Computer Science Division, University of
*>       California at Berkeley, USA \n
*>     Osni Marques, LBNL/NERSC, USA \n
*
*  =====================================================================
      SUBROUTINE dlalsdquad( UPLO, SMLSIZ, N, NRHS, D, E, B, LDB, RCOND,
     $                   RANK, WORK, IWORK, INFO )
*
*  -- LAPACK computational routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDB, N, NRHS, RANK, SMLSIZ
      real*16   RCOND
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      real*16   B( LDB, * ), D( * ), E( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      real*16   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            BX, BXST, C, DIFL, DIFR, GIVCOL, GIVNUM,
     $                   GIVPTR, I, ICMPQ1, ICMPQ2, IWK, J, K, NLVL,
     $                   NM1, NSIZE, NSUB, NWORK, PERM, POLES, S, SIZEI,
     $                   SMLSZP, SQRE, ST, ST1, U, VT, Z
      real*16   CS, EPS, ORGNRM, R, RCND, SN, TOL
*     ..
*     .. External Functions ..
      INTEGER            idamaxquad
      real*16   dlamchquad, dlanstquad
      EXTERNAL           idamaxquad, dlamchquad, dlanstquad
*     ..
*     .. External Subroutines ..
      EXTERNAL dcopyquad, dgemmquad, dlacpyquad, dlalsaquad, dlartgquad,
     $         dlasclquad, dlasdaquad, dlasdqquad, dlasetquad,
     $         dlasrtquad, drotquad, xerblaquad
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, real, INT, LOG, SIGN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.1 ) THEN
         INFO = -4
      ELSE IF( ( LDB.LT.1 ) .OR. ( LDB.LT.N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL xerblaquad( 'dlalsdquad', -INFO )
         RETURN
      END IF
*
      EPS = dlamchquad( 'Epsilon' )
*
*     Set up the tolerance.
*
      IF( ( RCOND.LE.ZERO ) .OR. ( RCOND.GE.ONE ) ) THEN
         RCND = EPS
      ELSE
         RCND = RCOND
      END IF
*
      RANK = 0
*
*     Quick return if possible.
*
      IF( N.EQ.0 ) THEN
         RETURN
      ELSE IF( N.EQ.1 ) THEN
         IF( D( 1 ).EQ.ZERO ) THEN
            CALL dlasetquad( 'A', 1, NRHS, ZERO, ZERO, B, LDB )
         ELSE
            RANK = 1
            CALL dlasclquad('G',0,0,D(1), ONE, 1, NRHS, B, LDB, INFO )
            D( 1 ) = ABS( D( 1 ) )
         END IF
         RETURN
      END IF
*
*     Rotate the matrix if it is lower bidiagonal.
*
      IF( UPLO.EQ.'L' ) THEN
         DO 10 I = 1, N - 1
            CALL dlartgquad( D( I ), E( I ), CS, SN, R )
            D( I ) = R
            E( I ) = SN*D( I+1 )
            D( I+1 ) = CS*D( I+1 )
            IF( NRHS.EQ.1 ) THEN
               CALL drotquad( 1, B( I, 1 ), 1, B( I+1, 1 ), 1, CS, SN )
            ELSE
               WORK( I*2-1 ) = CS
               WORK( I*2 ) = SN
            END IF
   10    CONTINUE
         IF( NRHS.GT.1 ) THEN
            DO 30 I = 1, NRHS
               DO 20 J = 1, N - 1
                  CS = WORK( J*2-1 )
                  SN = WORK( J*2 )
                  CALL drotquad(1,B(J, I ), 1, B( J+1, I ), 1, CS, SN )
   20          CONTINUE
   30       CONTINUE
         END IF
      END IF
*
*     Scale.
*
      NM1 = N - 1
      ORGNRM = dlanstquad( 'M', N, D, E )
      IF( ORGNRM.EQ.ZERO ) THEN
         CALL dlasetquad( 'A', N, NRHS, ZERO, ZERO, B, LDB )
         RETURN
      END IF
*
      CALL dlasclquad( 'G', 0, 0, ORGNRM, ONE, N, 1, D, N, INFO )
      CALL dlasclquad( 'G', 0, 0, ORGNRM, ONE, NM1, 1, E, NM1, INFO )
*
*     If N is smaller than the minimum divide size SMLSIZ, then solve
*     the problem with another solver.
*
      IF( N.LE.SMLSIZ ) THEN
         NWORK = 1 + N*N
         CALL dlasetquad( 'A', N, N, ZERO, ONE, WORK, N )
         CALL dlasdqquad('U',0,N,N, 0, NRHS, D, E, WORK, N, WORK, N, B,
     $                LDB, WORK( NWORK ), INFO )
         IF( INFO.NE.0 ) THEN
            RETURN
         END IF
         TOL = RCND*ABS( D( idamaxquad( N, D, 1 ) ) )
         DO 40 I = 1, N
            IF( D( I ).LE.TOL ) THEN
               CALL dlasetquad('A',1,NRHS, ZERO, ZERO, B( I, 1 ), LDB )
            ELSE
               CALL dlasclquad('G',0,0,D( I ), ONE, 1, NRHS, B( I, 1 ),
     $                      LDB, INFO )
               RANK = RANK + 1
            END IF
   40    CONTINUE
         CALL dgemmquad('T','N', N, NRHS, N, ONE, WORK, N, B, LDB, ZERO,
     $               WORK( NWORK ), N )
         CALL dlacpyquad( 'A', N, NRHS, WORK( NWORK ), N, B, LDB )
*
*        Unscale.
*
         CALL dlasclquad( 'G', 0, 0, ONE, ORGNRM, N, 1, D, N, INFO )
         CALL dlasrtquad( 'D', N, D, INFO )
         CALL dlasclquad('G', 0, 0, ORGNRM, ONE, N, NRHS, B, LDB, INFO )
*
         RETURN
      END IF
*
*     Book-keeping and setting up some constants.
*
      NLVL = INT( LOG( real( N ) / real( SMLSIZ+1 ) ) / LOG( TWO ) ) + 1
*
      SMLSZP = SMLSIZ + 1
*
      U = 1
      VT = 1 + SMLSIZ*N
      DIFL = VT + SMLSZP*N
      DIFR = DIFL + NLVL*N
      Z = DIFR + NLVL*N*2
      C = Z + NLVL*N
      S = C + N
      POLES = S + N
      GIVNUM = POLES + 2*NLVL*N
      BX = GIVNUM + 2*NLVL*N
      NWORK = BX + N*NRHS
*
      SIZEI = 1 + N
      K = SIZEI + N
      GIVPTR = K + N
      PERM = GIVPTR + N
      GIVCOL = PERM + NLVL*N
      IWK = GIVCOL + NLVL*N*2
*
      ST = 1
      SQRE = 0
      ICMPQ1 = 1
      ICMPQ2 = 0
      NSUB = 0
*
      DO 50 I = 1, N
         IF( ABS( D( I ) ).LT.EPS ) THEN
            D( I ) = SIGN( EPS, D( I ) )
         END IF
   50 CONTINUE
*
      DO 60 I = 1, NM1
         IF( ( ABS( E( I ) ).LT.EPS ) .OR. ( I.EQ.NM1 ) ) THEN
            NSUB = NSUB + 1
            IWORK( NSUB ) = ST
*
*           Subproblem found. First determine its size and then
*           apply divide and conquer on it.
*
            IF( I.LT.NM1 ) THEN
*
*              A subproblem with E(I) small for I < NM1.
*
               NSIZE = I - ST + 1
               IWORK( SIZEI+NSUB-1 ) = NSIZE
            ELSE IF( ABS( E( I ) ).GE.EPS ) THEN
*
*              A subproblem with E(NM1) not too small but I = NM1.
*
               NSIZE = N - ST + 1
               IWORK( SIZEI+NSUB-1 ) = NSIZE
            ELSE
*
*              A subproblem with E(NM1) small. This implies an
*              1-by-1 subproblem at D(N), which is not solved
*              explicitly.
*
               NSIZE = I - ST + 1
               IWORK( SIZEI+NSUB-1 ) = NSIZE
               NSUB = NSUB + 1
               IWORK( NSUB ) = N
               IWORK( SIZEI+NSUB-1 ) = 1
               CALL dcopyquad( NRHS, B( N, 1 ), LDB, WORK( BX+NM1 ), N )
            END IF
            ST1 = ST - 1
            IF( NSIZE.EQ.1 ) THEN
*
*              This is a 1-by-1 subproblem and is not solved
*              explicitly.
*
               CALL dcopyquad(NRHS, B( ST, 1 ), LDB, WORK( BX+ST1 ), N )
            ELSE IF( NSIZE.LE.SMLSIZ ) THEN
*
*              This is a small subproblem and is solved by dlasdqquad.
*
               CALL dlasetquad( 'A', NSIZE, NSIZE, ZERO, ONE,
     $                      WORK( VT+ST1 ), N )
               CALL dlasdqquad( 'U', 0, NSIZE, NSIZE, 0, NRHS, D( ST ),
     $                      E( ST ), WORK( VT+ST1 ), N, WORK( NWORK ),
     $                      N, B( ST, 1 ), LDB, WORK( NWORK ), INFO )
               IF( INFO.NE.0 ) THEN
                  RETURN
               END IF
               CALL dlacpyquad( 'A', NSIZE, NRHS, B( ST, 1 ), LDB,
     $                      WORK( BX+ST1 ), N )
            ELSE
*
*              A large problem. Solve it using divide and conquer.
*
               CALL dlasdaquad( ICMPQ1, SMLSIZ, NSIZE, SQRE, D( ST ),
     $                      E( ST ), WORK( U+ST1 ), N, WORK( VT+ST1 ),
     $                      IWORK( K+ST1 ), WORK( DIFL+ST1 ),
     $                      WORK( DIFR+ST1 ), WORK( Z+ST1 ),
     $                      WORK( POLES+ST1 ), IWORK( GIVPTR+ST1 ),
     $                      IWORK( GIVCOL+ST1 ), N, IWORK( PERM+ST1 ),
     $                      WORK( GIVNUM+ST1 ), WORK( C+ST1 ),
     $                      WORK( S+ST1 ), WORK( NWORK ), IWORK( IWK ),
     $                      INFO )
               IF( INFO.NE.0 ) THEN
                  RETURN
               END IF
               BXST = BX + ST1
               CALL dlalsaquad( ICMPQ2, SMLSIZ, NSIZE, NRHS, B( ST, 1 ),
     $                      LDB, WORK( BXST ), N, WORK( U+ST1 ), N,
     $                      WORK( VT+ST1 ), IWORK( K+ST1 ),
     $                      WORK( DIFL+ST1 ), WORK( DIFR+ST1 ),
     $                      WORK( Z+ST1 ), WORK( POLES+ST1 ),
     $                      IWORK( GIVPTR+ST1 ), IWORK( GIVCOL+ST1 ), N,
     $                      IWORK( PERM+ST1 ), WORK( GIVNUM+ST1 ),
     $                      WORK( C+ST1 ), WORK( S+ST1 ), WORK( NWORK ),
     $                      IWORK( IWK ), INFO )
               IF( INFO.NE.0 ) THEN
                  RETURN
               END IF
            END IF
            ST = I + 1
         END IF
   60 CONTINUE
*
*     Apply the singular values and treat the tiny ones as zero.
*
      TOL = RCND*ABS( D( idamaxquad( N, D, 1 ) ) )
*
      DO 70 I = 1, N
*
*        Some of the elements in D can be negative because 1-by-1
*        subproblems were not solved explicitly.
*
         IF( ABS( D( I ) ).LE.TOL ) THEN
            CALL dlasetquad('A',1, NRHS, ZERO, ZERO, WORK( BX+I-1 ), N )
         ELSE
            RANK = RANK + 1
            CALL dlasclquad( 'G', 0, 0, D( I ), ONE, 1, NRHS,
     $                   WORK( BX+I-1 ), N, INFO )
         END IF
         D( I ) = ABS( D( I ) )
   70 CONTINUE
*
*     Now apply back the right singular vectors.
*
      ICMPQ2 = 1
      DO 80 I = 1, NSUB
         ST = IWORK( I )
         ST1 = ST - 1
         NSIZE = IWORK( SIZEI+I-1 )
         BXST = BX + ST1
         IF( NSIZE.EQ.1 ) THEN
            CALL dcopyquad( NRHS, WORK( BXST ), N, B( ST, 1 ), LDB )
         ELSE IF( NSIZE.LE.SMLSIZ ) THEN
            CALL dgemmquad( 'T', 'N', NSIZE, NRHS, NSIZE, ONE,
     $                  WORK( VT+ST1 ), N, WORK( BXST ), N, ZERO,
     $                  B( ST, 1 ), LDB )
         ELSE
            CALL dlalsaquad(ICMPQ2,SMLSIZ, NSIZE, NRHS, WORK( BXST ), N,
     $                   B( ST, 1 ), LDB, WORK( U+ST1 ), N,
     $                   WORK( VT+ST1 ), IWORK( K+ST1 ),
     $                   WORK( DIFL+ST1 ), WORK( DIFR+ST1 ),
     $                   WORK( Z+ST1 ), WORK( POLES+ST1 ),
     $                   IWORK( GIVPTR+ST1 ), IWORK( GIVCOL+ST1 ), N,
     $                   IWORK( PERM+ST1 ), WORK( GIVNUM+ST1 ),
     $                   WORK( C+ST1 ), WORK( S+ST1 ), WORK( NWORK ),
     $                   IWORK( IWK ), INFO )
            IF( INFO.NE.0 ) THEN
               RETURN
            END IF
         END IF
   80 CONTINUE
*
*     Unscale and sort the singular values.
*
      CALL dlasclquad( 'G', 0, 0, ONE, ORGNRM, N, 1, D, N, INFO )
      CALL dlasrtquad( 'D', N, D, INFO )
      CALL dlasclquad( 'G', 0, 0, ORGNRM, ONE, N, NRHS, B, LDB, INFO )
*
      RETURN
*
*     End of dlalsdquad
*
      END
