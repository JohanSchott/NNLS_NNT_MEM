*> \brief \b dlanstquad returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a real symmetric tridiagonal matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download dlanstquad + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlanst.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlanst.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlanst.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       real*16 FUNCTION dlanstquad( NORM, N, D, E )
* 
*       .. Scalar Arguments ..
*       CHARACTER          NORM
*       INTEGER            N
*       ..
*       .. Array Arguments ..
*       real*16   D( * ), E( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> dlanstquad  returns the value of the one norm,  or the Frobenius norm, or
*> the  infinity norm,  or the  element of  largest absolute value  of a
*> real symmetric tridiagonal matrix A.
*> \endverbatim
*>
*> \return dlanstquad
*> \verbatim
*>
*>    dlanstquad = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*>             (
*>             ( norm1(A),         NORM = '1', 'O' or 'o'
*>             (
*>             ( normI(A),         NORM = 'I' or 'i'
*>             (
*>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*>
*> where  norm1  denotes the  one norm of a matrix (maximum column sum),
*> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*> normF  denotes the  Frobenius norm of a matrix (square root of sum of
*> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] NORM
*> \verbatim
*>          NORM is CHARACTER*1
*>          Specifies the value to be returned in dlanstquad as described
*>          above.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.  When N = 0, dlanstquad is
*>          set to zero.
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is real*16 array, dimension (N)
*>          The diagonal elements of A.
*> \endverbatim
*>
*> \param[in] E
*> \verbatim
*>          E is real*16 array, dimension (N-1)
*>          The (n-1) sub-diagonal or super-diagonal elements of A.
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
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      real*16 FUNCTION dlanstquad( NORM, N, D, E )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            N
*     ..
*     .. Array Arguments ..
      real*16   D( * ), E( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      real*16   ONE, ZERO
      PARAMETER          ( ONE = 1.0q0, ZERO = 0.0q0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      real*16   ANORM, SCALE, SUM
*     ..
*     .. External Functions ..
      LOGICAL            lsamequad, disnanquad
      EXTERNAL           lsamequad, disnanquad
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlassqquad
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 ) THEN
         ANORM = ZERO
      ELSE IF( lsamequad( NORM, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         ANORM = ABS( D( N ) )
         DO 10 I = 1, N - 1
            SUM = ABS( D( I ) )
            IF( ANORM .LT. SUM .OR. disnanquad( SUM ) ) ANORM = SUM
            SUM = ABS( E( I ) )
            IF( ANORM .LT. SUM .OR. disnanquad( SUM ) ) ANORM = SUM
   10    CONTINUE
      ELSE IF( lsamequad( NORM, 'O' ) .OR. NORM.EQ.'1' .OR.
     $         lsamequad( NORM, 'I' ) ) THEN
*
*        Find norm1(A).
*
         IF( N.EQ.1 ) THEN
            ANORM = ABS( D( 1 ) )
         ELSE
            ANORM = ABS( D( 1 ) )+ABS( E( 1 ) )
            SUM = ABS( E( N-1 ) )+ABS( D( N ) )
            IF( ANORM .LT. SUM .OR. disnanquad( SUM ) ) ANORM = SUM
            DO 20 I = 2, N - 1
               SUM = ABS( D( I ) )+ABS( E( I ) )+ABS( E( I-1 ) )
               IF( ANORM .LT. SUM .OR. disnanquad( SUM ) ) ANORM = SUM
   20       CONTINUE
         END IF
      ELSE IF ((lsamequad(NORM,'F')) .OR. (lsamequad(NORM,'E'))) THEN
*
*        Find normF(A).
*
         SCALE = ZERO
         SUM = ONE
         IF( N.GT.1 ) THEN
            CALL dlassqquad( N-1, E, 1, SCALE, SUM )
            SUM = 2*SUM
         END IF
         CALL dlassqquad( N, D, 1, SCALE, SUM )
         ANORM = SCALE*SQRT( SUM )
      END IF
*
      dlanstquad = ANORM
      RETURN
*
*     End of dlanstquad
*
      END
