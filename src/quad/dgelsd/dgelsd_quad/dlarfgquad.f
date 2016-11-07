*> \brief \b dlarfgquad generates an elementary reflector (Householder matrix).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download dlarfgquad + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarfg.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarfg.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarfg.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE dlarfgquad( N, ALPHA, X, INCX, TAU )
* 
*       .. Scalar Arguments ..
*       INTEGER            INCX, N
*       real*16   ALPHA, TAU
*       ..
*       .. Array Arguments ..
*       real*16   X( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> dlarfgquad generates a real elementary reflector H of order n, such
*> that
*>
*>       H * ( alpha ) = ( beta ),   H**T * H = I.
*>           (   x   )   (   0  )
*>
*> where alpha and beta are scalars, and x is an (n-1)-element real
*> vector. H is represented in the form
*>
*>       H = I - tau * ( 1 ) * ( 1 v**T ) ,
*>                     ( v )
*>
*> where tau is a real scalar and v is a real (n-1)-element
*> vector.
*>
*> If the elements of x are all zero, then tau = 0 and H is taken to be
*> the unit matrix.
*>
*> Otherwise  1 <= tau <= 2.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the elementary reflector.
*> \endverbatim
*>
*> \param[in,out] ALPHA
*> \verbatim
*>          ALPHA is real*16
*>          On entry, the value alpha.
*>          On exit, it is overwritten with the value beta.
*> \endverbatim
*>
*> \param[in,out] X
*> \verbatim
*>          X is real*16 array, dimension
*>                         (1+(N-2)*abs(INCX))
*>          On entry, the vector x.
*>          On exit, it is overwritten with the vector v.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>          The increment between elements of X. INCX > 0.
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is real*16
*>          The value tau.
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
*> \ingroup doubleOTHERauxiliary
*
*  =====================================================================
      SUBROUTINE dlarfgquad( N, ALPHA, X, INCX, TAU )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      real*16   ALPHA, TAU
*     ..
*     .. Array Arguments ..
      real*16   X( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      real*16   ONE, ZERO
      PARAMETER          ( ONE = 1.0q0, ZERO = 0.0q0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, KNT
      real*16   BETA, RSAFMN, SAFMIN, XNORM
*     ..
*     .. External Functions ..
      real*16   dlamchquad, dlapy2quad, dnrm2quad
      EXTERNAL           dlamchquad, dlapy2quad, dnrm2quad
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN
*     ..
*     .. External Subroutines ..
      EXTERNAL           dscalquad
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.1 ) THEN
         TAU = ZERO
         RETURN
      END IF
*
      XNORM = dnrm2quad( N-1, X, INCX )
*
      IF( XNORM.EQ.ZERO ) THEN
*
*        H  =  I
*
         TAU = ZERO
      ELSE
*
*        general case
*
         BETA = -SIGN( dlapy2quad( ALPHA, XNORM ), ALPHA )
         SAFMIN = dlamchquad( 'S' ) / dlamchquad( 'E' )
         KNT = 0
         IF( ABS( BETA ).LT.SAFMIN ) THEN
*
*           XNORM, BETA may be inaccurate; scale X and recompute them
*
            RSAFMN = ONE / SAFMIN
   10       CONTINUE
            KNT = KNT + 1
            CALL dscalquad( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHA = ALPHA*RSAFMN
            IF( ABS( BETA ).LT.SAFMIN )
     $         GO TO 10
*
*           New BETA is at most 1, at least SAFMIN
*
            XNORM = dnrm2quad( N-1, X, INCX )
            BETA = -SIGN( dlapy2quad( ALPHA, XNORM ), ALPHA )
         END IF
         TAU = ( BETA-ALPHA ) / BETA
         CALL dscalquad( N-1, ONE / ( ALPHA-BETA ), X, INCX )
*
*        If ALPHA is subnormal, it may lose relative accuracy
*
         DO 20 J = 1, KNT
            BETA = BETA*SAFMIN
 20      CONTINUE
         ALPHA = BETA
      END IF
*
      RETURN
*
*     End of dlarfgquad
*
      END
