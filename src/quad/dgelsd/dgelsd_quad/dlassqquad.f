*> \brief \b dlassqquad updates a sum of squares represented in scaled form.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download dlassqquad + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlassq.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlassq.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlassq.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE dlassqquad( N, X, INCX, SCALE, SUMSQ )
* 
*       .. Scalar Arguments ..
*       INTEGER            INCX, N
*       real*16   SCALE, SUMSQ
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
*> dlassqquad  returns the values  scl  and  smsq  such that
*>
*>    ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
*>
*> where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
*> assumed to be non-negative and  scl  returns the value
*>
*>    scl = max( scale, abs( x( i ) ) ).
*>
*> scale and sumsq must be supplied in SCALE and SUMSQ and
*> scl and smsq are overwritten on SCALE and SUMSQ respectively.
*>
*> The routine makes only one pass through the vector x.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of elements to be used from the vector X.
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is real*16 array, dimension (N)
*>          The vector for which a scaled sum of squares is computed.
*>             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>          The increment between successive values of the vector X.
*>          INCX > 0.
*> \endverbatim
*>
*> \param[in,out] SCALE
*> \verbatim
*>          SCALE is real*16
*>          On entry, the value  scale  in the equation above.
*>          On exit, SCALE is overwritten with  scl , the scaling factor
*>          for the sum of squares.
*> \endverbatim
*>
*> \param[in,out] SUMSQ
*> \verbatim
*>          SUMSQ is real*16
*>          On entry, the value  sumsq  in the equation above.
*>          On exit, SUMSQ is overwritten with  smsq , the basic sum of
*>          squares from which  scl  has been factored out.
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
      SUBROUTINE dlassqquad( N, X, INCX, SCALE, SUMSQ )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      real*16   SCALE, SUMSQ
*     ..
*     .. Array Arguments ..
      real*16   X( * )
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      real*16   ZERO
      PARAMETER          ( ZERO = 0.0q0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IX
      real*16   ABSXI
*     ..
*     .. External Functions ..
      LOGICAL            disnanquad
      EXTERNAL           disnanquad
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            ABSXI = ABS( X( IX ) )
            IF( ABSXI.GT.ZERO.OR.disnanquad( ABSXI ) ) THEN
               IF( SCALE.LT.ABSXI ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
*
*     End of dlassqquad
*
      END