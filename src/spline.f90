!-----------------------------------------------------------------------
! See LICENSE.TXT file for copying and redistribution conditions.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as
! published by the Free Software Foundation, either version 3 of the
! License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!-----------------------------------------------------------------------
!
! The following code was adapted from code the software package
! SPLINE for interpolation and approximation of data by John Burkhard,
! Florida State University, June 2008.
!
! http://people.sc.fsu.edu/~jburkardt/f_src/spline
!
! Modifications by Remko Scharroo, Altimetrics LLC, June 2012:
! - Reduced to only the subroutines least_set, least_val, spline_cubic_set,
!   spline_cubic_val, r83_np_fs, r8vec_bracket.
! - Added typesizes for portability.
! - Added intent() specifications and module wrapper for consistency checks.
! - Removed print statements.
! - Removed unnecessary checks that are implicit by or are already.
!   performed in the calling routines.
! - Combined last two arguments of r83_np_fs into one.
!-----------------------------------------------------------------------

module spline
use typesizes

contains

!*****************************************************************************80
!
!! LEAST_SET defines a least squares polynomial for given data.
!
subroutine least_set (point_num, x, f, w, nterms, b, c, d )
integer(fourbyteint), intent(in) :: point_num
real(eightbytereal), intent(in) :: x(point_num), f(point_num), w(point_num)
integer(fourbyteint), intent(in) :: nterms
real(eightbytereal), intent(out) :: b(nterms), c(nterms), d(nterms)
!
!  Discussion:
!
!    This routine is based on ORTPOL by Conte and deBoor.
!
!    The polynomial may be evaluated at any point X by calling LEAST_VAL.
!
!    Thanks to Andrew Telford for pointing out a mistake in the form of
!    the check that there are enough unique data points, 25 June 2008.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 June 2008
!
!  Author:
!
!    Original FORTRAN77 version by Samuel Conte, Carl deBoor.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Samuel Conte, Carl deBoor,
!    Elementary Numerical Analysis,
!    Second Edition,
!    McGraw Hill, 1972,
!    ISBN: 07-012446-4,
!    LC: QA297.C65.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of data values.
!
!    Input, real ( kind = 8 ) X(POINT_NUM), the abscissas of the data points.
!    At least NTERMS of the values in X must be distinct.
!
!    Input, real ( kind = 8 ) F(POINT_NUM), the data values at the points X(*).
!
!    Input, real ( kind = 8 ) W(POINT_NUM), the weights associated with
!    the data points.  Each entry of W should be positive.
!
!    Input, integer ( kind = 4 ) NTERMS, the number of terms to use in the
!    approximating polynomial.  NTERMS must be at least 1.
!    The degree of the polynomial is NTERMS-1.
!
!    Output, real ( kind = 8 ) B(NTERMS), C(NTERMS), D(NTERMS), are quantities
!    defining the least squares polynomial for the input data,
!    which will be needed to evaluate the polynomial.
!
!-----------------------------------------------------------------------
integer(fourbyteint) :: i, j
real(eightbytereal) :: p, pj(point_num), pjm1(point_num), s(nterms)
!
!  Set the values of P(-1,X) and P(0,X) at all data points.
!
pjm1 = 0d0
pj = 1d0
!
!  Now compute the value of P(J,X(I)) as
!
!    P(J,X(I)) = ( X(I) - B(J) ) * P(J-1,X(I)) - C(J) * P(J-2,X(I))
!
!  where
!
!    S(J) = < P(J,X), P(J,X) >
!    B(J) = < x*P(J,X), P(J,X) > / < P(J,X), P(J,X) >
!    C(J) = S(J) / S(J-1)
!
!  and the least squares coefficients are
!
!    D(J) = < F(X), P(J,X) > / < P(J,X), P(J,X) >
!
do j = 1, nterms

	d(j) = sum (w * f * pj)
	b(j) = sum (w * x * pj**2)
	s(j) = sum (w * pj**2)

	d(j) = d(j) / s(j)

	if (j == nterms) then
		c(j) = 0d0
		return
	endif

	b(j) = b(j) / s(j)

	if (j == 1) then
		c(j) = 0d0
	else
		c(j) = s(j) / s(j-1)
	endif

	do i = 1, point_num
		p = pj(i)
		pj(i) = (x(i) - b(j)) * pj(i) - c(j) * pjm1(i)
		pjm1(i) = p
	enddo
enddo
end subroutine least_set

!*****************************************************************************80
!
!! LEAST_VAL evaluates a least squares polynomial defined by LEAST_SET.
!
subroutine least_val (nterms, b, c, d, x, px)
integer(fourbyteint), intent(in) :: nterms
real(eightbytereal), intent(in) :: b(nterms), c(nterms), d(nterms), x
real(eightbytereal), intent(out) :: px
!
!  Discussion:
!
!    The least squares polynomial is assumed to be defined as a sum
!
!      P(X) = sum ( 1 <= I <= NTERMS ) D(I) * P(I-1,X)
!
!    where the orthogonal basis polynomials P(I,X) satisfy the following
!    three term recurrence:
!
!      P(-1,X) = 0
!      P(0,X) = 1
!      P(I,X) = ( X - B(I-1) ) * P(I-1,X) - C(I-1) * P(I-2,X)
!
!    Therefore, the least squares polynomial can be evaluated as follows:
!
!    If NTERMS is 1, then the value of P(X) is D(1) * P(0,X) = D(1).
!
!    Otherwise, P(X) is defined as the sum of NTERMS > 1 terms.  We can
!    reduce the number of terms by 1, because the polynomial P(NTERMS,X)
!    can be rewritten as a sum of polynomials;  Therefore, P(NTERMS,X)
!    can be eliminated from the sum, and its coefficient merged in with
!    those of other polynomials.  Repeat this process for P(NTERMS-1,X)
!    and so on until a single term remains.
!    P(NTERMS,X) of P(NTERMS-1,X)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Samuel Conte, Carl deBoor,
!    Elementary Numerical Analysis,
!    Second Edition,
!    McGraw Hill, 1972,
!    ISBN: 07-012446-4,
!    LC: QA297.C65.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTERMS, the number of terms in the least
!    squares polynomial.  NTERMS must be at least 1.  The input value of NTERMS
!    may be reduced from the value given to LEAST_SET.  This will
!    evaluate the least squares polynomial of the lower degree specified.
!
!    Input, real ( kind = 8 ) B(NTERMS), C(NTERMS), D(NTERMS), the information
!    computed by LEAST_SET.
!
!    Input, real ( kind = 8 ) X, the point at which the least squares polynomial
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) PX, the value of the least squares
!    polynomial at X.
!
!-----------------------------------------------------------------------
real(eightbytereal) :: prev, prev2
integer(fourbyteint) :: i

px = d(nterms)
prev = 0d0

do i = nterms-1, 1, -1
	prev2 = prev
	prev = px
	px = d(i) + (x - b(i)) * prev - c(i+1) * prev2
enddo
end subroutine least_val

!*****************************************************************************80
!
!! SPLINE_CUBIC_SET computes the second derivatives of a piecewise cubic spline.
!
subroutine spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp)
integer(fourbyteint), intent(in) :: n, ibcbeg, ibcend
real(eightbytereal), intent(in) :: t(n), y(n), ybcbeg, ybcend
real(eightbytereal), intent(out) :: ypp(n)
!
!  Discussion:
!
!    For data interpolation, the user must call SPLINE_CUBIC_SET to
!    determine the second derivative data, passing in the data to be
!    interpolated, and the desired boundary conditions.
!
!    The data to be interpolated, plus the SPLINE_CUBIC_SET output,
!    defines the spline.  The user may then call SPLINE_CUBIC_VAL to
!    evaluate the spline at any point.
!
!    The cubic spline is a piecewise cubic polynomial.  The intervals
!    are determined by the "knots" or abscissas of the data to be
!    interpolated.  The cubic spline has continous first and second
!    derivatives over the entire interval of interpolation.
!
!    For any point T in the interval T(IVAL), T(IVAL+1), the form of
!    the spline is
!
!      SPL(T) = A(IVAL)
!             + B(IVAL) * ( T - T(IVAL) )
!             + C(IVAL) * ( T - T(IVAL) )**2
!             + D(IVAL) * ( T - T(IVAL) )**3
!
!    If we assume that we know the values Y(*) and YPP(*), which represent
!    the values and second derivatives of the spline at each knot, then
!    the coefficients can be computed as:
!
!      A(IVAL) = Y(IVAL)
!      B(IVAL) = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
!        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
!      C(IVAL) = YPP(IVAL) / 2
!      D(IVAL) = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
!
!    Since the first derivative of the spline is
!
!      SPL'(T) =     B(IVAL)
!              + 2 * C(IVAL) * ( T - T(IVAL) )
!              + 3 * D(IVAL) * ( T - T(IVAL) )**2,
!
!    the requirement that the first derivative be continuous at interior
!    knot I results in a total of N-2 equations, of the form:
!
!      B(IVAL-1) + 2 C(IVAL-1) * (T(IVAL)-T(IVAL-1))
!      + 3 * D(IVAL-1) * (T(IVAL) - T(IVAL-1))**2 = B(IVAL)
!
!    or, setting H(IVAL) = T(IVAL+1) - T(IVAL)
!
!      ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
!      - ( YPP(IVAL) + 2 * YPP(IVAL-1) ) * H(IVAL-1) / 6
!      + YPP(IVAL-1) * H(IVAL-1)
!      + ( YPP(IVAL) - YPP(IVAL-1) ) * H(IVAL-1) / 2
!      =
!      ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * H(IVAL) / 6
!
!    or
!
!      YPP(IVAL-1) * H(IVAL-1) + 2 * YPP(IVAL) * ( H(IVAL-1) + H(IVAL) )
!      + YPP(IVAL) * H(IVAL)
!      =
!      6 * ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - 6 * ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
!
!    Boundary conditions must be applied at the first and last knots.
!    The resulting tridiagonal system can be solved for the YPP values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points; N must be
!    at least 2.
!
!    Input, real ( kind = 8 ) T(N), the points where data is specified.
!    The values should be distinct, and increasing.
!
!    Input, real ( kind = 8 ) Y(N), the data values to be interpolated.
!
!    Input, integer ( kind = 4 ) IBCBEG, the left boundary condition flag:
!
!      0: the spline should be a quadratic over the first interval;
!      1: the first derivative at the left endpoint should be YBCBEG;
!      2: the second derivative at the left endpoint should be YBCBEG.
!
!    Input, real ( kind = 8 ) YBCBEG, the left boundary value, if needed.
!
!    Input, integer ( kind = 4 ) IBCEND, the right boundary condition flag:
!
!      0: the spline should be a quadratic over the last interval;
!      1: the first derivative at the right endpoint should be YBCEND;
!      2: the second derivative at the right endpoint should be YBCEND.
!
!    Input, real ( kind = 8 ) YBCEND, the right boundary value, if needed.
!
!    Output, real ( kind = 8 ) YPP(N), the second derivatives of
!    the cubic spline.
!
!-----------------------------------------------------------------------
real(eightbytereal) :: a(3,n)
integer(fourbyteint) :: i

!
!  Set the first equation.
!
if (ibcbeg == 0) then
	ypp(1) = 0d0
	a(2,1) = 1d0
	a(1,2) = -1d0
else if (ibcbeg == 1) then
	ypp(1) = (y(2) - y(1)) / (t(2) - t(1)) - ybcbeg
	a(2,1) = (t(2) - t(1)) / 3d0
	a(1,2) = (t(2) - t(1)) / 6d0
else
	ypp(1) = ybcbeg
	a(2,1) = 1d0
	a(1,2) = 0d0
endif
!
!  Set the intermediate equations.
!
do i = 2, n-1
	ypp(i) = (y(i+1) - y(i)) / (t(i+1) - t(i)) - (y(i) - y(i-1)) / (t(i) - t(i-1))
	a(3,i-1) = (t(i) - t(i-1)) / 6d0
	a(2,i) = (t(i+1) - t(i-1)) / 3d0
	a(1,i+1) = (t(i+1) - t(i)) / 6d0
enddo
!
!  Set the last equation.
!
if (ibcend == 0) then
	ypp(n) = 0d0
	a(3,n-1) = -1d0
	a(2,n) = 1d0
else if (ibcend == 1) then
	ypp(n) = ybcend - (y(n) - y(n-1)) / (t(n) - t(n-1))
	a(3,n-1) = (t(n) - t(n-1)) / 6d0
	a(2,n) = (t(n) - t(n-1)) / 3d0
else
	ypp(n) = ybcend
	a(3,n-1) = 0d0
	a(2,n) = 1d0
endif
!
!  Special case:
!    N = 2, IBCBEG = IBCEND = 0.
!
if (n == 2 .and. ibcbeg == 0 .and. ibcend == 0) then
	ypp(1) = 0d0
	ypp(2) = 0d0
!
!  Solve the linear system.
!
else
	call r83_np_fs (n, a, ypp)
endif
end subroutine spline_cubic_set

!*****************************************************************************80
!
!! SPLINE_CUBIC_VAL evaluates a piecewise cubic spline at a point.
!
subroutine spline_cubic_val (n, t, y, ypp, tval, yval, ypval, yppval)
integer(fourbyteint), intent(in) :: n
real(eightbytereal), intent(in) :: t(n), y(n), ypp(n), tval
real(eightbytereal), intent(out) :: yval, ypval, yppval
!
!  Discussion:
!
!    SPLINE_CUBIC_SET must have already been called to define the
!    values of YPP.
!
!    For any point T in the interval T(IVAL), T(IVAL+1), the form of
!    the spline is
!
!      SPL(T) = A
!             + B * ( T - T(IVAL) )
!             + C * ( T - T(IVAL) )**2
!             + D * ( T - T(IVAL) )**3
!
!    Here:
!      A = Y(IVAL)
!      B = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
!        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
!      C = YPP(IVAL) / 2
!      D = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input, real ( kind = 8 ) T(N), the knot values.
!
!    Input, real ( kind = 8 ) Y(N), the data values at the knots.
!
!    Input, real ( kind = 8 ) YPP(N), the second derivatives of the
!    spline at the knots.
!
!    Input, real ( kind = 8 ) TVAL, a point, typically between T(1) and
!    T(N), at which the spline is to be evalulated.  If TVAL lies outside
!    this range, extrapolation is used.
!
!    Output, real ( kind = 8 ) YVAL, YPVAL, YPPVAL, the value of the spline, and
!    its first two derivatives at TVAL.
!
!-----------------------------------------------------------------------
real(eightbytereal) :: dt, h
integer(fourbyteint) :: left, right
!
!  Determine the interval [T(LEFT), T(RIGHT)] that contains TVAL.
!  Values below T(1) or above T(N) use extrapolation.
!
call r8vec_bracket ( n, t, tval, left, right )
!
!  Evaluate the polynomial.
!
dt = tval - t(left)
h = t(right) - t(left)

yval = y(left) &
	+ dt * ((y(right) - y(left)) / h - (ypp(right) / 6d0 + ypp(left) / 3d0) * h &
	+ dt * (0.5d0 * ypp(left) &
	+ dt * ((ypp(right) - ypp(left)) / (6d0 * h))))

ypval = (y(right) - y(left)) / h &
	- (ypp(right) / 6d0 + ypp(left) / 3d0) * h &
	+ dt * (ypp(left) &
	+ dt * (0.5d0 * (ypp(right) - ypp(left)) / h))

yppval = ypp(left) + dt * (ypp(right) - ypp(left)) / h
end subroutine spline_cubic_val

!*****************************************************************************80
!
!! R83_NP_FS factors and solves an R83 system.
!
subroutine r83_np_fs (n, a, x)
integer(fourbyteint), intent(in) :: n
real(eightbytereal), intent(inout) :: a(3,n), x(n)
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!    This algorithm requires that each diagonal entry be nonzero.
!    It does not use pivoting, and so can fail on systems that
!    are actually nonsingular.
!
!    Here is how an R83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input/output, real ( kind = 8 ) A(3,N).
!    On input, the tridiagonal matrix.
!    On output, the data in these vectors has been overwritten
!    by factorization information.
!
!    Input/output, real ( kind = 8 ) X(N).
!    On input, the right hand side of the linear system.
!    On output, the solution of the linear system.
!
!-----------------------------------------------------------------------
integer(fourbyteint) :: i
real(eightbytereal) :: xmult

do i = 2, n
	xmult = a(3,i-1) / a(2,i-1)
	a(2,i) = a(2,i) - xmult * a(1,i)
	x(i) = x(i) - xmult * x(i-1)
enddo

x(n) = x(n) / a(2,n)
do i = n-1, 1, -1
	x(i) = (x(i) - a(1,i+1) * x(i+1)) / a(2,i)
enddo
end subroutine r83_np_fs

!*****************************************************************************80
!
!! R8VEC_BRACKET searches a sorted R8VEC for successive brackets of a value.
!
subroutine r8vec_bracket (n, x, xval, left, right)
integer(fourbyteint), intent(in) :: n
real(eightbytereal), intent(in) :: x(n), xval
integer(fourbyteint), intent(out) :: left, right
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    If the values in the vector are thought of as defining intervals
!    on the real line, then this routine searches for the interval
!    nearest to or containing the given value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, length of input array.
!
!    Input, real ( kind = 8 ) X(N), an array sorted into ascending order.
!
!    Input, real ( kind = 8 ) XVAL, a value to be bracketed.
!
!    Output, integer ( kind = 4 ) LEFT, RIGHT, the results of the search.
!    Either:
!      XVAL < X(1), when LEFT = 1, RIGHT = 2;
!      X(N) < XVAL, when LEFT = N-1, RIGHT = N;
!    or
!      X(LEFT) <= XVAL <= X(RIGHT).
!
!-----------------------------------------------------------------------
integer(fourbyteint) :: i

do i = 2, n - 1
	if (xval < x(i)) then
		left = i - 1
		right = i
		return
	endif
enddo
left = n - 1
right = n
end subroutine r8vec_bracket

end module spline
