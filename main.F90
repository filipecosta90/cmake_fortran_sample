program main
  implicit none

  include "fftw3.f"

 integer ( kind = 4 ), parameter :: nx = 8
  integer ( kind = 4 ), parameter :: ny = 10

   real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) i
  complex ( kind = 8 ) in(nx,ny)
  complex ( kind = 8 ) in2(nx,ny)
  integer ( kind = 4 ) j
  complex ( kind = 8 ) out(nx,ny)
  integer ( kind = 8 ) plan_backward
  integer ( kind = 8 ) plan_forward
  integer ( kind = 4 ) seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Demonstrate FFTW3 on a 2D complex array'
  write ( *, '(a,i8)' ) '  NX = ', nx
  write ( *, '(a,i8)' ) '  NY = ', ny
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Transform data to FFT coefficients.'
  write ( *, '(a)' ) '  Backtransform FFT coefficients to recover'
  write ( *, '(a)' ) '  the data.'
  write ( *, '(a)' ) '  Compare recovered data to original data.'

!
!  Compute the data.
!
  do j = 1, ny
    do i = 1, nx
      a = r8_uniform_01 ( seed )
      b = r8_uniform_01 ( seed )
      in(i,j) = cmplx ( a, b, kind = 8 )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Input Data:'
  write ( *, '(a)' ) ' '

  do i = 1, nx
    do j = 1, ny
      write ( *, '(2x,i4,2x,i4,2x,2g14.6)' ) i, j, in(i,j)
    end do
  end do
!
!  Make a plan for the FFT, and forward transform the data.
!
  call dfftw_plan_dft_2d_ ( plan_forward, nx, ny, in, out, FFTW_FORWARD, &
    FFTW_ESTIMATE )

  call dfftw_execute_ ( plan_forward )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Output FFT Coefficients:'
  write ( *, '(a)' ) ' '

  do i = 1, nx
    do j = 1, ny
      write ( *, '(2x,i4,2x,i4,2x,2g14.6)' ) i, j, out(i,j)
    end do
  end do
!
!  Make a plan for the backward FFT, and recover the original data.
!
  call dfftw_plan_dft_2d_ ( plan_backward, nx, ny, out, in2, FFTW_BACKWARD, &
    FFTW_ESTIMATE )

  call dfftw_execute_ ( plan_backward )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Recovered input data divided by NX * NY:'
  write ( *, '(a)' ) ' '

  do i = 1, nx
    do j = 1, ny
      write ( *, '(2x,i4,2x,i4,2x,2g14.6)' ) &
        i, j, in2(i,j) / real ( nx * ny, kind = 8 )
    end do
  end do
!
!  Discard the information associated with the plans.
!
  call dfftw_destroy_plan_ ( plan_forward )
  call dfftw_destroy_plan_ ( plan_backward )

end program

function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r8_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
