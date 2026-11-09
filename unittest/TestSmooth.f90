!> @brief Unit tests of the well supply reduction curve in SMOOTH
!<
module TestSmooth
  use testdrive, only: check, error_type, new_unittest, unittest_type
  implicit none
  private
  public :: collect_smooth

  double precision, parameter :: tol = 1.0d-12

contains

  !> @brief Collect the SMOOTH tests
  !<
  subroutine collect_smooth(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("bounds", test_bounds), &
                new_unittest("interior", test_interior), &
                new_unittest("monotonic", test_monotonic) &
                ]
  end subroutine collect_smooth

  !> @brief Zero at and below 0, one at and above 1
  !<
  subroutine test_bounds(error)
    type(error_type), allocatable, intent(out) :: error
    double precision :: y

    call smooth(-1.0d0, y)
    call check(error, y, 0.0d0, thr=tol)
    if (allocated(error)) return
    call smooth(0.0d0, y)
    call check(error, y, 0.0d0, thr=tol)
    if (allocated(error)) return
    call smooth(1.0d0, y)
    call check(error, y, 1.0d0, thr=tol)
    if (allocated(error)) return
    call smooth(2.0d0, y)
    call check(error, y, 1.0d0, thr=tol)
  end subroutine test_bounds

  !> @brief Cubic 3x**2 - 2x**3 between 0 and 1
  !<
  subroutine test_interior(error)
    type(error_type), allocatable, intent(out) :: error
    double precision :: x, y
    integer :: i

    do i = 1, 9
      x = 0.1d0 * i
      call smooth(x, y)
      call check(error, y, 3.0d0 * x * x - 2.0d0 * x**3, thr=tol)
      if (allocated(error)) return
    end do
  end subroutine test_interior

  !> @brief Never decreases between 0 and 1
  !<
  subroutine test_monotonic(error)
    type(error_type), allocatable, intent(out) :: error
    double precision :: x, y, yold
    integer :: i

    yold = 0.0d0
    do i = 1, 100
      x = 0.01d0 * i
      call smooth(x, y)
      call check(error, y >= yold, "smooth decreased at x = 0.01*i")
      if (allocated(error)) return
      yold = y
    end do
  end subroutine test_monotonic

end module TestSmooth
