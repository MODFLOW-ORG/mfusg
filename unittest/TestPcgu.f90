!> @brief Unit tests of the PCGU incomplete LU preconditioner
!!
!! The matrix has a 1-2 and 1-3 connection but no 2-3 connection, so MILU
!! drops fill that changes the sign of the 2 and 3 pivots for any shift.
!<
module TestPcgu
  use testdrive, only: check, error_type, new_unittest, unittest_type
  use PCGUMODULE, only: SPCGU_PCCRS, SPCGU_PCU, SPCGU_PCILU0
  implicit none
  private
  public :: collect_pcgu

  integer, parameter :: neq = 3
  integer, parameter :: nnz = 7
  integer, parameter :: iout = 6
  double precision, parameter :: a = 0.95d0
  double precision, parameter :: dtiny = 1.0d-6
  double precision, parameter :: tol = 1.0d-12

contains

  !> @brief Collect the PCGU tests
  !<
  subroutine collect_pcgu(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("ilu0_diagonal", test_ilu0_diagonal), &
                new_unittest("milu0_return", test_milu0_return), &
                new_unittest("milu0_small_pivot", test_milu0_small_pivot), &
                new_unittest("milu0_shift_limit", test_milu0_shift_limit) &
                ]
  end subroutine collect_pcgu

  !> @brief Fill the test matrix and the preconditioner structure
  !<
  subroutine setup(ia, ja, ac, iapc, japc, iw)
    integer, intent(out) :: ia(neq + 1), ja(nnz)
    double precision, intent(out) :: ac(nnz)
    integer, intent(out) :: iapc(neq + 1), japc(nnz), iw(neq + 1)

    ia = [1, 4, 6, 8]
    ja = [1, 2, 3, 2, 1, 3, 1]
    ac = [1.0d0, a, a, 1.0d0, a, 1.0d0, a]
    call SPCGU_PCCRS(neq, nnz, ia, ja, iapc, japc, iw)
  end subroutine setup

  !> @brief ILU0 of a well-posed matrix has the exact first pivot
  !<
  subroutine test_ilu0_diagonal(error)
    type(error_type), allocatable, intent(out) :: error
    integer :: ia(neq + 1), ja(nnz), iapc(neq + 1), japc(nnz), iw(neq + 1)
    double precision :: ac(nnz), apc(nnz), w(neq)

    call setup(ia, ja, ac, iapc, japc, iw)
    call SPCGU_PCU(iout, nnz, neq, neq, nnz, nnz, 2, 0.0, &
                   ac, ia, ja, apc, iapc, japc, iw, w)
    call check(error, apc(1), 1.0d0, thr=tol)
    if (allocated(error)) return
    call check(error, apc(2), 1.0d0 / (1.0d0 - a * a), thr=tol)
  end subroutine test_ilu0_diagonal

  !> @brief MILU0 returns IZERO of 1 for a sign-changing pivot
  !<
  subroutine test_milu0_return(error)
    type(error_type), allocatable, intent(out) :: error
    integer :: ia(neq + 1), ja(nnz), iapc(neq + 1), japc(nnz), iw(neq + 1)
    double precision :: ac(nnz), apc(nnz), w(neq)
    integer :: izero

    call setup(ia, ja, ac, iapc, japc, iw)
    izero = 0
    call SPCGU_PCILU0(nnz, neq, ac, ia, ja, apc, iapc, japc, iw, w, &
                      1.0, izero, 0.5d0)
    call check(error, izero, 1)
  end subroutine test_milu0_return

  !> @brief MILU0 with IZERO of 2 replaces sign-changing pivots
  !<
  subroutine test_milu0_small_pivot(error)
    type(error_type), allocatable, intent(out) :: error
    integer :: ia(neq + 1), ja(nnz), iapc(neq + 1), japc(nnz), iw(neq + 1)
    double precision :: ac(nnz), apc(nnz), w(neq)
    integer :: izero

    call setup(ia, ja, ac, iapc, japc, iw)
    izero = 2
    call SPCGU_PCILU0(nnz, neq, ac, ia, ja, apc, iapc, japc, iw, w, &
                      1.0, izero, 0.5d0)
    call check(error, izero, 0)
    if (allocated(error)) return
    call check(error, apc(2), 1.0d0 / dtiny, rel=.true., thr=tol)
    if (allocated(error)) return
    call check(error, apc(3), 1.0d0 / dtiny, rel=.true., thr=tol)
  end subroutine test_milu0_small_pivot

  !> @brief MILU0 stops shifting at 0.5 and completes the factorization
  !<
  subroutine test_milu0_shift_limit(error)
    type(error_type), allocatable, intent(out) :: error
    integer :: ia(neq + 1), ja(nnz), iapc(neq + 1), japc(nnz), iw(neq + 1)
    double precision :: ac(nnz), apc(nnz), w(neq)

    call setup(ia, ja, ac, iapc, japc, iw)
    call SPCGU_PCU(iout, nnz, neq, neq, nnz, nnz, 3, 1.0, &
                   ac, ia, ja, apc, iapc, japc, iw, w)
    call check(error, apc(1), 1.0d0 / 1.5d0, thr=tol)
    if (allocated(error)) return
    call check(error, apc(2), 1.0d0 / dtiny, rel=.true., thr=tol)
  end subroutine test_milu0_shift_limit

end module TestPcgu
