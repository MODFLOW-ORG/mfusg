!> @brief Unit tests of the PCGU solver class
!!
!! The system is a 4 by 4 grid with a five-point stencil and a known
!! solution. Inactive nodes are identity rows whose right-hand side is the
!! known head, and the active rows keep their connections to them.
!<
module TestPcguSolver
  use testdrive, only: check, error_type, new_unittest, unittest_type
  use PCGUSOLVERMODULE, only: pcgusolver
  implicit none
  private
  public :: collect_pcgusolver

  integer, parameter :: nside = 4
  integer, parameter :: neq = nside * nside
  integer, parameter :: maxnnz = 5 * neq
  integer, parameter :: iout = 6
  real, parameter :: hclose = 1.0e-9
  real, parameter :: rclose = 1.0e-9
  double precision, parameter :: tol = 1.0d-6

contains

  !> @brief Collect the PCGU solver tests
  !<
  subroutine collect_pcgusolver(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("all_options", test_all_options), &
                new_unittest("l2_norm", test_l2_norm), &
                new_unittest("inactive", test_inactive), &
                new_unittest("inactive_ordering", test_inactive_ordering), &
                new_unittest("active_set_change", test_active_set_change), &
                new_unittest("constant_head", test_constant_head) &
                ]
  end subroutine collect_pcgusolver

  !> @brief Fill the grid structure, the coefficients, and the known heads
  !<
  subroutine setup(ia, ja, ac, nnz, htrue)
    integer, intent(out) :: ia(neq + 1), ja(maxnnz)
    double precision, intent(out) :: ac(maxnnz)
    integer, intent(out) :: nnz
    double precision, intent(out) :: htrue(neq)
    integer :: i, j, n, m, k
    integer :: di(4), dj(4)

    di = [-1, 0, 0, 1]
    dj = [0, -1, 1, 0]
    nnz = 0
    do i = 1, nside
      do j = 1, nside
        n = (i - 1) * nside + j
        ia(n) = nnz + 1
        nnz = nnz + 1
        ja(nnz) = n
        ac(nnz) = 4.1d0
        do k = 1, 4
          if (i + di(k) < 1 .or. i + di(k) > nside) cycle
          if (j + dj(k) < 1 .or. j + dj(k) > nside) cycle
          m = (i + di(k) - 1) * nside + j + dj(k)
          nnz = nnz + 1
          ja(nnz) = m
          ac(nnz) = -1.0d0
        end do
        htrue(n) = 10.0d0 + 0.5d0 * i - 0.25d0 * j
      end do
    end do
    ia(neq + 1) = nnz + 1
  end subroutine setup

  !> @brief Make the IBOUND .le. 0 rows identity rows and fill the rhs
  !<
  subroutine boundary(ia, ja, ac, ibound, htrue, rhs)
    integer, intent(in) :: ia(neq + 1), ja(maxnnz)
    double precision, intent(inout) :: ac(maxnnz)
    integer, intent(in) :: ibound(neq)
    double precision, intent(in) :: htrue(neq)
    double precision, intent(out) :: rhs(neq)
    integer :: n, j

    do n = 1, neq
      if (ibound(n) <= 0) then
        ac(ia(n)) = 1.0d0
        do j = ia(n) + 1, ia(n + 1) - 1
          ac(j) = 0.0d0
        end do
      end if
    end do
    call matvec(ia, ja, ac, htrue, rhs)
  end subroutine boundary

  !> @brief Multiply the coefficients by a vector
  !<
  subroutine matvec(ia, ja, ac, x, y)
    integer, intent(in) :: ia(neq + 1), ja(maxnnz)
    double precision, intent(in) :: ac(maxnnz)
    double precision, intent(in) :: x(neq)
    double precision, intent(out) :: y(neq)
    integer :: n, j

    do n = 1, neq
      y(n) = 0.0d0
      do j = ia(n), ia(n + 1) - 1
        y(n) = y(n) + ac(j) * x(ja(j))
      end do
    end do
  end subroutine matvec

  !> @brief Solve and return the largest error in the solution
  !<
  subroutine run(ilinmeth, ipc, iscl, iord, iaord, icnvgopt, ibound, &
                 errmax, icnvg)
    integer, intent(in) :: ilinmeth, ipc, iscl, iord, iaord, icnvgopt
    integer, intent(in) :: ibound(neq)
    double precision, intent(out) :: errmax
    integer, intent(out) :: icnvg
    type(pcgusolver) :: solver
    integer :: ia(neq + 1), ja(maxnnz), nnz, n, innerit
    double precision :: ac(maxnnz), htrue(neq), rhs(neq), hnew(neq)
    real :: relax

    relax = 0.0
    if (ipc == 3) relax = 0.97
    call setup(ia, ja, ac, nnz, htrue)
    call boundary(ia, ja, ac, ibound, htrue, rhs)
    do n = 1, neq
      hnew(n) = 0.0d0
      if (ibound(n) <= 0) hnew(n) = htrue(n)
    end do
    call solver%init(iout, 0, 500, ilinmeth, ipc, iscl, iord, iaord, &
                     icnvgopt, hclose, rclose, relax, neq, nnz, ia, ja)
    call solver%solve(2, 1, ibound, ac, hnew, rhs, icnvg, innerit)
    call solver%destroy()
    errmax = 0.0d0
    do n = 1, neq
      errmax = max(errmax, abs(hnew(n) - htrue(n)))
    end do
  end subroutine run

  !> @brief Every method, preconditioner, scaling, and ordering
  !<
  subroutine test_all_options(error)
    type(error_type), allocatable, intent(out) :: error
    integer :: ibound(neq), ilinmeth, ipc, iscl, iord, icnvg
    double precision :: errmax
    character(len=64) :: msg

    ibound = 1
    do ilinmeth = 1, 2
      do ipc = 0, 3
        do iscl = 0, 2
          do iord = 0, 2
            call run(ilinmeth, ipc, iscl, iord, 0, 0, ibound, errmax, icnvg)
            write (msg, '(a,4i2)') 'ilinmeth ipc iscl iord', &
              ilinmeth, ipc, iscl, iord
            call check(error, errmax < tol, trim(msg))
            if (allocated(error)) return
          end do
        end do
      end do
    end do
  end subroutine test_all_options

  !> @brief The L2 norm option accepts a residual reduced by EPFACT
  !<
  subroutine test_l2_norm(error)
    type(error_type), allocatable, intent(out) :: error
    type(pcgusolver) :: solver
    integer :: ia(neq + 1), ja(maxnnz), ibound(neq), nnz, n, ilinmeth
    integer :: icnvg, innerit
    double precision :: ac(maxnnz), htrue(neq), rhs(neq), hnew(neq), r(neq)
    double precision :: r0, r1

    ibound = 1
    do ilinmeth = 1, 2
      call setup(ia, ja, ac, nnz, htrue)
      call boundary(ia, ja, ac, ibound, htrue, rhs)
      do n = 1, neq
        hnew(n) = 0.0d0
      end do
      call solver%init(iout, 0, 500, ilinmeth, 2, 0, 0, 0, 1, hclose, &
                       1.0e-30, 0.0, neq, nnz, ia, ja)
      call solver%solve(2, 1, ibound, ac, hnew, rhs, icnvg, innerit)
      call solver%destroy()
      call matvec(ia, ja, ac, hnew, r)
      r0 = 0.0d0
      r1 = 0.0d0
      do n = 1, neq
        r0 = r0 + rhs(n) * rhs(n)
        r1 = r1 + (rhs(n) - r(n)) * (rhs(n) - r(n))
      end do
      r0 = sqrt(r0)
      r1 = sqrt(r1)
      call check(error, icnvg, 1)
      if (allocated(error)) return
      ! the second time step accepts a tenfold reduction
      call check(error, r1 <= 0.1d0 * r0, "l2 norm reduction")
      if (allocated(error)) return
    end do
  end subroutine test_l2_norm

  !> @brief Inactive nodes are not solved and keep their heads
  !<
  subroutine test_inactive(error)
    type(error_type), allocatable, intent(out) :: error
    integer :: ibound(neq), ilinmeth, ipc, icnvg
    double precision :: errmax

    ibound = 1
    ibound(1) = 0
    ibound(6) = 0
    ibound(16) = 0
    do ilinmeth = 1, 2
      do ipc = 0, 3
        call run(ilinmeth, ipc, 0, 0, 1, 0, ibound, errmax, icnvg)
        call check(error, errmax < tol, "inactive solution")
        if (allocated(error)) return
      end do
    end do
  end subroutine test_inactive

  !> @brief Inactive nodes with scaling and RCM or minimum degree ordering
  !<
  subroutine test_inactive_ordering(error)
    type(error_type), allocatable, intent(out) :: error
    integer :: ibound(neq), iscl, iord, icnvg
    double precision :: errmax

    ibound = 1
    ibound(2) = 0
    ibound(11) = 0
    do iscl = 1, 2
      do iord = 0, 2
        call run(2, 3, iscl, iord, 1, 0, ibound, errmax, icnvg)
        call check(error, icnvg, 1)
        if (allocated(error)) return
        call check(error, errmax < tol, "inactive ordering solution")
        if (allocated(error)) return
      end do
    end do
  end subroutine test_inactive_ordering

  !> @brief The active node ordering follows a change in IBOUND
  !<
  subroutine test_active_set_change(error)
    type(error_type), allocatable, intent(out) :: error
    type(pcgusolver) :: solver
    integer :: ia(neq + 1), ja(maxnnz), ibound(neq), nnz, n, k
    integer :: icnvg, innerit
    double precision :: ac(maxnnz), htrue(neq), rhs(neq), hnew(neq)
    double precision :: errmax

    call setup(ia, ja, ac, nnz, htrue)
    call solver%init(iout, 0, 500, 1, 2, 0, 0, 1, 0, hclose, rclose, 0.0, &
                     neq, nnz, ia, ja)
    do k = 1, 3
      ibound = 1
      select case (k)
      case (1)
        ibound(4) = 0
      case (2)
        ibound(4) = 0
        ibound(13) = 0
      case (3)
        ibound(7) = 0
      end select
      call setup(ia, ja, ac, nnz, htrue)
      call boundary(ia, ja, ac, ibound, htrue, rhs)
      do n = 1, neq
        hnew(n) = 0.0d0
        if (ibound(n) <= 0) hnew(n) = htrue(n)
      end do
      call solver%solve(2, 1, ibound, ac, hnew, rhs, icnvg, innerit)
      errmax = 0.0d0
      do n = 1, neq
        errmax = max(errmax, abs(hnew(n) - htrue(n)))
      end do
      call check(error, icnvg, 1)
      if (allocated(error)) exit
      call check(error, errmax < tol, "solution after an IBOUND change")
      if (allocated(error)) exit
    end do
    call solver%destroy()
  end subroutine test_active_set_change

  !> @brief Constant-head nodes stay in the active system
  !<
  subroutine test_constant_head(error)
    type(error_type), allocatable, intent(out) :: error
    integer :: ibound(neq), icnvg
    double precision :: errmax

    ibound = 1
    ibound(1) = -1
    ibound(8) = 0
    call run(2, 3, 0, 0, 1, 0, ibound, errmax, icnvg)
    call check(error, icnvg, 1)
    if (allocated(error)) return
    call check(error, errmax < tol, "constant head solution")
  end subroutine test_constant_head

end module TestPcguSolver
