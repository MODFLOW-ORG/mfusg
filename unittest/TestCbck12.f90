!> @brief Unit tests of the anisotropic conductance in CBCK12
!!
!! Two cells in layer 1 joined by connection 1, with harmonic averaging
!! (LAYAVG = 0) and horizontal anisotropy.
!<
module TestCbck12
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use testdrive, only: check, error_type, new_unittest, unittest_type
  use GLOBAL, only: NODES, ARAD, JAS
  use GWFBCFMODULE, only: LAYAVG, CHANI, HANI, IHANISO
  implicit none
  private
  public :: collect_cbck12

  double precision, parameter :: tol = 1.0d-5

contains

  !> @brief Collect the CBCK12 tests
  !<
  subroutine collect_cbck12(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("x_direction", test_x_direction), &
                new_unittest("rotated", test_rotated), &
                new_unittest("zero_conductivity", test_zero_conductivity), &
                new_unittest("zero_without_weight", test_zero_without_weight) &
                ]
  end subroutine collect_cbck12

  !> @brief Set up two cells with the anisotropy ratio and connection angle
  !<
  subroutine setup(chani1, angle)
    real, intent(in) :: chani1 !< anisotropy ratio, or < 0 to use HANI
    real, intent(in) :: angle !< connection angle in radians

    if (.not. associated(NODES)) allocate (NODES)
    if (.not. associated(IHANISO)) allocate (IHANISO)
    if (.not. allocated(ARAD)) allocate (ARAD(1))
    if (.not. allocated(JAS)) allocate (JAS(1))
    if (.not. allocated(LAYAVG)) allocate (LAYAVG(1))
    if (.not. allocated(CHANI)) allocate (CHANI(1))
    if (.not. allocated(HANI)) allocate (HANI(2))
    NODES = 2
    IHANISO = 1
    JAS(1) = 1
    ARAD(1) = angle
    LAYAVG(1) = 0
    CHANI(1) = chani1
    HANI(1) = 1.0
    HANI(2) = 1.0
  end subroutine setup

  !> @brief Harmonic mean of two transmissivities over the connection
  !<
  function hmean(t1, t2, el1, el2) result(t)
    double precision, intent(in) :: t1, t2
    real, intent(in) :: el1, el2
    double precision :: t
    t = t1 * t2 * (el1 + el2) / (t1 * el2 + t2 * el1)
  end function hmean

  !> @brief Angle 0 gives the x-direction mean
  !<
  subroutine test_x_direction(error)
    type(error_type), allocatable, intent(out) :: error
    real :: hk(2), el1, el2
    double precision :: anum, th1, th2, expected

    call setup(0.5, 0.0)
    hk(1) = 10.0
    hk(2) = 2.0
    el1 = 50.0
    el2 = 100.0
    th1 = 4.0d0
    th2 = 6.0d0
    call cbck12(1, 1, 2, 1, el1, el2, anum, hk, th1, th2)
    expected = hmean(dble(hk(1)) * th1, dble(hk(2)) * th2, el1, el2)
    call check(error, anum, expected, rel=.true., thr=tol)
  end subroutine test_x_direction

  !> @brief A 45 degree connection combines the x and y means
  !<
  subroutine test_rotated(error)
    type(error_type), allocatable, intent(out) :: error
    real :: hk(2), el1, el2
    double precision :: anum, th1, th2, tx, ty, expected

    call setup(0.5, real(atan(1.0d0)))
    hk(1) = 10.0
    hk(2) = 2.0
    el1 = 50.0
    el2 = 100.0
    th1 = 4.0d0
    th2 = 6.0d0
    call cbck12(1, 1, 2, 1, el1, el2, anum, hk, th1, th2)
    tx = hmean(dble(hk(1)) * th1, dble(hk(2)) * th2, el1, el2)
    ty = hmean(0.5d0 * hk(1) * th1, 0.5d0 * hk(2) * th2, el1, el2)
    expected = 1.0d0 / (0.5d0 / tx + 0.5d0 / ty)
    call check(error, anum, expected, rel=.true., thr=tol)
  end subroutine test_rotated

  !> @brief A cell with zero conductivity gives zero conductance
  !<
  subroutine test_zero_conductivity(error)
    type(error_type), allocatable, intent(out) :: error
    real :: hk(2), el1, el2
    double precision :: anum, th1, th2

    call setup(0.5, real(atan(1.0d0)))
    hk(1) = 0.0
    hk(2) = 2.0
    el1 = 50.0
    el2 = 100.0
    th1 = 4.0d0
    th2 = 6.0d0
    call cbck12(1, 1, 2, 1, el1, el2, anum, hk, th1, th2)
    call check(error, ieee_is_finite(anum), "conductance is not finite")
    if (allocated(error)) return
    call check(error, anum, 0.0d0)
  end subroutine test_zero_conductivity

  !> @brief A zero y-direction mean is ignored when the angle is 0
  !<
  subroutine test_zero_without_weight(error)
    type(error_type), allocatable, intent(out) :: error
    real :: hk(2), el1, el2
    double precision :: anum, th1, th2, expected

    ! negative CHANI uses HANI, which is zero in both cells
    call setup(-1.0, 0.0)
    HANI(1) = 0.0
    HANI(2) = 0.0
    hk(1) = 10.0
    hk(2) = 2.0
    el1 = 50.0
    el2 = 100.0
    th1 = 4.0d0
    th2 = 6.0d0
    call cbck12(1, 1, 2, 1, el1, el2, anum, hk, th1, th2)
    call check(error, ieee_is_finite(anum), "conductance is not finite")
    if (allocated(error)) return
    expected = hmean(dble(hk(1)) * th1, dble(hk(2)) * th2, el1, el2)
    call check(error, anum, expected, rel=.true., thr=tol)
  end subroutine test_zero_without_weight

end module TestCbck12
