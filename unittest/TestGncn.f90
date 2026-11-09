!> @brief Unit tests of GNCn array allocation and deallocation
!!
!! The input defines no ghost nodes and room for two parameter ghost nodes,
!! with the options that use SATCn, IRGNCn, and LGNCn turned off.
!<
module TestGncn
  use testdrive, only: check, error_type, new_unittest, unittest_type
  use GLOBAL, only: IOUT, IUNIT, INGNCn
  use PARAMMODULE, only: IPSUM
  use GNCnMODULE, only: GNCn, SATCn, IRGNCn, LGNCn, MXGNCn, MXADJn, &
                        ISYMGNCn, IFLALPHAn, I2Kn
  implicit none
  private
  public :: collect_gncn

  integer, parameter :: iugncn = 36 !< name file index of the GNC package

contains

  !> @brief Collect the GNCn tests
  !<
  subroutine collect_gncn(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("allocate_and_deallocate", test_allocation) &
                ]
  end subroutine collect_gncn

  !> @brief Allocate the optional arrays with their options off, then
  !! deallocate everything
  !<
  subroutine test_allocation(error)
    type(error_type), allocatable, intent(out) :: error
    integer :: inunit, outunit

    ! NPGNCn MXFBPn NGNCNPn MXADJn I2Kn ISYMGNCn IFLALPHAn
    open (newunit=inunit, status='scratch', action='readwrite')
    write (inunit, '(a)') '0 2 0 3 0 1 0'
    rewind (inunit)
    open (newunit=outunit, status='scratch', action='write')

    if (.not. associated(IOUT)) allocate (IOUT)
    if (.not. associated(INGNCn)) allocate (INGNCn)
    if (.not. associated(IPSUM)) allocate (IPSUM)
    if (.not. allocated(IUNIT)) allocate (IUNIT(iugncn))
    IOUT = outunit
    IPSUM = 0
    IUNIT(iugncn) = inunit

    call gncn2disu1ar(iugncn)
    close (inunit)
    close (outunit)

    call check(error, I2Kn == 0 .and. ISYMGNCn == 1 .and. IFLALPHAn == 0, &
               "optional array flags were not read as off")
    if (allocated(error)) return
    call check(error, allocated(SATCn), "SATCn is not allocated")
    if (allocated(error)) return
    call check(error, allocated(IRGNCn), "IRGNCn is not allocated")
    if (allocated(error)) return
    call check(error, allocated(LGNCn), "LGNCn is not allocated")
    if (allocated(error)) return
    call check(error, size(SATCn, 1) == MXADJn + 1 .and. &
               size(SATCn, 2) == MXGNCn, "SATCn has the wrong shape")
    if (allocated(error)) return
    call check(error, size(IRGNCn, 2) == MXADJn .and. &
               size(IRGNCn, 3) == MXGNCn, "IRGNCn has the wrong shape")
    if (allocated(error)) return
    call check(error, size(LGNCn) == MXGNCn, "LGNCn has the wrong size")
    if (allocated(error)) return

    call gncn2disu1da()
    call check(error, .not. allocated(GNCn), "GNCn is still allocated")
    if (allocated(error)) return
    call check(error, .not. allocated(SATCn), "SATCn is still allocated")
    if (allocated(error)) return
    call check(error, .not. allocated(IRGNCn), "IRGNCn is still allocated")
    if (allocated(error)) return
    call check(error, .not. allocated(LGNCn), "LGNCn is still allocated")
    if (allocated(error)) return
    call check(error, .not. associated(ISYMGNCn), &
               "ISYMGNCn is still allocated")
    if (allocated(error)) return
    call check(error, .not. associated(IFLALPHAn), &
               "IFLALPHAn is still allocated")
  end subroutine test_allocation

end module TestGncn
