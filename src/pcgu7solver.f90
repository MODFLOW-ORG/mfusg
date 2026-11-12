!> @brief Preconditioned conjugate gradient solver for unstructured grids
!!
!! CG or BiCGSTAB with Jacobi, ILU0, or MILU0 preconditioning, optional
!! scaling and RCM or minimum degree ordering, and an option to solve only
!! the active nodes.
!<
module PCGUSOLVERMODULE
  use PCGULIBMODULE, only: SPCGU_SCALE, SPCGU_PCU, SPCGU_JACA, &
                           SPCGU_ILU0A, SPCGU_TESTCNVG, SPCGU_PCCRS, &
                           SPCGU_MV, SPCGU_DP, SPCGU_CALC_ORDER, &
                           SPCGU_CALC_IAORDER, SPCGU_IAPERM, DPERM, VPERM
  implicit none
  private
  public :: pcgusolver

  double precision, parameter :: DZERO = 0.0d0

  type :: pcgusolver
    integer :: initialized = 0
    ! solver parameters
    integer :: iout = 0
    integer :: maxinner = 0
    integer :: ilinmeth = 0 !< 1 is CG, 2 is BiCGSTAB
    integer :: ipc = 0 !< preconditioner
    integer :: iscl = 0 !< scaling
    integer :: iord = 0 !< RCM or minimum degree ordering
    integer :: iaord = 0 !< solve only the active nodes
    integer :: icnvgopt = 0 !< 1 tests the L2 norm of the residual
    real :: hclose = 0.0
    real :: rclose = 0.0
    real :: relax = 0.0
    ! working variables
    integer :: iacpc = 1 !< recalculate the preconditioner structure
    integer :: iaord_init = 0
    integer :: niter = 0
    double precision :: epfact = DZERO
    double precision :: l2norm0 = DZERO
    ! full and active system sizes
    integer :: nia_all = 0
    integer :: nnz_all = 0
    integer :: nia = 0
    integer :: nnz = 0
    integer, allocatable, dimension(:) :: ia
    integer, allocatable, dimension(:) :: ja
    integer, allocatable, dimension(:) :: ibound !< active set of the ordering
    ! scaled full system
    double precision, allocatable, dimension(:) :: a
    double precision, allocatable, dimension(:) :: x
    double precision, allocatable, dimension(:) :: b
    ! preconditioner
    integer :: niapc = 0
    integer :: njapc = 0
    integer :: nnzapc = 0
    integer, allocatable, dimension(:) :: iapc
    integer, allocatable, dimension(:) :: japc
    double precision, allocatable, dimension(:) :: apc
    integer, allocatable, dimension(:) :: iw
    double precision, allocatable, dimension(:) :: w
    ! scaling
    double precision, allocatable, dimension(:) :: dscale
    double precision, allocatable, dimension(:) :: dscale2
    ! factors that return the residual and head change to model units,
    ! in the order of the system passed to the linear accelerators
    double precision, allocatable, dimension(:) :: rscale
    double precision, allocatable, dimension(:) :: hscale
    ! active node ordering
    integer, allocatable, dimension(:) :: lorder_act
    integer, allocatable, dimension(:) :: iorder_act
    integer, allocatable, dimension(:) :: ia0
    integer, allocatable, dimension(:) :: ja0
    double precision, allocatable, dimension(:) :: a0
    ! rcm and minimum degree ordering
    integer, allocatable, dimension(:) :: lorder
    integer, allocatable, dimension(:) :: iorder
    integer, allocatable, dimension(:) :: id
    ! system passed to the linear accelerators
    integer, allocatable, dimension(:) :: ia1
    integer, allocatable, dimension(:) :: ja1
    double precision, allocatable, dimension(:) :: a1
    ! cg working vectors
    double precision, allocatable, dimension(:) :: d
    double precision, allocatable, dimension(:) :: p
    double precision, allocatable, dimension(:) :: q
    double precision, allocatable, dimension(:) :: z
    ! additional bicgstab working vectors
    integer :: niabcgs = 0
    double precision, allocatable, dimension(:) :: t
    double precision, allocatable, dimension(:) :: v
    double precision, allocatable, dimension(:) :: dhat
    double precision, allocatable, dimension(:) :: phat
    double precision, allocatable, dimension(:) :: qhat
  contains
    procedure :: init
    procedure :: solve
    procedure :: destroy
    procedure, private :: reorder
    procedure, private :: set_cnvg_scale
    procedure, private :: cg
    procedure, private :: bcgs
  end type pcgusolver

contains

  !> @brief Set the solver parameters and allocate the solver storage
  !<
  subroutine init(this, iout, iprsms, maxinner, ilinmeth, ipc, iscl, iord, &
                  iaord, icnvgopt, hclose, rclose, relax, nia, nnz, ia, ja)
    class(pcgusolver), intent(inout) :: this
    integer, intent(in) :: iout
    integer, intent(in) :: iprsms !< 2 prints the node ordering
    integer, intent(in) :: maxinner
    integer, intent(in) :: ilinmeth, ipc, iscl, iord, iaord, icnvgopt
    real, intent(in) :: hclose, rclose, relax
    integer, intent(in) :: nia, nnz
    integer, dimension(nia + 1), intent(in) :: ia
    integer, dimension(nnz), intent(in) :: ja
    integer :: n, i, j, nscl, nreo, nact
    logical :: lreorder
2030 format(1X, A20, 1X, 6(I6, 1X))
2040 format(1X, 20('-'), 1X, 6(6('-'), 1X))
2050 format(1X, 62('-'), /)

    this%initialized = 1
    this%iout = iout
    this%maxinner = maxinner
    this%ilinmeth = ilinmeth
    this%ipc = ipc
    this%iscl = iscl
    this%iord = iord
    this%iaord = iaord
    this%icnvgopt = icnvgopt
    this%hclose = hclose
    this%rclose = rclose
    this%relax = relax
    this%iacpc = 1
    this%iaord_init = 0
    this%niter = 0
    this%nia_all = nia
    this%nnz_all = nnz
    this%nia = nia
    this%nnz = nnz
    lreorder = iaord /= 0 .or. iord /= 0

    ! full system structure
    allocate (this%ia(nia + 1))
    allocate (this%ja(nnz))
    do n = 1, nia + 1
      this%ia(n) = ia(n)
    end do
    do n = 1, nnz
      this%ja(n) = ja(n)
    end do

    ! scaled full system, only needed when the system is reordered
    nreo = 1
    if (lreorder) nreo = nnz
    allocate (this%a(nreo))
    allocate (this%x(nia))
    allocate (this%b(nia))

    ! scaling
    nscl = 1
    if (iscl /= 0) nscl = nia
    allocate (this%dscale(nscl))
    allocate (this%dscale2(nscl))
    allocate (this%rscale(nscl))
    allocate (this%hscale(nscl))
    do n = 1, nscl
      this%dscale(n) = 1.0d0
      this%dscale2(n) = 1.0d0
      this%rscale(n) = 1.0d0
      this%hscale(n) = 1.0d0
    end do

    ! preconditioner
    this%niapc = nia
    this%njapc = nnz
    this%nnzapc = nnz
    if (ipc == 0) then
      this%niapc = 1
      this%njapc = 1
      this%nnzapc = 1
    else if (ipc == 1) then
      this%niapc = 1
      this%njapc = 1
      this%nnzapc = nnz
    end if
    allocate (this%iapc(this%niapc + 1))
    allocate (this%japc(this%njapc))
    allocate (this%apc(this%nnzapc))
    allocate (this%iw(this%niapc + 1))
    allocate (this%w(this%niapc))
    do n = 1, this%nnzapc
      this%apc(n) = DZERO
    end do
    do n = 1, this%niapc
      this%iw(n) = 0
      this%w(n) = DZERO
    end do

    ! active node ordering
    nact = 1
    if (iaord /= 0) nact = nia
    allocate (this%ibound(nact))
    allocate (this%lorder_act(nact))
    allocate (this%iorder_act(nact))
    allocate (this%ia0(nact + 1))
    allocate (this%ja0(nreo))
    allocate (this%a0(nreo))

    ! rcm and minimum degree ordering
    n = 1
    if (iord /= 0) n = nia
    allocate (this%lorder(n))
    allocate (this%iorder(n))
    allocate (this%id(nia))

    ! system passed to the linear accelerators
    allocate (this%ia1(nia + 1))
    allocate (this%ja1(nnz))
    allocate (this%a1(nnz))
    do n = 1, nia + 1
      this%ia1(n) = ia(n)
    end do
    do n = 1, nnz
      this%ja1(n) = ja(n)
    end do

    ! working vectors
    allocate (this%d(nia))
    allocate (this%p(nia))
    allocate (this%q(nia))
    allocate (this%z(nia))
    this%niabcgs = 1
    if (ilinmeth == 2) this%niabcgs = nia
    allocate (this%t(this%niabcgs))
    allocate (this%v(this%niabcgs))
    allocate (this%dhat(this%niabcgs))
    allocate (this%phat(this%niabcgs))
    allocate (this%qhat(this%niabcgs))
    do n = 1, nia
      this%id(n) = 0
      this%d(n) = DZERO
      this%p(n) = DZERO
      this%q(n) = DZERO
      this%z(n) = DZERO
    end do
    do n = 1, this%niabcgs
      this%t(n) = DZERO
      this%v(n) = DZERO
      this%dhat(n) = DZERO
      this%phat(n) = DZERO
      this%qhat(n) = DZERO
    end do

    ! the full system ordering does not change, so calculate it once
    if (iord /= 0 .and. iaord == 0) then
      call SPCGU_CALC_ORDER(iord, nia, nnz, this%ia, this%ja, &
                            this%lorder, this%iorder)
      if (iprsms == 2) then
        do i = 1, nia, 6
          write (iout, 2030) 'ORIGINAL NODE      :', &
            (j, j=i, min(i + 5, nia))
          write (iout, 2040)
          write (iout, 2030) 'REORDERED INDEX    :', &
            (this%lorder(j), j=i, min(i + 5, nia))
          write (iout, 2030) 'REORDERED NODE     :', &
            (this%iorder(j), j=i, min(i + 5, nia))
          write (iout, 2050)
        end do
      end if
    end if
  end subroutine init

  !> @brief Solve the system of equations for HNEW
  !<
  subroutine solve(this, kstp, kiter, ibound, ac, hnew, rhs, icnvg, innerit)
    class(pcgusolver), intent(inout) :: this
    integer, intent(in) :: kstp
    integer, intent(in) :: kiter
    integer, dimension(this%nia_all), intent(in) :: ibound
    double precision, dimension(this%nnz_all), intent(in) :: ac
    double precision, dimension(this%nia_all), intent(inout) :: hnew
    double precision, dimension(this%nia_all), intent(in) :: rhs
    integer, intent(inout) :: icnvg
    integer, intent(inout) :: innerit
    integer :: n, itmax
    double precision :: tv, rmax

    icnvg = 0
    innerit = 0
    if (this%initialized == 0) return
    if (kstp == 1) then
      this%epfact = 0.01d0
    else
      this%epfact = 0.10d0
    end if
    if (kiter == 1) this%niter = 0

    ! fill and scale the system, in a1 when it is not reordered
    do n = 1, this%nia_all
      this%x(n) = hnew(n)
      this%b(n) = rhs(n)
    end do
    if (this%iaord == 0 .and. this%iord == 0) then
      do n = 1, this%nnz_all
        this%a1(n) = ac(n)
      end do
      if (this%iscl /= 0) then
        call SPCGU_SCALE(0, this%iscl, this%nia_all, this%nnz_all, &
                         this%ia, this%ja, this%a1, this%x, this%b, &
                         this%dscale, this%dscale2)
        call this%set_cnvg_scale()
      end if
    else
      do n = 1, this%nnz_all
        this%a(n) = ac(n)
      end do
      if (this%iscl /= 0) then
        call SPCGU_SCALE(0, this%iscl, this%nia_all, this%nnz_all, &
                         this%ia, this%ja, this%a, this%x, this%b, &
                         this%dscale, this%dscale2)
        call this%set_cnvg_scale()
      end if
      call this%reorder(ibound)
    end if

    ! preconditioner structure changes only with the ordering
    if (this%iacpc /= 0) then
      if (this%ipc == 2 .or. this%ipc == 3) then
        call SPCGU_PCCRS(this%nia, this%nnz, this%ia1, this%ja1, &
                         this%iapc, this%japc, this%iw)
      end if
      this%iacpc = 0
    end if
    call SPCGU_PCU(this%iout, this%nnz, this%nia, &
                   this%niapc, this%njapc, this%nnzapc, &
                   this%ipc, this%relax, &
                   this%a1, this%ia1, this%ja1, &
                   this%apc, this%iapc, this%japc, this%iw, this%w)

    ! initial residual
    do n = 1, this%nia
      this%d(n) = DZERO
      this%p(n) = DZERO
      this%q(n) = DZERO
      this%z(n) = DZERO
    end do
    call SPCGU_MV(this%nnz, this%nia, this%a1, this%x, this%d, &
                  this%ia1, this%ja1)
    rmax = DZERO
    this%l2norm0 = DZERO
    do n = 1, this%nia
      tv = this%b(n) - this%d(n)
      this%d(n) = tv
      if (abs(tv) > rmax) rmax = abs(tv)
      if (this%iscl /= 0) tv = tv * this%rscale(n)
      this%l2norm0 = this%l2norm0 + tv * tv
    end do
    this%l2norm0 = sqrt(this%l2norm0)
    itmax = this%maxinner
    if (rmax == DZERO) itmax = 0

    if (this%ilinmeth == 1) then
      call this%cg(icnvg, itmax, innerit)
    else if (this%ilinmeth == 2) then
      call this%bcgs(icnvg, itmax, innerit)
    end if

    ! back permute and unscale the solution
    if (this%iord /= 0) then
      call VPERM(this%nia, this%x, this%iorder)
    end if
    if (this%iaord /= 0) then
      call VPERM(this%nia_all, this%x, this%iorder_act)
    end if
    if (this%iscl /= 0) then
      do n = 1, this%nia_all
        this%x(n) = this%x(n) * this%dscale2(n)
      end do
    end if
    do n = 1, this%nia_all
      hnew(n) = this%x(n)
    end do
  end subroutine solve

  !> @brief Set the factors that return the residual and head change to
  !! model units, before the system is reordered
  !<
  subroutine set_cnvg_scale(this)
    class(pcgusolver), intent(inout) :: this
    integer :: n

    do n = 1, this%nia_all
      this%rscale(n) = 1.0d0 / this%dscale(n)
      this%hscale(n) = this%dscale2(n)
    end do
  end subroutine set_cnvg_scale

  !> @brief Reorder the scaled system into a1, ia1, and ja1
  !!
  !! The active node ordering is recalculated when IBOUND changes.
  !<
  subroutine reorder(this, ibound)
    class(pcgusolver), intent(inout) :: this
    integer, dimension(this%nia_all), intent(in) :: ibound
    integer :: n, j, jcol
    logical :: lchanged

    if (this%iaord == 0) then
      call DPERM(this%nia, this%a, this%ja, this%ia, this%a1, this%ja1, &
                 this%ia1, this%lorder, this%id, 1)
      call VPERM(this%nia, this%x, this%lorder)
      call VPERM(this%nia, this%b, this%lorder)
      if (this%iscl /= 0) then
        call VPERM(this%nia, this%rscale, this%lorder)
        call VPERM(this%nia, this%hscale, this%lorder)
      end if
      return
    end if

    lchanged = this%iaord_init == 0
    if (.not. lchanged) then
      do n = 1, this%nia_all
        if ((ibound(n) == 0) .neqv. (this%ibound(n) == 0)) then
          lchanged = .true.
          exit
        end if
      end do
    end if

    if (lchanged) then
      do n = 1, this%nia_all
        this%ibound(n) = ibound(n)
      end do
      call SPCGU_CALC_IAORDER(this%nia_all, this%nnz_all, this%ia, this%ja, &
                              this%a, this%ibound, this%nia, this%nnz, &
                              this%ia0, this%ja0, this%a0, &
                              this%lorder_act, this%iorder_act)
      if (this%iord /= 0) then
        call SPCGU_CALC_ORDER(this%iord, this%nia, this%nnz, this%ia0, &
                              this%ja0, this%lorder, this%iorder)
      else
        do n = 1, this%nia + 1
          this%ia1(n) = this%ia0(n)
        end do
        do n = 1, this%nnz
          this%ja1(n) = this%ja0(n)
        end do
      end if
      this%iaord_init = 1
      this%iacpc = 1
    else
      call SPCGU_IAPERM(this%nia_all, this%nnz_all, this%ia, this%ja, &
                        this%a, this%ibound, this%nia, this%a0, &
                        this%iorder_act)
    end if

    ! move the dropped inactive columns to the right-hand side
    do n = 1, this%nia_all
      if (this%ibound(n) == 0) cycle
      do j = this%ia(n) + 1, this%ia(n + 1) - 1
        jcol = this%ja(j)
        if (this%ibound(jcol) == 0) then
          this%b(n) = this%b(n) - this%a(j) * this%x(jcol)
        end if
      end do
    end do
    call VPERM(this%nia_all, this%x, this%lorder_act)
    call VPERM(this%nia_all, this%b, this%lorder_act)
    if (this%iscl /= 0) then
      call VPERM(this%nia_all, this%rscale, this%lorder_act)
      call VPERM(this%nia_all, this%hscale, this%lorder_act)
    end if

    if (this%iord /= 0) then
      call DPERM(this%nia, this%a0, this%ja0, this%ia0, this%a1, this%ja1, &
                 this%ia1, this%lorder, this%id, 1)
      call VPERM(this%nia, this%x, this%lorder)
      call VPERM(this%nia, this%b, this%lorder)
      if (this%iscl /= 0) then
        call VPERM(this%nia, this%rscale, this%lorder)
        call VPERM(this%nia, this%hscale, this%lorder)
      end if
    else
      do n = 1, this%nnz
        this%a1(n) = this%a0(n)
      end do
    end if
  end subroutine reorder

  !> @brief Preconditioned conjugate gradient
  !<
  subroutine cg(this, icnvg, itmax, innerit)
    class(pcgusolver), intent(inout) :: this
    integer, intent(inout) :: icnvg
    integer, intent(in) :: itmax
    integer, intent(inout) :: innerit
    integer :: n, iiter
    double precision :: dhclose, drclose
    double precision :: tv, deltax, rmax, l2norm, rcnvg
    double precision :: alpha, beta, rho, rho0, machprec

    innerit = 0
    rho0 = DZERO
    machprec = epsilon(DZERO)
    dhclose = dble(this%hclose)
    drclose = dble(this%rclose)
    inner: do iiter = 1, itmax
      innerit = innerit + 1
      this%niter = this%niter + 1
      select case (this%ipc)
      case (0)
        do n = 1, this%nia
          this%z(n) = this%d(n)
        end do
      case (1)
        call SPCGU_JACA(this%nia, this%apc, this%d, this%z)
      case (2, 3)
        call SPCGU_ILU0A(this%nnz, this%nia, this%apc, this%iapc, &
                         this%japc, this%d, this%z)
      end select
      rho = SPCGU_DP(this%nia, this%d, this%z)
      if (iiter == 1) then
        do n = 1, this%nia
          this%p(n) = this%z(n)
        end do
      else
        beta = rho / rho0
        do n = 1, this%nia
          this%p(n) = this%z(n) + beta * this%p(n)
        end do
      end if
      call SPCGU_MV(this%nnz, this%nia, this%a1, this%p, this%q, &
                    this%ia1, this%ja1)
      alpha = rho / SPCGU_DP(this%nia, this%p, this%q)
      deltax = DZERO
      rmax = DZERO
      l2norm = DZERO
      do n = 1, this%nia
        tv = alpha * this%p(n)
        this%x(n) = this%x(n) + tv
        if (this%iscl /= 0) tv = tv * this%hscale(n)
        deltax = max(abs(tv), deltax)
        tv = this%d(n) - alpha * this%q(n)
        this%d(n) = tv
        if (this%iscl /= 0) tv = tv * this%rscale(n)
        rmax = max(abs(tv), rmax)
        l2norm = l2norm + tv * tv
      end do
      l2norm = sqrt(l2norm)
      if (this%icnvgopt == 1) then
        rcnvg = l2norm
      else
        rcnvg = rmax
      end if
      call SPCGU_TESTCNVG(this%icnvgopt, icnvg, deltax, rcnvg, &
                          this%l2norm0, this%epfact, dhclose, drclose)
      if (rcnvg == DZERO) icnvg = 1
      if (icnvg == 1) exit inner
      if (abs(rho - rho0) < machprec) exit inner
      rho0 = rho
    end do inner
  end subroutine cg

  !> @brief Preconditioned bi-conjugate gradient stabilized
  !<
  subroutine bcgs(this, icnvg, itmax, innerit)
    class(pcgusolver), intent(inout) :: this
    integer, intent(inout) :: icnvg
    integer, intent(in) :: itmax
    integer, intent(inout) :: innerit
    integer :: n, iiter
    double precision :: dhclose, drclose
    double precision :: tv, deltax, rmax, l2norm, rcnvg
    double precision :: alpha, alpha0, beta, rho, rho0, omega, omega0
    double precision :: machprec, numer, denom

    innerit = 0
    machprec = epsilon(DZERO)
    dhclose = dble(this%hclose)
    drclose = dble(this%rclose)
    alpha = DZERO
    alpha0 = DZERO
    beta = DZERO
    rho = DZERO
    rho0 = DZERO
    omega = DZERO
    omega0 = DZERO
    do n = 1, this%nia
      this%dhat(n) = this%d(n)
    end do
    inner: do iiter = 1, itmax
      innerit = innerit + 1
      this%niter = this%niter + 1
      rho = SPCGU_DP(this%nia, this%dhat, this%d)
      if (iiter == 1) then
        do n = 1, this%nia
          this%p(n) = this%d(n)
        end do
      else
        beta = (rho / rho0) * (alpha0 / omega0)
        do n = 1, this%nia
          this%p(n) = this%d(n) + beta * (this%p(n) - omega * this%v(n))
        end do
      end if
      select case (this%ipc)
      case (0)
        do n = 1, this%nia
          this%phat(n) = this%p(n)
        end do
      case (1)
        call SPCGU_JACA(this%nia, this%apc, this%p, this%phat)
      case (2, 3)
        call SPCGU_ILU0A(this%nnz, this%nia, this%apc, this%iapc, &
                         this%japc, this%p, this%phat)
      end select
      call SPCGU_MV(this%nnz, this%nia, this%a1, this%phat, this%v, &
                    this%ia1, this%ja1)
      denom = SPCGU_DP(this%nia, this%dhat, this%v)
      alpha = rho / denom
      do n = 1, this%nia
        this%q(n) = this%d(n) - alpha * this%v(n)
      end do
      ! q below machine precision solves the system with the alpha step,
      ! and omega is undefined, so take that step and stop
      rmax = DZERO
      do n = 1, this%nia
        tv = this%q(n)
        if (this%iscl /= 0) tv = tv * this%rscale(n)
        if (abs(tv) > abs(rmax)) rmax = tv
      end do
      if (abs(rmax) <= machprec) then
        deltax = DZERO
        do n = 1, this%nia
          tv = alpha * this%phat(n)
          this%x(n) = this%x(n) + tv
          if (this%iscl /= 0) tv = tv * this%hscale(n)
          if (abs(tv) > abs(deltax)) deltax = tv
        end do
        call SPCGU_TESTCNVG(this%icnvgopt, icnvg, deltax, abs(rmax), &
                            this%l2norm0, this%epfact, dhclose, drclose)
        exit inner
      end if
      select case (this%ipc)
      case (0)
        do n = 1, this%nia
          this%qhat(n) = this%q(n)
        end do
      case (1)
        call SPCGU_JACA(this%nia, this%apc, this%q, this%qhat)
      case (2, 3)
        call SPCGU_ILU0A(this%nnz, this%nia, this%apc, this%iapc, &
                         this%japc, this%q, this%qhat)
      end select
      call SPCGU_MV(this%nnz, this%nia, this%a1, this%qhat, this%t, &
                    this%ia1, this%ja1)
      numer = SPCGU_DP(this%nia, this%t, this%q)
      denom = SPCGU_DP(this%nia, this%t, this%t)
      denom = denom + sign(machprec, denom)
      omega = numer / denom
      deltax = DZERO
      rmax = DZERO
      l2norm = DZERO
      do n = 1, this%nia
        tv = alpha * this%phat(n) + omega * this%qhat(n)
        this%x(n) = this%x(n) + tv
        if (this%iscl /= 0) tv = tv * this%hscale(n)
        if (abs(tv) > abs(deltax)) deltax = tv
        tv = this%q(n) - omega * this%t(n)
        this%d(n) = tv
        if (this%iscl /= 0) tv = tv * this%rscale(n)
        if (abs(tv) > abs(rmax)) rmax = tv
        l2norm = l2norm + tv * tv
      end do
      l2norm = sqrt(l2norm)
      if (this%icnvgopt == 1) then
        rcnvg = l2norm
      else
        rcnvg = rmax
      end if
      call SPCGU_TESTCNVG(this%icnvgopt, icnvg, deltax, rcnvg, &
                          this%l2norm0, this%epfact, dhclose, drclose)
      if (rcnvg == DZERO) icnvg = 1
      if (icnvg == 1) exit inner
      if (abs(rho - rho0) < machprec) exit inner
      if (abs(alpha - alpha0) < machprec) exit inner
      if (abs(omega - omega0) < machprec) exit inner
      rho0 = rho
      alpha0 = alpha
      omega0 = omega
    end do inner
  end subroutine bcgs

  !> @brief Deallocate the solver storage
  !<
  subroutine destroy(this)
    class(pcgusolver), intent(inout) :: this

    if (this%initialized == 0) return
    deallocate (this%ia, this%ja, this%ibound)
    deallocate (this%a, this%x, this%b)
    deallocate (this%iapc, this%japc, this%apc, this%iw, this%w)
    deallocate (this%dscale, this%dscale2, this%rscale, this%hscale)
    deallocate (this%lorder_act, this%iorder_act)
    deallocate (this%ia0, this%ja0, this%a0)
    deallocate (this%lorder, this%iorder, this%id)
    deallocate (this%ia1, this%ja1, this%a1)
    deallocate (this%d, this%p, this%q, this%z)
    deallocate (this%t, this%v, this%dhat, this%phat, this%qhat)
    this%initialized = 0
  end subroutine destroy

end module PCGUSOLVERMODULE
