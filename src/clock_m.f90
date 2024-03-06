module clock_m
  use, intrinsic :: iso_fortran_env
  implicit none
  private
  public :: clock
  integer(int32), parameter :: clock_max_state_limit = 50 ! 実装的に 50^5 が限界.
  real(real64), parameter :: pi = 4 * atan(1.0_real64)
  type :: clock
     private
     integer(int32) :: max_state_
     integer(int64) :: nx_, ny_, nall_
     real(real64) :: beta_
     real(real64) :: pi_state_inv_
     integer(int32), allocatable :: spins_(:)
     real(real64), allocatable :: energy_table_(:, :, :, :, :)
     real(real64), allocatable :: ws_(:, :, :, :, :, :)
   contains
     !> initializer.
     procedure, pass :: init => init_clock
     !> setter.
     procedure, pass :: set_ising_allup => set_ising_allup_clock
     procedure, pass :: set_ising_mixed_phase => set_ising_mixed_phase_clock
     procedure, pass :: set_ising_random => set_ising_random_clock
     procedure, pass :: set_kbt => set_kbt_clock
     procedure, pass :: set_beta => set_beta_clock
     !> generator.
     procedure, pass, private :: generate_random_spin => generate_random_spin_clock
     procedure, pass, private :: generate_other_random_spin => generate_other_random_spin_clock
     !> updater.
     procedure, pass :: update => update_clock
     procedure, pass, private :: update_onesite => update_onesite_clock
     procedure, pass, private :: update_norishiro => update_norishiro_clock
     procedure, pass, private :: update_energy_table => update_energy_table_clock
     procedure, pass, private :: update_ws => update_ws_clock
     !> calculator.
     procedure, pass :: calc_energy_summ => calc_total_energy_clock
     procedure, pass :: calc_magne_summ => calc_total_magne_clock
     !> getter.
     procedure, pass :: nx => nx_clock
     procedure, pass :: ny => ny_clock
     procedure, pass :: nall => nall_clock
     procedure, pass :: max_state => max_state_clock
     procedure, pass :: kbt => kbt_clock
     procedure, pass :: beta => beta_clock
     procedure, pass, private :: norishiro_begin => norishiro_begin_clock
     procedure, pass, private :: norishiro_end => norishiro_end_clock
  end type clock
contains
  !> init_clock: Initialize clock once.
  pure subroutine init_clock(this, nx, ny, kbt, state)
    class(clock), intent(inout) :: this
    integer(int64), intent(in) :: nx, ny
    real(real64), intent(in) :: kbt
    integer(int32), intent(in) :: state
    if (allocated(this%spins_)) return
    if (is_even(nx) .or. (.not. is_even(ny))) then
       error stop "The parity of size must be (x, y) == (odd, even)."
    end if
    if (state > clock_max_state_limit) then
       error stop "The state of Clock must be smaller than 50 for this object."
    end if
    this%nx_ = nx
    this%ny_ = ny
    this%nall_ = nx * ny
    this%max_state_ = state
    this%pi_state_inv_ = 2 * pi / state
    allocate(this%spins_(this%norishiro_begin() : this%norishiro_end()))
    call this%set_ising_allup()
    allocate(this%energy_table_(0:state-1, 0:state-1, 0:state-1, 0:state-1, 0:state-1))
    allocate(this%ws_(0:state-1, 0:state-1, 0:state-1, 0:state-1, 0:state-1, 0:state-1))
    call this%set_kbt(kbt)
  contains
    pure logical function is_even(v) result(res)
      integer(int64), intent(in) :: v
      res = iand(v, b'1') == 0_int64
    end function is_even
  end subroutine init_clock
  !> set_ising_allup_clock: Set spins `1`.
  pure subroutine set_ising_allup_clock(this)
    class(clock), intent(inout) :: this
    this%spins_(:) = 0_int32
  end subroutine set_ising_allup_clock
  !> set_ising_allup_clock: Set spin to `1` in the first half of the region; in the other half, set it randomly.
  impure subroutine set_ising_mixed_phase_clock(this)
    class(clock), intent(inout) :: this
    integer(int64) :: half
    real(real64), allocatable :: r(:)
    integer(int32) :: i
    !> ..... [16 ~ 20]
    !> ..... [11 ~ 15]
    !> ..... [ 6 ~ 10] 10 == (ny / 2) * nx
    !> ..... [ 1 ~  5]
    half = (this%ny_ / 2) * this%nx_
    this%spins_(1:half) = 0_int32
    allocate(r(half + 1:this%nall_))
    call random_number(r)
    do i = half + 1, this%nall_
       this%spins_(i) = this%generate_random_spin()
    end do
    call this%update_norishiro()
  end subroutine set_ising_mixed_phase_clock
  !> set_ising_random_clock: Set spin to an integer between `1` and `this%max_state_` randomly.
  impure subroutine set_ising_random_clock(this)
    class(clock), intent(inout) :: this
    real(real64), allocatable :: r(:)
    allocate(r(1:this%nall_))
    call random_number(r)
    this%spins_(1:this%nall_) = floor(r * this%max_state_) + 1
    call this%update_norishiro()
  end subroutine set_ising_random_clock
  !> set_kbt_clock: Set parameter `beta` as `1 / kbt`.
  pure subroutine set_kbt_clock(this, kbt)
    class(clock), intent(inout) :: this
    real(real64), intent(in) :: kbt
    call this%set_beta(1 / kbt)
  end subroutine set_kbt_clock
  !> set_beta_clock: Set parameter `beta` and update `this%exparr_`.
  pure subroutine set_beta_clock(this, beta)
    class(clock), intent(inout) :: this
    real(real64), intent(in) :: beta
    this%beta_ = beta
    call this%update_ws()
  end subroutine set_beta_clock
  pure subroutine update_energy_table_clock(this)
    class(clock), intent(inout) :: this
    integer(int32) :: center, i, j, k, l
    do center = 0, this%max_state_ - 1
       do l = 0, this%max_state_ - 1
          do k = 0, this%max_state_ - 1
             do j = 0, this%max_state_ - 1
                do i = 0, this%max_state_ - 1
                   this%energy_table_(i, j, k, l, center) = calc_local_energy(i, j, k, l, center)
                end do
             end do
          end do
       end do
    end do
  contains
    pure real(real64) function calc_local_energy(i, j, k, l, center) result(res)
      integer(int32), intent(in) :: i, j, k, l, center
      res = - sum(cos(this%pi_state_inv_ * [i - center, j - center, k - center, l - center]))
    end function calc_local_energy
  end subroutine update_energy_table_clock
  pure subroutine update_ws_clock(this)
    class(clock), intent(inout) :: this
    integer(int32) :: center_before, center_after, i, j, k, l
    call this%update_energy_table()
    do center_after = 0, this%max_state_ - 1
       do center_before = 0, this%max_state_ - 1
          do l = 0, this%max_state_ - 1
             do k = 0, this%max_state_ - 1
                do j = 0, this%max_state_ - 1
                   do i = 0, this%max_state_ - 1
                      associate(delta_e => &
                           & this%energy_table_(i, j, k, l, center_after) - this%energy_table_(i, j, k, l, center_before))
                        if (delta_e <= 0.0_real64) then
                           this%ws_(i, j, k, l, center_before, center_after) = 1.0_real64
                        else
                           this%ws_(i, j, k, l, center_before, center_after) = exp(- this%beta_ * delta_e)
                        end if
                      end associate
                   end do
                end do
             end do
          end do
       end do
    end do
  end subroutine update_ws_clock
  !> generate_random_spin_clock: Generate the spin of clock model at random.
  impure integer(int32) function generate_random_spin_clock(this) result(res)
    class(clock), intent(in) :: this
    real(real64) :: r
    call random_number(r)
    res = floor(r * this%max_state_)
  end function generate_random_spin_clock
  !> generate_other_random_spin_clock: Generate the spin of clock model which is not equal `state` at random.
  impure integer(int32) function generate_other_random_spin_clock(this, state) result(res)
    class(clock), intent(in) :: this
    integer(int32), intent(in) :: state
    real(real64) :: r
    call random_number(r)
    res = floor(r * (this%max_state_ - 1)) + state + 1
    if (res >= this%max_state_) &
         & res = res - this%max_state_
  end function  generate_other_random_spin_clock
  !> update_clock: Update the system by Metropolis method.
  impure subroutine update_clock(this)
    class(clock), intent(inout) :: this
    real(real64), allocatable :: r(:)
    integer(int64) :: i, j
    allocate(r(this%nall_))
    call random_number(r)
    do j = 1, 2
       do i = j, this%nall_, 2
          call this%update_onesite(i, r(i), this%generate_other_random_spin(this%spins_(i)))
       end do
       call this%update_norishiro()
    end do
  end subroutine update_clock
  !> update_onesite_clock: Update a spin of the system.
  pure subroutine update_onesite_clock(this, idx, r, next_state)
    class(clock), intent(inout) :: this
    integer(int64), intent(in) :: idx
    real(real64), intent(in) :: r
    integer(int32), intent(in) :: next_state
    if (r < this%ws_( this%spins_(idx + 1), this%spins_(idx - 1), &
                    & this%spins_(idx + this%nx_), this%spins_(idx - this%nx_), &
                    & this%spins_(idx), next_state )) then
       this%spins_(idx) = next_state
    end if
  end subroutine update_onesite_clock
  !> update_norishiro_clock: Update norishiro.
  pure subroutine update_norishiro_clock(this)
    class(clock), intent(inout) :: this
    integer(int64) :: i
    do i = 1_int64, this%nx_
       this%spins_(this%norishiro_begin() + i - 1) = this%spins_(this%nall_ - this%nx_ + i)
       this%spins_(this%norishiro_end() - this%nx_ + i) = this%spins_(i)
    end do
  end subroutine update_norishiro_clock

  !> calc_total_energy_clock: Calculate the total energy.
  pure integer(int64) function calc_total_energy_clock(this) result(res)
    class(clock), intent(in) :: this
    integer(int64) :: i
    res = 0_int64
    do i = 1_int64, this%nall_
       res = res - this%energy_table_(this%spins_(i - 1), this%spins_(i + 1), &
            & this%spins_(i - this%nx_), this%spins_(i + this%nx_), this%spins_(i))
    end do
    res = res / 2
  end function calc_total_energy_clock
  !> calc_total_magne_clock: Calculate the total magne.
  pure integer(int64) function calc_total_magne_clock(this) result(res)
    class(clock), intent(in) :: this
    integer(int64) :: i
    res = 0_int64
    do i = 1_int64, this%nall_
       res = res + cos(this%pi_state_inv_ * this%spins_(i))
    end do
  end function calc_total_magne_clock

  !> nx_clock: Return size of `x` of the system.
  pure integer(int64) function nx_clock(this) result(res)
    class(clock), intent(in) :: this
    res = this%nx_
  end function nx_clock
  !> ny_clock: Return size of `y` of the system.
  pure integer(int64) function ny_clock(this) result(res)
    class(clock), intent(in) :: this
    res = this%ny_
  end function ny_clock
  !> nall_clock: Return size of the system.
  pure integer(int64) function nall_clock(this) result(res)
    class(clock), intent(in) :: this
    res = this%nall_
  end function nall_clock
  !> max_state_clock: Return maximum state of the system.
  pure integer(int64) function max_state_clock(this) result(res)
    class(clock), intent(in) :: this
    res = this%max_state_
  end function max_state_clock
  !> kbt_clock: Return temperature of the system.
  pure real(real64) function kbt_clock(this) result(res)
    class(clock), intent(in) :: this
    res = 1 / this%beta_
  end function kbt_clock
  !> beta_clock: Return inverse temperature of the system.
  pure real(real64) function beta_clock(this) result(res)
    class(clock), intent(in) :: this
    res = this%beta_
  end function beta_clock
  !> norishiro_begin_clock: Return start index of `this%spins_(:)`.
  pure integer(int64) function norishiro_begin_clock(this) result(res)
    class(clock), intent(in) :: this
    res = 1 - this%nx_
  end function norishiro_begin_clock
  !> norishiro_end_clock: Return end index of `this%spins_(:)`.
  pure integer(int64) function norishiro_end_clock(this) result(res)
    class(clock), intent(in) :: this
    res = this%nall_ + this%nx_
  end function norishiro_end_clock
end module clock_m
