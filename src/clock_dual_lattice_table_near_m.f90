module sixclock_periodic_dual_lattice_table_near_m
  use, intrinsic :: iso_fortran_env
  use msmt19937
  implicit none
  private
  real(real64), parameter :: pi = 4 * atan(1.0d0)

  integer(int32), parameter :: mstate = 6
  real(real64), parameter :: pi_state_inv = 2 * pi / mstate

  integer(int64), parameter :: nx = 100_int64, ny = nx, nall = nx * ny
  real(real64), parameter :: nall_inv = 1d0 / nall
  real(real64), parameter :: kbt = 0.9d0, beta = 1 / kbt

  real(real64) :: magne_table(0:mstate-1)
  real(real64) :: energy_table_ru(0:mstate-1, 0:mstate-1, 0:mstate-1)
  real(real64) :: energy_table(0:mstate-1, 0:mstate-1, 0:mstate-1, 0:mstate-1, 0:mstate-1)
  real(real64) :: energy_diff_table(0:mstate-1, 0:mstate-1, 0:mstate-1, 0:mstate-1, 0:mstate-1, 0:mstate-1)
  real(real64) :: probability_table(0:mstate-1, 0:mstate-1, 0:mstate-1, 0:mstate-1, 0:mstate-1, 0:mstate-1)

  integer(int32), allocatable :: sixclock_even(:, :), sixclock_odd(:, :)

  real(real64), allocatable :: rnds(:, :, :)

  integer(int32), parameter :: nd = 4, up = 1, right = 2, down = 3, left = 4 !> 1: up, 2: right, 3: down, 4: left.
  integer(int64), parameter :: dx(nd, 0:1) = &
       & reshape([[0_int64, 1_int64, 0_int64, 0_int64], [0_int64, 0_int64, 0_int64, -1_int64]], shape = [nd, 2])

  public :: mstate, nx, ny, nall, kbt, beta
  public :: init_sixclock, init_sixclock_order, update_metropolis, calc_energy, calc_magne

contains
  !> init_sixclock: Initialize 6-state clock model.
  impure subroutine init_sixclock()
    allocate(sixclock_even(nx / 2, ny))
    allocate(sixclock_odd(nx / 2, ny))
    allocate(rnds(1:2, nx, ny))
    call init_tables()
  end subroutine init_sixclock
  !> init_tables: Initialize `energy_table` and `probability_table`.
  impure subroutine init_tables()
    real(real64) :: delta_e
    integer(int32) :: u, r, d, l, c, nc
    do c = 0, mstate - 1
       magne_table(c) = cos(c * pi_state_inv)
    end do
    do c = 0, mstate - 1
       do l = 0, mstate - 1
          do d = 0, mstate - 1
             energy_table_ru(l, d, c) = calc_local_energy_local(l, l, d, d, c) / 2
             do r = 0, mstate - 1
                do u = 0, mstate - 1
                   energy_table(u, r, d, l, c) = calc_local_energy_local(u, r, d, l, c)
                end do
             end do
          end do
       end do
    end do
    ! write(error_unit, *) energy_table
    do nc = 0, mstate - 1
       do c = 0, mstate - 1
          do l = 0, mstate - 1
             do d = 0, mstate - 1
                do r = 0, mstate - 1
                   do u = 0, mstate - 1
                      delta_e = energy_table(u, r, d, l, nc) - energy_table(u, r, d, l, c)
                      energy_diff_table(u, r, d, l, c, nc) = delta_e
                      if (delta_e <= 0.0d0) then
                         probability_table(u, r, d, l, c, nc) = 2d0
                      else
                         probability_table(u, r, d, l, c, nc) = exp(- beta * delta_e)
                      end if
                   end do
                end do
             end do
          end do
       end do
    end do
  contains
    pure real(real64) function calc_local_energy_local(u, r, d, l, c) result(res)
      integer(int32), intent(in) :: u, r, d, l, c
      res = - sum(cos((c - [u, r, d, l]) * pi_state_inv))
    end function calc_local_energy_local
  end subroutine init_tables
  !> init_sixclock_order: Set spins with the all-alinged state.
  impure subroutine init_sixclock_order()
    sixclock_even(:, :) = 0_int32
    sixclock_odd(:, :) = 0_int32
  end subroutine init_sixclock_order
  !> update_metropolis: Update the lattice with 1MCS.
  impure subroutine update_metropolis()
    integer(int64) :: x, y
    integer(int64) :: uy, dy
    do y = 1, ny
       do x = 1, nx
          rnds(1, x, y) = grnd()
          rnds(2, x, y) = grnd()
       end do
    end do
    !> even.
    do y = 1, ny
       uy = y + 1
       if (uy > ny) uy = 1_int64
       dy = y - 1
       if (dy < 1) dy = ny
       do x = 1, nx / 2
          call local_flip_even(x, y, uy, dy)
       end do
    end do
    !> odd.
    do y = 1, ny
       uy = y + 1
       if (uy > ny) uy = 1_int64
       dy = y - 1
       if (dy < 1) dy = ny
       do x = 1, nx / 2
          call local_flip_odd(x, y, uy, dy)
       end do
    end do
  end subroutine update_metropolis
  !> local_flip_even: Flip a spin of (x, y) position of lattice in the even checkerboard.
  !> @param x A x-position of the spin.
  !> @param y A y-position of the spin.
  impure subroutine local_flip_even(x, y, uy, dy)
    integer(int64), intent(in) :: x, y, uy, dy
    integer(int32) :: candidate_states
    integer(int64) :: lx, rx
    real(real64) :: prob
    rx = x + dx(right, iand(y, b'1'))
    if (rx > nx / 2) rx = 1_int64
    lx = x + dx(left, iand(y, b'1'))
    if (lx < 1) lx = nx / 2
    candidate_states = sixclock_even(x, y) + (1 + floor(mstate * rnds(1, 2 * x - 1 + iand(y + 1, b'1'), y)))
    if (candidate_states >= mstate) candidate_states = candidate_states - mstate
    prob = probability_table(&
         & sixclock_odd(x, uy), sixclock_odd(x, dy), sixclock_odd(lx, y), sixclock_odd(rx, y), &
         & sixclock_even(x, y), &
         & candidate_states)
    if (prob >= 1d0) then
       sixclock_even(x, y) = candidate_states
       return
    end if
    if (rnds(2, 2 * x - 1 + iand(y + 1, b'1'), y) < prob) then
       sixclock_even(x, y) = candidate_states
    end if
  end subroutine local_flip_even
  !> local_flip_odd: Flip a spin of (x, y) position of lattice in the odd checkerboard.
  !> @param x A x-position of the spin.
  !> @param y A y-position of the spin.
  impure subroutine local_flip_odd(x, y, uy, dy)
    integer(int64), intent(in) :: x, y, uy, dy
    integer(int32) :: candidate_states
    integer(int64) :: lx, rx
    real(real64) :: prob
    rx = x + dx(right, iand(y + 1, b'1'))
    if (rx > nx / 2) rx = 1_int64
    lx = x + dx(left, iand(y + 1, b'1'))
    if (lx < 1) lx = nx / 2
    candidate_states = sixclock_odd(x, y) + (1 + floor(mstate * rnds(1, 2 * x - 1 + iand(y, b'1'), y)))
    if (candidate_states >= mstate) candidate_states = candidate_states - mstate
    prob = probability_table(&
         & sixclock_even(x, uy), sixclock_even(x, dy), sixclock_even(lx, y), sixclock_even(rx, y), &
         & sixclock_odd(x, y), &
         & candidate_states)
    if (prob >= 1d0) then
       sixclock_odd(x, y) = candidate_states
       return
    end if
    if (rnds(2, 2 * x - 1 + iand(y, b'1'), y) < prob) then
       sixclock_odd(x, y) = candidate_states
    end if
  end subroutine local_flip_odd
  !> calc_energy: Calculate the energy density.
  pure real(real64) function calc_energy() result(res)
    integer(int64) :: x, y
    integer(int64) :: rx, uy
    res = 0d0
    do y = 1, ny
       uy = y + 1
       if (uy > ny) uy = 1_int64
       do x = 1, nx / 2
          rx = x + dx(right, iand(y, b'1'))
          if (rx > nx / 2) rx = 1_int64
          res = res + energy_table_ru(sixclock_odd(x, uy), sixclock_odd(rx, y), sixclock_even(x, y))
       end do
       do x = 1, nx / 2
          rx = x + dx(right, iand(y + 1, b'1'))
          if (rx > nx / 2) rx = 1_int64
          res = res + energy_table_ru(sixclock_even(x, uy), sixclock_even(rx, y), sixclock_odd(x, y))
       end do
    end do
    res = res * nall_inv
  end function calc_energy
  !> calc_magne: Calculate the magnetism density.
  pure real(real64) function calc_magne() result(res)
    integer(int64) :: x, y
    res = 0d0
    do y = 1, ny
       do x = 1, nx / 2
          res = res + magne_table(sixclock_even(x, y))
          res = res + magne_table(sixclock_odd(x, y))
       end do
    end do
    res = res * nall_inv
  end function calc_magne
end module sixclock_periodic_dual_lattice_table_near_m