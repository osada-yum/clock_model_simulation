module sixclock_periodic_table_m
  use, intrinsic :: iso_fortran_env
  use msmt19937
  implicit none
  private
  character(len=*), parameter :: version = "table"
  public :: print_version

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

  integer(int32), allocatable :: sixclock(:, :)

  real(real64), allocatable :: rnds(:, :, :)

  integer(int64), parameter :: nd = 4
  integer(int64), parameter :: dy(nd) = [0_int64, 1_int64, 0_int64, -1_int64]
  integer(int64), parameter :: dx(nd) = [1_int64, 0_int64, -1_int64, 0_int64]

  public :: mstate, nx, ny, nall, kbt, beta
  public :: init_sixclock, init_sixclock_order, update_metropolis, calc_energy, calc_magne

contains
  impure subroutine print_version()
    write(output_unit, '(a)') "#"//version
    write(error_unit, '(a)') "#"//version
  end subroutine print_version
  !> init_sixclock: Initialize 6-state clock model.
  impure subroutine init_sixclock()
    allocate(sixclock(nx, ny))
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
    sixclock(:, :) = 0_int32
  end subroutine init_sixclock_order
  !> update_metropolis: Update the lattice with 1MCS.
  impure subroutine update_metropolis()
    integer(int64) :: x, y
    do y = 1, ny
       do x = 1, nx
          rnds(1, x, y) = grnd()
          rnds(2, x, y) = grnd()
       end do
    end do
    !> even.
    do y = 1, ny
       do x = 1 + iand(y - 1, b'1'), nx, 2
          call local_flip(x, y)
       end do
    end do
    !> odd.
    do y = 1, ny
       do x = 1 + iand(y, b'1'), nx, 2
          call local_flip(x, y)
       end do
    end do
  end subroutine update_metropolis
  !> local_flip: Flip a spin of (x, y) position of lattice.
  !> @param x A x-position of the spin.
  !> @param y A y-position of the spin.
  impure subroutine local_flip(x, y)
    integer(int64), intent(in) :: x, y
    integer(int32) :: candidate_states
    integer(int32) :: nearest_states(nd)
    integer(int64) :: near_x, near_y
    real(real64) :: prob
    integer(int32) :: d
    do d = 1, nd
       near_x = x + dx(d)
       if (near_x > nx) then
          near_x = 1_int64
       else if (near_x < 1_int64) then
          near_x = nx
       end if
       near_y = y + dy(d)
       if (near_y > ny) then
          near_y = 1_int64
       else if (near_y < 1_int64) then
          near_y = ny
       end if
       nearest_states(d) = sixclock(near_x, near_y)
    end do
    candidate_states = sixclock(x, y) + (1 + floor((mstate - 1) * rnds(1, x, y)))
    if (candidate_states >= mstate) candidate_states = candidate_states - mstate
    prob = probability_table(&
         & nearest_states(1), &
         & nearest_states(2), &
         & nearest_states(3), &
         & nearest_states(4), &
         & sixclock(x, y), &
         & candidate_states)
    if (prob >= 1d0) then
       sixclock(x, y) = candidate_states
       return
    end if
    if (rnds(2, x, y) < prob) then
       sixclock(x, y) = candidate_states
    end if
  end subroutine local_flip
  !> calc_energy: Calculate the energy density.
  pure real(real64) function calc_energy() result(res)
    integer(int64) :: x, y
    integer(int64) :: rx, uy
    res = 0d0
    do y = 1, ny
       uy = y + 1
       if (uy > ny) uy = 1_int64
       do x = 1, nx
          rx = x + 1
          if (rx > nx) rx = 1_int64
          res = res + energy_table_ru(sixclock(x, uy), sixclock(rx, y), sixclock(x, y))
       end do
    end do
    res = res * nall_inv
  end function calc_energy
  !> calc_magne: Calculate the magnetism density.
  pure real(real64) function calc_magne() result(res)
    integer(int64) :: x, y
    res = 0d0
    do y = 1, ny
       do x = 1, nx
          res = res + magne_table(sixclock(x, y))
       end do
    end do
    res = res * nall_inv
  end function calc_magne
end module sixclock_periodic_table_m
