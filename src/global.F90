!! Copyright (C) 2018 M. Moreira
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 3, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!!

module global

  implicit none

  private

  !> some mathematical constants
  REAL(8), public, parameter :: M_Pi          = 3.1415926535897932384626433832795029
  REAL(8), public, parameter :: M_E           = 2.7182818284590452353602874713526625
  REAL(8), public, parameter :: M_TWO         = 2.0
  REAL(8), public, parameter :: M_THREE       = 3.0
  REAL(8), public, parameter :: M_FOUR        = 4.0
  REAL(8), public, parameter :: M_FIVE        = 5.0
  REAL(8), public, parameter :: M_SIX         = 6.0
  REAL(8), public, parameter :: M_SEVEN       = 7.0
  REAL(8), public, parameter :: M_EIGHT       = 8.0
  REAL(8), public, parameter :: M_NINE        = 9.0
  REAL(8), public, parameter :: M_TEN         = 10.0
  REAL(8), public, parameter :: M_HALF        = 0.5

  !> fixed model parameters
  REAL, public :: dt                       = 0.005
  REAL, public :: radius                   = 6.0
  REAL, public :: alpha_v                  = 1.0
  REAL, public :: epsilon                  = 1.0

  !> some strings
  CHARACTER(len=:), allocatable, public :: string, dir_name, sim_id, format_string
  CHARACTER(len=10), public :: arg_iseed, arg_density, arg_chi, arg_gamma, arg_eta, arg_simid, file_name, file_id

  INTEGER, allocatable, public :: sphere(:,:)
  INTEGER, public :: ip, s(3), nstep, tstep, iseed, output_counter, output_period, np_sphere
  REAL, public :: cm(3), cm_old(3), dr(3), volume, density, eta, gamma, chi, gradient(3), time_init, time_end, ctime, volume_target

  

end module global
