module global_m

  implicit none

  private

  type mesh_t
  
    real :: phi, mu, lapl_phi, lapl_mu, lapl_h
    integer :: itype
  end type mesh_t

  public :: mesh_t

  !> some mathematical constants
  REAL(8), public, parameter :: M_Pi        = 3.1415926535897932384626433832795029
  REAL(8), public, parameter :: M_E         = 2.7182818284590452353602874713526625
  REAL(8), public, parameter :: M_TWO       = 2.0
  REAL(8), public, parameter :: M_THREE     = 3.0
  REAL(8), public, parameter :: M_FOUR      = 4.0
  REAL(8), public, parameter :: M_FIVE      = 5.0
  REAL(8), public, parameter :: M_SIX       = 6.0
  REAL(8), public, parameter :: M_SEVEN     = 7.0
  REAL(8), public, parameter :: M_EIGHT     = 8.0
  REAL(8), public, parameter :: M_NINE      = 9.0
  REAL(8), public, parameter :: M_TEN       = 10.0
  REAL(8), public, parameter :: M_HALF      = 0.5

  !  begin parameters
  real, public :: cell_radius, interface_width,  dt, eta, chi, gamma, density
  integer, public :: tstep, iseed, temp, ip_part2, ip_part
  character(len=10), public :: buffer
  character(len=:), allocatable, public :: string
  character(len=:), allocatable, public :: dir_name
  character(len=:), allocatable, public :: sim_id
  character(len=10), public :: arg_iseed
  character(len=10), public :: arg_density
  character(len=10), public :: arg_chi
  character(len=10), public :: arg_gamma
  character(len=10), public :: arg_eta
  character(len=10), public :: arg_simid

  character(len=10), public :: file_name, file_id
  character(len=:), allocatable, public :: format_string
  real, public ::  depletion_weight, phi_total
  ! mesh variables
  integer, allocatable, public :: lxyz(:,:), lxyz_inv(:,:,:),  grid_cell_domain(:),ncell(:)
  integer, public :: Lsize(1:3), np_bndry
  integer, public :: np, np_tt, ip
  real, allocatable, public :: gg(:,:), r_cm(:,:), r_cm_part(:,:), r_cmg(:,:)
  type(mesh_t), allocatable, public :: cell(:,:), aux(:,:)
  ! chemical and substrate variables
  real, allocatable, public :: chem(:), gchem(:,:), s(:), shfield(:), shfield_lapl(:)
  real, public :: chemresponse, metcoef, scoef, adhs
  ! velocity
  real, public, allocatable :: velocity(:), path(:)
  real, public :: rt(3), vt(3)
  integer, public :: vcounter
  ! new
  real, public :: self_int, sum_int, wgamma(10)
  integer, public :: jcell, icell, tcell, i, j, k
  integer, allocatable, public :: r(:)
  real, allocatable, public :: f(:), mu(:), lapl(:)
  ! misc
  integer, public :: nstep, counter = 0, ndim, output_period, extra_steps
  real, public :: hs, time_init, time_end, ctime, phi_max, hphi, fnu, porosity
  logical, public :: periodic

  integer, public :: np_sphere, sphere(1000,1:3),  porous(1000,1:3), np_porous, nParticle
  ! local variables
  integer,allocatable, public :: lxyz_part(:,:), lxyz_inv_part(:,:,:), gammaw(:)
  integer, public :: Lsize_part(3), np_part, np_part_tt, ntype, nleap, dri(3), drf(3), dr(3),&
  itype, np_part_bndry, ip2, ip_global

  ! coupling constants between the cell_radius
  real, allocatable, public :: volume(:)
  real, allocatable, public :: adhesion(:,:), hfield(:,:), hfield_lapl(:,:)
  real, public :: adh1, adh2, v1, v2, volume_target, vol_lagrangian
  integer, public :: cm_calc_counter


end module global_m
