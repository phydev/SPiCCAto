program main

  use run_cells_m
  use global_m

  implicit none

  call get_command_argument(1,arg_simid) ! simulation id
  call get_command_argument(2,arg_iseed) ! seed for the random generator
  call get_command_argument(3,arg_density) ! density of fibres
  call get_command_argument(4,arg_eta) ! adhesion coefficient
  call get_command_argument(5,arg_chi) ! chemotactic response
  call get_command_argument(6,arg_gamma) ! depletion force between the cell and fibres

  sim_id = trim(arg_simid)
  read (arg_iseed,*) iseed
  read (arg_density,*) density
  read (arg_eta,*) eta
  read (arg_chi,*) chi
  read (arg_gamma,*) gamma
  dir_name = trim(sim_id)

  tstep = 10000
  ntype = 1
  adh1 = 0.0
  adh2 = 0.0
  cell_radius = 5.0
  interface_width = 1.0
  dt = 0.005
  vol_lagrangian = 0.005/dt
  Lsize = (/ 30, 30, 30 /)
  dr = (/ 1, 1, 1 /)
  np_bndry = 100
  depletion_weight = 1.0

  periodic = .true.
  output_period = 1000
  write(*,*) sim_id, iseed, density, eta, chi, gamma, dir_name

  call system('mkdir '//trim(dir_name))

  call run_cells()



end program main
