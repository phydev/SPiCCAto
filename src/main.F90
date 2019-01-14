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

program main
 
    use class_mesh
    use global
    use routines 

    implicit none

    type(cell), allocatable :: phi(:)
    type(cell) :: phi_old
    type(mesh) :: sub, aux

    call get_command_argument(1,arg_simid) ! simulation id
    call get_command_argument(2,arg_iseed) ! seed for the random generator
    call get_command_argument(3,arg_density) ! density of fibres
    call get_command_argument(4,arg_eta) ! adhesion coefficient 
    call get_command_argument(5,arg_chi) ! chemotactic response
    call get_command_argument(6,arg_gamma) ! depletion force between the cell and fibres
    call get_command_argument(7,arg_number_of_cells) ! number of cells


    sim_id = trim(arg_simid)
    read (arg_iseed,*) iseed
    read (arg_density,*) density
    read (arg_eta,*) eta
    read (arg_chi,*) chi
    read (arg_gamma,*) gamma
    read (arg_number_of_cells,*) nc 
    dir_name = trim(sim_id)


    ! ic - index used to identify a given cell ic = 1 ... nc

    ALLOCATE(box_position(1:nc,3))
    ALLOCATE(phi(nc))

    call system('mkdir '//trim(dir_name))

    ! Variables that may be changed by the user
    L = (/ 30, 30, 30 /) ! size of the simulation box
    box_length = (/30, 30, 30 /) ! lengths of the small box around the cell - if high deformation is required, please increase the size
                                 ! for small deformations 20x20x20 works just fine!
    tsteps = 50000  ! total number of iterations
    output_period = 500 ! the program will give outputs at every output_period iterations
    output_counter = 0  ! must be set to zero - do not change
    collision = .True. ! if you want to run the test simulation where a cell collide with another, set it to .True. and use nc = 2
    volume_target = (M_FOUR/M_THREE)*M_PI*radius**3 ! target volume of each cell used in the Lagrange multiplier
    
    call print_header(sim_id, iseed, density, eta, chi, gamma, nc, L, tsteps)
 

    !! example : how to input a grid from file:
    !! call phi(1)%initialize(box_length(1),box_length(2),box_length(3),'Neumann',cm,radius) ! creates the grid 
    !! file_name = 'ad3/phii'  ! do not write the file extension here, just the folder and file name
    !! call phi(1)%input(file_name) 
    !! you must write a file with the box_position too and use that as input file

    if(collision) then
        box_position(1,:) = (/15,14,14 /) ! position of the first cell
        box_position(2,:) = (/25,7,7 /) ! position of the second cell
        phi(1)%gcom = box_position(1,:)
        phi(2)%gcom = box_position(2,:)
    else
        ic = 1
        do i=0,L(1)-10,10
            do j=0,L(2)-10,10
                box_position(ic,:) = (/i, j, 10/) 
                phi(ic)%gcom = real(box_position(ic,:))
                ic = ic + 1
                if(ic.gt.nc) EXIT
            end do
            if(ic.gt.nc) EXIT
        end do
    end if


    do ic=1,nc
        cm =  box_position(ic,:)
        call phi(ic)%initialize(box_length(1),box_length(2),box_length(3),'Neumann',cm,radius)
    end do

    call sub%initialize(L(1),L(2),L(3),'periodic') ! substrate initialization 
    call aux%initialize(L(1),L(2),L(3),'periodic') ! auxiliar field that contains all cells of the same type

    call gen_cell_points(2.0,sphere,np_sphere)

    print*, "Initializing substrate . . ."
    !call substrate_init(sub, density, sphere, np_sphere, iseed, cm, radius)
    print*, "Smoothing interfaces . . ."
    !call sub%smoothing

    call phi(1)%smoothing ! integrating the cell 1 in order to smooth the interfaces
    do ic=2,nc
        call phi(ic)%copy(phi(1)) ! copying the integrated field to the other cells
    end do
    file_name = sim_id//'/auxi'
    call aux%output(file_name)
    file_name = sim_id//'/phii'
    call sub%output(file_name,phi(1))
    

    nstep = 0
    nprint = 0
    print*, "Starting simulation . . ."
    do while (nstep.le.tsteps)

        call CPU_TIME(time_init)

        call mesh_calculate_auxiliar_field(aux,phi) ! calculates auxiliar field

        do ic=1, nc  ! running the loop over all cells

            call phi(ic)%com ! calculates center of mass of each cell

            call phi_old%copy(phi(ic)) ! copy phi(ic) to phi_old

            do ip=0, phi(ic)%nodes-1 ! loop over each element of the cells' box

                s_local = box_position(ic,:) - phi(ic)%L/2 ! subtracting half the size of the box, in order to obtain the position of the box corner (0,0,0)

                call vec_local2global(s, int(anint(s_local)), phi(ic)%position(ip),L) ! calculates the global vector relative to the cell position
              
                ip_global = aux%ip(s) ! calculates the ip global from the global vector
                
                gradient = 0.0 !phi(nc)%gradient(ip)
                if(ic.eq.1) gradient = phi(1)%gradient(ip) ! giving a velocity only to the first cell
                !if(ic.eq.2) gradient = -phi(2)%gradient(ip)

                phi(ic)%grid(ip) = phi_old%gt(ip) + dt*(-chi*gradient(1)   + phi_old%laplacian(ip) + &
                    epsilon*phi_old%gt(ip)*(1.0-phi_old%gt(ip))*(phi_old%gt(ip) -0.5 &
                 + alpha_v*(volume_target-phi(ic)%volume) - gamma*h(sub%gt(ip_global)) - gamma*( h(aux%gt(ip_global) ) &
                    - gamma*h(phi_old%gt(ip)) )  + eta*aux%laplacian(ip_global)   ) )   ! integrating

            end do
        
            call phi(ic)%com  ! calculating the new center of mass
      
            dr = int(phi(ic)%lcom-(/ phi(ic)%L(1)/2, phi(ic)%L(2)/2, phi(ic)%L(3)/2 /) )
     
            phi(ic)%gcom = phi(ic)%gcom + anint(dr) ! integrating the center of mass
            box_position(ic,:) = box_position(ic,:) + anint(dr) 

            ! verifying boundary conditions
            do i=1,3
             k = box_position(ic,i)
             j = phi(ic)%gcom(i) 
             call check_boundary(k,L(i),sub%b)
             call check_boundary(j,L(i),sub%b)
             box_position(ic,i) = k
             phi(ic)%gcom(i) = j 
            end do
            ! re-centering the cell inside the small box
            do ip=0, phi_old%nodes-1
               s =  phi(ic)%position(ip) + dr(1:3)
               phi_old%grid(ip) = phi(ic)%gt(phi(ic)%ip(s) ) ! the re-centered cell is being stored in phi_old
            end do
            call phi(ic)%copy(phi_old) ! copying the re-centered cell from phi_old to phi

        end do

        if(output_counter.ge.output_period) then
            output_counter = 0
            write(*,*) nstep
            call format_this(nstep,format_string)   ! formating the string 
            write(file_name,format_string) nstep ! transforming an integer to a string
            file_name = sim_id//"/phi1_"//file_name
            call aux%output(file_name,phi(1)) !! output the cell field phi inside the simulation box aux   
            write(file_name,format_string) nstep
            file_name = sim_id//"/phi2_"//file_name
            call aux%output(file_name,phi(2)) !! output the cell field phi inside the simulation box of  aux 
            !write(file_name,format_string) nstep
            !file_name = sim_id//"/sub"//file_name  !! uncomment if you're using a substrate
            !call sub%output(file_name)
            !call phi%output(file_name)
        end if


        call CPU_TIME(time_end)

        ctime = ctime + (time_end - time_init)

        if(nstep.eq.500) then
          ctime = (ctime*(tsteps-nstep) )/6000.d0

          if( ctime>60.d0) then
            write(*,'(A,F10.2)') " Estimated time (hour): ",ctime/60.d0
          else
            write(*,'(A,F10.2)') " Estimated time (min): ",ctime
          end if
        end if
   

        output_counter = output_counter + 1
        nstep = nstep + 1
        nprint = nprint + 1
    end do


end program main




