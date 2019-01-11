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
    !call get_command_argument(7,arg_number_of_cells) ! number of cells

    sim_id = trim(arg_simid)
    read (arg_iseed,*) iseed
    read (arg_density,*) density
    read (arg_eta,*) eta
    read (arg_chi,*) chi
    read (arg_gamma,*) gamma
    !read (arg_number_of_cells,*) arg_number_of_cells
    dir_name = trim(sim_id)

    ! nc - number of cells
    ! ic - index used to identify a given cell ic = 1 ... nc
    nc = 9
    ALLOCATE(box_position(1:nc,3))
    ALLOCATE(phi(nc))

    call system('mkdir '//trim(dir_name))

    write(*,'(A,A)') " Simulation ID : ", sim_id
    write(*,'(A,I10)') " Random seed   : ", iseed
    write(*,'(A,F10.2)') " Density       : ", density
    write(*,'(A,F10.2)') " Adhesion      : ", eta
    write(*,'(A,F10.2)') " Chemotaxis    : ", chi
    write(*,'(A,F10.2)') " Repulsion     : ", gamma
 
    box_length = (/20, 20, 20 /)
    L = (/ 30, 30, 20 /)
    
    do ic=1,nc
        call phi(ic)%initialize(box_length(1),box_length(2),box_length(3),'Neumann',cm,radius)
    end do

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

    cm =  box_position(1,:) !(/15, 10, 10 /) 
    volume_target = (M_FOUR/M_THREE)*M_PI*radius**3


    call sub%initialize(L(1),L(2),L(3),'periodic')
    call aux%initialize(L(1),L(2),L(3),'periodic')
    ! do ic=1, nc
    !     call format_this(ic,format_string)   
    !     write(file_name,format_string) ic
    !     file_name = sim_id//"/phi"//file_name
    !     print*, phi(ic)%gcom
    !     call sub%output(file_name,phi(ic)) 
    ! end do
 
    call mesh_calculate_auxiliar_field(aux,phi)
    file_name = 'auxiliar'
    call aux%output(file_name)

    tstep = 50000
    output_period = 500
    output_counter = 0

    call gen_cell_points(2.0,sphere,np_sphere)

    print*, "Initializing substrate . . ."
    !call substrate_init(sub, density, sphere, np_sphere, iseed, cm, radius)
    print*, "Smoothing interfaces . . ."
    !call sub%smoothing
    call phi(1)%smoothing
    do ic=2,nc
        call phi(ic)%copy(phi(1))
    end do
    file_name = sim_id//'/auxi'
    call aux%output(file_name)
    file_name = sim_id//'/phii'
    call sub%output(file_name,phi(1))
    

    nstep = 0
    nprint = 0
    print*, "Starting simulation . . ."
    do while (nstep.le.tstep)

        call CPU_TIME(time_init)

        call mesh_calculate_auxiliar_field(aux,phi)

        do ic=1, nc 

            call phi(ic)%com

            call phi_old%copy(phi(ic))

            do ip=0, phi(ic)%nodes-1

                s_local = box_position(ic,:) - phi(ic)%L/2

                call vec_local2global(s, int(anint(s_local)), phi(ic)%position(ip),L)
              
                ip_global = aux%ip(s)
                
                gradient = 0.0 !phi(nc)%gradient(ip)
                if(ic.eq.1) gradient = phi(1)%gradient(ip)
                phi(ic)%grid(ip) = phi_old%gt(ip) + dt*(-chi*sum(gradient(1:2))   + phi_old%laplacian(ip) + &
                    epsilon*phi_old%gt(ip)*(1.0-phi_old%gt(ip))*(phi_old%gt(ip) -0.5 &
                 + alpha_v*(volume_target-phi(ic)%volume) - gamma*h(sub%gt(ip_global)) - gamma*( h(aux%gt(ip_global) ) &
                    - gamma*h(phi_old%gt(ip)) ) ) )  

            end do
        
            call phi(ic)%com 
      
            dr = int(phi(ic)%lcom-(/ phi(ic)%L(1)/2, phi(ic)%L(2)/2, phi(ic)%L(3)/2 /) )
            !dr = (/ img(dr(1),L(1)),img(dr(2),L(2)),img(dr(3),L(3)) /) ! I think I can just comment this. . .
     
            phi(ic)%gcom = phi(ic)%gcom + anint(dr)
            box_position(ic,:) = box_position(ic,:) + anint(dr) 


            do i=1,3
             k = box_position(ic,i)
             j = phi(ic)%gcom(i) 
             call check_boundary(k,L(i),sub%b)
             call check_boundary(j,L(i),sub%b)
             box_position(ic,i) = k
             phi(ic)%gcom(i) = j 
            end do

            do ip=0, phi_old%nodes-1
               s =  phi(ic)%position(ip) + dr(1:3)
               phi_old%grid(ip) = phi(ic)%gt(phi(ic)%ip(s) )
            end do
            call phi(ic)%copy(phi_old)

        end do

        if(output_counter.ge.output_period) then
            output_counter = 0
            write(*,*) nstep
            call format_this(nstep,format_string)   
            write(file_name,format_string) nstep
            file_name = sim_id//"/phi"//file_name
            call aux%output(file_name,phi(1)) !! output the cell field phi inside the simulation box of sub.   
            write(file_name,format_string) nstep
            file_name = sim_id//"/aux"//file_name
            call aux%output(file_name)
            !call phi%output(file_name)
        end if


        call CPU_TIME(time_end)

        ctime = ctime + (time_end - time_init)

        if(nstep.eq.500) then
          ctime = (ctime*(tstep-nstep) )/6000.d0

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




