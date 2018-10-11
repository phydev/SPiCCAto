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

    type(cell) :: phi, phi_old
    type(mesh) :: sub

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

    call system('mkdir '//trim(dir_name))

    write(*,'(A,A)') " Simulation ID : ", sim_id
    write(*,'(A,I10)') " Random seed   : ", iseed
    write(*,'(A,F10.2)') " Density       : ", density
    write(*,'(A,F10.2)') " Adhesion      : ", eta
    write(*,'(A,F10.2)') " Chemotaxis    : ", chi
    write(*,'(A,F10.2)') " Repulsion     : ", gamma
 

    cm = (/M_TEN, M_TEN, M_TEN/)

    volume_target = (M_FOUR/M_THREE)*M_PI*radius**3
    call phi%initialize(20,20,20,'Neumann',cm,radius)
    call sub%initialize(200,30,30,'periodic')

    tstep = 200000
    output_period = 1000
    output_counter = 0

    call phi_old%copy(phi)
    call gen_cell_points(2.0,sphere,np_sphere)

    print*, "Initializing substrate . . ."
    call substrate_init(sub, density, sphere, np_sphere, iseed, cm, radius)
    print*, "Smoothing interfaces . . ."
    call sub%smoothing
    call phi%smoothing

    file_name = 's'
    call sub%output(file_name)

    nstep = 0
    print*, "Starting simulation . . ."
    do while (nstep.le.tstep)

        call CPU_TIME(time_init)

        call phi%com
      

!        write(*,'(F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,I3)'), phi%gcom, phi%lcom, nstep

        do ip=0, phi%nodes

            call vec_local2global(s, int( anint(phi%gcom)), phi%position(ip))
            !print*, s
            ip_global = sub%ip(s)
            gradient = phi%gradient(ip)
            phi%grid(ip) = phi_old%gt(ip) + dt*(-chi*gradient(1) + phi_old%laplacian(ip) + &
                epsilon*phi_old%gt(ip)*(1.0-phi_old%gt(ip))*(phi_old%gt(ip) -0.5 &
             + alpha_v*(volume_target-phi%volume) - gamma*h(sub%gt(ip_global)) ) )  

        end do

        call phi%com 
      
        dr = anint(phi%lcom-(/ phi%L(1)/2, phi%L(1)/2, phi%L(1)/2 /))

         if(abs(dr(1)).ge.1.or.abs(dr(2)).ge.1.or.abs(dr(3)).ge.1. ) then
            phi%gcom = phi%gcom + dr
         !   print*, dr
            do ip=0, phi_old%nodes
               s =  phi%position(ip) + dr(1:3)
               phi_old%grid(ip) = phi%gt(phi%ip(s) )
           end do
           call phi%copy(phi_old)

         else
           call phi_old%copy(phi)

         end if

        if(output_counter.ge.output_period) then
            output_counter = 0
            write(*,*) nstep
            call format_this(nstep,format_string)   
            write(file_name,format_string) nstep
            file_name = sim_id//"/phi"//file_name
            call sub%output(file_name,phi) !! output the cell field phi inside the simulation box of sub.
           !!call phi%output(file_name)
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
    end do

end program main



