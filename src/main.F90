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
 
    box_length = (/30, 30, 30 /)
    L = (/ 60, 60, 60 /)
    
    box_position = (/0, 30, 30 /)

    cm =  box_position !(/15, 10, 10 /) 
    volume_target = (M_FOUR/M_THREE)*M_PI*radius**3
    call phi%initialize(box_length(1),box_length(2),box_length(3),'Neumann',cm,radius)
    call sub%initialize(L(1),L(2),L(3),'periodic')
    phi%gcom = cm
    tstep = 20000
    output_period = 1000
    output_counter = 0

    call phi_old%copy(phi)
    call gen_cell_points(2.0,sphere,np_sphere)

    print*, "Initializing substrate . . ."
    call substrate_init(sub, density, sphere, np_sphere, iseed, cm, radius)
    print*, "Smoothing interfaces . . ."
    call sub%smoothing
    call phi%smoothing

    file_name = sim_id//'/s'
    call sub%output(file_name)
    file_name = sim_id//'/phii'
    call sub%output(file_name,phi)
    file_name = sim_id//'/phil'
    call phi%output(file_name)



    ! the position of the CoM along the time will be registered in the following file
    open(unit=1001001, file=sim_id//"/v.out" )

    nstep = 0
    nprint = 0
    print*, "Starting simulation . . ."
    do while (nstep.le.tstep)

        call CPU_TIME(time_init)

        call phi%com
      


        do ip=0, phi%nodes-1

            call vec_local2global(s, int( anint(box_position - phi%L)), phi%position(ip),L)
            !print*, s
            ip_global = sub%ip(s)
            gradient = phi%gradient(ip)
            phi%grid(ip) = phi_old%gt(ip) + dt*(-chi*gradient(1) + phi_old%laplacian(ip) + &
                epsilon*phi_old%gt(ip)*(1.0-phi_old%gt(ip))*(phi_old%gt(ip) -0.5 &
             + alpha_v*(volume_target-phi%volume) - gamma*h(sub%gt(ip_global)) ) )  

        end do

        call phi%com 
      
        dr = int(phi%lcom-(/ phi%L(1)/2, phi%L(2)/2, phi%L(3)/2 /) )
        dr = (/ img(dr(1),L(1)),img(dr(2),L(2)),img(dr(3),L(3)) /)
   !     if(abs(dr(1)).ge.1.or.abs(dr(2)).ge.1.or.abs(dr(3)).ge.1. ) then
             phi%gcom = phi%gcom + anint(dr)
             box_position = box_position + anint(dr) 


             do i=1,3
              k = box_position(i)
              j = phi%gcom(i) 
              call check_boundary(k,L(i),sub%b)
              call check_boundary(j,L(i),sub%b)
              box_position(i) = k
              phi%gcom(i) = j 
             end do

            do ip=0, phi_old%nodes-1
               s =  phi%position(ip) + dr(1:3)
               phi_old%grid(ip) = phi%gt(phi%ip(s) )
            end do
            call phi%copy(phi_old)

  !        else

!           call phi_old%copy(phi)

!        end if



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
   
        if(nprint.eq.200) then !! the interval of 200 iterations was careful estimated - do not change!
            nprint = 0
            write(1001001,'(F10.2,F10.2,F10.2,F10.2)')  nstep*dt, phi%gcom(1), phi%gcom(2), phi%gcom(3)
        end if

        output_counter = output_counter + 1
        nstep = nstep + 1
        nprint = nprint + 1

        if(phi%gcom(1).ge.190) EXIT !! if the cell reaches the end of the box, the simulation is terminated

    end do

    call integral_path(sim_id//"v.out")

    close(1001001)
end program main




