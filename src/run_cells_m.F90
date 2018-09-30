!! Copyright (C) 2016 M. Moreira
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
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

module run_cells_m

  use global_m
  use derivatives_m
  use sim_init_m
  use misc_m

  implicit none

  private

  public :: run_cells


  contains

    subroutine run_cells()

      implicit none

!      real, intent(in) :: porosity
 !     integer, intent(in) :: iseed
 !     character(3), intent(in) :: sim_id
     

      

      dir_name = sim_id

      ! initializing parameters
      call  write_parameters(cell_radius, density, eta, gamma, chi, interface_width, tstep, dt, Lsize,  dir_name, iseed,  periodic)

      ! number of points in the mesh
      np = 8*Lsize(1)*Lsize(2)*Lsize(3) ! number of points
      np_tt = 8*(Lsize(1)+2*np_bndry)*(Lsize(2)+2*np_bndry)*(Lsize(3)+2*np_bndry) ! number of points plus boundary points


      ntype = 1
      ALLOCATE(ncell(1:ntype))     
      ncell(1) = 1 ! (np*density(1:ntype))/(4*cell_radius**2)
      tcell = 1!sum(ncell(1:ntype))

      ! partition mesh
      Lsize_part(1:3) = (/20, 20, 20/)
      np_part = 8*Lsize_part(1)*Lsize_part(2)*Lsize_part(3)
      np_part_bndry = 2*Lsize(1)
      np_part_tt = 8*(Lsize_part(1) + 2*np_part_bndry)*(Lsize_part(2) + 2*np_part_bndry)*(Lsize_part(3) + 2*np_part_bndry)
      ! allocating matrices and vectors

      ALLOCATE(lxyz(np_tt,1:3))
      ALLOCATE(lxyz_inv(-Lsize(1)-np_bndry:Lsize(1)+np_bndry, &
                        -Lsize(2)-np_bndry:Lsize(2)+np_bndry,&
                        -Lsize(3)-np_bndry:Lsize(3)+np_bndry))
      ALLOCATE(lxyz_part(np_part_tt,1:3))
      ALLOCATE(lxyz_inv_part(-Lsize_part(1)-np_part_bndry:Lsize_part(1)+np_part_bndry, &
                        -Lsize_part(2)-np_part_bndry:Lsize_part(2)+np_part_bndry, &
                        -Lsize_part(3)-np_part_bndry:Lsize_part(3)+np_part_bndry))

      ALLOCATE(cell(0:np_part,tcell))
      ALLOCATE(aux(np,ntype))
      ALLOCATE(adhesion(0:np_part,tcell))
      ALLOCATE(hfield(0:np_part,tcell))
      ALLOCATE(hfield_lapl(0:np_part,tcell))
      ALLOCATE(gg(1:np,1:3))
      ALLOCATE(s(1:np))
      ALLOCATE(shfield(0:np))
      ALLOCATE(shfield_lapl(0:np))
      ALLOCATE(chem(1:np))
      ALLOCATE(gchem(1:np,1:3))
      ALLOCATE(r(1:tcell))
      ALLOCATE(velocity(1:tstep))
      ALLOCATE(r_cm(1:tcell,1:3))
      ALLOCATE(r_cm_part(1:tcell,1:3))
      ALLOCATE(r_cmg(1:tcell,1:3))
      ALLOCATE(path(1:np))
      ALLOCATE(volume(1:tcell))
      !call system('rm 001/*')

!      ALLOCATE(gammaw(ntype))
!      gammaw(1:3) = (/ 1.0, 2.0 /)
      cell(0,:)%phi = 0.d0
      cell(0,:)%mu = 0.d0
      cell(0,:)%lapl_mu = 0.d0
      cell(0,:)%lapl_phi = 0.d0
      cell(0,:)%lapl_h = 0.d0
      ! initializing space matrices
      call space_init(Lsize, lxyz, lxyz_inv, np_bndry, np, periodic)  ! auxiliar fields
      call space_init(Lsize_part, lxyz_part, lxyz_inv_part, np_part_bndry, np_part, .false.) 
      ! single cell fields Dirichlet boundary condition
  
      ! generating sphere points
      call gen_cell_points(cell_radius,sphere,np_sphere)

      call single_cell_init(cell, tcell, ncell, lxyz_part, lxyz_inv_part, sphere, np_sphere, np_part, first=.true.)

      call cahn_hilliard(cell(0:np_part,1), 100, np_part, 1.0, lxyz_part, lxyz_inv_part)

      if(tcell.gt.1) then
         call single_cell_init(cell, tcell, ncell, lxyz_part, lxyz_inv_part, sphere, np_sphere, np_part, first=.false.)
      end if


      volume_target = (4.d0/3.d0)*M_PI*cell_radius**3
      ! initializing simulation box

      if(tcell.eq.2) then
        r(1) = lxyz_inv(0,-7,0) !(-6,10)
        r(2) = lxyz_inv(0,7,0) !(8,-10)
      elseif(tcell.eq.1) then
        r(1) = lxyz_inv(-15,-5,0)
      else
        dr = (/ int(2.0*cell_radius),  int(2.0*cell_radius), ntype*int(anint(2.0*cell_radius)) /)
        dri = (/ 2*int(cell_radius), 2*int(cell_radius), 2*int(cell_radius) /)

        do itype=1, ntype
          drf = (/ 2*cell_radius, 2*cell_radius, cell_radius + (ntype-itype)*(int(anint(2*cell_radius))) /)
          dr = (/ int(2.0*cell_radius), int(2.0*cell_radius), ntype*int(anint(2.0*cell_radius)) /)
          nleap =  tcell-ncell(itype)
          temp = ncell(itype) + (itype-1)*ncell(itype-1)
          call cell_pos_init(r, icell, nleap,  dr, dri, drf, temp, cell_radius, &
                                          lxyz, lxyz_inv, Lsize, np, sphere, np_sphere, iseed)
          dri(3) = dri(3) +   2*int(cell_radius)


        end do
      end if


    !  call hfield_calc(cell, aux, r, lxyz, lxyz_inv, lxyz_part, lxyz_inv_part, ncell, tcell, ntype, np_part )


      ! printing header
      call print_header(Lsize, tcell, ntype, ncell, dir_name, periodic)


      call gen_cell_points(2.0,porous,np_porous)
      call substrate_init(density, s, np, porous, np_porous, Lsize, lxyz, lxyz_inv, iseed)

      do i=1, np_sphere
        s(lxyz_inv(lxyz(r(1),1)+sphere(i,1),lxyz(r(1),2)+sphere(i,2),lxyz(r(1),3)+sphere(i,3) ) ) = 0.d0
      end do

      ! calculating the h(s) function
      do ip=1, np
        shfield(ip) = hfunc(s(ip))
      end do
      ! calculating the laplacian of the h(s) field substrate
      call dderivatives_lapl(shfield(0:np), shfield_lapl(1:np), np, dr, lxyz, lxyz_inv)


      call output(s, 300, 7, 'subs0000',dir_name, 1, ntype, np, lxyz)


      cell(0,:)%phi = 0.d0
      cell(0,:)%mu = 0.d0

      nstep = 0
      counter = 0
      cm_calc_counter = 0
      vol_lagrangian = 1.d0
      adh2 = 1.49*adh1
      metcoef = 0.d0
      vcounter = 100
      path(:) = 0.d0
      icell = 1
      open(UNIT=500, FILE=dir_name//'/vt.dat')
      write(*,'(A)') "Initiating the core program... "



      !*!*!*!*!*!*!*!*!*!*!*!
      !!     MAIN LOOP     !!
      !*!*!*!*!*!*!*!*!*!*!*!
      do while(nstep<=tstep)
         nstep = nstep + 1

         call CPU_TIME(time_init)

         ! calculating gradient of vegf
         !call dderivatives_grad(cell, gg, np, lxyz, lxyz_inv, dr)

         !call hfield_calc(cell, aux, r, lxyz, lxyz_inv, lxyz_part, lxyz_inv_part, ncell, tcell, ntype, np_part )

   

         ! calculating laplacian of phi
        call dderivatives_lapl(cell(0:np,icell)%phi, cell(1:np_part,icell)%lapl_phi ,&
                 np_part, dr, lxyz_part, lxyz_inv_part)


        call volume_calc(volume, cell, tcell, np, np_part, lxyz, lxyz_inv, lxyz_inv_part)

        adhesion(:,:) = 0.d0

        ! adhesion correction in order to avoid divergence
        call dderivatives_lapl(hfield(0:np_part,icell), hfield_lapl(1:np_part,icell), &
        np_part, dr, lxyz_part, lxyz_inv_part)
        
        do ip = 1, np_part

          call vec_local2global(ip_global, r(icell), ip, lxyz, lxyz_inv, lxyz_part)
       
          !! chemical response only the x-axis
          chemresponse = cell(lxyz_inv_part(lxyz_part(ip,1)+1,lxyz_part(ip,2),lxyz_part(ip,3)),icell)%phi - &
            cell(lxyz_inv_part(lxyz_part(ip,1)-1,lxyz_part(ip,2),lxyz_part(ip,3)),icell)%phi 

          cell(ip,icell)%mu = - chi*chemresponse/2.d0 + interface_width*cell(ip,icell)%lapl_phi +&
            cell(ip,icell)%phi*(1.d0-cell(ip,icell)%phi)*(cell(ip,icell)%phi - 0.50 + &
            vol_lagrangian*(volume_target-volume(icell))  - gamma*hfunc(s(ip_global)) + &
            shfield_lapl(ip_global) + eta*hfield_lapl(ip,icell)  ) 
        
        end do


        ! calculating laplacian(Gamma_l - hfunc(phi_m)*delta_k)

 

        cell(1:np_part, 1:tcell)%phi = cell(1:np_part,1:tcell)%phi + dt*(cell(1:np_part,1:tcell)%mu)

        !cm_calc_counter = cm_calc_counter + 1
        !if(cm_calc_counter.eq.10) then

        rt(1:3) = r_cmg(1,1:3)
        cm_calc_counter = 0
        call cm_calc_local(r_cm_part, cell, tcell, np_part, r, lxyz_part, lxyz_inv_part)

        ip = lxyz_inv_part(int(r_cm_part(icell,1)),int(r_cm_part(icell,2)),int(r_cm_part(icell,3)))
        call vec_local2global(ip_global, r(icell), ip, lxyz, lxyz_inv, lxyz_part)
        r_cm(icell,1:3) = real(lxyz(ip_global,1:3)) + FRACTION(r_cm_part(icell,1:3))

        path(ip_global) = 1.0

        call move(cell, r_cm, r, np, np_part, tcell, lxyz, lxyz_part, lxyz_inv, lxyz_inv_part, Lsize)



        call CPU_TIME(time_end)

        ctime = ctime + (time_end - time_init)

        if(nstep.eq.100) then
          ctime = (ctime*(tstep-nstep) )/6000.d0

          if( ctime>60.d0) then
            write(*,'(A,F10.2)') "Estimated time (hour): ",ctime/60.d0
          else
            write(*,'(A,F10.2)') "Estimated time (min): ",ctime
          end if
        end if



        ! output

        counter = counter + 1


        if(counter.eq.output_period) then
          write(*,*) nstep

          call format_this(nstep,format_string)
       
          write(file_name,format_string) nstep

          open(UNIT=101010, FILE='phi'//trim(file_name))
          !!do ip=1, np_part
          !!  call vec_local2global(ip_global, r(1), ip, lxyz, lxyz_inv, lxyz_part)
          !!  write(101010,'(I10,I10,I10,F10.2)') lxyz(ip_global,1:3), cell(ip,1)%phi
          !!end do
          call output_vti(file_name)
          counter = 0

        end if
         ! end of the output


      end do
      close(500)



    velocity(1) = (sum(path(1:np))*1.25)/(nstep*dt*0.26)
    write(*,*) velocity(1)

    end subroutine run_cells


    subroutine cm_calc_local(r_cm_part, f, tcell, np_part, r, lxyz_part, lxyz_inv_part)

      implicit none

      integer, intent(in) :: tcell, np_part
      integer, allocatable, intent(in) :: lxyz_part(:,:), lxyz_inv_part(:,:,:)
      type(mesh_t), allocatable, intent(inout) :: f(:,:)
      real, allocatable, intent(inout) :: r_cm_part(:,:)
      integer, allocatable, intent(inout) :: r(:)

      integer :: ip, icell, ncell, ip_part, ip_global_m
      real :: volume(1:tcell)

      volume(:) = 0.d0
      r_cm_part(:,:) = 0.d0
      do ip_part=1, np_part

         do icell=1, tcell
            !ip_global_m = r(icell)

            !call vec_global2local(ip_part, ip_global_m, ip, lxyz, lxyz_inv, lxyz_inv_part)

            if(f(ip_part,icell)%phi.gt.0.d0) then
               volume(icell) = volume(icell) + f(ip_part,icell)%phi
               r_cm_part(icell,1:3) = r_cm_part(icell,1:3) + f(ip_part,icell)%phi*lxyz_part(ip_part,1:3)
            end if
         end do
      end do
      ! Volume = sum phi_i
      ! (sum phi_i r_i )/Volume
      r_cm_part(1:tcell,1) = r_cm_part(1:tcell,1)/volume(1:tcell)
      r_cm_part(1:tcell,2) = r_cm_part(1:tcell,2)/volume(1:tcell)
      r_cm_part(1:tcell,3) = r_cm_part(1:tcell,3)/volume(1:tcell)


    end subroutine cm_calc_local


    subroutine volume_calc(volume, f, tcell, np, np_part, lxyz, lxyz_inv, lxyz_inv_part)

      implicit none

      real, allocatable, intent(inout) :: volume(:)
      integer, intent(in) :: np, tcell, np_part
      integer, allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:,:), lxyz_inv_part(:,:,:)
      type(mesh_t), allocatable, intent(inout) :: f(:,:)
      integer :: ip, icell, ncell, ip_part, ip_global_m

      volume(:) = 0.d0
      do ip=1, np

         do icell=1, tcell
            ip_global_m = r(icell)

            call vec_global2local(ip_part, ip_global_m, ip, lxyz, lxyz_inv, lxyz_inv_part)

            if(f(ip_part,icell)%phi.gt.0.d0) then
               volume(icell) = volume(icell) + f(ip_part,icell)%phi
            end if
         end do
      end do
      ! Volume = sum phi_i


    end subroutine volume_calc

    subroutine cm_calc(r_cm, f, tcell, np, np_part, r, lxyz, lxyz_inv, lxyz_inv_part)

      implicit none

      integer, intent(in) :: np, tcell, np_part
      integer, allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:,:), lxyz_inv_part(:,:,:)
      type(mesh_t), allocatable, intent(inout) :: f(:,:)
      real, allocatable, intent(inout) :: r_cm(:,:)
      integer, allocatable, intent(inout) :: r(:)

      integer :: ip, icell, ncell, ip_part, ip_global_m
      real :: volume(1:tcell)

      volume = 0.d0
      r_cm(tcell,1:3) = 0.d0
      do ip=1, np

         do icell=1, tcell
            ip_global_m = r(icell)

            call vec_global2local(ip_part, ip_global_m, ip, lxyz, lxyz_inv, lxyz_inv_part)

            if(f(ip_part,icell)%phi.gt.0.d0) then
               volume(icell) = volume(icell) + f(ip_part,icell)%phi
               r_cm(icell,1:3) = r_cm(icell,1:3) + f(ip_part,icell)%phi*lxyz(ip,1:3)
            end if
         end do
      end do
      ! Volume = sum phi_i
      ! (sum phi_i r_i )/Volume
      r_cm(1:tcell,1) = r_cm(1:tcell,1)/volume(1:tcell)
      r_cm(1:tcell,2) = r_cm(1:tcell,2)/volume(1:tcell)
      r_cm(1:tcell,3) = r_cm(1:tcell,3)/volume(1:tcell)

    end subroutine cm_calc

    subroutine move(f, r_cm, r, np, np_part, tcell, lxyz, lxyz_part, lxyz_inv, lxyz_inv_part,Lsize)

      type(mesh_t), allocatable, intent(inout) :: f(:,:)
      real, allocatable, intent(in) :: r_cm(:,:)
      integer, allocatable, intent(inout) :: r(:)
      integer, intent(in) :: np, np_part, tcell, Lsize(3)
      integer, allocatable, intent(in) :: lxyz(:,:), lxyz_part(:,:), lxyz_inv(:,:,:), lxyz_inv_part(:,:,:)
      integer :: ip, icell, ncell, ip_part, ip_global_m
      real ::  delta_r(1:tcell,1:3), ftemp(1:np_part), dimg(1:tcell,1:3)

      do icell=1, tcell
        ! the geringonca is working, so ignore my previous comments Jan 12 2016
        ! I think the problem is the round error in the delta_r
        delta_r(icell,1:3) = r_cm(icell,1:3)-lxyz(r(icell),1:3)

        !dimg(icell,1) = min(abs(delta_r(icell,1)), 2*Lsize(1)-1-abs(delta_r(icell,1)))
        !dimg(icell,2) = min(abs(delta_r(icell,2)), 2*Lsize(2)-1-abs(delta_r(icell,2)))

        !if( sqrt(delta_r(icell,1)*delta_r(icell,1)+delta_r(icell,2)*delta_r(icell,2)) .ge. 1.d0 ) then

          do ip_part=1, np_part
              ! min image method is needed here! maybe.. i don't know
              ftemp(ip_part) = f(lxyz_inv_part( lxyz_part(ip_part,1) + int(anint(delta_r(icell,1))), &
                                                lxyz_part(ip_part,2) + int(anint(delta_r(icell,2))), &
                                                lxyz_part(ip_part,3) + int(anint(delta_r(icell,3)))), icell)%phi
          end do

          f(1:np_part,icell)%phi = ftemp(1:np_part)
          r(icell) = lxyz_inv(int(anint(r_cm(icell,1))),&
                              int(anint(r_cm(icell,2))),&
                              int(anint(r_cm(icell,3))) )
        !end if
      end do

    end subroutine move



end module run_cells_m
