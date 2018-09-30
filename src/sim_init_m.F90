module sim_init_m

  use global_m
  use misc_m

  implicit none

  private

  public :: write_parameters, substrate_init, single_cell_init, cell_pos_init, space_init, print_header

  contains

    subroutine substrate_init(density_target, s, np, sphere, np_sphere, Lsize, lxyz, lxyz_inv, idum)

      real, allocatable, intent(inout) :: s(:)
      real, intent(in) :: density_target
      integer, allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:,:)
      integer, intent(in) :: np, Lsize(3), sphere(1000,1:3), np_sphere, idum
      integer :: iParticle, iseed, l
      real :: density, s_volume, t_volume

      t_volume = 8.d0*Lsize(1)*Lsize(2)*Lsize(3)

      s(:) = 0.d0
      porosity = 1.d0

      do while (density.le.density_target)
        ip = ran2(idum)*np

        do l=1, np_sphere
          s(lxyz_inv(lxyz(ip,1)+sphere(l,1),lxyz(ip,2)+sphere(l,2),lxyz(ip,3)+sphere(l,3) ) ) = 1.d0
        end do
        s_volume = 0.d0
        do ip=1, np
          if(s(ip).gt.0) s_volume = s_volume +1.d0
        end do

        density = s_volume/t_volume

      end do

    end subroutine substrate_init



   subroutine  single_cell_init(cell, tcell, ncell, lxyz_part, lxyz_inv_part, sphere, np_sphere, np_part, first)

      implicit none

      ! input variables
      type(mesh_t), allocatable, intent(inout) :: cell(:,:)
      integer, allocatable, intent(in) :: lxyz_part(:,:), lxyz_inv_part(:,:,:), ncell(:)
      integer, intent(in) ::  tcell, sphere(:,:), np_sphere, np_part
      logical, intent(in) :: first

      ! internal only
      integer :: ip, i, j, l, icell, itype

      itype = 1
      i = 1
      if(first) then
         cell(0:np_part,1)%phi = 0.d0
         cell(0:np_part,1)%itype = 1.d0
         do l = 1, np_sphere
            cell(lxyz_inv_part(sphere(l,1),sphere(l,2),sphere(l,3)),1)%phi = 1.d0
            cell(lxyz_inv_part(sphere(l,1),sphere(l,2),sphere(l,3)),1)%itype = 1.d0
         end do
      else
         do icell=2,tcell
            i = i +1
            cell(0:np_part,icell)%phi = cell(0:np_part,1)%phi
            if(i.gt.ncell(itype)) then
               itype = itype + 1
               i = 1
            end if
            cell(0:np_part,icell)%itype = itype

         end do

      end if

    end subroutine single_cell_init



    subroutine cell_pos_init(r, icell, nleap,  dr, dri, drf, ncell, cell_radius, lxyz, lxyz_inv,&
         Lsize, np, sphere, np_sphere, iseed)

      implicit none

      ! input variables
      integer, allocatable, intent(inout) :: r(:)
      integer, intent(inout) ::  dr(3), icell
      integer, allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:,:)
      integer, intent(in) ::  np, iseed, Lsize(3),  ncell, sphere(:,:),&
           np_sphere, drf(3), dri(3), nleap
      real, intent(in) :: cell_radius
      ! internal only
      integer :: ip, i, j, k, l


      icell = 1
      i=-Lsize(1)+dri(1)
      do while  (i<= Lsize(1)-drf(1))
          j= -Lsize(2)+dri(2)
         do while (j<= Lsize(2)-drf(2))
            k = -Lsize(3)+dri(3)
            do while (k<= Lsize(3) - drf(3))
                ip = lxyz_inv(i,j,k)

                r(icell) = ip

                icell = icell + 1
                k = k + dr(3)

               !if(icell.ge.nleap) dr(2) = 2*cell_radius
            if (icell.gt.ncell) EXIT
          end do
          j = j + dr(2)
        end do
         i = i + dr(1)
         if (icell.gt.ncell) EXIT
      end do



    end subroutine cell_pos_init

    subroutine space_init(Lsize, lxyz, lxyz_inv, np_bndry, np, periodic)

      implicit none

      ! input/output variables
      integer, intent(in) ::  Lsize(1:3), np_bndry
      integer, allocatable, intent(inout) :: lxyz(:,:), lxyz_inv(:,:,:)
      ! internal variables
      integer :: i, j, k, l, m, n, ip, ip_part, np, ki, kf
      real :: hs(1:3)
      logical :: boundary, periodic

      ip = 0
      ip_part = np

      ! bulk points
      ! allocated from 1 to np

      do k=-Lsize(3), Lsize(3)-1
         do j=-Lsize(2), Lsize(2)-1
           do i=-Lsize(1), Lsize(1)-1

                ip = ip + 1

                lxyz(ip,1) = i
                lxyz(ip,2) = j
                lxyz(ip,3) = k

                lxyz_inv(i,j,k) = ip

             end do
         end do
      end do

      ! boundary points
      ! allocated from np to np_part

      do k=-Lsize(3)-np_bndry, Lsize(3)-1+np_bndry
         do j=-Lsize(2)-np_bndry, Lsize(2)-1+np_bndry
           do i=-Lsize(1)-np_bndry, Lsize(1)-1+np_bndry
               l = i
               m = j
               n = k

               boundary = .false.

               hs(1) = heaviside(real(i))
               hs(2) = heaviside(real(j))
               hs(3) = heaviside(real(k))
               if( abs(i)>Lsize(1)-hs(1)) then
                  boundary = .true.
                  if(periodic) then
                     l = i - SIGN(1,i)*(2*Lsize(1))
                  else
                     l = i - SIGN(1,i)
                  end if
               end if

               if( abs(j)>Lsize(2)-hs(2)) then
                  boundary = .true.

                  if(periodic) then
                     m = j  - SIGN(1,j)*(2*Lsize(2))
                  else
                     m = j - SIGN(1,j)
                  end if
               end if

               if( abs(k)>Lsize(3)-hs(3)) then
                  boundary = .true.

                  if(periodic) then
                     n = k  - SIGN(1,k)*(2*Lsize(3))
                  else
                     n = k - SIGN(1,k)
                  end if
               end if


               if(boundary) then
                  if(periodic) then
                     ip_part = ip_part + 1

                     lxyz(ip_part,1) = l
                     lxyz(ip_part,2) = m
                     lxyz(ip_part,3) = n

                     lxyz_inv(i,j,k) = lxyz_inv(l,m,n)
                  else
                     !ip_part = ip_part + 1
                     !lxyz(ip_part,1) = l
                     !lxyz(ip_part,2) = m

                     lxyz_inv(i,j,k) =  lxyz_inv(l,m,n)

                  end if
               end if



               ! the updating of boundaries should be
               ! cell(ip_part) = cell( lxyz_inv( lxyz(np_part,1), lxyz(np_part,2), lxyz(np_part,3) ) )
               ! because:
               ! lxyz(np_part,1), lxyz(np_part,2), lxyz(np_part,3) ->(l, m, n) in the bulk
               ! then lxyz_inv(m,n,l) will give the updated value in the bulk
               ! for the respective boundary position cell(ip_part)
             end do
         end do
      end do

    end subroutine space_init


    subroutine write_parameters(cell_radius, density, eta, gamma, chi, interface_width,&
                                 tstep, dt, Lsize,  dir_name, iseed,  periodic)

      implicit none

      real, intent(in) :: density, gamma, eta, chi, cell_radius,  interface_width, dt
      integer, intent(in) :: tstep, Lsize(3)
      integer, intent(in) :: iseed
      character(*), intent(in) :: dir_name
      character(len=255) :: temp
      logical :: periodic


      OPEN (UNIT=2,FILE=trim(dir_name)//'/parameters'//trim(dir_name)//'.init')
      write(2,'(F10.2,A)') cell_radius, " cell_radius" ! R_c - Cell Radius
      write(2,'(F10.2,A)') gamma, "repulsion" ! repulsive force coefficient
      write(2,'(F10.2,A)') eta, " adhesion eta" ! adhesion eta
      write(2,'(F10.2,A)') chi, " chemotaxis chi" ! chemotaxis chi
      write(2,'(F10.2,A)') density, " density" ! density
      write(2,'(F10.2,A)') interface_width, " interface_width" ! Eps - Interface Witdh
      write(2,'(I10,A)')   tstep, " tstep" ! tstep - Total time step
      write(2,'(F10.5,A)') dt, " dt" ! dt - Time increment
      write(2,'(I3, I3, I3,  A)') Lsize(1:3),  " box_length_xyz" !Box Length - x,y,z
      write(2,'(A,A)') dir_name, " dir_name" ! Simulation name
      write(2,'(I10,A)') iseed, " iseed" ! Initial Seed for RAN2
      write(2,'(L1,A) ') periodic, " periodic" ! boundary conditions
      CLOSE(2)
    end subroutine write_parameters



    subroutine print_header(Lsize, tcell, ntypes, ncell, dir_name, periodic)

      implicit none
      integer, intent(in) :: Lsize(1:3), ntypes, ncell(:), tcell
      character(len=255) :: cwd, hostname
      character(len=32) :: username
      character(8)  :: date
      character(10) :: time
      character(5)  :: zone
      character(3)  :: dir_name
      integer, dimension(8) :: values
      logical :: periodic

      call date_and_time(date,time,zone,values)
      call date_and_time(DATE=date,ZONE=zone)
      call date_and_time(TIME=time)

      call hostnm(hostname)
      call getcwd(cwd)
      call getlog(username)
      write(*,'(A)') "                                Running Cell3D"
      write(*,'(A)') "       "
      write(*,'(A)') "Version        :       1.0.s (December 14, 2016)"
      write(*,'(A,A)') "Locate         :       ", trim(cwd)
      write(*,'(A,A)') "User           :       ", trim(username)
      write(*,'(A)') "Developer      :       Moreira, M."
      write(*,'(A)') "       "
      write(*,'(A,A)') "                      The code is running in ", trim(hostname)
      write(*,'(A)') "       "
      write(*,'(A)') "       "
      write(*,'(A,2X,A,2X,A)') "             Calculation started on", date(7:8)//"/"//date(5:6)//"/"//date(1:4),&
           time(1:2)//":"//time(3:4)//":"//time(5:6)
      write(*,'(A)') "       "
      write(*,'(A)') "************************************ Grid *************************************"
      write(*,'(A)') "Simulation Box:"
      write(*,'(A,I3,A,I3,A,I3,A,I3,A)') "  Lengths = (",Lsize(1),",", Lsize(2), ",",Lsize(3), ")"
      write(*,'(A)') "  the code will run in 3 dimension(s)."
      if(periodic) then
         write(*,'(A)') "  the code will treat the system as periodic in 3 dimension(s)."
      else
         write(*,'(A)') "  the code will treat the system as non-periodic."
      end if
      write(*,'(A)') "*******************************************************************************"
      write(*,'(A,I3)') "  Number of cells:", tcell
      write(*,'(A,I3)') "  Cell types:", ntypes
      write(*,*) " Cells per type:", ncell(1:ntypes)
      write(*,'(A,A)') "  Simulation ID: ", dir_name
    end subroutine print_header

end module sim_init_m
