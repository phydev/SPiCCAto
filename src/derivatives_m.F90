module derivatives_m

  use global_m
  use misc_m
  implicit none

  private

  public :: cahn_hilliard, dderivatives_grad, dderivatives_lapl

  contains

    subroutine cahn_hilliard(f, nstep, np, interface_width, lxyz, lxyz_inv)

      implicit none
      integer, intent(in) :: nstep, np
      type(mesh_t), intent(inout) :: f(0:np)
      integer,allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:,:)

      real, intent(in) :: interface_width
      ! internal variables
      integer :: tstep, ip, dr(3)
      real :: dt
      dr(1:2) = 1
      dt = 0.0001

      do tstep=1, nstep
         f(0)%phi = 0.d0
         f(0)%mu = 0.d0
         f(0)%lapl_mu = 0.d0
         ! calculating laplacian of phi
         call dderivatives_lapl(f(0:np)%phi,f(1:np)%lapl_phi, np, dr, lxyz, lxyz_inv)

         ! Chemical potential:  phi**3 - phi - epsilon*laplacian(phi )
         do ip = 1, np
            f(ip)%mu = (2.0*f(ip)%phi-1.0)**3 - 2.0*f(ip)%phi + 1.0 - 2.0*interface_width*f(ip)%lapl_phi
         end do

         ! Calculating laplacian of mu
         call dderivatives_lapl(f(0:np)%mu, f(1:np)%lapl_mu, np, dr, lxyz, lxyz_inv)

         do ip=1, np
            f(ip)%phi = f(ip)%phi + dt*f(ip)%lapl_mu
         end do

      end do

    end subroutine cahn_hilliard

    subroutine dderivatives_grad(f, gradient, np, lxyz, lxyz_inv, dr)

      implicit none

      real, allocatable, intent(in) :: f(:)
      real, allocatable, intent(inout) :: gradient(:,:)
      integer, allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:,:)
      integer, intent(in) :: np,  dr(3)
      integer :: ip, i, j, k

      do ip=1, np

         i = lxyz(ip,1)
         j = lxyz(ip,2)
         k = lxyz(ip,3)

         gradient(ip,1) = (f( lxyz_inv(i+dr(1),j,k) ) -&
              f( lxyz_inv(i-dr(1),j,k) ) )/(2.d0*dr(1))

         gradient(ip,2) = (f( lxyz_inv(i,j+dr(2),k) ) -&
              f( lxyz_inv(i,j-dr(2),k) ) )/(2.d0*dr(2))

        gradient(ip,3) = (f( lxyz_inv(i,j,k+dr(3)) ) -&
                   f( lxyz_inv(i,j,k-dr(3)) ) )/(2.d0*dr(3))

      end do

    end subroutine dderivatives_grad


    subroutine dderivatives_lapl(f, lapl, np, dr, lxyz, lxyz_inv)

      implicit none
      integer, intent(in) :: np, dr(2)
      real,  intent(in) :: f(0:np)
      integer, allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:,:)
      real,  intent(inout) :: lapl(:)
      integer :: i, j, k, ip, a, b, c

      do ip=1, np

         i = lxyz(ip,1)
         j = lxyz(ip,2)
         k = lxyz(ip,3)
         lapl(ip) = f(lxyz_inv(i+1,j,k)) + f(lxyz_inv(i-1,j,k)) + &
                    f(lxyz_inv(i,j+1,k)) + f(lxyz_inv(i,j-1,k)) + &
                    f(lxyz_inv(i,j,k+1)) + f(lxyz_inv(i,j,k-1))-6.d0*f(ip)

      end do


    end subroutine dderivatives_lapl

  end module derivatives_m
