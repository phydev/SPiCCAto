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

module routines
    
    implicit none

    private

    public :: check_boundary, &
              vec_local2global, &
              vec_global2local, &
              img, &
              format_this, &
              gen_cell_points, &
              h, &
              integral_path,&
              count_lines, &
              print_header, &
              ran2 


contains

  subroutine check_boundary(x, L, b)
    implicit none

    integer, intent(inout) :: x
    character(len=:),allocatable, intent(in) :: b
    integer, intent(in) :: L

    if(x<0) then
      if(b.eq.'periodic') then
        x = L + x 
      else if(b.eq.'Neumann') then
        x = 0 
      end if
    else if(x>=L) then
      if(b.eq.'periodic') then
        x = x - L
      else if(b.eq.'Neumann') then
        x = L - 1 
      end if
    end if
  end subroutine check_boundary

  subroutine vec_local2global(s_global, s_box, s_local, box_length)

    implicit none

    integer, intent(in) :: s_box(3), s_local(3), box_length(3)
    integer, intent(out) :: s_global(3)
    character(len=:),allocatable :: b
    b = 'periodic'
    s_global =  s_box  + s_local 

    call check_boundary(s_global(1),box_length(1),b)
    call check_boundary(s_global(2),box_length(2),b)
    call check_boundary(s_global(3),box_length(3),b)
  end subroutine vec_local2global


  subroutine vec_global2local(s_local, s_global, s_box, box_length)

    implicit none

    integer, intent(in) :: s_global(3), s_box(3), box_length(3)
    integer, intent(out) :: s_local(3)
    character(len=:),allocatable :: b
    b = 'periodic'

    s_local =  s_global - s_box

    call check_boundary(s_local(1),box_length(1),b)
    call check_boundary(s_local(2),box_length(2),b)
    call check_boundary(s_local(3),box_length(3),b)
  end subroutine vec_global2local

  function img(dx, x_size) result(x)

    implicit none

    real, intent(in) :: dx 
    integer, intent(in) :: x_size
    integer :: x 
    x = dx - anint(real(dx) / real(x_size)) * x_size
  end function img

  subroutine format_this(number,format_string)

    implicit none

    integer, intent(in) :: number
    character(len=:), allocatable, intent(inout) :: format_string    

    if (number .lt. 10) then
      format_string = "(I1)"
    else if(number.ge.10 .and. number.lt.100) then
      format_string = "(I2)"
    else if(number.ge.100 .and. number.lt.1000) then
      format_string = "(I3)"
    else if(number.ge.1000 .and. number.lt.10000) then
      format_string = "(I4)"
    else if(number.ge.10000 .and. number.lt.100000) then
      format_string = "(I5)"
    else if(number.ge.100000 .and. number.lt.1000000) then
      format_string = "(I6)"
    endif
  end subroutine format_this

  subroutine gen_cell_points(Rc,R,ndim)

    implicit none

    real, intent(in) :: Rc
    integer, intent(inout) :: ndim
    integer, allocatable, intent(inout) :: R(:,:)
    real :: dx, i, j, k, M_Pi
    integer :: counter, np_estimative
    M_Pi = 3.14159265359
    counter = 0
    np_estimative = anint((4.0/3.0)*3.14*Rc**3)
    ALLOCATE(R(1:1000,1:3))
    i = -Rc
    do while(i<=Rc)
       j = -Rc
       do while (j<=Rc)
         k = -Rc
         do while(k<=Rc)
             dx = sqrt(i**2 + j**2 + k**2)
             if(dx <= Rc) then
                counter = counter + 1
                R(counter,1) = int(anint(i))
                R(counter,2) = int(anint(j))
                R(counter,3) = int(anint(k))
             end if
             k = k + 1.d0
           end do
          j = j +1.d0
       end do
       i = i + 1.d0
    end do

    ndim = counter
  end subroutine gen_cell_points

  function h(x) result(hx)
    real :: hx
    real, intent(in) :: x
    hx = x**2*(3.0-2.0*x) 
  end function h

!! calculate the integral path
  subroutine integral_path(filename)

    character(len=*), intent(in) :: filename
    real, dimension(:), allocatable ::  x, y, z, t, dx, dy, dz
    integer :: number_of_lines, i
    real :: length

    number_of_lines = count_lines(filename)
    length = 0.0

    ALLOCATE(x(0:number_of_lines))
    ALLOCATE(y(0:number_of_lines))
    ALLOCATE(z(0:number_of_lines))
    ALLOCATE(t(0:number_of_lines))

    ALLOCATE(dx(0:number_of_lines-1))
    ALLOCATE(dy(0:number_of_lines-1))
    ALLOCATE(dz(0:number_of_lines-1))


    do i=0, number_of_lines
      read(101010,'(F10.2,F10.2,F10.2,F10.2)') t(i), x(i),  y(i), z(i)
    end do
 
    close(101010)

    do i=0, number_of_lines-1
       dx(i) = x(i+1)-x(i)
       dy(i) = y(i+1)-y(i)
       dz(i) = z(i+1)-z(i)

       length = length + sqrt(dx(i)**2 + dy(i)**2 + dz(i)**2)
    end do

    print*, "Path Length: ", length
  end subroutine integral_path

  function count_lines(filename) result(nlines)
    !! function by Alexander Vogt && user31436
    !! https://stackoverflow.com/questions/32682309/number-of-lines-of-a-text-file
    implicit none
    character(len=*)    :: filename
    integer             :: nlines 
    integer             :: io

    open(10,file=filename, iostat=io, status='old')
    if (io/=0) stop 'Cannot open file! '

    nlines = 0
    do
      read(10,*,iostat=io)
      if (io/=0) exit
      nlines = nlines + 1
    end do
    close(10)
  end function count_lines

  subroutine print_header(sim_id, iseed, density, eta, chi, gamma, nc, L, nsteps)
    
    implicit none
    character(len=*), intent(in) :: sim_id
    integer, intent(in) :: iseed, nc, nsteps, L(3)
    real, intent(in) :: density, eta, chi, gamma

    write(*,'(A)') "                            Running SPiCCAto - Software PaCkage for Cell reseArch - v.3.0 "
    write(*,'(A)') " Developer       : mms@uc.pt "
    write(*,'(A)') " Repository      : https://gitlab.com/phydev/SPiCCAto "
    write(*,'(A,I3,A,I3,A,I3,A)') " Box Size        :  (",L(1)," ",L(2)," ", L(3),")"
    write(*,'(A,A)') " Simulation ID   :        ", sim_id
    write(*,'(A,I10)') " Random seed     : ", iseed
    write(*,'(A,F10.2)') " Density         : ", density
    write(*,'(A,F10.2)') " Adhesion        : ", eta
    write(*,'(A,F10.2)') " Chemotaxis      : ", chi
    write(*,'(A,F10.2)') " Repulsion       : ", gamma
    write(*,'(A,I10)') " Number of cells : ", nc
    write(*,'(A,I10)') " Iterations      : ", nsteps

  end subroutine print_header

  function ran2(idum) !! Numerical Recipes in Fortrans
    integer :: idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
    real :: ran2,AM,EPS,RNMX
    parameter (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
         IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
         NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
    integer :: idum2,j,k,iv(NTAB),iy
    save iv,iy,idum2
    data idum2/123456789/, iv/NTAB*0/, iy/0/
    if (idum.le.0) then
       idum=max(-idum,1)
       idum2=idum
       do 11 j=NTAB+8,1,-1

          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11          continue
          iy=iv(1)
       endif
       k=idum/IQ1
       idum=IA1*(idum-k*IQ1)-k*IR1
       if (idum.lt.0) idum=idum+IM1
       k=idum2/IQ2
       idum2=IA2*(idum2-k*IQ2)-k*IR2
       if (idum2.lt.0) idum2=idum2+IM2
       j=1+iy/NDIV
       iy=iv(j)-idum2
       iv(j)=idum
       if(iy.lt.1)iy=iy+IMM1
       ran2=min(AM*iy,RNMX)
       return
  end function ran2
end module routines
