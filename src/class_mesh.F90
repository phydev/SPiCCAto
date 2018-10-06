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

module class_mesh

  use routines

  implicit none

  private

  type, public :: mesh

    real, allocatable :: grid(:)
    integer :: nodes, L(3)
    character(len=:), allocatable :: b

    contains
      procedure :: position => mesh_position
      procedure :: initialize => mesh_initialize  
      procedure :: gt => mesh_get_item
      procedure :: ip => mesh_get_index
      procedure :: laplacian => grid_laplacian
      procedure :: gradient => grid_gradient
      procedure :: copy => mesh_copy
      procedure :: output => grid_output
      procedure :: smoothing => grid_smoothing
  end type mesh

  type, extends(mesh), public :: cell
    real :: rcom(3), volume

    contains         
      procedure :: com => calculate_com
      procedure :: initialize => cell_initialize
  end type cell

  public :: substrate_init

contains

  subroutine mesh_initialize(this, x1, x2, x3, b, s, radius)
    class(mesh) :: this
    integer, intent(in) :: x1, x2, x3
    character(*), intent(in) :: b
    real, optional :: s(3), radius

    this%L = (/x1, x2, x3/)
    this%nodes = x1*x2*x3
    this%b = b

    ALLOCATE(this%grid(0:this%nodes-1))
  end subroutine mesh_initialize

  subroutine cell_initialize(this, x1, x2, x3, b, s, radius)
    class(cell) :: this
    integer, intent(in) :: x1, x2, x3
    character(*), intent(in) :: b
    real,  optional :: s(3), radius
    integer :: ip
    real :: r(3)

    this%L = (/x1, x2, x3/)
    this%nodes = x1*x2*x3
    this%b = b
    ALLOCATE(this%grid(0:this%nodes-1))

    if(present(s)) this%rcom = s

    do ip=0, this%nodes
      r = this%position(ip)

      if( sqrt(sum((s-r)*(s-r))) <= radius) then
       this%grid(ip) = 1.0
      end if
    end do
    
  end subroutine cell_initialize

  function mesh_position(this, index) result(position)
    class(mesh) :: this
    integer, intent(in) :: index
    integer :: position(3)

    if(index<=this%nodes) then
        position(3) = floor( real(index)/real((this%L(1)*this%L(2))) )      
        position(2) = floor(real(index - position(3)*this%L(1)*this%L(2))/real(this%L(1)))
        position(1) = index - this%L(1)*position(2) - position(3)*this%L(1)*this%L(2)
    else
        STOP "class_mesh.F90 -> mesh_position % Error: index out of bounds!"
    end if
  end function mesh_position

  function mesh_get_index(this, s) result(index)
    class(mesh) :: this
    integer, intent(in) :: s(3)
    integer :: index, x(3)
    x = s
    call check_boundary(x(1),this%L(1),this%b)
    call check_boundary(x(2),this%L(2),this%b)
    call check_boundary(x(3),this%L(3),this%b)

    index = x(1) + this%L(1)*x(2) + this%L(1)*this%L(2)*x(3)
  end function mesh_get_index

  function mesh_get_item(this, index) result(item)
    class(mesh) :: this
    integer, intent(in) :: index
    real :: item 

    if(index<=this%nodes) then
        item = this%grid(index)
    else
        STOP "class_mesh.F90 -> mesh_get_item % Error: index out of bounds!"
    end if
  end function mesh_get_item

  function grid_laplacian(this,index) result(laplacian)
    class(mesh) :: this
    integer, intent(in) :: index
    real :: laplacian, y, yh, yl, weight
    integer :: s(3), i

    s = this%position(index)
    y = this%gt(index)
    weight = 1.0
    laplacian = 0.0

    do i=1, 3
        s(i) = s(i) + 1
        yh = this%gt(this%ip(s))
        s(i) = s(i) - 2
        yl = this%gt(this%ip(s))
        s(i) = s(i) + 1         
        !weight = 1.0/(phi.dr*phi.dr)
        laplacian = laplacian + (yh + yl - 2.0*y)*weight
    end do
  end function grid_laplacian

  function grid_gradient(this,index) result(gradient)
    class(mesh) :: this
    integer, intent(in) :: index
    real ::  gradient(3),  yh, yl, weight
    integer :: s(3), i

    s = this%position(index)
    weight = 0.5
    gradient = (/ 0.0, 0.0, 0.0 /)

    do i=1, 3
        s(i) = s(i) + 1
        yh = this%gt(this%ip(s))
        s(i) = s(i) - 2
        yl = this%gt(this%ip(s))
        s(i) = s(i) + 1         
        !weight = 1.0/(phi.dr*phi.dr)
        gradient(i) = (yh - yl) * weight
    end do
  end function grid_gradient

  subroutine grid_smoothing(this)
    class(mesh) :: this
    integer :: cycles, ip
    real :: volume_target, volume
    real, allocatable :: laplacian(:)

    ALLOCATE(laplacian(0:this%nodes))
    volume = 0.0
    volume_target = 0.0

    do ip=0, this%nodes
      volume_target = volume_target + this%gt(ip)
    end do

    do cycles = 1, 100
      volume = 0.0
      do ip=0, this%nodes
        laplacian(ip) = this%laplacian(ip)
        volume = volume + this%gt(ip)
      end do

      do ip=0, this%nodes
        this%grid(ip) = this%grid(ip) + 0.0001*( laplacian(ip) + this%gt(ip)*(1.0-this%gt(ip))*(this%gt(ip) -0.5 &
         + (volume_target-volume) ) ) 
      end do
    end do

    DEALLOCATE(laplacian)
  end subroutine grid_smoothing

  subroutine mesh_copy(this,other)
    class(mesh), intent(inout) :: this
    class(mesh), intent(in) :: other

    if(.not.ALLOCATED(this%grid) ) then
      ALLOCATE(this%grid(0:other%nodes))
    else if(this%nodes .ne. other%nodes) then
      DEALLOCATE(this%grid)
      ALLOCATE(this%grid(0:other%nodes))
    end if

    this%nodes = other%nodes
    this%grid = other%grid
    this%L = other%L
    this%b = other%b
  end subroutine mesh_copy

  subroutine calculate_com(this)

    class(cell) :: this
    real :: com(3), volume
    integer :: icom(3)
    integer :: ip, s(3)

    volume = 0.0
    com = this%rcom
    icom = int(this%rcom)

    do ip=0, this%nodes
      volume = volume + this%gt(ip)

      call vec_local2global(s, icom, this%position(ip))

      com = com + this%gt(ip)*s
    end do

    this%rcom = com/volume
    this%volume = volume
  end subroutine

  
  subroutine grid_output(this, filename, other, rcom)

    class(mesh) :: this
    class(mesh), optional, intent(in) :: other
    integer, optional, intent(in) :: rcom(3)
    character(len=10), intent(in) :: filename
    integer :: tag, ip, L(3), nodes, s(3), icom(3), ip_new

    tag = 33423
    if(present(other)) then
      L = other%L
      nodes = other%nodes
    else
      L = this%L
      nodes = this%nodes
    end if

    OPEN(UNIT=tag, FILE=trim(filename)//".vti" )
    write(tag,'(A)')'<?xml version="1.0"?>' 
    write(tag,'(A)')'<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">'
    write(tag,'(A,I1,A,I3,A,I1,A,I3,A,I1,A,I3,A)')'  <ImageData WholeExtent="',0,' ',L(1),&
    ' ',0,' ',L(2),' ',0,' ',L(3),'"   Origin="0 0 0" Spacing="1 1 1">'
    write(tag,'(A,I1,A,I3,A,I1,A,I3,A,I1,A,I3,A)')'    <Piece Extent="',0,' ',L(1),' ',0,' ',L(2),' ',0,' ',L(3),'">'
    write(tag,*)'      <CellData> '
    write(tag,*)'        <DataArray Name="scalar_data" type="Float64" format="ascii">'
    do ip=0, nodes
      if(present(other)) then 
        ip_new = this%ip(other%position(ip))
      else
        ip_new = ip
      end if
      write(tag,'(F10.2)', ADVANCE='no') this%gt(ip_new)
    end do
    write(tag,'(A)') "" 
    write(tag,'(A)')"         </DataArray>"  
    write(tag,'(A)')"      </CellData>"
    write(tag,'(A)')"    </Piece>"
    write(tag,'(A)')"</ImageData>"
    write(tag,'(A)')"</VTKFile>"
    CLOSE(tag)

  end subroutine grid_output


  subroutine substrate_init(this, density_target, sphere, np_sphere, idum, com, radius)

    type(mesh), intent(inout) :: this
    real, intent(in) :: density_target
    integer, allocatable, intent(in) :: sphere(:,:)
    integer, intent(in) :: np_sphere, idum
    real, optional, intent(in) :: com(3), radius
    integer :: l, s(3), ip, r(3)
    real :: density, s_volume, t_volume

    t_volume = this%L(1)*this%L(2)*this%L(3)
    density = 0.0

    do while (density.le.density_target)
      ip = ran2(idum)*this%nodes

      do l=1, np_sphere
        s = this%position(ip) + sphere(l,1:3)
        this%grid(this%ip(s)) = 1.d0
      end do
      s_volume = 0.d0

      do ip=0, this%nodes
        if(this%grid(ip).gt.0.0) s_volume = s_volume +1.d0
      end do

      density = s_volume/t_volume

    end do

    if(present(com).and.present(radius)) then
      do ip=0, this%nodes
       r = this%position(ip)
       if( sqrt(sum((com-r)*(com-r))) <= radius) then
        this%grid(ip) = 0.0
        end if
      end do
    end if

  end subroutine substrate_init
end module class_mesh
