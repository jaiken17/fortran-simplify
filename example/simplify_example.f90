!   Copyright (c) 2023, Joshua Aiken
!   All rights reserved.
!
!   This source code is licensed under the BSD-style license found in the
!   LICENSE file in the root directory of this source tree. 


program simplify_example
    use precision
    use simplify
    use io
    implicit none

    interface write_vector
        procedure :: write_vectorNoUnit
        procedure :: write_vector_unit
    end interface write_vector




    real(dp),dimension(:),allocatable :: x
    real(dp),dimension(10,2) :: y
    real(dp),dimension(:,:),allocatable :: simpleY
    integer,dimension(:),allocatable :: simple_y_indices
    character(len=20),dimension(:),allocatable :: headers
    integer :: i

    call chdir("example")   ! FPM always runs in 'root' workspace directory

    y(:,1) = 1._dp*[ 0, 20, 45, 50, 100, 130, 165, 180, 200, 198 ]
    y(:,2) = 1._dp*[ 0, 110, 105, 110, 100, -20, -10, -20, 100, 102 ]
    
    allocate(x(10))

    do i=1,10
        x(i) = 1.0_dp*i
    end do
    
    write(*,'(A15)',advance='no') "x="
    call write_vector(x)

    write(*,'(A15)',advance='no') "3th of x=" 
    call write_vector(nth_point(x,3))
    write(*,'(A15)',advance='no') "within 2 of x="
    call write_vector(radial_distance(x,1._dp))

    headers = (/'#x', ' y'/)
    call output_matrix_with_headers(y,headers,"curve.data",20)

    simpleY = perpendicular_distance(y,20._dp)
    simple_y_indices = perpendicular_distance_indices(y,20._dp)

    call output_matrix_with_headers(simpleY,headers,"perp_simple_curve.data",21)
    call output_matrix_with_headers(y(simple_y_indices,:),headers,"perp_simple_curve_indices.data",23)
    

    simpleY = reumann_witkam(y,20._dp)
    simple_y_indices = reumann_witkam_indices(y,20._dp)

    call output_matrix_with_headers(simpleY,headers,"reumann_witkam_simple_curve.data",22)
    call output_matrix_with_headers(y(simple_y_indices,:),headers,"reumann_witkam_simple_curve_indices.data",24)

contains

    subroutine write_vectorNoUnit(vector,format)
        real(dp),dimension(:),intent(in) :: vector
        character(len=*),intent(in),optional :: format

        integer :: i, length
        character(len=:),allocatable :: fmt

        if (present(format)) then
            fmt = format
        else
            fmt = '(F20.8)'
        end if

        length = size(vector)

        do i=1,length
            write(unit=*,fmt=fmt,advance='no') vector(i)
        end do
        write(*,*)  ! insert newline

    end subroutine write_vectorNoUnit

    subroutine write_vector_unit(vector,unit,format)
        real(dp),dimension(:),intent(in) :: vector
        integer,intent(in) :: unit
        character(len=*),intent(in),optional :: format

        integer :: i, length
        character(len=:),allocatable :: fmt

        if (present(format)) then
            fmt = format
        else
            fmt = '(F15.8)'
        end if

        length = size(vector)
        
        do i=1,length
            write(unit=unit,fmt=fmt,advance='no') vector(i)
        end do
        write(*,*)  ! insert newline

    end subroutine write_vector_unit


end program simplify_example