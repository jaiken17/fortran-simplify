program PolyTest
    use Precision
    use Simplify
    use IO
    implicit none

    interface writeVector
        procedure :: writeVectorNoUnit
        procedure :: writeVectorUnit
    end interface writeVector




    real(dp),dimension(:),allocatable :: x
    real(dp),dimension(9,2) :: y = 1._dp*(/ (/ 0, 20, 45, 50, 100, 130, 165, 180, 200 /)   ,      &
                                       (/ 0, 110, 105, 110, 100, -20, -10, -20, 100 /) /)
    real(dp),dimension(:,:),allocatable :: simpleY
    character(len=20),dimension(:),allocatable :: headers
    integer :: i

    allocate(x(10))

    do i=1,10
        x(i) = 1.0_dp*i
    end do
    
    write(*,'(A15)',advance='no') "x="
    call writeVector(x)

    write(*,'(A15)',advance='no') "3th of x=" 
    call writeVector(nthPoint(x,3))
    write(*,'(A15)',advance='no') "within 2 of x="
    call writeVector(radialDistance(x,1._dp))

    headers = (/"#x", "y"/)
    call outputMatrixWithHeaders(y,headers,"curve.data",20)

    simpleY = perpendicularDistance(y,2._dp)

    call outputMatrixWithHeaders(simpleY,headers,"perp_simple_curve.data",21)

contains

    subroutine writeVectorNoUnit(vector,format)
        real(dp),dimension(:),intent(in) :: vector
        character(len=*),intent(in),optional :: format

        integer :: i, length
        integer :: unitOut
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

    end subroutine writeVectorNoUnit

    subroutine writeVectorUnit(vector,unit,format)
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
        write(*,*)  ! inser newline

    end subroutine writeVectorUnit


end program PolyTest