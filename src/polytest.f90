program PolyTest
    use Precision
    use Simplify
    implicit none

    interface writeVector
        procedure :: writeVectorNoUnit
        procedure :: writeVectorUnit
    end interface writeVector




    real(dp),dimension(:),allocatable :: x
    real(dp),dimension(10,2) :: y = (/ (/ (i, i=1,10,1) /)   ,      &
                                       (/ (i**2, i=1,10,1) /) /)
    real(dp),dimension(:,:),allocatable :: simpleY
    integer :: i

    allocate(x(10))

    do i=1,10
        x(i) = 1.0_dp*i
    end do
    
    write(*,'(A15)',advance='no') "x="
    call writeVector(x)
    write(*,*)

    write(*,'(A15)',advance='no') "3th of x=" 
    call writeVector(nthPoint(x,3))
    write(*,*)
    write(*,*) "within 2 of x=", radialDistance(x,1._dp)

    write(*,*) "y(:)=", y(:,1)
    write(*,*) "     ", y(:,2)
    simpleY = perpendicularDistance(y,2._dp)
    write(*,*) "perp of y, tol=2:", simpleY(:,1)
    write(*,*) "                 ", simpleY(:,2)

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

    end subroutine writeVectorUnit


end program PolyTest