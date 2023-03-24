program PolyTest
    use Precision
    use Simplify
    implicit none

    real(dp),dimension(:),allocatable :: x
    integer :: i

    allocate(x(10))

    do i=1,10
        x(i) = 1.0_dp*i
    end do
    
    write(*,*) "x=", x

    write(*,*) "3th of x=", nthPoint(x,3)
    write(*,*) "within 2 of x=", radialDistance(x,1._dp)

end program PolyTest