module PolySimplify
    use Precision
    implicit none


    private
    public nthPoint




    interface nthPoint
        procedure :: nthPointMulti
        procedure :: nthPointSingle
    end interface nthPoint



contains


! ~~~~~~~ Nth Point ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    function nthPointMulti(curve,n) result(simpleCurve)
        ! Function implements nth point algorithm of polyline simplification.
        ! Function acts on curve where shape is assumed to be (i,j) where j 
        ! is cooridinate index and i is parametric coordinate.

        real(dp),dimension(:,:),allocatable :: simpleCurve
        real(dp),dimension(:,:),intent(in) :: curve
        integer,intent(in) :: n

        integer :: i,length
        integer :: counter,newLength

        length = size(curve,dim=1)


        allocate(simpleCurve(length,size(curve,dim=2)))

        simpleCurve(1,:) = curve(1,:)
        counter = 1; newLength = 1;
        do i=2,length-1
            counter = counter + 1
            if (counter == n) then
                newLength = newLength + 1
                simpleCurve(newLength,:) = curve(i,:)
                counter = 0 ! reset counter
            end if
        end do
        newLength = newLength + 1
        simpleCurve(newLength,:) = curve(length,:)
        simpleCurve = simpleCurve(1:newLength,:)


    end function nthPointMulti

    function nthPointSingle(curve,n) result(simpleCurve)
        ! Function implements nth point algorithm of polyline simplification
        ! on vector.

        real(dp),dimension(:),allocatable :: simpleCurve
        real(dp),dimension(:),intent(in) :: curve
        integer,intent(in) :: n

        integer :: i,length
        integer :: counter,newLength

        length = size(curve)


        allocate(simpleCurve(length))

        simpleCurve(1) = curve(1)
        counter = 1; newLength = 1;
        do i=2,length-1
            counter = counter + 1
            if (counter == n) then
                newLength = newLength + 1
                simpleCurve(newLength) = curve(i)
                counter = 0 ! reset counter
            end if
        end do
        newLength = newLength + 1
        simpleCurve(newLength) = curve(length)
        simpleCurve = simpleCurve(1:newLength)


    end function nthPointSingle

! ~~~~~~~ End Nth Point ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


end module PolySimplify