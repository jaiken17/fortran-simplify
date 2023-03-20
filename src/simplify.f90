module Simplify
    use Precision
    implicit none


    private
    public nthPoint


    interface nthPoint
        procedure :: nthPointMulti
        procedure :: nthPointSingle
    end interface nthPoint


    interface euclidDistance
        procedure :: euclidDistanceMulti
        procedure :: euclidDistance1D
    end interface euclidDistance


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

        if (length < 2) then
            stop 'curve is too short'
        end if


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



! ~~~~~~~ Radial Distance ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    function radialDistanceMulti(curve,tolerance) result(simpleCurve)
        real(dp),dimension(:,:),allocatable :: simpleCurve
        real(dp),dimension(:,:),intent(in) :: curve
        real(dp),intent(in) :: tolerance

        integer :: i, length
        integer :: newLength
        logical :: notLast
        real(dp) :: distance
        real(dp),dimension(:),allocatable :: current

        if (tolerance <= 0._dp) then
            stop 'tolerance must be >0'
        end if

        if (size(curve,dim=1) < 2) then
            stop 'curve is too short'
        end if

        length = size(curve,dim=1)

        allocate(simpleCurve(length,size(curve,dim=2)))


        simpleCurve(1,:) = curve(1,:)
        newLength = 1
        i = 1
        notLast = .true.
        do while (notLast)
            current = simpleCurve(newLength,:)
            i = i + 1
            if (i==length) then
                exit
            end if
            distance = euclidDistance(current,curve(i,:))
            if (distance > tolerance) then
                newLength = newLength + 1
                simpleCurve(newLength,:) = curve(i,:)
            end if
        end do

        newLength = newLength + 1
        simpleCurve(newLength,:) = curve(length,:)

        simpleCurve = simpleCurve(1:newLength,:)

    end function radialDistanceMulti

    function radialDistanceSingle(curve,tolerance) result(simpleCurve)
        real(dp),dimension(:),allocatable :: simpleCurve
        real(dp),dimension(:),intent(in) :: curve
        real(dp),intent(in) :: tolerance

        integer :: i, length
        integer :: newLength
        logical :: notLast
        real(dp) :: distance
        real(dp) :: current

        if (tolerance <= 0._dp) then
            stop 'tolerance must be >0'
        end if

        if (size(curve,dim=1) < 2) then
            stop 'curve is too short'
        end if

        length = size(curve,dim=1)

        allocate(simpleCurve(length))


        simpleCurve(1) = curve(1)
        newLength = 1
        i = 1
        notLast = .true.
        do while (notLast)
            current = simpleCurve(newLength)
            i = i + 1
            if (i==length) then
                exit
            end if
            distance = euclidDistance(current,curve(i))
            if (distance > tolerance) then
                newLength = newLength + 1
                simpleCurve(newLength) = curve(i)
            end if
        end do

        newLength = newLength + 1
        simpleCurve(newLength) = curve(length)

        simpleCurve = simpleCurve(1:newLength)

    end function radialDistanceSingle

    function euclidDistanceMulti(point1,point2) result(distance)
        real(dp) :: distance
        real(dp),dimension(:),intent(in) :: point1,point2

        real(dp) :: squareSum
        integer :: i,length

        length = size(point1,dim=1)
        if (length /= size(point2,dim=1)) then
            stop 'length of vectors not equal'
        end if

        squareSum = 0._dp
        do i=1,length
            squareSum = squareSum + (point1(i) - point2(i))**2
        end do

        distance = sqrt(squareSum)

    end function euclidDistanceMulti

    function euclidDistance1D(point1,point2) result(distance)
        real(dp) :: distance
        real(dp),intent(in) :: point1,point2

        distance = abs(point1-point2)

    end function euclidDistance1D


end module Simplify