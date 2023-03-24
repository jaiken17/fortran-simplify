module Simplify
    use Precision
    implicit none


    private
    public nthPoint, radialDistance


    interface nthPoint
        procedure :: nthPointMulti
        procedure :: nthPointSingle
    end interface nthPoint

    interface radialDistance
        procedure :: radialDistanceMulti
        procedure :: radialDistanceSingle
    end interface radialDistance




    interface squareEuclidDistance
        procedure :: squareEuclidDistanceMulti
        procedure :: squareEuclidDistance1D
    end interface squareEuclidDistance


contains


! ~~~~~~~ Square Euclid Distnace ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! squareEuclidDistance is a helper function that calculates the square 
    ! euclidean distance between any two, n-dimensional points (including 1D)

    function squareEuclidDistanceMulti(point1,point2) result(distance)
        real(dp) :: distance
        real(dp),dimension(:),intent(in) :: point1,point2

        real(dp) :: squareSum
        integer :: i,length

        length = size(point1,dim=1)
        if (length /= size(point2,dim=1)) then
            stop 'length of vectors not equal'
        end if

        distance = dot_product(point1,point2)

    end function squareEuclidDistanceMulti

    function squareEuclidDistance1D(point1,point2) result(distance)
        real(dp) :: distance
        real(dp),intent(in) :: point1,point2
        
        real(dp),dimension(1) :: point1Vec, point2Vec
        
        point1Vec(1) = point1
        point2Vec(1) = point2

        distance = dot_product(point1Vec,point2Vec)

    end function squareEuclidDistance1D


! ~~~~~~~ End Square Euclid Distnace ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



! ~~~~~~~ Distance From Line to Point ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! d2LineToPoint is a helper function that computes the perpendicular 
    ! square distance from a line (defined by two points) to a point. This function
    ! works on n-dimensional points (including the trivial case of 1D).

    function d2LinetoPoint(linePoint1,linePoint2,point) result(distance2)
        real(dp) :: distance2
        real(dp),dimension(:),intent(in) :: linePoint1, linePoint2, point

        real(dp),dimension(size(linePoint1,dim=1)) :: v,xMinusA,bMinusA,vMinusX
        real(dp) :: t
        integer :: length

        length = size(linePoint1,dim=1)

        if ((length /= size(linePoint2,dim=1)) .or. (length /= size(point,dim=1))) then
            stop 'linePoint1, linePoint2, point not all the same length'
        end if


        xMinusA = point - linePoint1
        bMinusA = linePoint2 - linePoint1

        t = dot_product(xMinusA,bMinusA)/dot_product(bMinusA,bMinusA)
        v = linePoint1 + t*bMinusA

        distance2 = dot_product(vMinusX,vMinusX)

    end function d2LinetoPoint

! ~~~~~~~ End Distance From Line to Point ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! ~~~~~~~ Nth Point ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! nthPoint is a function that implements the nth point polyline 
    ! simplification algorithm. It works on n-dimensional curves (including 1D)
    ! and takes parameter "n" which is interpreted as every nth element to be 
    ! kept.

    function nthPointMulti(curve,n) result(simpleCurve)

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

    ! radialDistance is a function that implements the radial distance 
    ! polyline simplification algorithm. It works on n-dimensional curves
    ! (including 1D) and takes a parameter "tolerance" which is the minimum 
    ! distance allowed between any two points on the simplified curve (excluding
    ! the last two points).

    function radialDistanceMulti(curve,tolerance) result(simpleCurve)
        real(dp),dimension(:,:),allocatable :: simpleCurve
        real(dp),dimension(:,:),intent(in) :: curve
        real(dp),intent(in) :: tolerance
        
        real(dp) :: squareTolerance
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

        squareTolerance = tolerance*tolerance

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
            distance = squareEuclidDistance(current,curve(i,:))
            if (distance > squareTolerance) then
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
        
        real(dp) :: squareTolerance
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

        squareTolerance = tolerance*tolerance

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
            distance = squareEuclidDistance(current,curve(i))
            if (distance > squareTolerance) then
                newLength = newLength + 1
                simpleCurve(newLength) = curve(i)
            end if
        end do

        newLength = newLength + 1
        simpleCurve(newLength) = curve(length)

        simpleCurve = simpleCurve(1:newLength)

    end function radialDistanceSingle

! ~~~~~~~ End Radial Distance ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




! ~~~~~~~ Perpindicular Distance ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! perpindicularDistance is a function that implement the perpendicular
    ! distance polyline simplification algorithm. It works on n>1-dimensional
    ! curves and takes a parameter "tolerance" which is the minimum 
    ! perpendicular distance a middle point can be from the line connecting 
    ! the two points on either side of that middle point. It also takes 
    ! a parameter "repeat" which is how many times to repeat the algorithm,
    ! since at most 50% of the curves points can be removed on a single pass. 


    function perpDistance(curve,tolerance) result(simpleCurve)
        real(dp),dimension(:,:),allocatable :: simpleCurve
        real(dp),dimension(:,:),intent(in) :: curve
        real(dp),intent(in) :: tolerance

        real(dp) :: squareTolerance
        integer :: i, length, newLength
        real(dp) :: distance
        




    end function perpDistance


! ~~~~~~~ End Perpindicular Distan~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


end module Simplify