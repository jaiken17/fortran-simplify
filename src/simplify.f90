module Simplify
    use Precision
    implicit none


    private
    public nthPoint, radialDistance, perpendicularDistance, reumannWitkam


    interface nthPoint
        procedure :: nthPointMulti
        procedure :: nthPointSingle
    end interface nthPoint

    interface radialDistance
        procedure :: radialDistanceMulti
        procedure :: radialDistanceSingle
    end interface radialDistance

    interface perpendicularDistance
        procedure :: perpDistance
        procedure :: perpendicularDistanceRepeat
    end interface perpendicularDistance


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
        
        vMinusX = v - point

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
        
        ! Make necessary checks
        if (length < 2) then
            stop 'curve is too short'
        end if

        if (n < 1) then
            stop 'n is too small (<1)'
        end if

        ! Allocate maximum memory needed, will be trimmed at the end
        allocate(simpleCurve(length,size(curve,dim=2)))

        simpleCurve(1,:) = curve(1,:)   ! First element is always a key
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
        simpleCurve(newLength,:) = curve(length,:)  ! Last element is always a key
        
        simpleCurve = simpleCurve(1:newLength,:)    ! Resize final curve


    end function nthPointMulti

    function nthPointSingle(curve,n) result(simpleCurve)

        real(dp),dimension(:),allocatable :: simpleCurve
        real(dp),dimension(:),intent(in) :: curve
        integer,intent(in) :: n

        integer :: i,length
        integer :: counter,newLength

        length = size(curve)

        ! Make necessary checks
        if (length < 2) then
            stop 'curve is too short'
        end if

        if (n < 1) then
            stop 'n is too small (<1)'
        end if

        ! Allocate maximum memory needed, will be trimmed at the end
        allocate(simpleCurve(length))

        simpleCurve(1) = curve(1)   ! First element is always a key
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
        simpleCurve(newLength) = curve(length)  ! Last element is always a key
        
        simpleCurve = simpleCurve(1:newLength)  ! Resize final curve


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

        ! Make necessary checks
        if (tolerance <= 0._dp) then
            stop 'tolerance must be >0'
        end if

        if (size(curve,dim=1) < 2) then
            stop 'curve is too short'
        end if

        length = size(curve,dim=1)

        ! Allocate absolute maximum memory needed, will be trimmed at end
        allocate(simpleCurve(length,size(curve,dim=2)))

        squareTolerance = tolerance*tolerance

        simpleCurve(1,:) = curve(1,:)   ! First element is always a key
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
        simpleCurve(newLength,:) = curve(length,:)  ! Last element is always a key

        simpleCurve = simpleCurve(1:newLength,:)    ! Resize final curve

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

        ! Make necessary checks

        if (tolerance <= 0._dp) then
            stop 'tolerance must be >0'
        end if

        if (size(curve,dim=1) < 2) then
            stop 'curve is too short'
        end if

        length = size(curve,dim=1)


        ! Allocate absolute maximum amount of memory needed, will be trimmed at end

        allocate(simpleCurve(length))

        squareTolerance = tolerance*tolerance

        simpleCurve(1) = curve(1)   ! First element is always a key
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
        simpleCurve(newLength) = curve(length)  ! Last element is always a key

        simpleCurve = simpleCurve(1:newLength)  ! Resize final curve

    end function radialDistanceSingle

! ~~~~~~~ End Radial Distance ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




! ~~~~~~~ Perpendicular Distance ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! perpendicularDistance is a function that implement the perpendicular
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

        integer :: i, length, newLength
        logical :: last_is_key = .false., skip_iter = .false.
        real(dp) :: squareDistance, squareTolerance
        

        length = size(curve,dim=1)

        ! Make necessary checks

        if (tolerance <= 0._dp) then
            stop 'tolerance must be >0'
        end if

        if (size(curve,dim=1) < 2) then
            stop 'curve is too short'
        end if

        squareTolerance = tolerance*tolerance

        ! Allocate absolute maximum amount of memory, will get trimmed at end

        allocate(simpleCurve(length,size(curve,dim=2)))

        simpleCurve(1,:) = curve(1,:)   ! First element is always a key
        newLength = 1

        do i=2,length-1
            if (skip_iter) then
                skip_iter = .false.
                cycle
            end if

            squareDistance = d2LineToPoint(simpleCurve(newLength,:),curve(i+1,:),curve(i,:))
            if (squareDistance >= squareTolerance) then
                newLength = newLength + 1
                simpleCurve(newLength,:) = curve(i,:)
            else
                newLength = newLength + 1
                simpleCurve(newLength,:) = curve(i+1,:)
                skip_iter = .true.
                if (i+1 == length) last_is_key = .true.
            end if
        end do

        if (.not. last_is_key) then   ! last element may already be set as key
            newLength = newLength + 1
            simpleCurve(newLength,:) = curve(length,:)  ! Last element is always a key
        end if

        simpleCurve = simpleCurve(1:newLength,:)    ! Resize simpleCurve to only needed elements

    end function perpDistance


    function perpendicularDistanceRepeat(curve,tolerance,repeat) result(simpleCurve)
        real(dp),dimension(:,:),allocatable :: simpleCurve
        real(dp),dimension(:,:),intent(in) :: curve
        real(dp),intent(in) :: tolerance
        integer,intent(in) :: repeat

        integer :: i
        
        if (repeat < 1) then
            stop 'number of repeats too small'
        end if

        simpleCurve = curve

        do i=1,repeat
            simpleCurve = perpDistance(simpleCurve,tolerance)
        end do


    end function perpendicularDistanceRepeat

! ~~~~~~~ End Perpendicular Distan~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





! ~~~~~~~ Reumann-Witkam ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    function reumannWitkam(curve,tolerance) result(simpleCurve)
        real(dp),dimension(:,:),allocatable :: simpleCurve
        real(dp),dimension(:,:),intent(in) :: curve
        real(dp),intent(in) :: tolerance

        integer :: i, j, length, newLength
        logical :: not_last,not_key
        real(dp) :: squareTolerance, squareDistance
        real(dp),dimension(:),allocatable :: current, next_point, test_point, next_potential_key


        length = size(curve,dim=1)

        ! Make necessary checks

        if (tolerance <= 0._dp) then
            stop 'tolerance must be >0'
        end if

        if (length < 2) then
            stop 'curve is too short'
        end if

        squareTolerance = tolerance*tolerance

        ! allocate maximum amount of memory
        allocate(simpleCurve(length,size(curve,dim=2)))

        simpleCurve(1,:) = curve(1,:)   ! First element is always a key
        newLength = 1

        i = 1
        not_last = .true.

        main_loop: do while(not_last)
            i = i + 1
            if (i==length) then ! last point is always a key
                newLength = newLength + 1
                simpleCurve(newLength,:) = curve(i,:)
                exit main_loop
            end if
            current = simpleCurve(newLength,:)  ! current is most recent key
            next_point = curve(i,:)             ! next point is point next to most recent key
            next_potential_key = next_point
            not_key = .true.
            j = i + 1
            do while(not_key)   ! loop over points until one is outside of tolerance. key is the last within tolerance
                test_point = curve(j,:)
                squareDistance = d2LineToPoint(current,next_point,test_point)
                if (squareDistance >= squareTolerance) then     ! test point is outside of tolerance so previous point is key
                    newLength = newLength + 1
                    simpleCurve(newLength,:) = next_potential_key
                    not_key = .false.
                else if (j == length) then  ! last point is always a key
                    newLength = newLength + 1
                    simpleCurve(newLength,:) = curve(i,:)
                    exit main_loop
                else            ! test point still within tolerance so is upgraded to next potential key
                    next_potential_key = test_point
                    j = j + 1
                end if
            end do
            i = j-1

        end do main_loop


        simpleCurve = simpleCurve(1:newLength,:)    ! resize array

    end function reumannWitkam



! ~~~~~~~ End Reumann-Witkam ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



end module Simplify