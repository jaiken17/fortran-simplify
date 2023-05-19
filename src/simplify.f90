!   Copyright (c) 2023, Joshua Aiken
!   All rights reserved.
!
!   This source code is licensed under the BSD-style license found in the
!   LICENSE file in the root directory of this source tree. 


module simplify
    use precision
    implicit none


    private
    public nth_point, nth_point_indices, radial_distance, radial_distance_indices, perpendicular_distance,  &
        perpendicular_distance_indices, reumann_witkam, reumann_witkam_indices


    interface nth_point
        procedure :: nth_point_multi
        procedure :: nth_point_single
    end interface nth_point

    interface nth_point_indices
        procedure :: nth_point_multi_indices
        procedure :: nth_point_single_indices
    end interface nth_point_indices

    interface radial_distance
        procedure :: radial_distance_multi
        procedure :: radial_distance_single
    end interface radial_distance

    interface radial_distance_indices
        procedure :: radial_distance_multi_indices
        procedure :: radial_distance_single_indices
    end interface radial_distance_indices

    interface perpendicular_distance
        procedure :: perp_distance
        procedure :: prependicular_distance_repeat
    end interface perpendicular_distance

    interface perpendicular_distance_indices
        procedure :: perp_distance_indices
        !procedure :: perp_distance_repeat_indices
    end interface perpendicular_distance_indices


    interface square_euclid_distance
        procedure :: square_euclid_distance_multi
        procedure :: square_euclid_distance_1d
    end interface square_euclid_distance


contains


! ~~~~~~~ Square Euclid Distnace ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! square_euclid_distance is a helper function that calculates the square 
    ! euclidean distance between any two, n-dimensional points (including 1D)

    function square_euclid_distance_multi(point1,point2) result(distance)
        real(dp) :: distance
        real(dp),dimension(:),intent(in) :: point1,point2

        integer :: length

        length = size(point1,dim=1)
        if (length /= size(point2,dim=1)) then
            stop 'length of vectors not equal'
        end if

        distance = dot_product(point1,point2)

    end function square_euclid_distance_multi

    function square_euclid_distance_1d(point1,point2) result(distance)
        real(dp) :: distance
        real(dp),intent(in) :: point1,point2
        
        real(dp),dimension(1) :: point1_vec, point2_vec
        
        point1_vec(1) = point1
        point2_vec(1) = point2

        distance = dot_product(point1_vec,point2_vec)

    end function square_euclid_distance_1d


! ~~~~~~~ End Square Euclid Distnace ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



! ~~~~~~~ Distance From Line to Point ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! d2LineToPoint is a helper function that computes the perpendicular 
    ! square distance from a line (defined by two points) to a point. This function
    ! works on n-dimensional points (including the trivial case of 1D).

    function d2_line_to_point(line_point1,line_point2,point) result(distance2)
        real(dp) :: distance2
        real(dp),dimension(:),intent(in) :: line_point1, line_point2, point

        real(dp),dimension(size(line_point1,dim=1)) :: v,x_minus_a,b_minus_a,v_minus_x
        real(dp) :: t
        integer :: length

        length = size(line_point1,dim=1)

        if ((length /= size(line_point2,dim=1)) .or. (length /= size(point,dim=1))) then
            stop 'line_point1, line_point2, point not all the same length'
        end if

        x_minus_a = point - line_point1
        b_minus_a = line_point2 - line_point1

        t = dot_product(x_minus_a,b_minus_a)/dot_product(b_minus_a,b_minus_a)
        v = line_point1 + t*b_minus_a
        
        v_minus_x = v - point

        distance2 = dot_product(v_minus_x,v_minus_x)

    end function d2_line_to_point

! ~~~~~~~ End Distance From Line to Point ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



! ~~~~~~~ Nth Point ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! nthPoint is a function that implements the nth point polyline 
    ! simplification algorithm. It works on n-dimensional curves (including 1D)
    ! and takes parameter "n" which is interpreted as every nth element to be 
    ! kept.

    function nth_point_multi(curve,n) result(simple_curve)

        real(dp),dimension(:,:),allocatable :: simple_curve
        real(dp),dimension(:,:),intent(in) :: curve
        integer,intent(in) :: n

        integer :: i,length
        integer :: counter,new_length

        length = size(curve,dim=1)
        
        ! Make necessary checks
        if (length < 2) then
            stop 'curve is too short'
        end if

        if (n < 1) then
            stop 'n is too small (<1)'
        end if

        ! Allocate maximum memory needed, will be trimmed at the end
        allocate(simple_curve(length,size(curve,dim=2)))

        simple_curve(1,:) = curve(1,:)   ! First element is always a key
        counter = 1; new_length = 1;
        do i=2,length-1
            counter = counter + 1
            if (counter == n) then
                new_length = new_length + 1
                simple_curve(new_length,:) = curve(i,:)
                counter = 0 ! reset counter
            end if
        end do
        new_length = new_length + 1
        simple_curve(new_length,:) = curve(length,:)  ! Last element is always a key
        
        simple_curve = simple_curve(1:new_length,:)    ! Resize final curve


    end function nth_point_multi

    function nth_point_single(curve,n) result(simple_curve)

        real(dp),dimension(:),allocatable :: simple_curve
        real(dp),dimension(:),intent(in) :: curve
        integer,intent(in) :: n

        integer :: i,length
        integer :: counter,new_length

        length = size(curve)

        ! Make necessary checks
        if (length < 2) then
            stop 'curve is too short'
        end if

        if (n < 1) then
            stop 'n is too small (<1)'
        end if

        ! Allocate maximum memory needed, will be trimmed at the end
        allocate(simple_curve(length))

        simple_curve(1) = curve(1)   ! First element is always a key
        counter = 1; new_length = 1;
        do i=2,length-1
            counter = counter + 1
            if (counter == n) then
                new_length = new_length + 1
                simple_curve(new_length) = curve(i)
                counter = 0 ! reset counter
            end if
        end do

        new_length = new_length + 1
        simple_curve(new_length) = curve(length)  ! Last element is always a key
        
        simple_curve = simple_curve(1:new_length)  ! Resize final curve


    end function nth_point_single



    ! Index returning versions

    function nth_point_multi_indices(curve,n) result(simple_curve_indices)

        integer,dimension(:),allocatable :: simple_curve_indices
        real(dp),dimension(:,:),intent(in) :: curve
        integer,intent(in) :: n

        integer :: i,length
        integer :: counter,new_length

        length = size(curve,dim=1)
        
        ! Make necessary checks
        if (length < 2) then
            stop 'curve is too short'
        end if

        if (n < 1) then
            stop 'n is too small (<1)'
        end if

        ! Allocate maximum memory needed, will be trimmed at the end
        allocate(simple_curve_indices(length))

        simple_curve_indices(1) = 1   ! First element is always a key
        counter = 1; new_length = 1;
        do i=2,length-1
            counter = counter + 1
            if (counter == n) then
                new_length = new_length + 1
                simple_curve_indices(new_length) = i    ! add element as key
                counter = 0 ! reset counter
            end if
        end do
        new_length = new_length + 1
        simple_curve_indices(new_length) = length  ! Last element is always a key
        
        simple_curve_indices = simple_curve_indices(1:new_length)    ! Resize final curve


    end function nth_point_multi_indices

    function nth_point_single_indices(curve,n) result(simple_curve_indices)

        real(dp),dimension(:),allocatable :: simple_curve_indices
        real(dp),dimension(:),intent(in) :: curve
        integer,intent(in) :: n

        integer :: i,length
        integer :: counter,new_length

        length = size(curve)

        ! Make necessary checks
        if (length < 2) then
            stop 'curve is too short'
        end if

        if (n < 1) then
            stop 'n is too small (<1)'
        end if

        ! Allocate maximum memory needed, will be trimmed at the end
        allocate(simple_curve_indices(length))

        simple_curve_indices(1) = 1   ! First element is always a key
        counter = 1; new_length = 1;
        do i=2,length-1
            counter = counter + 1
            if (counter == n) then
                new_length = new_length + 1
                simple_curve_indices(new_length) = i
                counter = 0 ! reset counter
            end if
        end do

        new_length = new_length + 1
        simple_curve_indices(new_length) = length  ! Last element is always a key
        
        simple_curve_indices = simple_curve_indices(1:new_length)  ! Resize final curve


    end function nth_point_single_indices



! ~~~~~~~ End Nth Point ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



! ~~~~~~~ Radial Distance ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! radial_distance is a function that implements the radial distance 
    ! polyline simplification algorithm. It works on n-dimensional curves
    ! (including 1D) and takes a parameter "tolerance" which is the minimum 
    ! distance allowed between any two points on the simplified curve (excluding
    ! the last two points).

    function radial_distance_multi(curve,tolerance) result(simple_curve)
        real(dp),dimension(:,:),allocatable :: simple_curve
        real(dp),dimension(:,:),intent(in) :: curve
        real(dp),intent(in) :: tolerance
        
        real(dp) :: square_tolerance
        integer :: i, length
        integer :: new_length
        logical :: not_last
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
        allocate(simple_curve(length,size(curve,dim=2)))

        square_tolerance = tolerance*tolerance

        simple_curve(1,:) = curve(1,:)   ! First element is always a key
        new_length = 1
        i = 1
        not_last = .true.
        do while (not_last)
            current = simple_curve(new_length,:)
            i = i + 1
            if (i==length) then
                exit
            end if
            distance = square_euclid_distance(current,curve(i,:))
            if (distance > square_tolerance) then
                new_length = new_length + 1
                simple_curve(new_length,:) = curve(i,:)
            end if
        end do

        new_length = new_length + 1
        simple_curve(new_length,:) = curve(length,:)  ! Last element is always a key

        simple_curve = simple_curve(1:new_length,:)    ! Resize final curve

    end function radial_distance_multi

    function radial_distance_single(curve,tolerance) result(simple_curve)
        real(dp),dimension(:),allocatable :: simple_curve
        real(dp),dimension(:),intent(in) :: curve
        real(dp),intent(in) :: tolerance
        
        real(dp) :: square_tolerance
        integer :: i, length
        integer :: new_length
        logical :: not_last
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

        allocate(simple_curve(length))

        square_tolerance = tolerance*tolerance

        simple_curve(1) = curve(1)   ! First element is always a key
        new_length = 1
        i = 1
        not_last = .true.
        do while (not_last)
            current = simple_curve(new_length)
            i = i + 1
            if (i==length) then
                exit
            end if
            distance = square_euclid_distance(current,curve(i))
            if (distance > square_tolerance) then
                new_length = new_length + 1
                simple_curve(new_length) = curve(i)
            end if
        end do

        new_length = new_length + 1
        simple_curve(new_length) = curve(length)  ! Last element is always a key

        simple_curve = simple_curve(1:new_length)  ! Resize final curve

    end function radial_distance_single



    ! Index returning versions


    function radial_distance_multi_indices(curve,tolerance) result(simple_curve_indices)
        integer,dimension(:),allocatable :: simple_curve_indices
        real(dp),dimension(:,:),intent(in) :: curve
        real(dp),intent(in) :: tolerance
        
        real(dp) :: square_tolerance
        integer :: i, length
        integer :: new_length
        logical :: not_last
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
        allocate(simple_curve_indices(length))

        square_tolerance = tolerance*tolerance

        simple_curve_indices(1) = 1     ! First element is always a key
        new_length = 1
        i = 1
        not_last = .true.
        do while (not_last)
            current = curve(simple_curve_indices(new_length),:)
            i = i + 1
            if (i==length) then
                exit
            end if
            distance = square_euclid_distance(current,curve(i,:))
            if (distance > square_tolerance) then
                new_length = new_length + 1
                simple_curve_indices(new_length) = i
            end if
        end do

        new_length = new_length + 1
        simple_curve_indices(new_length) = length  ! Last element is always a key

        simple_curve_indices = simple_curve_indices(1:new_length)    ! Resize final curve

    end function radial_distance_multi_indices

    function radial_distance_single_indices(curve,tolerance) result(simple_curve_indices)
        integer,dimension(:),allocatable :: simple_curve_indices
        real(dp),dimension(:),intent(in) :: curve
        real(dp),intent(in) :: tolerance
        
        real(dp) :: square_tolerance
        integer :: i, length
        integer :: new_length
        logical :: not_last
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

        allocate(simple_curve_indices(length))

        square_tolerance = tolerance*tolerance

        simple_curve_indices(1) = 1   ! First element is always a key
        new_length = 1
        i = 1
        not_last = .true.
        do while (not_last)
            current = curve(simple_curve_indices(new_length))
            i = i + 1
            if (i==length) then
                exit
            end if
            distance = square_euclid_distance(current,curve(i))
            if (distance > square_tolerance) then
                new_length = new_length + 1
                simple_curve_indices(new_length) = i
            end if
        end do

        new_length = new_length + 1
        simple_curve_indices(new_length) = length  ! Last element is always a key

        simple_curve_indices = simple_curve_indices(1:new_length)  ! Resize final curve

    end function radial_distance_single_indices


! ~~~~~~~ End Radial Distance ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




! ~~~~~~~ Perpendicular Distance ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! perpendicularDistance is a function that implement the perpendicular
    ! distance polyline simplification algorithm. It works on n>1-dimensional
    ! curves and takes a parameter "tolerance" which is the minimum 
    ! perpendicular distance a middle point can be from the line connecting 
    ! the two points on either side of that middle point. It also takes 
    ! a parameter "repeat" which is how many times to repeat the algorithm,
    ! since at most 50% of the curves points can be removed on a single pass. 


    function perp_distance(curve,tolerance) result(simple_curve)
        real(dp),dimension(:,:),allocatable :: simple_curve
        real(dp),dimension(:,:),intent(in) :: curve
        real(dp),intent(in) :: tolerance

        integer :: i, length, new_length
        logical :: last_is_key = .false., skip_iter = .false.
        real(dp) :: square_distance, square_tolerance
        

        length = size(curve,dim=1)

        ! Make necessary checks

        if (tolerance <= 0._dp) then
            stop 'tolerance must be >0'
        end if

        if (size(curve,dim=1) < 2) then
            stop 'curve is too short'
        end if

        square_tolerance = tolerance*tolerance

        ! Allocate absolute maximum amount of memory, will get trimmed at end

        allocate(simple_curve(length,size(curve,dim=2)))

        simple_curve(1,:) = curve(1,:)   ! First element is always a key
        new_length = 1

        do i=2,length-1
            if (skip_iter) then
                skip_iter = .false.
                cycle
            end if

            square_distance = d2_line_to_point(simple_curve(new_length,:),curve(i+1,:),curve(i,:))
            if (square_distance >= square_tolerance) then
                new_length = new_length + 1
                simple_curve(new_length,:) = curve(i,:)
            else
                new_length = new_length + 1
                simple_curve(new_length,:) = curve(i+1,:)
                skip_iter = .true.
                if (i+1 == length) last_is_key = .true.
            end if
        end do

        if (.not. last_is_key) then   ! last element may already be set as key
            new_length = new_length + 1
            simple_curve(new_length,:) = curve(length,:)  ! Last element is always a key
        end if

        simple_curve = simple_curve(1:new_length,:)    ! Resize simple_curve to only needed elements

    end function perp_distance


    function prependicular_distance_repeat(curve,tolerance,repeat) result(simple_curve)
        real(dp),dimension(:,:),allocatable :: simple_curve
        real(dp),dimension(:,:),intent(in) :: curve
        real(dp),intent(in) :: tolerance
        integer,intent(in) :: repeat

        integer :: i
        
        if (repeat < 1) then
            stop 'number of repeats too small'
        end if

        simple_curve = curve

        do i=1,repeat
            simple_curve = perp_distance(simple_curve,tolerance)
        end do


    end function prependicular_distance_repeat



    ! Index returning versions


function perp_distance_indices(curve,tolerance) result(simple_curve_indices)
        integer,dimension(:),allocatable :: simple_curve_indices
        real(dp),dimension(:,:),intent(in) :: curve
        real(dp),intent(in) :: tolerance

        integer :: i, length, new_length
        logical :: last_is_key = .false., skip_iter = .false.
        real(dp) :: square_distance, square_tolerance
        

        length = size(curve,dim=1)

        ! Make necessary checks

        if (tolerance <= 0._dp) then
            stop 'tolerance must be >0'
        end if

        if (size(curve,dim=1) < 2) then
            stop 'curve is too short'
        end if

        square_tolerance = tolerance*tolerance

        ! Allocate absolute maximum amount of memory, will get trimmed at end

        allocate(simple_curve_indices(length))

        simple_curve_indices(1) = 1   ! First element is always a key
        new_length = 1

        do i=2,length-1
            if (skip_iter) then
                skip_iter = .false.
                cycle
            end if

            square_distance = d2_line_to_point(curve(simple_curve_indices(new_length),:),curve(i+1,:),curve(i,:))
            if (square_distance >= square_tolerance) then
                new_length = new_length + 1
                simple_curve_indices(new_length) = i
            else
                new_length = new_length + 1
                simple_curve_indices(new_length) = i+1
                skip_iter = .true.
                if (i+1 == length) last_is_key = .true.
            end if
        end do

        if (.not. last_is_key) then   ! last element may already be set as key
            new_length = new_length + 1
            simple_curve_indices(new_length) = length  ! Last element is always a key
        end if

        simple_curve_indices = simple_curve_indices(1:new_length)    ! Resize simple_curve to only needed elements

    end function perp_distance_indices


    ! Converting this procedure to return indices from original curve is non-trivial
    ! For now, just call perp_distance_indices multiple times.

    ! function prependicular_distance_repeat_indices(curve,tolerance,repeat) result(simple_curve_indices)
    !     real(dp),dimension(:,:),allocatable :: simple_curve_indices
    !     real(dp),dimension(:,:),intent(in) :: curve
    !     real(dp),intent(in) :: tolerance
    !     integer,intent(in) :: repeat

    !     integer :: i
        
    !     if (repeat < 1) then
    !         stop 'number of repeats too small'
    !     end if

    !     simple_curve = curve

    !     do i=1,repeat
    !         simple_curve = perp_distance(simple_curve,tolerance)
    !     end do


    ! end function prependicular_distance_repeat_indices




! ~~~~~~~ End Perpendicular Distan~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





! ~~~~~~~ Reumann-Witkam ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! rumann_witkam is a function that implements the Reumann-Witkam polyline 
    ! simplification algorithm. It works on n-dimensional curves (n>2)
    ! and takes parameter "tolerance" which is interpreted as how close points
    ! can be to a line defined by the current key and next element of the
    ! curve. Smaller values of tolerance result in more points being included as
    ! keys in the simplified curve.


    function reumann_witkam(curve,tolerance) result(simple_curve)
        real(dp),dimension(:,:),allocatable :: simple_curve
        real(dp),dimension(:,:),intent(in) :: curve
        real(dp),intent(in) :: tolerance

        integer :: i, j, length, new_length
        logical :: not_last,not_key
        real(dp) :: square_tolerance, square_distance
        real(dp),dimension(:),allocatable :: current, next_point, test_point, next_potential_key


        length = size(curve,dim=1)

        ! Make necessary checks

        if (tolerance <= 0._dp) then
            stop 'tolerance must be >0'
        end if

        if (length < 2) then
            stop 'curve is too short'
        end if

        square_tolerance = tolerance*tolerance

        ! allocate maximum amount of memory
        allocate(simple_curve(length,size(curve,dim=2)))

        simple_curve(1,:) = curve(1,:)   ! First element is always a key
        new_length = 1

        i = 1
        not_last = .true.

        main_loop: do while(not_last)
            i = i + 1
            if (i==length) then ! last point is always a key
                new_length = new_length + 1
                simple_curve(new_length,:) = curve(i,:)
                exit main_loop
            end if
            current = simple_curve(new_length,:)  ! current is most recent key
            next_point = curve(i,:)             ! next point is point next to most recent key
            next_potential_key = next_point
            not_key = .true.
            j = i + 1
            do while(not_key)   ! loop over points until one is outside of tolerance. key is the last within tolerance
                test_point = curve(j,:)
                square_distance = d2_line_to_point(current,next_point,test_point)
                if (square_distance >= square_tolerance) then     ! test point is outside of tolerance so previous point is key
                    new_length = new_length + 1
                    simple_curve(new_length,:) = next_potential_key
                    not_key = .false.
                else if (j == length) then  ! last point is always a key
                    new_length = new_length + 1
                    simple_curve(new_length,:) = curve(j,:)
                    exit main_loop
                else            ! test point still within tolerance so is upgraded to next potential key
                    next_potential_key = test_point
                    j = j + 1
                end if
            end do
            i = j-1

        end do main_loop


        simple_curve = simple_curve(1:new_length,:)    ! resize array

    end function reumann_witkam



    ! Index returning version

    function reumann_witkam_indices(curve,tolerance) result(simple_curve_indices)
        integer,dimension(:),allocatable :: simple_curve_indices
        real(dp),dimension(:,:),intent(in) :: curve
        real(dp),intent(in) :: tolerance

        integer :: i, j, length, new_length
        logical :: not_last,not_key
        real(dp) :: square_tolerance, square_distance
        real(dp),dimension(:),allocatable :: current, next_point, test_point, next_potential_key
        integer :: next_potential_key_index


        length = size(curve,dim=1)

        ! Make necessary checks

        if (tolerance <= 0._dp) then
            stop 'tolerance must be >0'
        end if

        if (length < 2) then
            stop 'curve is too short'
        end if

        square_tolerance = tolerance*tolerance

        ! allocate maximum amount of memory
        allocate(simple_curve_indices(length))

        simple_curve_indices(1) = 1   ! First element is always a key
        new_length = 1

        i = 1
        not_last = .true.

        main_loop: do while(not_last)
            i = i + 1
            if (i==length) then ! last point is always a key
                new_length = new_length + 1
                simple_curve_indices(new_length) = i
                exit main_loop
            end if
            current = curve(simple_curve_indices(new_length),:)  ! current is most recent key
            next_point = curve(i,:)             ! next point is point next to most recent key
            next_potential_key = next_point
            next_potential_key_index = i
            not_key = .true.
            j = i + 1
            do while(not_key)   ! loop over points until one is outside of tolerance. key is the last within tolerance
                test_point = curve(j,:)
                square_distance = d2_line_to_point(current,next_point,test_point)
                if (square_distance >= square_tolerance) then     ! test point is outside of tolerance so previous point is key
                    new_length = new_length + 1
                    simple_curve_indices(new_length) = next_potential_key_index
                    not_key = .false.
                else if (j == length) then  ! last point is always a key
                    new_length = new_length + 1
                    simple_curve_indices(new_length) = j
                    exit main_loop
                else            ! test point still within tolerance so is upgraded to next potential key
                    next_potential_key = test_point
                    next_potential_key_index = j
                    j = j + 1
                end if
            end do
            i = j-1

        end do main_loop


        simple_curve_indices = simple_curve_indices(1:new_length)    ! resize array

    end function reumann_witkam_indices



! ~~~~~~~ End Reumann-Witkam ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



end module simplify