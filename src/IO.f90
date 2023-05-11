module IO
    use Precision
    implicit none

    ! overload subroutine
    interface outputMatrixWithHeaders
        module procedure :: outputMatrixWithHeadersTerminal
        module procedure :: outputMatrixWithHeadersFile
        module procedure :: outputVectorWithHeadersTerminal
        module procedure :: outputVectorWithHeadersFile
    end interface outputMatrixWithHeaders


contains


    subroutine readDataFile(filename, unit, num_rows, num_cols, data)
        character(*),intent(in) :: filename  ! name of data file
        integer,intent(in) :: unit, num_cols, num_rows
        real(dp),dimension(:, :),allocatable,intent(out) :: data  ! matrix to store data
    
        integer :: i, j
        real(dp) :: temp

        open(unit=unit, file=filename, status='old')
    
        ! allocate memory for data matrix
        allocate(data(num_rows, num_cols))
    
        ! read data into matrix
        do i = 1, num_rows
            read(unit,*) data(i,1:num_cols)
            ! do j = 1, num_cols
                ! read(unit, *) temp
                ! data(i, j) = temp
            ! end do
        end do
    
        close(unit)

    end subroutine readDataFile


    subroutine outputMatrix(matrix,reverse_rows,reverse_cols)
        real(dp), dimension(:, :), intent(in) :: matrix  ! matrix to output
        logical,intent(in),optional :: reverse_rows, reverse_cols
        integer :: num_rows, num_cols         ! number of rows and columns in matrix
        integer :: i, j, order_row, order_col
        integer :: row_start, row_end, col_start, col_end
      
        num_rows = size(matrix,dim=1)
        num_cols = size(matrix,dim=2)

        ! initialize default output order
        order_row = 1
        row_start = 1
        row_end = num_rows

        order_col = 1
        col_start = 1
        col_end = num_cols

        ! check if output order needs reversing
        if (present(reverse_rows)) then
            if (reverse_rows) then
                order_row = -1
                row_start = num_rows
                row_end = 1
            end if 
        end if

        if (present(reverse_cols)) then
            if (reverse_cols) then
                order_col = -1
                col_start = num_cols
                col_end = 1
            end if
        end if


        ! print matrix rows with row numbers
        do i = row_start, row_end, order_row
            write(*, '(i5, 1x)', advance='no') i
            do j = col_start, col_end, order_col
                write(*, '(f8.3, 1x)', advance='no') matrix(i, j)
            end do
            write(*, *) ! print newline after each row
        end do

        ! print column numbers
        write(*,'(a5,1x)',advance='no') "     "
        do j = col_start, col_end, order_col
            write(*, '(i8,1x)', advance='no') j
        end do
        write(*,*) ! newline

    end subroutine outputMatrix
      






! OUTPUT MATRIX WITH HEADERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    subroutine outputMatrixWithHeadersTerminal(matrix, headers)
        implicit none
        real(dp), dimension(:, :), intent(in) :: matrix     ! matrix to output
        character(len=*), dimension(:), intent(in) :: headers ! column headers
      
        
        integer :: num_rows, num_cols            ! number of rows and columns in matrix
        integer :: i, j
      
        num_rows = size(matrix,dim=1)
        num_cols = size(matrix,dim=2)

        if (num_cols /= size(headers)) then
            stop 'number of columns does not match number of headers'
        end if

        ! print column headers
        print *, "   "
        do j = 1, num_cols
            write(*, '(a19,1x)', advance='no') trim(headers(j))
        end do
        print *
      
        ! print matrix rows
        do i = 1, num_rows
            do j = 1, num_cols
                write(*, '(f19.10, 1x)', advance='no') matrix(i, j)
            end do
            write(*, *) ! print newline after each row
        end do
      end subroutine outputMatrixWithHeadersTerminal
      
      subroutine outputMatrixWithHeadersFile(matrix, headers, filename, unit)
        implicit none
        real(dp), dimension(:, :), intent(in) :: matrix     ! matrix to output
        character(len=*), dimension(:), intent(in) :: headers ! column headers
        character(*),intent(in) :: filename
        integer,intent(in) :: unit


        integer :: num_rows, num_cols            ! number of rows and columns in matrix
        integer :: i, j
      

        num_rows = size(matrix,dim=1)
        num_cols = size(matrix,dim=2)

        if (num_cols /= size(headers)) then
            stop 'number of columns does not match number of headers'
        end if

        open(unit=unit, file=filename, status='replace')


        ! print column headers
        do j = 1, num_cols
            write(unit, '(a19,1x)', advance='no') trim(headers(j))
        end do
        write(unit,*) ! newline
      
        ! print matrix rows
        do i = 1, num_rows
            do j = 1, num_cols
                write(unit, '(f19.10, 1x)', advance='no') matrix(i, j)
            end do
            write(unit, *) ! print newline after each row
        end do

        close(unit)

      end subroutine outputMatrixWithHeadersFile


      subroutine outputVectorWithHeadersTerminal(vector, headers)
        implicit none
        real(dp), dimension(:), intent(in) :: vector     ! vector to output
        character(len=*), intent(in) :: headers ! column headers
      
        
        integer :: num_rows            ! number of rows in vector
        integer :: i
      
        num_rows = size(vector,dim=1)


        ! print column headers
        print *, "   "
        write(*, '(a19,1x)') trim(headers)
      
        ! print vector rows
        do i = 1, num_rows
            write(*, '(f19.10, 1x)') vector(i)
        end do
      end subroutine outputVectorWithHeadersTerminal
      
      subroutine outputVectorWithHeadersFile(vector, headers, filename, unit)
        implicit none
        real(dp), dimension(:), intent(in) :: vector     ! vector to output
        character(len=*), intent(in) :: headers ! column headers
        character(*),intent(in) :: filename
        integer,intent(in) :: unit


        integer :: num_rows            ! number of rows in vector
        integer :: i
      
        open(unit=unit, file=filename, status='replace')


        num_rows = size(vector,dim=1)

        ! print column headers
        write(unit, '(a19,1x)') trim(headers)
      
        ! print vector rows
        do i = 1, num_rows
            write(unit, '(f19.10, 1x)') vector(i)
        end do

        close(unit)

      end subroutine outputVectorWithHeadersFile




end module IO