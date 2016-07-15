module util

contains 

    real* 8 function mean(arr)
        real *8 arr(:)

        mean = sum(arr)/(max(1,size(arr)))
    end function

    subroutine cumsum(x, y)
      real(8), intent(in)  :: x(:)
      real(8), intent(out) :: y(:)

      integer i
      y(1) = x(1)
      do i = 2,size(x)
         y(i) = y(i-1) + x(i)
      end do
    end subroutine cumsum

    integer function searchsorted(x, a)
      real(8) :: x(:), a

      integer i

      do i=1,size(x)
         if (a <= x(i)) then
            searchsorted = i
            return
         end if
      end do

      ! this only is evaluated if a >= then all elements of x
      searchsorted = -1
    end function searchsorted
end module
