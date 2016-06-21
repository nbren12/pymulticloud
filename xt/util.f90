module util

contains 

    real* 8 function mean(arr)
        real *8 arr(:)

        mean = sum(arr)/(max(1,size(arr)))
    end function
end module
