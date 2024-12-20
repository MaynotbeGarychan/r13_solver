      function array_sum(num,array) result(val)
        !============================================================
        ! sum the array
        !------------------------------------------------------------
        ! input: 
        ! array(num)      - array
        ! num             - length of the array
        ! output:
        ! val             - sum of the array
        !------------------------------------------------------------
        implicit none
        integer i
        integer num
        double precision array(num)
        double precision val

        val=0.
        do i=1,num
          val=val+array(i)
        enddo
        
      end function array_sum

      function array_maxidx(arr,num) result(maxidx)
      !============================================================
      ! Return the index of a maximum value inside an array
      !------------------------------------------------------------
      ! input: 
      ! arr(num) - a given array
      ! num      - dimension of an array
      ! output:
      ! maxidx   - integer of the index of maximum value
      !------------------------------------------------------------
      implicit none
      integer i,num,maxidx
      double precision arr(num)
      double precision maxval

      maxval=arr(1)
      maxidx=1

      do i=2,num
          if(arr(i)>maxval)then
              maxval=arr(i)
              maxidx=i
          endif
      enddo

      end function array_maxidx

      subroutine array_linear(bgn,end,num,arr)
      !============================================================
      ! Return the index of a maximum value inside an array
      !------------------------------------------------------------
      ! input: 
      ! bgn, end - begin and end of the array
      ! num      - dimension of an array
      ! output:
      ! arr(num)   - array that returned
      !------------------------------------------------------------

        implicit none
        double precision bgn,end
        integer num
        double precision itv
        integer i
        double precision arr(num)

        itv=(end-bgn)/(num-1)

        do i=1,num
          arr(i)=bgn+(i-1)*itv
        enddo

      end subroutine array_linear

      subroutine array_fill(val,num,arr)
      !============================================================
      ! return an array filled with specific value
      !------------------------------------------------------------
      ! input: 
      ! val      - filled value
      ! num      - dimension of an array
      ! output:
      ! arr(num)   - array that returned
      !------------------------------------------------------------
        implicit none
        integer num
        double precision val
        double precision arr(num)
        integer i

        do i=1,num
          arr(i)=val
        enddo

      end subroutine array_fill
