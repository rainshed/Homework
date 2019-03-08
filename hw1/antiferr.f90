!operation with array
!-----------------------------------------------------------------------------------------------
module array
	!real,external :: vec
	integer,parameter :: dp = selected_real_kind(8)
	contains

!create a array with given initial value, end value and interval
!---------------------------------------------------------------------
	function vec(init,ed,inter) result(arr)
		implicit none
		integer :: i
		real(dp) :: init,ed,inter
		integer :: num
		real(dp) :: arr(int((ed-init)/inter+1))
		num = (ed-init)/inter+1
		do i=1,num
			arr(i) = init+(i-1)*inter
		end do
	end function
!-------------------------------------------------------------------------

!create a array with given initial value,end value and length of array
!-------------------------------------------------------------------------
	function vecn(init,ed,num) result(arr)
		implicit none
		integer :: i,num
		real(dp) :: init,ed,inter,arr(num)
		inter = (ed-init)/(real(num)-1d0)
		do i = 1,num
		  arr(i) = init + (i-1)*inter
		end do
	end function vecn
end module array
!---------------------------------------------------------------------------------------------



!------------------------------------------------------------------------------------------------

program main
	use array
	write(*,*) vecn(1d0,10d0,10)
end 

	
