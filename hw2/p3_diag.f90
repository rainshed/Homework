!this program is to solve Hamitonian equation by exact diagonalization
program main
	use array
	implicit none
	integer,parameter :: x_num=201,base_num=100
	real(8),parameter :: x_min=-1, x_max=1,a=1 !a = x_max-x_min
	real(8),parameter :: pi = acos(-1d0)
	real(8) :: x(x_num),v(x_num),psi(x_num,base_num),base(base_num,x_num),H(base_num,base_num)
	real(8) :: dx = (x_max-x_min)/(x_num-1),V0
	integer :: xnum,n,i,j

	!-----parameters for DSYEV--------------------
	real(8) :: eigenvalue(base_num),Work(1000)
	integer :: Lwork = 1000,info



	x = vecn(x_min,x_max,x_num)
	v = 0d0
	do xnum=1,x_num
		!define v(x)
		if (abs(x(xnum))<0.5d0) then
			v(xnum)=10 !which is dependent on different question
		end if
		!define base function
		do n = 1,base_num
			if (mod(n,2)==0) then
				base(n,xnum) = sqrt(1d0/a)*sin(n*pi*x(xnum)/(2*a))
			else
				base(n,xnum) = sqrt(1d0/a)*cos(n*pi*x(xnum)/(2*a))
			end if
		end do
	end do

	H =0	
	do i=1,base_num
		do j=1,base_num
			call Simpson(x_num,base(i,:),base(j,:),v,dx,V0)
			H(i,j) = V0
		!	write(*,*) V0
			if (i==j) then
				H(i,j) = (dble(i)*pi)**2/(8d0*a**2) + V0
			!	write(*,*) H(i,j)
			end if
			H(j,i)=H(i,j)
		end do
	end do
!	write(*,*) H

	call DSYEV('V','U',base_num,H,base_num,eigenvalue,Work,Lwork,info)
	if (info /= 0) then
		write(*,*) 'DSYEV error'
		stop
	end if

	write(*,*) eigenvalue(1)
	open(unit=11,file='data/p3_diag_eigenvalue.dat',status='unknown')
	write(11,*) eigenvalue
	close(11)
	
	psi = matmul(transpose(base),H)
	open(unit=10,file='data/p3_diag.dat',status='unknown')
	do i = 1,x_num
		write(10,*) x(i),(psi(i,1))**2
	end do
	close(10)

	
end program

subroutine Simpson(N,u,v,Vpot,dx,V0)
	implicit none
	integer :: N
	real(8),intent(in) :: u(N),v(N),Vpot(N),dx
	real(8),intent(out) :: V0

	!--- local vars --------
	integer :: i
	real(8) :: Work(N) 

	if (mod(N,2) == 0) then
		write(*,*) 'Array with even elements, Simpson does not know how to work'
		stop
	end if 

	do i=1,N
		Work(i) = Vpot(i)*u(i)*v(i)
	end do

	V0 = Work(1)+Work(n)
	do i=2,N-1,2
		V0 = V0+4.d0*Work(i)
	end do

	do i=3,N-2,2
		V0 = V0 + 2.d0*Work(i)
	end do

	V0 = V0*dx/3d0

end subroutine
