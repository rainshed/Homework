!the program is to calculate question 1 of exercise 1 of computational physics 
program main
	use array	
	implicit none	
	integer,parameter :: x_num = 300
	integer :: xnum,i,j=11
	real(8),parameter :: hbar=1,m=1
	real(8),parameter :: x_min=0 , x_max=2  ! which is need to set
	real(8) :: x(x_num),v(x_num),g(x_num),f(x_num),y(x_num),y_old(x_num)
	real(8) :: E=0, ytotal
	real(8) :: dx = (x_max-x_min)/(x_num-1),dE=0.0001

	open(unit=10,file='./data/p1_energy_value.dat',status='unknown')
	open(unit=11,file='./data/p1_wave_function1.dat',status='unknown')
	open(unit=12,file='./data/p1_wave_function2.dat',status='unknown')
	open(unit=13,file='./data/p1_wave_function3.dat',status='unknown')

	do i = 1,150000 

		x = vecn(x_min,x_max,x_num)
		y = 0
		y(1) = 0
		y(2) = 0.01
		ytotal = y(2)*dx 
		do xnum=1,x_num
			v(xnum)=0d0  !which can be set
			g(xnum)=(2*m/hbar**2)*(E-v(xnum))
			f(xnum)=1d0+1d0/12d0*g(xnum)*dx**2
			if (xnum < 3) cycle
			y(xnum)=((12-10*f(xnum-1))*y(xnum-1)-f(xnum-2)*y(xnum-2))/f(xnum)
			ytotal = ytotal + y(xnum)**2*dx
		end do
		
		do xnum=1,x_num
			y(xnum) = y(xnum)/sqrt(ytotal)
		end do
	!	y_old(1) = 0
	!	do xnum=1,x_num
	!		y_old(1) = y_old(1)+y(xnum)*dx
	!	end do
	!	write(*,*) y_old(1)
		
		if (y_old(x_num)*y(x_num)<0.0d0) then
			write(*,*) E,'is the',j-11,'th energy' 
			write(10,*) E,'is the energy'
			do xnum=1,x_num
				write(j,*) x(xnum),y(xnum)
			end do
			E = E +0.2
			j = j+1
		end if
	
		E = E+dE
		y_old = y

	end do
	close(10)
	close(11)


end program
