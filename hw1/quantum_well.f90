program main
	use array
	implicit none
	integer,parameter :: Nj = 512
	real(8),parameter :: pi = acos(-1d0)
	real(8),parameter :: x0 = 5.0d0,k0=7.0d0
	complex(8),parameter :: xi = complex(0,1)
	real(8) :: Xmax = 10.d0, Delta_t = 0.02d0, Delta_x
	complex(8) :: M1(Nj,Nj),M2(Nj,NJ),inv(Nj,Nj)
	integer :: i,j
	complex(8) :: x(NJ), psi_old(Nj),psi_new(Nj)
	real(8) :: y(Nj),norm
	
	Delta_x = Xmax/Nj

	x = vecn(0d0,Xmax,Nj)
	norm = 0
!	open(unit = 10 ,file = 'psi.dat',status='unknown')
	do i =1,Nj
		psi_old(i) = sqrt(1/(2*pi)) * exp(xi*k0*x(i)-((x(i)-x0)**2)/2)
!		y(i) = sqrt(real(psi_old(i))**2 + aimag(psi_old(i))**2) 
!		norm = norm + y(i)*Delta_x 
!		write(10,*) real(x(i)),y(i)
	end do
!	write(*,*) norm
!
!	stop


	M1 = complex(0,0)
	M2 = complex(0,0)

	do i = 1,Nj
		do j = 1,Nj
			if (i==j) then
				M1(i,j) = 1d0 + xi*Delta_t/Delta_x**2
				M2(i,j) = 1d0 - xi*Delta_t/Delta_x**2
			end if

			if (i==j+1 .or. i==j-1) then 
				M1(i,j) = -xi * Delta_t/(2d0*Delta_x**2) 
				M2(i,j) = -M1(i,j)
			end if
			
		end do
	end do

	
	
	inv = inverse(M1,Nj)
	M1 = matmul(inv,M2)
	do j=1,10
	psi_old = matmul(M1,psi_old)
	
	open(unit=11, file='final_results.dat',status = 'unknown')
	do i = 1,Nj
		y(i) = sqrt(real(psi_old(i))**2 + aimag(psi_old(i))**2) 
		write(11,*) dble(x(i)),y(i)
	end do
	end do 
	

end
