p4n:
	gfortran module/array_mod.f90 p4_numerov.f90 -llapack
	./a.out
	gnuplot --persist plot/wavefun_p4.gp
	rm a.out
p1:
	gfortran module/array_mod.f90 numerov.f90 -llapack
	./a.out
	gnuplot --persist plot/wavefun.gp
	rm a.out
p3d:
	gfortran module/array_mod.f90 p3_diag.f90 -llapack
	./a.out
	gnuplot --persist plot/p3_diag.gp
	rm a.out
p3n:
	gfortran module/array_mod.f90 p3_numerov.f90 -llapack
	./a.out
	gnuplot --persist plot/p3_numerov.gp
	rm a.out
p4d:
	gfortran module/array_mod.f90 p4_diag.f90 -llapack
	./a.out
	gnuplot --persist plot/p4_diag.gp
	rm a.out
potential:
	gfortran module/array_mod.f90 p4_numerov.f90 -llapack
	./a.out
	gnuplot --persist plot/potential.gp
	rm a.out
