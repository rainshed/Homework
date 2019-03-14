p4:
	gfortran module/array_mod.o p4_numerov.f90 -llapack
	./a.out
	gnuplot plot/wavefun_p4.gp
	rm a.out
p1:
	gfortran module/array_mod.o numerov.f90 -llapack
	./a.out
	gnuplot plot/wavefun.gp
	rm a.out
p3:
	gfortran module/array_mod.o p3_diag.f90 -llapack
	./a.out
	gnuplot --persist plot/p3_diag.gp
	rm a.out
