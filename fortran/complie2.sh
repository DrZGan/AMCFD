./clean.sh
gfortran -c -mcmodel=large mod_const.f90 mod_param.f90 mod_geom.f90 mod_init.f90 mod_laser.f90 mod_dimen.f90 mod_bound.f90 mod_discret.f90 mod_entot.f90 mod_sour.f90 mod_flux.f90 mod_prop.f90 mod_resid.f90 mod_revise.f90 mod_solve.f90 mod_print.f90 mod_converge.f90 mod_toolpath.f90 
gfortran -c -mcmodel=large main.f90

gfortran mod_const.o mod_param.o mod_geom.o mod_init.o mod_laser.o mod_dimen.o mod_bound.o mod_discret.o mod_entot.o mod_sour.o mod_flux.o mod_prop.o mod_resid.o mod_revise.o mod_solve.o mod_print.o mod_converge.o mod_toolpath.o main.o -o  cluster_main




