#CC				= gcc
#FC				= gfortran
#LD				= g++
CC				= /opt/intel/bin/icc
FC				= /opt/intel/bin/ifort -nofor_main
IFC				= /opt/intel/bin/ifort
LD				= /opt/intel/bin/icpc
CFLAGS				= -g -O2
BINDIR				= ../bin
#MINUITINC			= -I/usr/local/c-minuit/c-minuit-CVS-20031201
#MINUITLIB			= -L/usr/local/c-minuit/c-minuit-CVS-20031201/minuit/code -lminuit
MINUITINC			= -I/usr/local/c-minuit/c-minuit-CVS-20031201-intel
MINUITLIB			= -L/usr/local/c-minuit/c-minuit-CVS-20031201-intel/minuit/code -lminuit
MINUITDEF			= -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DPACKAGE=\"c-minuit\" -DVERSION=\"CVS-20031201\" -DHAVE_LIBM=1 -Df2cFortran=1 -DSTDC_HEADERS=1 -DTIME_WITH_SYS_TIME=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_UNISTD_H=1 -DHAVE_GETOPT_H=1 -DHAVE_GETLINE=1
ROOTCLIB0S			:= $(shell root-config --cflags)
ROOTLIBS			:= $(shell root-config --libs) -lMinuit
ROOTGLIBS			:= $(shell root-config --glibs)
ROOTINC				:= $(shell root-config --incdir)
INCDIR				= -I$(shell root-config --incdir)
LIBDIR				= -L$(shell root-config --libdir)
INTEL_LIBDIR			= -Wl,-rpath-link,/opt/intel/composer_xe_2015.3.187/compiler/lib/intel64
#F2CLIB				= -L/usr/local/f2c/lib -lf2c -u MAIN__
#F2CLIB				= /usr/lib/libf2c.a
F2CLIB				= -lf2c -u MAIN__
F2CINC				= -I/usr/local/f2c/include
LIB0				= -lm
LIB1				= -lcfitsio
LIB2				= -ljpeg
LIB3				= -lgsl -lgslcblas
LIB4				= -lifcore
WARNING				= -Wall
FC_WARNING			= -warn all
OBJ1				= $(BINDIR)/strutil.o
OBJ2				= $(BINDIR)/numutil.o
OBJ3				= $(BINDIR)/tutil.o
OBJ4				= $(BINDIR)/mie2new.o
OBJ5				= $(BINDIR)/strutil.o $(BINDIR)/MIEV0.o $(BINDIR)/ErrPack.o $(BINDIR)/aerprf.o $(BINDIR)/prfdta.o
OBJ6				= $(BINDIR)/aeros_bb.o $(BINDIR)/airmie.o $(BINDIR)/append_f.o $(BINDIR)/findpar.o $(BINDIR)/legendre_coeff.o $(BINDIR)/r_min_max_comput.o $(BINDIR)/size.o
OBJ7				= $(BINDIR)/tmx.o $(BINDIR)/lpd.o
OBJ8				= $(BINDIR)/tmq.o $(BINDIR)/lpq.o
OBJ8L				= $(BINDIR)/tmq_large.o $(BINDIR)/lpq.o
OBJ9				= $(BINDIR)/newrad.o
OBJS				= $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6) $(OBJ7) $(OBJ8) $(OBJ8L) $(OBJ9)
OBJM				= $(BINDIR)/MIEV0.o $(BINDIR)/ErrPack.o
PROGS				= common aeros_db aeros_dls aeros_fit aeros_mix aeros_opt aeros_param aeros_mie2new aeros_miev0 aeros_mixcomp aeros_simul aeros_table aeros_tmx aeros_tmq mie2new.exe dls mie mie_size_dist mie_eff2mod ms720_correct_fov ms720_correct_tmp ms720_fov ms720_skyrad phase_integ read_mie read_tmq resolution_effect tmx_size_dist db
#-----------------------------------------------------------------------------------------
all				: clean $(PROGS)
				rm -f $(OBJS)
clean				:
				rm -f $(OBJS)
common				: common.c $(OBJ1)
				$(CC) $@.c $(OBJ1) -o $(BINDIR)/$@ $(LIB0) $(WARNING)
aeros_db			: aeros_common.c aeros_db.c $(OBJ1)
				$(CC) $@.c $(OBJ1) -o $(BINDIR)/$@ $(LIB0) $(WARNING)
aeros_dls			: aeros_common.c aeros_dls.c $(OBJ1)
				$(CC) $@.c $(OBJ1) -o $(BINDIR)/$@ $(LIB0) $(WARNING)
aeros_mixcomp			: aeros_common.c aeros_mixcomp.c $(OBJ1)
				$(CC) $@.c $(OBJ1) -o $(BINDIR)/$@ $(LIB0) $(WARNING)
aeros_mie2new			: aeros_common.c aeros_mie2new.c $(OBJ1) $(OBJ9)
				$(CC) $@.c $(OBJ1) $(OBJ9) -o $(BINDIR)/$@ $(LIB0) $(WARNING)
aeros_miev0			: aeros_common.c aeros_miev0.c $(OBJ5)
				$(CC) -c $@.c -o $(BINDIR)/$@.o $(WARNING)
				$(FC) $(BINDIR)/$@.o $(OBJ5) -o $(BINDIR)/$@ $(FC_WARNING)
				rm -f $(BINDIR)/$@.o
test				: test.c
				$(LD) $@.c -o $(BINDIR)/$@ $(INCDIR) -I$(ROOTINC) $(LIBDIR) $(ROOTLIBS) $(LIB0) $(OPT) $(WARNING)
#aeros_fit			: aeros_fit.c $(OBJ5) $(OBJ7)
#				$(CC) -c $@.c -o $(BINDIR)/$@.o $(MINUITINC) $(WARNING)
#				$(FC) $(MINUITDEF) $(BINDIR)/$@.o $(OBJ5) $(OBJ7) $(MINUITLIB) -o $(BINDIR)/$@ $(MINUITINC) $(F2CINC) $(LIBDIR) $(F2CLIB) $(LIB0) $(CFLAGS) $(FC_WARNING) -mcmodel=medium
#				rm -f $(BINDIR)/$@.o
aeros_fit			: aeros_fit.c $(OBJ5) $(OBJ7)
				$(CC) $(MINUITDEF) $@.c $(OBJ5) $(OBJ7) $(MINUITLIB) -o $(BINDIR)/$@ $(MINUITINC) $(F2CINC) $(LIBDIR) $(F2CLIB) $(LIB0) $(LIB4) $(CFLAGS) $(WARNING) -mcmodel=medium
aeros_mie			: aeros_mie.c $(OBJ5) $(OBJ7)
				$(CC) -c $@.c -o $(BINDIR)/$@.o $(MINUITINC) $(WARNING)
				$(FC) $(BINDIR)/$@.o $(OBJ5) $(OBJ7) -o $(BINDIR)/$@ $(LIBDIR) $(LIB0) $(CFLAGS) $(FC_WARNING)
				rm -f $(BINDIR)/$@.o
#aeros_mix			: aeros_mix.c $(OBJ5) $(OBJ7)
#				$(CC) -c $@.c -o $(BINDIR)/$@.o $(MINUITINC) $(WARNING)
#				$(FC) $(MINUITDEF) $(BINDIR)/$@.o $(OBJ5) $(OBJ7) $(MINUITLIB) -o $(BINDIR)/$@ $(MINUITINC) $(F2CINC) $(LIBDIR) $(F2CLIB) $(LIB0) $(LIB3) $(CFLAGS) $(FC_WARNING)
#				rm -f $(BINDIR)/$@.o
aeros_mix			: aeros_mix.c $(OBJ5) $(OBJ7)
				$(CC) $(MINUITDEF) $@.c $(OBJ5) $(OBJ7) $(MINUITLIB) -o $(BINDIR)/$@ $(MINUITINC) $(F2CINC) $(LIBDIR) $(F2CLIB) $(LIB0) $(LIB3) $(LIB4) $(INTEL_LIBDIR) $(CFLAGS) $(WARNING)
#aeros_opt			: aeros_opt.c $(OBJ5) $(OBJ7)
#				$(CC) -c $@.c -o $(BINDIR)/$@.o $(MINUITINC) $(WARNING)
#				$(FC) $(MINUITDEF) $(BINDIR)/$@.o $(OBJ5) $(OBJ7) $(MINUITLIB) -o $(BINDIR)/$@ $(MINUITINC) $(F2CINC) $(LIBDIR) $(F2CLIB) $(LIB0) $(CFLAGS) $(FC_WARNING)
#				rm -f $(BINDIR)/$@.o
aeros_opt			: aeros_opt.c $(OBJ5) $(OBJ7)
				$(CC) $(MINUITDEF) $@.c $(OBJ5) $(OBJ7) $(MINUITLIB) -o $(BINDIR)/$@ $(MINUITINC) $(F2CINC) $(LIBDIR) $(F2CLIB) $(LIB0) $(LIB4) $(INTEL_LIBDIR) $(CFLAGS) $(WARNING)
#aeros_param			: aeros_param.c $(OBJ5) $(OBJ7)
#				$(CC) -c $@.c -o $(BINDIR)/$@.o $(MINUITINC) $(WARNING)
#				$(FC) $(MINUITDEF) $(BINDIR)/$@.o $(OBJ5) $(OBJ7) $(MINUITLIB) -o $(BINDIR)/$@ $(MINUITINC) $(F2CINC) $(LIBDIR) $(F2CLIB) $(LIB0) $(LIB3) $(CFLAGS) $(FC_WARNING)
#				rm -f $(BINDIR)/$@.o
aeros_param			: aeros_param.c $(OBJ5) $(OBJ7)
				$(CC) $(MINUITDEF) $@.c $(OBJ5) $(OBJ7) $(MINUITLIB) -o $(BINDIR)/$@ $(MINUITINC) $(F2CINC) $(LIBDIR) $(F2CLIB) $(LIB0) $(LIB3) $(LIB4) $(INTEL_LIBDIR) $(CFLAGS) $(WARNING)
aeros_read_profile		: aeros_read_profile.c
				$(CC) $@.c -o $(BINDIR)/$@ $(LIB0) $(WARNING)
aeros_simul			: aeros_simul.c $(OBJ5) $(OBJ7)
				$(CC) -c $@.c -o $(BINDIR)/$@.o $(MINUITINC) $(WARNING)
				$(FC) $(BINDIR)/$@.o $(OBJ5) $(OBJ7) -o $(BINDIR)/$@ $(LIBDIR) $(LIB0) $(CFLAGS) $(FC_WARNING)
				rm -f $(BINDIR)/$@.o
#aeros_test			: aeros_test.c $(OBJ5) $(OBJ7)
#				$(CC) -c $@.c -o $(BINDIR)/$@.o $(MINUITINC) $(WARNING)
#				$(FC) $(MINUITDEF) $(BINDIR)/$@.o $(OBJ5) $(OBJ7) $(MINUITLIB) -o $(BINDIR)/$@ $(MINUITINC) $(F2CINC) $(LIBDIR) $(F2CLIB) $(LIB0) $(CFLAGS) $(FC_WARNING)
#				rm -f $(BINDIR)/$@.o
aeros_test			: aeros_test.c $(OBJ5) $(OBJ7)
				$(CC) $(MINUITDEF) $@.c $(OBJ5) $(OBJ7) $(MINUITLIB) -o $(BINDIR)/$@ $(MINUITINC) $(F2CINC) $(LIBDIR) $(F2CLIB) $(LIB0) $(LIB4) $(CFLAGS) $(WARNING)
				rm -f $(BINDIR)/$@.o
aeros_table			: aeros_common.c aeros_table.c $(OBJ5)
				$(CC) -c $@.c -o $(BINDIR)/$@.o $(WARNING)
				$(FC) $(BINDIR)/$@.o $(OBJ5) $(LIB3) -o $(BINDIR)/$@ $(FC_WARNING)
				rm -f $(BINDIR)/$@.o
aod_vs_vis			: aod_vs_vis.c $(OBJ5) $(OBJ7)
				$(CC) -c $@.c -o $(BINDIR)/$@.o $(MINUITINC) $(WARNING)
				$(FC) $(BINDIR)/$@.o $(OBJ5) $(OBJ7) -o $(BINDIR)/$@ $(LIBDIR) $(LIB0) $(CFLAGS) $(FC_WARNING)
				rm -f $(BINDIR)/$@.o
dls				: matrix_optchr.f matrix_intrpl_vS.f matrix_fixget_vS.f sizedstr.f spline.f
				$(FC) matrix_optchr.f matrix_intrpl_vS.f matrix_fixget_vS.f sizedstr.f spline.f -o $(BINDIR)/$@.run
				rm -f alloc.mod alloc1.mod
mie2new.exe			: mie2new.f $(OBJ4)
				$(FC) $(OBJ4) -o $(BINDIR)/$@ -O
ms720_correct_fov		: ms720_correct_fov.c $(OBJ5) $(OBJ7)
				$(CC) -c $@.c -o $(BINDIR)/$@.o $(MINUITINC) $(WARNING)
				$(FC) $(BINDIR)/$@.o $(OBJ5) $(OBJ7) -o $(BINDIR)/$@ $(LIBDIR) $(LIB0) $(CFLAGS) $(FC_WARNING)
				rm -f $(BINDIR)/$@.o
ms720_correct_tmp		: ms720_correct_tmp.c $(OBJ5) $(OBJ7)
				$(CC) -c $@.c -o $(BINDIR)/$@.o $(MINUITINC) $(WARNING)
				$(FC) $(BINDIR)/$@.o $(OBJ5) $(OBJ7) -o $(BINDIR)/$@ $(LIBDIR) $(LIB0) $(CFLAGS) $(FC_WARNING)
				rm -f $(BINDIR)/$@.o
ms720_fov			: ms720_fov.c $(OBJ1)
				$(LD) $@.c $(OBJ1) -o $(BINDIR)/$@ $(INCDIR) -I$(ROOTINC) $(LIBDIR) $(ROOTLIBS) $(OPT) $(WARNING)
ms720_skyrad			: ms720_skyrad.c $(OBJ1)
				$(CC) $@.c $(OBJ1) -o $(BINDIR)/$@ $(LIB0) $(WARNING)
$(BINDIR)/mie2new.o		: mie2new.f
				$(FC) -c mie2new.f -o $@ -O
mie				: mie.c $(OBJ5)
				$(CC) -c $@.c -o $(BINDIR)/$@.o $(WARNING)
				$(FC) $(BINDIR)/$@.o $(OBJ5) -o $(BINDIR)/$@ $(FC_WARNING)
				rm -f $(BINDIR)/$@.o
mie_size_dist			: mie_size_dist.c $(OBJ5)
				$(CC) -c $@.c -o $(BINDIR)/$@.o $(WARNING)
				$(FC) $(BINDIR)/$@.o $(OBJ5) -o $(BINDIR)/$@ $(FC_WARNING)
				rm -f $(BINDIR)/$@.o
mie_eff2mod			: mie_eff2mod.c $(OBJ5) $(OBJ7)
				$(CC) -c $@.c -o $(BINDIR)/$@.o $(MINUITINC) $(WARNING)
				$(FC) $(BINDIR)/$@.o $(OBJ5) $(OBJ7) -o $(BINDIR)/$@ $(LIBDIR) $(LIB0) $(CFLAGS) $(FC_WARNING)
				rm -f $(BINDIR)/$@.o
phase_integ			: phase_integ.c $(OBJ1)
				$(CC) $@.c $(OBJ1) -o $(BINDIR)/$@ $(LIB0) $(WARNING)
resolution_effect		: resolution_effect.c
				$(CC) $@.c -o $(BINDIR)/$@ $(LIB0) $(WARNING)
$(BINDIR)/MIEV0.o		: MIEV0.f
				$(FC) -c $(*F).f -o $@
$(BINDIR)/ErrPack.o		: ErrPack.f
				$(FC) -c $(*F).f -o $@
db				: main.f $(OBJ6)
				$(FC) main.f $(OBJ6) -o $(BINDIR)/$@
$(BINDIR)/aeros_bb.o		: aeros_bb.f
				$(FC) -c $(*F).f -o $@
$(BINDIR)/airmie.o		: airmie.f
				$(FC) -c $(*F).f -o $@
$(BINDIR)/append_f.o		: append_f.f
				$(FC) -c $(*F).f -o $@
$(BINDIR)/findpar.o		: findpar.f
				$(FC) -c $(*F).f -o $@
$(BINDIR)/legendre_coeff.o	: legendre_coeff.f
				$(FC) -c $(*F).f -o $@
$(BINDIR)/r_min_max_comput.o	: r_min_max_comput.f
				$(FC) -c $(*F).f -o $@
$(BINDIR)/size.o		: size.f
				$(FC) -c $(*F).f -o $@
$(BINDIR)/aerprf.o		: aerprf.f
				$(FC) -c $(*F).f -o $(BINDIR)/$@ $(FC_WARNING)
				rm -f $(*F)__genmod.f90
				rm -f $(*F)__genmod.mod
$(BINDIR)/prfdta.o		: prfdta.f
				$(FC) -c $(*F).f -o $(BINDIR)/$@ $(FC_WARNING)
$(BINDIR)/strutil.o		: strutil.c
				$(CC) -c $(*F).c -o $@
tmx_size_dist			: tmx_size_dist.c $(OBJ1) $(OBJ7)
				$(CC) -c $@.c -o $(BINDIR)/$@.o $(WARNING)
				$(FC) -o $(BINDIR)/$@ $(BINDIR)/$@.o $(OBJ1) $(OBJ7) -lm
				rm -f $(BINDIR)/$@.o
tmx_test			: tmx_test.c $(OBJ1) $(OBJ7)
				$(CC) -c $@.c -o $(BINDIR)/$@.o $(WARNING)
				$(FC) -o $(BINDIR)/$@ $(BINDIR)/$@.o $(OBJ1) $(OBJ7) -lm
				rm -f $(BINDIR)/$@.o
aeros_tmx			: aeros_tmx.c $(OBJ1) $(OBJ7)
				$(CC) -c $@.c -o $(BINDIR)/$@.o $(WARNING)
				$(FC) -o $(BINDIR)/$@ $(BINDIR)/$@.o $(OBJ1) $(OBJ7) -lm
				rm -f $(BINDIR)/$@.o
aeros_tmq			: aeros_tmq.c $(OBJ1) $(OBJ8)
				$(CC) -c $@.c -o $(BINDIR)/$@.o $(WARNING)
				$(FC) -o $(BINDIR)/$@ $(BINDIR)/$@.o $(OBJ1) $(OBJ8) -lm
				rm -f $(BINDIR)/$@.o
tmq_calc			: tmq_calc.c $(OBJ1) $(OBJ8)
				$(CC) -c $@.c -o $(BINDIR)/$@.o $(WARNING)
				$(FC) -o $(BINDIR)/$@ $(BINDIR)/$@.o $(OBJ1) $(OBJ8) -lm
				rm -f $(BINDIR)/$@.o
tmq_large			: tmq_calc.c $(OBJ1) $(OBJ8L)
				$(CC) -c tmq_calc.c -o $(BINDIR)/tmq_calc.o $(WARNING)
				$(FC) -o $(BINDIR)/$@ $(BINDIR)/tmq_calc.o $(OBJ1) $(OBJ8L) -lm
				rm -f $(BINDIR)/tmq_calc.o
mie_calc			: mie_calc.c $(OBJ1) $(OBJM)
				$(CC) -c $@.c -o $(BINDIR)/$@.o $(WARNING)
				$(FC) $(BINDIR)/$@.o $(OBJ1) $(OBJM) -o $(BINDIR)/$@
				rm -f $(BINDIR)/$@.o
read_mie			: read_mie.c $(OBJ1)
				$(CC) -o $(BINDIR)/$@ $@.c $(OBJ1) $(LIB0) $(WARNING)
read_tmq			: read_tmq.c $(OBJ1)
				$(CC) -o $(BINDIR)/$@ $@.c $(OBJ1) $(LIB0) $(WARNING)
$(BINDIR)/tmx.o			: tmx.f tmd.par.f
				$(FC) -c $(*F).f -o $@ -mcmodel=medium
$(BINDIR)/lpd.o			: lpd.f
				$(FC) -c $(*F).f -o $@
$(BINDIR)/tmq.o			: tmq.f tmq.par.f
				$(FC) -c $(*F).f -o $@ -mcmodel=medium
$(BINDIR)/tmq_large.o		: tmq_large.f tmq.par.f
				$(FC) -c $(*F).f -o $@ -mcmodel=medium
$(BINDIR)/lpq.o			: lpq.f
				$(FC) -c $(*F).f -o $@
$(BINDIR)/newrad.o		: newrad.f
				$(FC) -c $(*F).f -o $@ -O -r8
