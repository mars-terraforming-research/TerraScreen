
#     GCM v23   2010
#     Ames Mars GCM group
#     Jeffery Hollingsworth, PI
#     NASA Ames Research Center

F90_COMP=gfortran
#F_OPTS=-O3 -finit-local-zero -frecord-marker=4
F_OPTS=-O3 -frecord-marker=4

all: modules.mod TerraScreen

TerraScreen:  TerraScreen.o radsetup.o setrad.o setspi.o setspv.o\
         initinterp.o lagrange.o ini_optdst.o ini_optcld.o\
         settozero.o fillpt.o filltaucum.o jsrchgt.o dustprofile.o\
         optcv.o optci.o tpindex.o sfluxv.o gfluxv.o\
         sfluxi.o gfluxi.o getdetau.o dsolver.o dtridgl.o cldprofile.o\
         modules.o boxinterp.o convect.o

	$(F90_COMP) -o TerraScreen\
               TerraScreen.o radsetup.o setrad.o setspi.o setspv.o\
               initinterp.o lagrange.o ini_optdst.o ini_optcld.o\
               settozero.o fillpt.o filltaucum.o jsrchgt.o dustprofile.o\
               optcv.o optci.o tpindex.o sfluxv.o gfluxv.o cldprofile.o\
               sfluxi.o gfluxi.o getdetau.o dsolver.o dtridgl.o\
               modules.o boxinterp.o convect.o

TerraScreen.o: TerraScreen.f
	$(F90_COMP) -c $(F_OPTS) TerraScreen.f
radsetup.o:  radsetup.f90
	$(F90_COMP) -c $(F_OPTS) radsetup.f90
setrad.o: setrad.f90
	$(F90_COMP) -c $(F_OPTS) setrad.f90
initinterp.o: initinterp.f90
	$(F90_COMP) -c $(F_OPTS) initinterp.f90
boxinterp.o: initinterp.f90
	$(F90_COMP) -c $(F_OPTS) boxinterp.f90
lagrange.o: lagrange.f90
	$(F90_COMP) -c $(F_OPTS) lagrange.f90
setspi.o: setspi.f90
	$(F90_COMP) -c $(F_OPTS) setspi.f90
setspv.o: setspv.f90
	$(F90_COMP) -c $(F_OPTS) setspv.f90
ini_optdst.o: ini_optdst.f
	$(F90_COMP) -c $(F_OPTS) ini_optdst.f
ini_optcld.o: ini_optcld.f
	$(F90_COMP) -c $(F_OPTS) ini_optcld.f
settozero.o: settozero.f
	$(F90_COMP) -c $(F_OPTS) settozero.f
fillpt.o: fillpt.f90
	$(F90_COMP) -c $(F_OPTS) fillpt.f90
filltaucum.o: filltaucum.f
	$(F90_COMP) -c $(F_OPTS) filltaucum.f
jsrchgt.o: jsrchgt.f90
	$(F90_COMP) -c $(F_OPTS) jsrchgt.f90
dustprofile.o: dustprofile.f
	$(F90_COMP) -c $(F_OPTS) dustprofile.f
cldprofile.o: cldprofile.f
	$(F90_COMP) -c $(F_OPTS) cldprofile.f
optcv.o: optcv.f90
	$(F90_COMP) -c $(F_OPTS) optcv.f90
optci.o: optci.f90
	$(F90_COMP) -c $(F_OPTS) optci.f90
tpindex.o: tpindex.f90
	$(F90_COMP) -c $(F_OPTS) tpindex.f90
sfluxv.o: sfluxv.f90
	$(F90_COMP) -c $(F_OPTS) sfluxv.f90
gfluxv.o: gfluxv.f90
	$(F90_COMP) -c $(F_OPTS) gfluxv.f90
sfluxi.o: sfluxi.f90
	$(F90_COMP) -c $(F_OPTS) sfluxi.f90
gfluxi.o: gfluxi.f90
	$(F90_COMP) -c $(F_OPTS) gfluxi.f90
getdetau.o: getdetau.f90
	$(F90_COMP) -c $(F_OPTS) getdetau.f90
dsolver.o: dsolver.f90
	$(F90_COMP) -c $(F_OPTS) dsolver.f90
dtridgl.o: dtridgl.f90
	$(F90_COMP) -c $(F_OPTS) dtridgl.f90
convect.o: convect.f
	$(F90_COMP) -c $(F_OPTS) convect.f	
modules.mod: modules.o
modules.o: modules.f90
	$(F90_COMP) -c $(F_OPTS) modules.f90
clean:
	rm -f *.o *.mod *.out TerraScreen
