.POSIX:
#
# VH1 requires the NetCDF library:
#
#   - Apt:     `apt-get install libnetcdf-dev`
#   - Pacman:  `pacman -S netcdf`
#
# ----------------------------------------------------------------------------

# compilers
FC=gfortran
MPIFC=mpif90

# compiler flags
INCLUDE_DIRS=-I/usr/include #-I/usr/local/include -I/opt/local/include
LIB_DIRS=-L/usr/lib #-L/usr/local/lib -L/opt/local/lib
FFLAGS=$(INCLUDE_DIRS) -O3 -Wall -march=native -Wa,-q
LNETCDF=$(LIB_DIRS) -lnetcdf -lnetcdff

# number of processors
NP=16
# note: NP is not equal to P :)

# VH1 remote origin
VH1_URL=http://astro.physics.ncsu.edu/pub/VH-1/VH1.tar

# VH1 local directory
VH1_PAR=dist/deps
VH1_DIR=$(VH1_PAR)/VH1

# dummy file used to check if VH1 has been downloaded
VH1=$(VH1_PAR)/.VH1-cookie

# ----------------------------------------------------------------------------

vh1-serial: $(VH1_DIR)/vh1-serial

all: vh1-mpi vh1-serial vh1-starter

clean:
	@rm -fr dist

vh1-mpi: $(VH1_DIR)/vh1-mpi

vh1-starter: $(VH1_DIR)/vh1-starter

run: $(VH1_DIR)/vh1-serial $(VH1_DIR)/indat
	@rm -fr $(VH1_DIR)/output/*
	@cd $(VH1_DIR) && ./vh1-serial

run-mpi: $(VH1_DIR)/vh1-mpi $(VH1_DIR)/indat
	@rm -fr $(VH1_DIR)/output/*
	@cd $(VH1_DIR) && mpirun -n $(NP) ./vh1-mpi

run-starter: $(VH1_DIR)/vh1-starter $(VH1_DIR)/indat
	@rm -fr $(VH1_DIR)/output/*
	@cd $(VH1_DIR) && ./vh1-starter

$(VH1_DIR)/vh1-serial: \
    $(VH1) \
    $(VH1_DIR)/src/Serial/dtcon.f90 \
    $(VH1_DIR)/src/Serial/dump.f90 \
    $(VH1_DIR)/src/Serial/images.f90 \
    $(VH1_DIR)/src/Serial/init.f90 \
    $(VH1_DIR)/src/Serial/prin.f90 \
    $(VH1_DIR)/src/Serial/sweepx.f90 \
    $(VH1_DIR)/src/Serial/sweepy.f90 \
    $(VH1_DIR)/src/Serial/sweepz.f90 \
    $(VH1_DIR)/src/Serial/vhone.f90 \
    $(VH1_DIR)/src/Serial/zonemod.f90 \
    $(VH1_DIR)/src/PPMLR/boundary.f90 \
    $(VH1_DIR)/src/PPMLR/evolve.f90 \
    $(VH1_DIR)/src/PPMLR/flatten.f90 \
    $(VH1_DIR)/src/PPMLR/forces.f90 \
    $(VH1_DIR)/src/PPMLR/parabola.f90 \
    $(VH1_DIR)/src/PPMLR/ppmlr.f90 \
    $(VH1_DIR)/src/PPMLR/remap.f90 \
    $(VH1_DIR)/src/PPMLR/riemann.f90 \
    $(VH1_DIR)/src/PPMLR/states.f90 \
    $(VH1_DIR)/src/PPMLR/vh1mods.f90 \
    $(VH1_DIR)/src/PPMLR/volume.f90
	@cd $(VH1_DIR)/src/Serial && \
	$(MAKE) F90="$(FC)" \
	        FFLAGS="$(FFLAGS) -c" \
	        LDR="$(FC)" \
	        LDRFLAGS="$(LDFLAGS)" \
	        LIBS="$(LNETCDF)"

$(VH1_DIR)/vh1-mpi: \
    $(VH1) \
    $(VH1_DIR)/src/Parallel/dtcon.f90 \
    $(VH1_DIR)/src/Parallel/dump.f90 \
    $(VH1_DIR)/src/Parallel/images.f90 \
    $(VH1_DIR)/src/Parallel/init.f90 \
    $(VH1_DIR)/src/Parallel/mergeslabs.f90 \
    $(VH1_DIR)/src/Parallel/prin.f90 \
    $(VH1_DIR)/src/Parallel/sweepx1.f90 \
    $(VH1_DIR)/src/Parallel/sweepx2.f90 \
    $(VH1_DIR)/src/Parallel/sweepy.f90 \
    $(VH1_DIR)/src/Parallel/sweepz.f90 \
    $(VH1_DIR)/src/Parallel/vhone.f90 \
    $(VH1_DIR)/src/Parallel/zonemod.f90 \
    $(VH1_DIR)/src/PPMLR/boundary.f90 \
    $(VH1_DIR)/src/PPMLR/evolve.f90 \
    $(VH1_DIR)/src/PPMLR/flatten.f90 \
    $(VH1_DIR)/src/PPMLR/forces.f90 \
    $(VH1_DIR)/src/PPMLR/parabola.f90 \
    $(VH1_DIR)/src/PPMLR/ppmlr.f90 \
    $(VH1_DIR)/src/PPMLR/remap.f90 \
    $(VH1_DIR)/src/PPMLR/riemann.f90 \
    $(VH1_DIR)/src/PPMLR/states.f90 \
    $(VH1_DIR)/src/PPMLR/vh1mods.f90 \
    $(VH1_DIR)/src/PPMLR/volume.f90
	@cd $(VH1_DIR)/src/Parallel && \
	$(MAKE) F90="$(MPIFC)" \
	        FFLAGS="$(FFLAGS) -c" \
	        LDR="$(MPIFC)" \
	        LDRFLAGS="$(LDFLAGS)" \
	        LIBS="$(LNETCDF)"

$(VH1_DIR)/vh1-starter: \
    $(VH1) \
    $(VH1_DIR)/src/Starter/vhone.f90 \
    $(VH1_DIR)/src/PPMLR/boundary.f90 \
    $(VH1_DIR)/src/PPMLR/evolve.f90 \
    $(VH1_DIR)/src/PPMLR/flatten.f90 \
    $(VH1_DIR)/src/PPMLR/forces.f90 \
    $(VH1_DIR)/src/PPMLR/parabola.f90 \
    $(VH1_DIR)/src/PPMLR/ppmlr.f90 \
    $(VH1_DIR)/src/PPMLR/remap.f90 \
    $(VH1_DIR)/src/PPMLR/riemann.f90 \
    $(VH1_DIR)/src/PPMLR/states.f90 \
    $(VH1_DIR)/src/PPMLR/vh1mods.f90 \
    $(VH1_DIR)/src/PPMLR/volume.f90
	@cd $(VH1_DIR)/src/Starter && \
	$(MAKE) F90="$(FC)" \
	        FFLAGS="$(FFLAGS) -c" \
	        LDR="$(FC)" \
	        LDRFLAGS="$(LDFLAGS)"

$(VH1_DIR)/src/PPMLR/boundary.f90: $(VH1)
$(VH1_DIR)/src/PPMLR/evolve.f90: $(VH1)
$(VH1_DIR)/src/PPMLR/flatten.f90: $(VH1)
$(VH1_DIR)/src/PPMLR/parabola.f90: $(VH1)
$(VH1_DIR)/src/PPMLR/ppmlr.f90: $(VH1)
$(VH1_DIR)/src/PPMLR/remap.f90: $(VH1)
$(VH1_DIR)/src/PPMLR/riemann.f90: $(VH1)
$(VH1_DIR)/src/PPMLR/states.f90: $(VH1)
$(VH1_DIR)/src/PPMLR/vh1mods.f90: $(VH1)
$(VH1_DIR)/src/PPMLR/volume.f90: $(VH1)
$(VH1_DIR)/src/Parallel/dtcon.f90: $(VH1)
$(VH1_DIR)/src/Parallel/dump.f90: $(VH1)
$(VH1_DIR)/src/Parallel/images.f90: $(VH1)
$(VH1_DIR)/src/Parallel/init.f90: $(VH1)
$(VH1_DIR)/src/Parallel/mergeslabs.f90: $(VH1)
$(VH1_DIR)/src/Parallel/prin.f90: $(VH1)
$(VH1_DIR)/src/Parallel/sweepx1.f90: $(VH1)
$(VH1_DIR)/src/Parallel/sweepx2.f90: $(VH1)
$(VH1_DIR)/src/Parallel/sweepy.f90: $(VH1)
$(VH1_DIR)/src/Parallel/sweepz.f90: $(VH1)
$(VH1_DIR)/src/Parallel/vhone.f90: $(VH1)
$(VH1_DIR)/src/Parallel/zonemod.f90: $(VH1)
$(VH1_DIR)/src/Serial/dtcon.f90: $(VH1)
$(VH1_DIR)/src/Serial/dump.f90: $(VH1)
$(VH1_DIR)/src/Serial/images.f90: $(VH1)
$(VH1_DIR)/src/Serial/prin.f90: $(VH1)
$(VH1_DIR)/src/Serial/sweepx.f90: $(VH1)
$(VH1_DIR)/src/Serial/sweepy.f90: $(VH1)
$(VH1_DIR)/src/Serial/sweepz.f90: $(VH1)
$(VH1_DIR)/src/Starter/vhone.f90: $(VH1)

$(VH1_DIR)/src/Serial/init.f90: init.f90 $(VH1)
	@cp init.f90 $@

$(VH1_DIR)/src/Serial/vhone.f90: vhone.f90 $(VH1)
	@cp vhone.f90 $@

$(VH1_DIR)/src/Serial/zonemod.f90: zonemod.f90 $(VH1)
	@cp zonemod.f90 $@

$(VH1_DIR)/src/PPMLR/forces.f90: forces.f90 $(VH1)
	@cp forces.f90 $@

$(VH1_DIR)/indat: indat $(VH1)
	@cp indat $@

$(VH1):
	@mkdir -p $(VH1_PAR)
	@cd $(VH1_PAR) && \
	download() { \
	    if command >/dev/null 2>&1 -v curl; \
	    then curl -fL -- "$$1"; \
	    else \
	        if command >/dev/null 2>&1 -v wget; \
	        then wget -O- -- "$$1"; \
	        else \
	            echo >&2 "error: missing 'curl' or 'wget'" && \
	            return 127; \
	        fi; \
	    fi; \
	} && \
	download $(VH1_URL) | tar xf -
	@touch $@
