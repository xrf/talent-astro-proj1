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
FFLAGS=-O3 -march=native -I/usr/include
LNETCDF=-lnetcdf -lnetcdff

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

all: vh1-mpi vh1-serial vh1-starter

clean:
	@rm -fr dist

vh1-mpi: $(VH1_DIR)/vh1-mpi

vh1-serial: $(VH1_DIR)/vh1-serial

vh1-starter: $(VH1_DIR)/vh1-starter

run: $(VH1_DIR)/vh1-serial
	cd $(VH1_DIR) && ./vh1-serial

run-mpi: $(VH1_DIR)/vh1-mpi
	cd $(VH1_DIR) && mpirun -n $(NP) ./vh1-mpi

run-starter: $(VH1_DIR)/vh1-starter
	cd $(VH1_DIR) && ./vh1-starter

$(VH1_DIR)/vh1-serial: $(VH1) $(VH1_DIR)/src/Serial/vhone.f90
	@cd $(VH1_DIR)/src/Serial && \
	$(MAKE) F90="$(FC)" \
	        FFLAGS="$(FFLAGS) -c" \
	        LDR="$(FC)" \
	        LDRFLAGS="$(LDFLAGS)" \
	        LIBS="$(LNETCDF)"

$(VH1_DIR)/vh1-mpi: $(VH1)
	@cd $(VH1_DIR)/src/Parallel && \
	$(MAKE) F90="$(MPIFC)" \
	        FFLAGS="$(FFLAGS) -c" \
	        LDR="$(MPIFC)" \
	        LDRFLAGS="$(LDFLAGS)" \
	        LIBS="$(LNETCDF)"

$(VH1_DIR)/vh1-starter: $(VH1)
	@cd $(VH1_DIR)/src/Starter && \
	$(MAKE) F90="$(FC)" \
	        FFLAGS="$(FFLAGS) -c" \
	        LDR="$(FC)" \
	        LDRFLAGS="$(LDFLAGS)"

$(VH1_DIR)/src/Serial/vhone.f90: vhone.f90 $(VH1)
	@cp vhone.f90 $@

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