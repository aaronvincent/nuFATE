# Compiler
#CC=clang
#CXX=clang++
#AR=ar
#LD=clang++


DYN_SUFFIX=.dylib
DYN_OPT=-dynamiclib -compatibility_version $(VERSION) -current_version $(VERSION)

VERSION=1.0.0

ifeq (${PREFIX},)
PREFIX=/usr/local
endif

PATH_nuFATE=$(shell pwd)

SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(patsubst src/%.cpp,build/%.o,$(SOURCES))

EXAMPLES := examples/example

CXXFLAGS= -std=c++11 -g

# Directories

GSL_CFLAGS=-I/usr/local/Cellar/gsl/2.4/include
GSL_LDFLAGS=-L/usr/local/Cellar/gsl/2.4/lib -lgsl -lgslcblas -lm
HDF5_CFLAGS=-I/usr/local/Cellar/hdf5/1.10.1_1/lib/../include
HDF5_LDFLAGS=-L/usr/local/Cellar/hdf5/1.10.1_1/lib -lhdf5_hl -lhdf5 -L/usr/local/opt/szip/lib -lsz -lz -ldl -lm

INCnuFATE=$(PATH_nuFATE)/include/nuFATE
LIBnuFATE=$(PATH_nuFATE)/lib

# FLAGS
CFLAGS= -O3 -fPIC -I$(INCnuFATE) $(GSL_CFLAGS) $(HDF5_CFLAGS)
LDFLAGS= -Wl,-rpath -Wl,$(LIBnuFATE) -L$(LIBnuFATE)
LDFLAGS+= $(SQUIDS_LDFLAGS) $(GSL_LDFLAGS) $(HDF5_LDFLAGS)
EXMAPLES_FLAGS=-I$(INCnuFATE) $(CXXFLAGS) $(CFLAGS)

# Project files
NAME=nuFATE
STAT_PRODUCT=lib/lib$(NAME).a
DYN_PRODUCT=lib/lib$(NAME)$(DYN_SUFFIX)

# Compilation rules
all: $(STAT_PRODUCT) $(DYN_PRODUCT)

examples : $(EXAMPLES)

$(DYN_PRODUCT) : $(OBJECTS)
	@echo Linking $(DYN_PRODUCT)
	@$(CXX) $(DYN_OPT)  $(LDFLAGS) -o $(DYN_PRODUCT) $(OBJECTS)

$(STAT_PRODUCT) : $(OBJECTS)
	@echo Linking $(STAT_PRODUCT)
	@$(AR) -rcs $(STAT_PRODUCT) $(OBJECTS)

build/%.o : src/%.cpp
	@echo Compiling $< to $@
	@$(CXX) $(CXXFLAGS) -c $(CFLAGS) $< -o $@

examples/example: $(DYN_PRODUCT) examples/example.cpp
	@echo Compatible example
	@$(CXX) $(EXMAPLES_FLAGS) examples/example.cpp -lnuFATE $(LDFLAGS) -o $@

.PHONY: install uninstall clean test docs
clean:
	@echo Erasing generated files
	@rm -f $(PATH_nuSQUIDS)/build/*.o
	@rm -f $(PATH_nuSQUIDS)/$(STAT_PRODUCT) $(PATH_nuSQUIDS)/$(DYN_PRODUCT)

doxygen:
	@doxygen src/doxyfile
docs:
	@doxygen src/doxyfile

test: $(DYN_PRODUCT) $(STAT_PRODUCT)
	@cd ./test ; ./run_tests

install: $(DYN_PRODUCT) $(STAT_PRODUCT)
	@echo Installing headers in $(PREFIX)/include/nuFATE
	@mkdir -p $(PREFIX)/include/nuFATE
	@cp $(INCnuSQUIDS)/nuFATE/*.h $(PREFIX)/include/nuFATE
	@echo Installing libraries in $(PREFIX)/lib
	@mkdir -p $(PREFIX)/lib
	@cp $(DYN_PRODUCT) $(STAT_PRODUCT) $(PREFIX)/lib
	@echo Installing config information in $(PREFIX)/lib/pkgconfig
	@mkdir -p $(PREFIX)/lib/pkgconfig

uninstall:
	@echo Removing headers from $(PREFIX)/include/nuFATE
	@rm -rf $(PREFIX)/include/nuFATE
	@echo Removing libraries from $(PREFIX)/lib
	@rm -f $(PREFIX)/lib/$(DYN_PRODUCT)
	@rm -f $(PREFIX)/lib/$(STAT_PRODUCT)

