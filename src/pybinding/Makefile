
# Compiler
CC=gcc
CXX=g++
AR=ar
LD=ld

DYN_SUFFIX=.so
DYN_OPT=-shared -Wl,-soname,$(shell basename $(DYN_PRODUCT))

# add python and boost_python
PYTHON_VERSION= $(shell echo `python -c 'import sys; print(str(sys.version_info.major)+"."+str(sys.version_info.minor))'`) 
CVMFS_DIR= $(shell echo `which python | sed -e 's/\/bin\/python//'`)

#PYTHON_LIBPATH=$(shell echo `python -c 'import sys; import re; print([ y for y in sys.path if re.search("\/lib\/python2.7",y)!=None and re.search("\/cvmfs",y)!=None ][0]);'`)
PYTHON_LIBPATH=$(CVMFS_DIR)/lib/python$(PYTHON_VERSION)
CVMFS_LIBPATH= $(CVMFS_DIR)/lib
PYTHON_INCPATH= $(CVMFS_DIR)/include/python$(PYTHON_VERSION)

#.PHONY: debug
#debug:
#	@echo CVMFS_DIR is $(CVMFS_DIR)

LDFLAGS+= -L$(PYTHON_LIBPATH) -L$(CVMFS_LIBPATH)
LDFLAGS+= -lpython$(PYTHON_VERSION)
LDFLAGS+= -Wl,-rpath -Wl,$(CVMFS_LIBPATH) -L$(CVMFS_LIBPATH) -lboost_python

INCCFLAGS+= -I$(PYTHON_INCPATH) -Wno-abi -I$(CVMFS_LIBPATH)/../include -I../include/

PATH_nuFATE= $(shell pwd)/../../
PATH_nuFATEpy= $(shell pwd)

SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

# Directories
LIBnuFATE=$(PATH_nuFATE)/lib
INCnuFATE=$(PATH_nuFATE)/include

LDFLAGS+= -L$(LIBnuFATE) -lnuFATE
INCCFLAGS+= -I$(INCnuFATE) -I$(PATH_nuFATEpy)
CXXFLAGS= -O3 -fPIC -std=c++11 $(INCCFLAGS)

# Project files
NAME=nuFATEpy
STAT_PRODUCT=$(PATH_nuFATEpy)/$(NAME).a
DYN_PRODUCT=$(PATH_nuFATEpy)/$(NAME)$(DYN_SUFFIX)

OS_NAME=$(shell uname -s)
ifeq ($(OS_NAME),Linux)
	DYN_SUFFIX=.so
	DYN_OPT=-shared -Wl,-soname,$(DYN_PRODUCT)
endif
ifeq ($(OS_NAME),Darwin)
	DYN_SUFFIX=.so
	DYN_OPT=-dynamiclib -install_name $(PATH_nuSQUIDSpy)/$(DYN_PRODUCT)
endif

# Compilation rules
all: $(STAT_PRODUCT) $(DYN_PRODUCT)

$(DYN_PRODUCT) : $(OBJECTS)
	@echo Linking $(DYN_PRODUCT)
	@$(CXX) $(DYN_OPT)  $(LDFLAGS) -o $(DYN_PRODUCT) $(OBJECTS)

$(STAT_PRODUCT) : $(OBJECTS)
	@echo Linking $(STAT_PRODUCT)
	@$(AR) -rcs $(STAT_PRODUCT) $(OBJECTS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f *.o *.so *.a

