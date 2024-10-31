##########################################
# This file is distributed under the GNU #
# LPGL-2.1-or-later Open Source License. #
# See LICENSE file for details.          #
#                                        #
# Copyright (c) 2023 Andreas F. Tillack  #
#           Forli Lab @ Scripps Research #
##########################################

CXX=g++
UNAME_CMD := $(shell uname -s)
# for OpenMP support on macOS use clang++ (i.e. from Brew)
ifeq ($(UNAME_CMD),Darwin)
	CXX=clang++
endif
GIT_VERSION := $(shell git describe --abbrev=7 --dirty --always --tags | sed 's/dirty/mod/g')
CXXFLAGS=-DCXK_VERSION=\"$(GIT_VERSION)\" -march=native -Wextra -std=c++11 -O3 -fopenmp
LDFLAGS=
EXECUTABLE=cryoXkit

SRCS=ScalarMat.cpp
OBJS=$(SRCS:.cpp=.o)

all: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(EXECUTABLE) cryoXkit.cpp $(OBJS)

DEPS := $(patsubst %.o,%.d,$(OBJS))
-include $(DEPS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

.PHONY: $(EXECUTABLE) python

$(EXECUTABLE): $(OBJS)
	$(CXX) $(LDFLAGS) -o $(EXECUTABLE) cryoXkit.cpp $(OBJS)

python:
	rm -rf python/build/
	rm -rf python/src/
	cd python ; python setup.py build install

clean:
	rm -rf python/build/
	rm -rf python/src/
	rm -rf python/cryoXkit.egg-info/
	rm -rf python/cryoXkit/__pycache__/
	rm -rf python/cryoXkit/__init__.py
	rm -f python/cryoXkit/cryoXkit_wrap.cpp
	rm -f python/cryoXkit/cryoXkit_wrapper.py
	rm -f $(OBJS) *.d $(EXECUTABLE) example/*.map example/*.grid.mrc
