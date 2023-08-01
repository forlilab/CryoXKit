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
CXXFLAGS=-DC2G_VERSION=\"$(GIT_VERSION)\" -march=native -Wextra -std=c++11 -O3 -fopenmp
LDFLAGS=
EXECUTABLE=cryo2grid

SRCS=
OBJS=$(SRCS:.cpp=.o)

all: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(EXECUTABLE) cryo2grid.cpp $(OBJS)

DEPS := $(patsubst %.o,%.d,$(OBJS))
-include $(DEPS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

.PHONY: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(CXX) $(LDFLAGS) -o $(EXECUTABLE) cryo2grid.cpp $(OBJS)

clean:
	rm -f $(OBJS) *.d $(EXECUTABLE) example/*.map example/*.grid.mrc
