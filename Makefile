##########################################
# This file is distributed under the GNU #
# LPGL-2.1-or-later Open Source License. #
# See LICENSE file for details.          #
#                                        #
# Copyright (c) 2023 Andreas F. Tillack  #
#           Forli Lab @ Scripps Research #
##########################################

CXX=g++
CXXFLAGS=-Wall -march=native -Wextra -std=c++11 -O3
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
