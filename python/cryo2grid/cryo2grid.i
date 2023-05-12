%define DOCSTRING
"############################################\n\
# This file is distributed under the GNU   #\n\
# LPGL-2.1-or-later Open Source License.   #\n\
# See LICENSE file in top directory for    #\n\
# details.                                 #\n\
#                                          #\n\
# Copyright (c) 2020 ADsandbox developers: #\n\
#               Andreas F. Tillack         #\n\
#               Althea A. Hansel           #\n\
#               Matthew Holcomb            #\n\
#           Forli Lab @ Scripps Research   #\n\
############################################\n"
%enddef

%module(docstring=DOCSTRING, package="cryo2grid") cryo2grid_wrapper

%begin %{
#define SWIG_PYTHON_2_UNICODE
//#define SWIG_FILE_WITH_INIT
%}

// Add standard C++ library
%include "std_array.i"
%include "std_list.i"
%include "std_map.i"
%include "std_vector.i"
%include "std_string.i"

// Help SWIG to understand some special types, like list of strings
namespace std {
    %template(IntVector)          vector<int>;
    %template(FloatVector)        vector<float>;
    %template(FloatVectorVector)  vector<vector<float>>;
    %template(DoubleVector)       vector<double>;
    %template(DoubleVectorVector) vector<vector<double>>;
    %template(StringVector)       vector<string>;
    %template(ConstCharVector)    vector<const char*>;
}

%include "c2g.i"

