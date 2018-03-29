#!/bin/bash

# ./grading-script.sh <test-dir> <grading-scheme.txt>

test_dir="$1"

# Set up scratch space for grading
dir="grading"
mkdir -p $dir
cp *.h *.cpp SConstruct Makefile $dir

# Check for cheating
token=`mktemp XXXXXXXXXXXXXXXXXXXXXXXX`
rm $token

g++ -Wall -g -O3 -std=c++11 -M -x c++ - < minigl.cpp | tr ' \\' '\n' | grep . | sort -u > $dir/actual_include.txt
g++ -Wall -g -O3 -std=c++11 -M -x c++ - <<EOF        | tr ' \\' '\n' | grep . | sort -u > $dir/max_include.txt
#include "minigl.h"
#include "vec.h"
#include "mat.h"
#include <cstdlib>
#include <csignal>
#include <csetjmp>
#include <cstdarg>
#include <typeinfo>
#include <bitset>
#include <functional>
#include <utility>
#include <ctime>
#include <cstddef>
#include <new>
#include <memory>
#include <climits>
#include <cfloat>
#include <limits>
#include <exception>
#include <stdexcept>
#include <cassert>
#include <cerrno>
#include <cctype>
#include <cwctype>
#include <cstring>
#include <cwchar>
#include <string>
#include <vector>
#include <deque>
#include <list>
#include <set>
#include <map>
#include <stack>
#include <queue>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <complex>
#include <valarray>
#include <numeric>
#include <iosfwd>
#include <ios>
#include <istream>
#include <ostream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <strstream>
#include <iomanip>
#include <streambuf>
#include <cstdio>
#include <locale>
#include <clocale>
#include <ciso646>
#include <assert.h>
#include <complex.h>
#include <ctype.h>
#include <errno.h>
#include <fenv.h>
#include <float.h>
#include <inttypes.h>
#include <iso646.h>
#include <limits.h>
#include <locale.h>
#include <math.h>
#include <setjmp.h>
#include <signal.h>
#include <stdalign.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>
#include <time.h>
#include <uchar.h>
#include <wchar.h>
#include <wctype.h>
EOF

if comm -23 $dir/actual_include.txt $dir/max_include.txt | grep . ; then
    echo "FAIL: forbidden include"
    exit 1
fi

# Compile project
{
    pushd $dir
    scons 2> /dev/null
    if [ $? -ne 0 ] ; then
        echo "SCONS FAIL: Did not compile, trying compiling with make ..."
	make clean; make;
	if [ $? -ne 0 ] ; then
		echo "FAIL: Did not compile with make, exiting!"
		exit 1;
	fi
    fi
    echo "COMPILED!"
    popd >/dev/null
}

# Check results against grading criteria
./grading-helper.pl $dir "$test_dir" $token
