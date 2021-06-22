CC = g++
#CC = /usr/local/Cellar/gcc/5.3.0/bin/g++-5
#CC = /usr/local/Cellar/gcc/6.3.0_1/bin/g++-6
#CFLAGS = -std=c++14 -stdlib=libc++ -O1 -march=native -fno-pic -Werror -Wall -Wextra -o
CFLAGS = -std=c++14 -stdlib=libc++ -O3 -march=native -fno-pic -Werror -Wall -Wextra -o
#CFLAGS = -std=gnu++14 -O3 -march=native -fno-pic -Werror -Wall -Wextra -o
CFLAGSDEBUG = -std=c++11 -stdlib=libc++ -Wall -g -o
CFLAGS2 = -L/usr/local/Cellar/boost/1.59.0 -L/usr/local/Cellar/mpfr/3.1.3 -lmpfr -std=c++14 -stdlib=libc++ -Werror -Wall -Wextra -O3 -march=native -o
CFLAGS3 = -L/usr/local/Cellar/mpfr/4.0.2 -lmpfr -L/usr/local/Cellar/mpfi/1.5.3 -lmpfi -std=c++14 -stdlib=libc++ -Werror -Wall -Wextra -O3 -march=native -o
CFLAGSTESTING = -L/usr/local/Cellar/mpfr/4.0.2 -lmpfr -L/usr/local/Cellar/mpfi/1.5.3 -lmpfi -std=c++14 -stdlib=libc++ -o

trace_curve: WBtemplate.hpp DSTools.hpp multiprec.hpp Interval.hpp ValidatedKAM.hpp trace_curve.cpp
	$(CC) $(CFLAGS3) trace_curve trace_curve.cpp

trace_curve_mpfr: WBtemplate.hpp DSTools.hpp multiprec.hpp Interval.hpp TransformTools.hpp trace_curve_mpfr.cpp
	$(CC) $(CFLAGS3) trace_curve_mpfr trace_curve_mpfr.cpp

interval_test: LinAlgTools.hpp Interval.hpp multiprec.hpp interval_test.cpp
	$(CC) $(CFLAGSTESTING) interval_test interval_test.cpp

fft_test: Interval.hpp multiprec.hpp ValidatedKAM.hpp fft_test.cpp
	$(CC) $(CFLAGS3) fft_test fft_test.cpp

validate_curve: Interval.hpp multiprec.hpp TransformTools.hpp validate_curve.cpp
	$(CC) $(CFLAGS3) validate_curve validate_curve.cpp

validate_curves: Interval.hpp multiprec.hpp TransformTools.hpp validate_curves.cpp
	$(CC) $(CFLAGS3) validate_curves validate_curves.cpp

qft_test: Interval.hpp multiprec.hpp ValidatedKAM.hpp qft_test.cpp
	$(CC) $(CFLAGSTESTING) qft_test qft_test.cpp

transform_demo: LinAlgTools.hpp Interval.hpp multiprec.hpp TransformTools.hpp transform_demo.cpp
	$(CC) $(CFLAGS3) transform_demo transform_demo.cpp

debugging: LinAlgTools.hpp Interval.hpp multiprec.hpp TransformTools.hpp debugging.cpp
	$(CC) $(CFLAGSTESTING) debugging debugging.cpp


qft_debug: LinAlgTools.hpp Interval.hpp multiprec.hpp TransformTools.hpp qft_debug.cpp
	$(CC) $(CFLAGS3) qft_debug qft_debug.cpp

