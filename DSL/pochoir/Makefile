#!/bin/bash
default: apop heat_1D_NP
CC = pochoir
ICC = icpc
#OPT_FLAGS = -O3 -DNDEBUG -Wall -Werror -funroll-loops -xHOST -fno-alias -fno-fnalias -fp-model precise -std=c++11
#OPT_FLAGS = -O3 -DNDEBUG -Wall -Werror -funroll-loops -xHOST -fno-alias -fno-fnalias -fp-model precise -std=c++0x
OPT_FLAGS = -O3 -g -Wall -Werror -funroll-loops -xHOST -fno-alias -fno-fnalias -fp-model precise -std=c++0x
#POCHOIR_DEBUG_FLAGS = -O0 -g3 -DDEBUG -debug -Wall -Werror -std=c++11
POCHOIR_DEBUG_FLAGS = -O0 -g3 -DDEBUG -debug -Wall -Werror -std=c++0x
#ICC_DEBUG_FLAGS = -O0 -g3 -DDEBUG -Wall -Werror -std=c++0x -I${INTEL_CILK_HEADER} -include cilk_stub.h  
ICC_DEBUG_FLAGS = -O3 -Wall -Werror -std=c++0x -I${INTEL_CILK_HEADER} 
#ICC_DEBUG_FLAGS = -O3 -Wall -Werror -std=c++0x -I${INTEL_CILK_HEADER} -include cilk_stub.h  

ICC_IGNORE = -wd1292 -wd780 -wd488 -wd161 -wd3286 -wd501 -wd14 -wd3699 -wd3946 -wd3947 

heat1d_2oa_exp : heat1d_2oa_exp.cpp
	${CC} -o qheat_1D_NP ${OPT_FLAGS} heat1d_2oa_exp.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o xheat1d_2oa_exp ${POCHOIR_DEBUG_FLAGS} heat1d_2oa_exp.cpp
#   Phase-II compilation
#	${CC} -o qheat_1D_NP ${OPT_FLAGS} heat1d_2oa_exp.cpp
#	${ICC} -o qheat_1D_NP_gdb ${ICC_DEBUG_FLAGS} ${ICC_IGNORE} heat1d_2oa_exp_pochoir.cpp
#	Phase-I compilation with debugging aid
#  ${CC} -o xheat1d_2oa_exp ${POCHOIR_DEBUG_FLAGS} heat1d_2oa_exp.cpp
#

heat2d_2oa_exp : heat2d_2oa_exp.cpp
	${CC} -o xheat_2d_exp ${OPT_FLAGS} heat2d_2oa_exp.cpp

heat3d_2oa_exp : heat3d_2oa_exp.cpp
	${CC} -o xheat_3d_exp ${OPT_FLAGS} heat3d_2oa_exp.cpp

clean: 
	rm -f *.o *.i *_pochoir *_gdb *_pochoir.cpp *.out
