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

prepare:
	mkdir -p target

heat_2oa_exp: heat2d_2oa_exp heat3d_2oa_exp
heat_4oa_exp: heat2d_4oa_exp heat3d_4oa_exp

wave_2oa_exp: wave2d_2oa_exp wave3d_2oa_exp
wave_4oa_exp: wave2d_4oa_exp wave3d_4oa_exp

lapla_2oa_exp: lapla2d_2oa_exp lapla3d_2oa_exp
lapla_4oa_exp: lapla2d_4oa_exp lapla3d_4oa_exp

2oa_exp: heat_2oa_exp wave_2oa_exp lapla_2oa_exp
4oa_exp: heat_4oa_exp wave_4oa_exp lapla_4oa_exp

all: 2oa_exp 4oa_exp

clean: 
	rm -rf target
	rm -f *.o *.i *_pochoir *_gdb *_pochoir.cpp *.out

#####################################
##### 2nd Order Accurate Builds
#####################################

heat2d_2oa_exp : heat2d_2oa_exp.cpp prepare
	rm -rf target/xheat2d_2oa_exp
	${CC} -o target/xheat2d_2oa_exp ${OPT_FLAGS} heat2d_2oa_exp.cpp

heat3d_2oa_exp : heat3d_2oa_exp.cpp prepare
	rm -rf target/xheat3d_2oa_exp
	${CC} -o target/xheat3d_2oa_exp ${OPT_FLAGS} heat3d_2oa_exp.cpp

wave2d_2oa_exp : wave2d_2oa_exp.cpp prepare
	rm -rf target/xwave2d_2oa_exp
	${CC} -o target/xwave2d_2oa_exp ${OPT_FLAGS} wave2d_2oa_exp.cpp

wave3d_2oa_exp : wave3d_2oa_exp.cpp prepare
	rm -rf target/xwave3d_2oa_exp
	${CC} -o target/xwave3d_2oa_exp ${OPT_FLAGS} wave3d_2oa_exp.cpp

lapla2d_2oa_exp : lapla2d_2oa_exp.cpp prepare
	rm -rf target/xlapla2d_2oa_exp
	${CC} -o target/xlapla2d_2oa_exp ${OPT_FLAGS} lapla2d_2oa_exp.cpp

lapla3d_2oa_exp : lapla3d_2oa_exp.cpp prepare
	rm -rf target/xlapla3d_2oa_exp
	${CC} -o target/xlapla3d_2oa_exp ${OPT_FLAGS} lapla3d_2oa_exp.cpp



#####################################
##### 4th Order Accurate Builds
#####################################

heat2d_4oa_exp : heat2d_4oa_exp.cpp prepare
	rm -rf target/xheat2d_4oa_exp
	${CC} -o target/xheat2d_4oa_exp ${OPT_FLAGS} heat2d_4oa_exp.cpp

heat3d_4oa_exp : heat3d_4oa_exp.cpp prepare
	rm -rf target/xheat3d_4oa_exp
	${CC} -o target/xheat3d_4oa_exp ${OPT_FLAGS} heat3d_4oa_exp.cpp

wave2d_4oa_exp : wave2d_4oa_exp.cpp prepare
	rm -rf target/xwave2d_4oa_exp
	${CC} -o target/xwave2d_4oa_exp ${OPT_FLAGS} wave2d_4oa_exp.cpp

wave3d_4oa_exp : wave3d_4oa_exp.cpp prepare
	rm -rf target/xwave3d_4oa_exp
	${CC} -o target/xwave3d_4oa_exp ${OPT_FLAGS} wave3d_4oa_exp.cpp

lapla2d_4oa_exp : lapla2d_4oa_exp.cpp prepare
	rm -rf target/xlapla2d_4oa_exp
	${CC} -o target/xlapla2d_4oa_exp ${OPT_FLAGS} lapla2d_4oa_exp.cpp

lapla3d_4oa_exp : lapla3d_4oa_exp.cpp prepare
	rm -rf target/xlapla3d_4oa_exp
	${CC} -o target/xlapla3d_4oa_exp ${OPT_FLAGS} lapla3d_4oa_exp.cpp


