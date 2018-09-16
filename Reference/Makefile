CC=icc
CFLAGS=-O3
#LDFLAGS=-lm -fopenmp
LDFLAGS=-lm -qopenmp
#LDFLAGS=-lm 
REF_TEST_DIR=#!date +%Y_%m_%d:%H_%M


heat_2oa_exp: heat1d_2oa_exp heat2d_2oa_exp heat3d_2oa_exp
wave_2oa_exp: wave1d_2oa_exp wave2d_2oa_exp wave3d_2oa_exp
lapla_2oa_exp: lapla1d_2oa_exp lapla2d_2oa_exp lapla3d_2oa_exp

heat_4oa_exp: heat1d_4oa_exp heat2d_4oa_exp heat3d_4oa_exp
wave_4oa_exp: wave1d_4oa_exp wave2d_4oa_exp wave3d_4oa_exp

2oa_exp: heat_2oa_exp wave_2oa_exp lapla_2oa_exp
4oa_exp: heat_4oa_exp wave_4oa_exp

all: 2oa_exp 4oa_exp

prepare:
	mkdir -p target

prep_test:
	echo ${REF_TEST_DIR}
	#mkdir -p target/$(REF_TEST_DIR)

###########################################
### 2nd order accurate, explicit methods###
###########################################
heat1d_2oa_exp: heat1d_2oa_exp.c prepare
	rm -rf target/xheat1d_2oa_exp
	$(CC) $(CFLAGS) heat1d_2oa_exp.c $(LDFLAGS) -o target/xheat1d_2oa_exp

heat2d_2oa_exp: heat2d_2oa_exp.c prepare
	rm -rf target/xheat2d_2oa_exp
	$(CC) $(CFLAGS) heat2d_2oa_exp.c $(LDFLAGS)  -o target/xheat2d_2oa_exp

heat3d_2oa_exp: heat3d_2oa_exp.c prepare
	rm -rf target/xheat3d_2oa_exp
	$(CC) $(CFLAGS) heat3d_2oa_exp.c $(LDFLAGS)  -o target/xheat3d_2oa_exp

wave1d_2oa_exp: wave1d_2oa_exp.c prepare
	rm -rf target/xwave1d_2oa_exp
	$(CC) $(CFLAGS) wave1d_2oa_exp.c $(LDFLAGS)  -o target/xwave1d_2oa_exp

wave2d_2oa_exp: wave2d_2oa_exp.c prepare
	rm -rf target/xwave2d_2oa_exp
	$(CC) $(CFLAGS) wave2d_2oa_exp.c $(LDFLAGS)  -o target/xwave2d_2oa_exp

wave3d_2oa_exp: wave3d_2oa_exp.c prepare
	rm -rf target/xwave3d_2oa_exp
	$(CC) $(CFLAGS) wave3d_2oa_exp.c $(LDFLAGS)  -o target/xwave3d_2oa_exp

lapla1d_2oa_exp: lapla1d_2oa_exp.c prepare
	rm -rf target/xlapla1d_2oa_exp
	$(CC) $(CFLAGS) lapla1d_2oa_exp.c $(LDFLAGS)  -o target/xlapla1d_2oa_exp

lapla2d_2oa_exp: lapla2d_2oa_exp.c prepare
	rm -rf target/xlapla2d_2oa_exp
	$(CC) $(CFLAGS) lapla2d_2oa_exp.c $(LDFLAGS)  -o target/xlapla2d_2oa_exp

lapla3d_2oa_exp: lapla3d_2oa_exp.c prepare
	rm -rf target/xlapla3d_2oa_exp
	$(CC) $(CFLAGS) lapla3d_2oa_exp.c $(LDFLAGS)  -o target/xlapla3d_2oa_exp

###########################################
### 4th order accurate, explicit methods###
###########################################
heat1d_4oa_exp: heat1d_4oa_exp.c prepare
	rm -rf target/xheat1d_4oa_exp
	$(CC) $(CFLAGS) heat1d_4oa_exp.c $(LDFLAGS) -o target/xheat1d_4oa_exp

heat2d_4oa_exp: heat2d_4oa_exp.c prepare
	rm -rf target/xheat2d_4oa_exp
	$(CC) $(CFLAGS) heat2d_4oa_exp.c $(LDFLAGS) -o target/xheat2d_4oa_exp

heat3d_4oa_exp: heat3d_4oa_exp.c prepare
	rm -rf target/xheat3d_4oa_exp
	$(CC) $(CFLAGS) heat3d_4oa_exp.c $(LDFLAGS) -o target/xheat3d_4oa_exp

wave1d_4oa_exp: wave1d_4oa_exp.c prepare
	rm -rf target/xwave1d_4oa_exp
	$(CC) $(CFLAGS) wave1d_4oa_exp.c $(LDFLAGS) -o target/xwave1d_4oa_exp

wave2d_4oa_exp: wave2d_4oa_exp.c prepare
	rm -rf target/xwave2d_4oa_exp
	$(CC) $(CFLAGS) wave2d_4oa_exp.c $(LDFLAGS) -o target/xwave2d_4oa_exp

wave3d_4oa_exp: wave3d_4oa_exp.c prepare
	rm -rf target/xwave3d_4oa_exp
	$(CC) $(CFLAGS) wave3d_4oa_exp.c $(LDFLAGS) -o target/xwave3d_4oa_exp



clean:
	rm -rf target

