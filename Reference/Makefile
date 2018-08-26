CC=gcc
CFLAGS=-O3
LDFLAGS=-lm


heat_2oa_exp: heat1d_2oa_exp heat2d_2oa_exp heat3d_2oa_exp
wave_2oa_exp: wave1d_2oa_exp wave2d_2oa_exp wave3d_2oa_exp
lapla_2oa_exp: lapla1d_2oa_exp lapla2d_2oa_exp lapla3d_2oa_exp

all: heat_2oa_exp wave_2oa_exp lapla_2oa_exp

prepare:
	mkdir -p target

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

clean:
	rm -rf target

