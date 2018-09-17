CC=gcc
CFLAGS=-O3
LDFLAGS=-lm -fopenmp
#LDFLAGS=-lm -qopenmp
REF_TEST_DIR=#!date +%Y_%m_%d:%H_%M


heat_2oa_exp: heat1d_2oa_exp heat2d_2oa_exp heat3d_2oa_exp
wave_2oa_exp: wave1d_2oa_exp wave2d_2oa_exp wave3d_2oa_exp
lapla_2oa_exp: lapla1d_2oa_exp lapla2d_2oa_exp lapla3d_2oa_exp

heat_4oa_exp: heat1d_4oa_exp heat2d_4oa_exp heat3d_4oa_exp
wave_4oa_exp: wave1d_4oa_exp wave2d_4oa_exp wave3d_4oa_exp
lapla_4oa_exp: lapla1d_4oa_exp lapla2d_4oa_exp lapla3d_4oa_exp

2oa_exp: heat_2oa_exp wave_2oa_exp lapla_2oa_exp
4oa_exp: heat_4oa_exp wave_4oa_exp lapla_4oa_exp

all: 2oa_exp 4oa_exp

prepare:
	mkdir -p target

clean:
	rm -rf target


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

lapla1d_4oa_exp: lapla1d_4oa_exp.c prepare
	rm -rf target/xlapla1d_4oa_exp
	$(CC) $(CFLAGS) lapla1d_4oa_exp.c $(LDFLAGS) -o target/xlapla1d_4oa_exp

lapla2d_4oa_exp: lapla2d_4oa_exp.c prepare
	rm -rf target/xlapla2d_4oa_exp
	$(CC) $(CFLAGS) lapla2d_4oa_exp.c $(LDFLAGS) -o target/xlapla2d_4oa_exp

lapla3d_4oa_exp: lapla3d_4oa_exp.c prepare
	rm -rf target/xlapla3d_4oa_exp
	$(CC) $(CFLAGS) lapla3d_4oa_exp.c $(LDFLAGS) -o target/xlapla3d_4oa_exp



###########################################
### Testing rules #########################
###########################################

T_MAX=100
1d_params=10000000 0 0 $(T_MAX)
2d_params=10000 10000 0 $(T_MAX)
3d_params=256 256 256 $(T_MAX)

test_all: test_heat test_wave test_lapla

test_heat_2oa: test_heat1d_2oa_exp test_heat2d_2oa_exp test_heat3d_2oa_exp
test_heat_4oa: test_heat1d_4oa_exp test_heat2d_4oa_exp test_heat3d_4oa_exp

test_wave_2oa: test_wave1d_2oa_exp test_wave2d_2oa_exp test_wave3d_2oa_exp
test_wave_4oa: test_wave1d_4oa_exp test_wave2d_4oa_exp test_wave3d_4oa_exp

test_lapla_2oa: test_lapla1d_2oa_exp test_lapla2d_2oa_exp test_lapla3d_2oa_exp
test_lapla_4oa: test_lapla1d_4oa_exp test_lapla2d_4oa_exp test_lapla3d_4oa_exp

test_heat: test_heat_2oa test_heat_4oa
test_wave: test_wave_2oa test_wave_4oa
test_lapla: test_lapla_2oa test_lapla_4oa

test_heat1d_2oa_exp: heat1d_2oa_exp
	  ./target/xheat1d_2oa_exp $(1d_params)
test_heat2d_2oa_exp: heat2d_2oa_exp
	  ./target/xheat2d_2oa_exp $(2d_params)
test_heat3d_2oa_exp: heat3d_2oa_exp
	  ./target/xheat3d_2oa_exp $(3d_params)
test_heat1d_4oa_exp: heat1d_4oa_exp
	  ./target/xheat1d_4oa_exp $(1d_params)
test_heat2d_4oa_exp: heat2d_4oa_exp
	  ./target/xheat2d_4oa_exp $(2d_params)
test_heat3d_4oa_exp: heat3d_4oa_exp
	  ./target/xheat3d_4oa_exp $(3d_params)
test_wave1d_2oa_exp: wave1d_2oa_exp
	  ./target/xwave1d_2oa_exp $(1d_params)
test_wave2d_2oa_exp: wave2d_2oa_exp
	  ./target/xwave2d_2oa_exp $(2d_params)
test_wave3d_2oa_exp: wave3d_2oa_exp
	  ./target/xwave3d_2oa_exp $(3d_params)
test_wave1d_4oa_exp: wave1d_4oa_exp
	  ./target/xwave1d_4oa_exp $(1d_params)
test_wave2d_4oa_exp: wave2d_4oa_exp
	  ./target/xwave2d_4oa_exp $(2d_params)
test_wave3d_4oa_exp: wave3d_4oa_exp
	  ./target/xwave3d_4oa_exp $(3d_params)
test_lapla1d_2oa_exp: lapla1d_2oa_exp
	  ./target/xlapla1d_2oa_exp $(1d_params)
test_lapla2d_2oa_exp: lapla2d_2oa_exp
	  ./target/xlapla2d_2oa_exp $(2d_params)
test_lapla3d_2oa_exp: lapla3d_2oa_exp
	  ./target/xlapla3d_2oa_exp $(3d_params)
test_lapla1d_4oa_exp: lapla1d_4oa_exp
	  ./target/xlapla1d_4oa_exp $(1d_params)
test_lapla2d_4oa_exp: lapla2d_4oa_exp
	  ./target/xlapla2d_4oa_exp $(2d_params)
test_lapla3d_4oa_exp: lapla3d_4oa_exp
	  ./target/xlapla3d_4oa_exp $(3d_params)
