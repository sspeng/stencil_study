CC=gcc
CFLAGS=-O3

all: heat1d_2oa_exp heat2d_2oa_exp

prepare:
	mkdir -p target

heat1d_2oa_exp: heat1d_2oa_exp.c prepare
	rm -rf target/xheat1d_2oa_exp
	$(CC) $(FLAGS) heat1d_2oa_exp.c -o target/xheat1d_2oa_exp

heat2d_2oa_exp: heat2d_2oa_exp.c prepare
	rm -rf target/xheat2d_2oa_exp
	$(CC) $(CFLAGS) heat2d_2oa_exp.c -o target/xheat2d_2oa_exp



clean:
	rm -rf target

