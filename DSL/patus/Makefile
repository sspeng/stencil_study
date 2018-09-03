
prepare:
	mkdir -p target


heat1d_2oa_exp: heat1d_2oa_exp.stc prepare
	rm -rf target/heat1d_2oa_exp
	patus --outdir=target/heat1d_2oa_exp --architecture=x86_64 --strategy=configs/simple.stg --create-validation=no heat1d_2oa_exp.stc

heat2d_2oa_exp: heat2d_2oa_exp.stc prepare
	rm -rf target/heat2d_2oa_exp
	patus --outdir=target/heat2d_2oa_exp --architecture=x86_64  --create-validation=no heat2d_2oa_exp.stc

heat3d_2oa_exp: heat3d_2oa_exp.stc prepare
	rm -rf target/heat3d_2oa_exp
	patus --outdir=target/heat3d_2oa_exp --architecture=x86_64 --create-validation=no heat3d_2oa_exp.stc

wave1d_2oa_exp: wave1d_2oa_exp.stc prepare
	rm -rf target/wave1d_2oa_exp
	patus --outdir=target/wave1d_2oa_exp --architecture=x86_64 --strategy=configs/simple.stg --create-validation=no wave1d_2oa_exp.stc

wave2d_2oa_exp: wave2d_2oa_exp.stc prepare
	rm -rf target/wave2d_2oa_exp
	patus --outdir=target/wave2d_2oa_exp --architecture=x86_64 --create-validation=no wave2d_2oa_exp.stc

wave3d_2oa_exp: wave3d_2oa_exp.stc prepare
	rm -rf target/wave3d_2oa_exp
	patus --outdir=target/wave3d_2oa_exp --architecture=x86_64 --create-validation=no wave3d_2oa_exp.stc

lapla1d_2oa_exp: lapla1d_2oa_exp.stc prepare
	rm -rf target/lapla1d_2oa_exp
	patus --outdir=target/lapla1d_2oa_exp --architecture=x86_64 --strategy=configs/simple.stg --create-validation=no lapla1d_2oa_exp.stc

lapla2d_2oa_exp: lapla2d_2oa_exp.stc prepare
	rm -rf target/lapla2d_2oa_exp
	patus --outdir=target/lapla2d_2oa_exp --architecture=x86_64 --create-validation=no lapla2d_2oa_exp.stc

lapla3d_2oa_exp: lapla3d_2oa_exp.stc prepare
	rm -rf target/lapla3d_2oa_exp
	patus --outdir=target/lapla3d_2oa_exp --architecture=x86_64 --create-validation=no lapla3d_2oa_exp.stc

heat2oa_exp: heat1d_2oa_exp heat2d_2oa_exp heat3d_2oa_exp
wave2oa_exp: wave1d_2oa_exp wave2d_2oa_exp wave3d_2oa_exp
lapla2oa_exp: lapla1d_2oa_exp lapla2d_2oa_exp lapla3d_2oa_exp

2oa_exp: heat2oa_exp wave2oa_exp lapla2oa_exp

all: 2oa_exp

clean:
	rm -rf target/*

