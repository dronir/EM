
BINPATH=../../build/bin

all: medium_test_uniform.nc medium_test_constant.nc

medium_test_uniform.nc:
	$(BINPATH)/medgen test_uniform.in

medium_test_constant.nc:
	$(BINPATH)/medgen test_constant.in

clean:
	rm -rf medium_test_constant.nc
	rm -rf medium_test_uniform.nc
