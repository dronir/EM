
export BINPATH=$(realpath ../test_build/bin)

.PHONY: default
default: all

.PHONY: all
all: test-medgen test-vScatter #test-fit_2hg test-geomScatter

.PHONY: test-medgen
test-medgen:
	$(MAKE) -e -C medgen

test-vScatter:
	$(MAKE) -e -C vScatter

clean:
	$(MAKE) -e -C medgen clean
	$(MAKE) -e -C vScatter clean
