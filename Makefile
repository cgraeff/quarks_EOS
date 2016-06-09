all:
	cd src; make
run:
	mkdir -p tests/data
	./eos $(ARGS)
graph:
	cd data; make
	cd tests; make
clean:
	-rm -f eos
	cd src; make clean
	cd data; make clean
	cd tests; make clean
