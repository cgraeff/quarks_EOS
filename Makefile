SHELL := /bin/bash # Use bash as shell

all:
	cd src; make
run:
	./eos $(ARGS)
graph:
	for dir in `echo output/*/`; do cd "$$dir" && gnuplot gnuplot.gpi && cd ../..; done
tests:
	./eos -t $(ARGS)
tgraph:
	for dir in `echo tests/*/`; do cd "$$dir" && gnuplot gnuplot.gpi && cd ../..; done
clean:
	-rm -f eos
	cd src; make clean
	find output -name "*.dat" -type f -delete
	find output -name "*.log" -type f -delete
	find output -name "*.png" -type f -delete
	find output -name "*.tex" -type f -delete
	find tests -name "*.dat" -type f -delete
	find tests -name "*.log" -type f -delete
	find tests -name "*.png" -type f -delete
	find tests -name "*.tex" -type f -delete
