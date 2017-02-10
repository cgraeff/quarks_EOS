SHELL := /bin/bash # Use bash as shell
TARGET = qeos

# List set for multirun
MULTIRUN_SETS = Buballa_1 Buballa_2 Buballa_3 \
		BuballaR_2 \
		BuballaR_2_GV BuballaR_2_GV_0.25 BuballaR_2_GV_0.35 BuballaR_2_GV_0.45 \
		BuballaR_2_GV_0.55 BuballaR_2_GV_0.65 BuballaR_2_GV_0.75 \
		PCP-0.0 PCP-0.1 PCP-0.2 PCP-0.3 PCP-0.4 PCP-0.5

.PHONY: all run graph tests tgraph clean

all:
	@echo "[Compiling...]"
	@cd src; make
	@echo "[done.]"
run:
	@./$(TARGET) -d $(ARGS)
graph:
	@echo "[Plotting...]"
	@cd output; \
	for dir in `echo */`; do \
		echo "$$dir"; \
		cd "$$dir"; \
		gnuplot gnuplot.gpi; \
		cd ..; \
	done
	@echo "[done.]"
multirun:
	@echo "[Running for multiple parameterizations...]"
	@for key in $(MULTIRUN_SETS); do \
		./$(TARGET) -d -p "$$key" $(ARGS); \
		if [ -d multioutput/"$$key" ]; then rm -r multioutput/"$$key"; fi; \
		cp -r output multioutput/"$$key"; \
	done
	@echo "[done.]"
mgraph:
	@echo "[Plotting for multiple parameterizations...]"
	@for dir in $(MULTIRUN_SETS); do \
		echo "$$dir"; \
		cd "multioutput/$$dir"; \
		for subdir in `echo */`; do \
			cd "$$subdir"; \
			gnuplot gnuplot.gpi; \
			cd ..; \
		done; \
		cd ../..; \
	done; \
	cd multioutput; gnuplot gnuplot.gpi
	@echo "[done.]"
tests:
	@echo "[Running tests...]"
	@./$(TARGET) -a $(ARGS)
	@echo "[done.]"
tgraph:
	@echo "[Plotting tests ...]"
	@cd tests/ ; \
	for dir in `echo */`; do \
		cd "$$dir"; \
		echo "$$dir"; \
		gnuplot gnuplot.gpi; \
		cd ../; \
	done;
	@echo "[done.]"
clean:
	@echo "[Cleaning...]"
	@-rm -f $(TARGET)
	@cd src; make clean
	@find . -name "*.dat" -type f -delete
	@find . -name "*.log" -type f -delete
	@find . -name "*.png" -type f -delete
	@find . -name "*.tex" -type f -delete
	@cd multioutput; rm -rf $(MULTIRUN_SETS)
	@echo "[done.]"
