.PHONY: all manifest verify build selftest one summary clean

all: manifest verify build

manifest:
	./rebuild_manifest.sh

verify:
	./verify_campaign.sh

build:
	./build.sh

selftest:
	./build.sh

one:
	@echo "Use: sbatch --export=ALL,TASK_ID=<0..95> run_one_case.slurm"

summary:
	python3 tools/collect_summary.py

clean:
	rm -rf build results logs RMT_BOOKFIX_V3_campaign_summary.csv
