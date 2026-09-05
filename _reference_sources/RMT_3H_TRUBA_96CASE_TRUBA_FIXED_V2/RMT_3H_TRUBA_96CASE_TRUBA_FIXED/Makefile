.PHONY: all build manifest verify smoke status summary scaling clean

all: manifest build verify

manifest:
	./rebuild_manifest.sh

build:
	./build.sh

verify:
	./verify_campaign.sh

smoke:
	./smoke_test.sh

status:
	./check_campaign_status.sh

summary:
	python3 collect_summary.py

scaling:
	python3 analyze_scaling.py

clean:
	./clean_results.sh
