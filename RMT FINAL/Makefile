.PHONY: all build smoke first-pressure manifest clean
all: build
build:
	./build.sh
smoke: build
	./smoke_local.sh
first-pressure: build
	./tests/first_pressure_matrix.sh
manifest:
	./rebuild_manifest.sh
clean:
	rm -rf build smoke_output first_pressure_logs results campaign_manifest.csv RMT_FINAL_campaign_summary.csv RMT_FINAL_scaling_summary.csv
