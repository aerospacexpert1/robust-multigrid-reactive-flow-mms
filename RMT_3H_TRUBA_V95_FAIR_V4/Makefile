.PHONY: all build smoke manifest clean
all: build
build:
	./build.sh
smoke: build
	./smoke_local.sh
manifest:
	./rebuild_manifest.sh
clean:
	rm -rf build smoke_output results campaign_manifest.csv
