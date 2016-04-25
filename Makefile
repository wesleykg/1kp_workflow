# all:

clean: 
	rm -rf wanted_assembly_URLs.txt

%.fa: data/wanted_species.txt
	python scripts/build_url.py data/wanted_species.txt
	scripts/download.sh 

.PHONY: all clean
.DELETE_ON_ERROR: