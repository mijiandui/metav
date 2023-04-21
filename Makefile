all:build runtime

build:
	./build.sh


runtime:
	conda create -n vs2 virsorter=2 checkv
	conda create --name dvf python=3.6 numpy theano=1.0.3 keras=2.2.4 scikit-learn Biopython h5py


.PNONY: clean
clean:
	rm -rf third_party_tools/DeepVirFinder third_party_tools/bwa/ third_party_tools/fastp/ third_party_tools/htslib/ \
	       third_party_tools/megahit-1.2.9/ third_party_tools/samtools/ third_party_tools/SPAdes-3.15.5/	
