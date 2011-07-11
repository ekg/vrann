all: vcftrain

vcflib/Variant.o:
	cd vcflib && $(MAKE)

vcftrain: vcflib/Variant.o vcftrain.cpp
	g++ -I/usr/include/fann/ vcflib/Variant.o vcflib/split.o vcflib/smithwaterman/SmithWatermanGotoh.o vcflib/tabixpp/tabix.o vcflib/tabixpp/bgzf.o vcftrain.cpp -o vcftrain -lfann -lm -lz -Lvcflib/ -Lvcflib/tabixpp/ -ltabix
