all: vcfneuraltrain vcfneuralapply

clean:
	rm vcfneuraltrain
	cd vcflib && $(MAKE) clean

vcflib/Variant.o:
	cd vcflib && $(MAKE)

vcfneuraltrain: vcflib/Variant.o vcfneuraltrain.cpp
	g++ -I/usr/include/fann/ vcflib/Variant.o vcflib/split.o vcflib/smithwaterman/SmithWatermanGotoh.o vcflib/tabixpp/tabix.o vcflib/tabixpp/bgzf.o vcfneuraltrain.cpp -o vcfneuraltrain -lfann -lm -lz -Lvcflib/ -Lvcflib/tabixpp/ -ltabix

vcfneuralapply: vcflib/Variant.o vcfneuralapply.cpp
	g++ -I/usr/include/fann/ vcflib/Variant.o vcflib/split.o vcflib/smithwaterman/SmithWatermanGotoh.o vcflib/tabixpp/tabix.o vcflib/tabixpp/bgzf.o vcfneuralapply.cpp -o vcfneuralapply -lfann -lm -lz -Lvcflib/ -Lvcflib/tabixpp/ -ltabix
