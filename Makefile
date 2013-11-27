all: vcfneuraltrain vcfneuralapply

FANNDIR=fann
FANNOBJ=$(FANNDIR)/build/src/CMakeFiles/doublefann.dir/doublefann.c.o
CFLAGS=-O3

clean:
	rm -f vcfneuraltrain vcfneuralapply
	cd vcflib && $(MAKE) clean
	cd $(FANNDIR)/build && $(MAKE) clean

vcflib/Variant.o:
	cd vcflib && $(MAKE)

$(FANNOBJ):
	cd $(FANNDIR) && mkdir -p build && cd build && cmake .. && $(MAKE)

vcfneuraltrain: vcflib/Variant.o vcfneuraltrain.cpp $(FANNOBJ)
	g++ -I$(FANNDIR)/src/include $(FANNOBJ) vcflib/Variant.o vcflib/split.o vcflib/smithwaterman/SmithWatermanGotoh.o vcflib/smithwaterman/LeftAlign.o vcflib/smithwaterman/disorder.c vcflib/smithwaterman/Repeats.o vcflib/smithwaterman/IndelAllele.o vcflib/tabixpp/tabix.o vcflib/tabixpp/bgzf.o vcfneuraltrain.cpp -o vcfneuraltrain -lm -lz -Lvcflib/ -Lvcflib/tabixpp/ -ltabix $(CFLAGS)

vcfneuralapply: vcflib/Variant.o vcfneuralapply.cpp $(FANNOBJ)
	g++ -I$(FANNDIR)/src/include $(FANNOBJ) vcflib/Variant.o vcflib/split.o vcflib/smithwaterman/SmithWatermanGotoh.o vcflib/smithwaterman/LeftAlign.o vcflib/smithwaterman/disorder.c vcflib/smithwaterman/Repeats.o vcflib/smithwaterman/IndelAllele.o vcflib/tabixpp/tabix.o vcflib/tabixpp/bgzf.o vcfneuralapply.cpp -o vcfneuralapply -lm -lz -Lvcflib/ -Lvcflib/tabixpp/ -ltabix $(CFLAGS)
