all: vcfneuraltrain vcfneuralapply

FANNDIR=fann-2.1.0

clean:
	rm vcfneuraltrain
	cd vcflib && $(MAKE) clean

vcflib/Variant.o:
	cd vcflib && $(MAKE)

$(FANNDIR)/src/doublefann.o: $(FANNDIR)/src/doublefann.c
	cd $(FANNDIR) && configure && $(MAKE)

vcfneuraltrain: vcflib/Variant.o vcfneuraltrain.cpp $(FANNDIR)/src/doublefann.o
	g++ -I$(FANNDIR)/src/include $(FANNDIR)/src/doublefann.o vcflib/Variant.o vcflib/split.o vcflib/smithwaterman/SmithWatermanGotoh.o vcflib/tabixpp/tabix.o vcflib/tabixpp/bgzf.o vcfneuraltrain.cpp -o vcfneuraltrain -lm -lz -Lvcflib/ -Lvcflib/tabixpp/ -ltabix

vcfneuralapply: vcflib/Variant.o vcfneuralapply.cpp $(FANNDIR)/src/doublefann.o
	g++ -I$(FANNDIR)/src/include $(FANNDIR)/src/doublefann.o vcflib/Variant.o vcflib/split.o vcflib/smithwaterman/SmithWatermanGotoh.o vcflib/tabixpp/tabix.o vcflib/tabixpp/bgzf.o vcfneuralapply.cpp -o vcfneuralapply -lm -lz -Lvcflib/ -Lvcflib/tabixpp/ -ltabix
