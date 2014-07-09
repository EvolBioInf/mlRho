#!/bin/make

# Global Makefile

.PHONY: all
all:
	cd src && $(MAKE) all;

.PHONY: clean
clean:
	cd src && $(MAKE) clean;

# FIXME
doc:
	cd ../Doc; make clean; make pdf;

# FIXME
tarfile:
	cd ../Doc; make clean; make pdf; cd ../$(DIRECTORY)_$(VERSION)
	mkdir $(DIRECTORY)_$(VERSION)
	cp -rf $(SRCFILES) valgrind.sh *.h  Scripts Makefile README COPYING \
	../Doc/mlRhoDoc.pdf $(DIRECTORY)_$(VERSION)
	tar cvzfh $(EXECFILE)_$(VERSION).tgz $(DIRECTORY)_$(VERSION)
	mv $(EXECFILE)_$(VERSION).tgz ../
	/bin/rm -rf $(DIRECTORY)_$(VERSION)

