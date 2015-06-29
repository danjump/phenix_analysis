AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = -O2 -Werror -m32

BUILT_SOURCES = \
  testexternals.C

INCLUDES = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I/opt/phenix/include \
  -I$(ROOTSYS)/include 

lib_LTLIBRARIES = \
  lib_bunch_shuffling.la

AM_LDFLAGS = -L$(libdir) -L$(OFFLINE_MAIN)/lib `root-config --evelibs` -lRooFitCore -lMinuit -lFoam

lib_bunch_shuffling_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lphool \
  -L/opt/phenix/lib \
  -luspin \
  -lodbc \
  -lodbc++

noinst_HEADERS = \
  W2eGetPol.h \
  spin_event.h

lib_bunch_shuffling_la_SOURCES = \
  bunch_shuffling.C \
  bunch_shuffling_Dict.C \
  W2eGetPol.C \
  spin_event.C

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

bunch_shuffling_Dict.C: bunch_shuffling.h W2eGetPol.h spin_event.h bunch_shufflingLinkDef.h
	rootcint -f $@ -c $(CINTFLAGS) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict.*
