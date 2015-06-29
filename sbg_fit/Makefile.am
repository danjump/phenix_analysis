AUTOMAKE_OPTIONS = foreign

INCLUDES = \
  -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include -I/opt/phenix/include -I/opt/phenix/include/odbc++

EXLIBS = -L/opt/phenix/lib -lodbc++ -L$(OFFLINE_MAIN)/lib -lphool

LDADD =  $(EXLIBS)

lib_LTLIBRARIES = \
  libsbg_fit.la

AM_LDFLAGS = -L$(libdir) -L$(OFFLINE_MAIN)/lib `root-config --evelibs` -lRooFitCore -lMinuit -lFoam

AM_CXXFLAGS = `root-config --cflags`

libsbg_fit_la_LDFLAGS = \
   -L$(OFFLINE_MAIN)/lib -lFROG \
   `root-config --libs`

libsbg_fit_la_SOURCES = \
  sbg_fit.C \
  sbg_fit_Dict.C 

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.C
testexternals_LDADD = libsbg_fit.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

# rules to make root dictionaries
sbg_fit_Dict.C: sbg_fit.h sbg_fitLinkDef.h
	rootcint -f $@ -c $(CINTFLAGS) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict* $(BUILT_SOURCES)

