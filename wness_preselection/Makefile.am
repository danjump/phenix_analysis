AUTOMAKE_OPTIONS = foreign

INCLUDES = \
  -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include -I/opt/phenix/include -I/opt/phenix/include/odbc++

EXLIBS = -L/opt/phenix/lib -lodbc++ -L$(OFFLINE_MAIN)/lib -lphool

LDADD =  $(EXLIBS)

lib_LTLIBRARIES = \
  libwness_preselection.la

AM_LDFLAGS = -L$(libdir) -L$(OFFLINE_MAIN)/lib `root-config --evelibs` -lRooFitCore -lMinuit -lFoam

AM_CXXFLAGS = `root-config --cflags`

libwness_preselection_la_LDFLAGS = \
   -L$(OFFLINE_MAIN)/lib -lFROG \
   `root-config --libs`

libwness_preselection_la_SOURCES = \
  wness_preselection.C \
  wness_preselection_Dict.C 

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.C
testexternals_LDADD = libwness_preselection.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

# rules to make root dictionaries
wness_preselection_Dict.C: wness_preselection.h wness_preselectionLinkDef.h
	rootcint -f $@ -c $(CINTFLAGS) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict* $(BUILT_SOURCES)

