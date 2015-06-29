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
  libasymmetry_calculation.la

AM_LDFLAGS = -L$(libdir) -L$(OFFLINE_MAIN)/lib `root-config --evelibs` -lRooFitCore -lMinuit -lFoam

libasymmetry_calculation_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lphool \
  -L/opt/phenix/lib \
  -luspin \
  -lodbc \
  -lodbc++

noinst_HEADERS = \
  W2eGetPol.h


libasymmetry_calculation_la_SOURCES = \
  asymmetry_calculation.C \
  asymmetry_calculation_Dict.C \
  W2eGetPol.C 

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

asymmetry_calculation_Dict.C: asymmetry_calculation.h W2eGetPol.h asymmetry_calculationLinkDef.h
	rootcint -f $@ -c $(CINTFLAGS) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict.*
