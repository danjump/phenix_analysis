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
  libW2eGetPol.la

libW2eGetPol_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L/opt/phenix/lib \
  -luspin \
  -lodbc \
  -lodbc++

noinst_HEADERS = \
  W2eGetPol.h


libW2eGetPol_la_SOURCES = \
  W2eGetPol_Dict.C \
  W2eGetPol.C 

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

W2eGetPol_Dict.C: W2eGetPol.h W2eGetPolLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict.*
