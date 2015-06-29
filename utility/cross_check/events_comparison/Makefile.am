AUTOMAKE_OPTIONS = foreign

# these are flags sent to the C++ compiler
# -Wall : I don't know what this does
# -Werror: this will treat undefined and unused variables as errors
AM_CXXFLAGS = -Wall -Werror

# This flag will build 'testexternals.C' as a test for the build environment
BUILT_SOURCES = \
  testexternals.C

# This is where we set the directories where we put our neccessary header files
AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libSimpleClass.la

noinst_HEADERS = \
	dontIncludeMe.h

libSimpleClass_la_LDFLAGS = \
  -L$(ROOTSYS)/lib `root-config --libs`

AM_LDFLAGS = \
  -L$(OFFLINE_MAIN)/lib \
  -L$(ROOTSYS)/lib \
  -L$(libdir)

# These sources are compiled into the final library
libSimpleClass_la_SOURCES = \
  SimpleClass.C \
  events_comparison.C \
  runevtobj.C \
  DictOutput.cxx

noinst_PROGRAMS = \
  testexternals

testexternals_LDADD = \
  libSimpleClass.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

# This is where the dictionary file is generated
DictOutput.cxx: SimpleClass.h events_comparison.h runevtobj.h SimpleClassLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f *Dict.*
