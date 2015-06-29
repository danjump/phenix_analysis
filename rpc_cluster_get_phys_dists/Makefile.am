AUTOMAKE_OPTIONS = foreign

INCLUDES = \
  -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include -I/opt/phenix/include -I/opt/phenix/include/odbc++

EXLIBS = -L/opt/phenix/lib -lodbc++ -L$(OFFLINE_MAIN)/lib -lphool

LDADD =  $(EXLIBS)

lib_LTLIBRARIES = \
  librpc_cluster_get_phys_dists.la

AM_LDFLAGS = -L$(libdir) -L$(OFFLINE_MAIN)/lib

librpc_cluster_get_phys_dists_la_LDFLAGS = \
   -L$(OFFLINE_MAIN)/lib -lFROG \
   `root-config --libs`

librpc_cluster_get_phys_dists_la_SOURCES = \
  rpc_cluster_get_phys_dists.C \
  rpc_cluster_get_phys_dists_Dict.C 

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.C
testexternals_LDADD = librpc_cluster_get_phys_dists.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

# rules to make root dictionaries
rpc_cluster_get_phys_dists_Dict.C: rpc_cluster_get_phys_dists.h rpc_cluster_get_phys_distsLinkDef.h
	rootcint -f $@ -c $(CINTFLAGS) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict* $(BUILT_SOURCES)

