#! /bin/csh
source /opt/phenix/bin/phenix_setup.csh;
rm -r /direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/asymmetry_calculation/spin_db_access/build
rm -r autom4te.cache/ clocal.m4  config.guess config.sub configure depcomp install-sh missing temp/ ltmain.sh

mkdir -p /direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/asymmetry_calculation/spin_db_access/build
mkdir -p /direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/asymmetry_calculation/spin_db_access/install

cd /direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/asymmetry_calculation/spin_db_access/build;
/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/asymmetry_calculation/spin_db_access/autogen.sh \
    --prefix=/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/asymmetry_calculation/spin_db_access/install;
make;
make install;

