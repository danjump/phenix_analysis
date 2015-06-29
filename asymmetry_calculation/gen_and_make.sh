#! /bin/bash
echo "run this code only from the code directory!"

rm -r /direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/asymmetry_calculation/build
rm -r /direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/asymmetry_calculation/install/lib
rm -r autom4te.cache/ aclocal.m4  config.guess config.sub configure depcomp install-sh Makefile.in missing temp/ ltmain.sh

mkdir -p /direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/asymmetry_calculation/build
mkdir -p /direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/asymmetry_calculation/install

base_dir=/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/asymmetry_calculation
cd $base_dir/build;
$base_dir/autogen.sh \
    --prefix=$base_dir/install;
make;
make install;

