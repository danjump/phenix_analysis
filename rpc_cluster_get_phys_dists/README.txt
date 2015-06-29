This directory contains code for the rpc QA fast analysis. 
More specifically, the process_pdst code which runs over pdst's to produce histograms. 
It also has the scrips/config files for compiling the process_pdst code into a library.  
Lastly, in the execute_manually directory are a couple of scripts and macros
to run over data by hand.


Instructions for compiling:

first make build and install directories:
mkdir build; mkdir install

next, move to the build directory and run autogen.sh as follows (substituting
this full direcory path for <source path>):
cd build
<source path>/autogen.sh --prefix=<source path>/install

finally, once autogen.sh runs sucessfully, while still in the build directory,
run make install:
make install

The code will now be compiled into a library in <source
path>/install/lib/librpc_qa.so

