first you must create a list of the files you want to run over.
this is just a text file containing full paths to the files you want,
one on each line. See "filelists/makefilelist.sh" for examples
of generating a file list.

submit jobs to condor by executing condor_pp510.sh

condor_pp510.sh loops over the file list and submits one job for each file.
note, this script has output paths hard-coded which must be modified,
including condor log directories which it generates.



Individual jobs are defined by condor_pp510.job and condor_pp510.cmd 

condor_pp510.job defines general condor parameters

condor_pp510.cmd actually calls a root macro. This scripts has parameters passed
to it by condor_pp510.sh



once jobs are submitted, you can check your status with this command:

condor_status -submitters | grep <username>


jobs can be canceled with this command:

condor_rm <username>

