
void events_comparison_execute() {
  
//  char * code_dir = "/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/wness_preselection";
  char lib_file[300];
	sprintf(lib_file,"../install/lib/libSimpleClass.so");
	gSystem->Load(lib_file);
        events_comparison();
}
