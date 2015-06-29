

sig_file=/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run12/wness_tree/run12_simulation_367593_combined_hists.root
bkg_file=/phenix/spin/phnxsp01/danielj/w_asymmetry_analysis_files/run13/wness_tree/run13_data_hists_combined.root

root -b -q plot_wness_pdf_hists.C+\(\"${sig_file}\",\"${bkg_file}\"\)
