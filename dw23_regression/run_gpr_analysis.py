import sys
import numpy as np

#contains main class for executing the analysis:
from gpr_analysis_obj import gpr_analysis_obj

def main(args):
    args_dict=dict()

    args_dict['do_method_raw'   ] = True
    args_dict['do_method_scaled'] = False
    args_dict['do_method_log'   ] = False
    args_dict['do_threshold'    ] = np.array([[.9,True],[.8,True],[.7,True]])
    args_dict['input_file'      ] = 'input/dw23_vs_wness_list.txt'
    
    gpr_obj = gpr_analysis_obj(args_dict)
    gpr_obj.get_data()
    gpr_obj.run_gpr()

if __name__=='__main__':
        main(sys.argv[1:])

