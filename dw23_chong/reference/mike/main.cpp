#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>

struct gpr_fit_values {
  int index, arm, charge;
  float dw23_bin_center, wness_bin_center, value, uncertainty;
  std::string GetDataString() {
    std::stringstream ss;
    ss << index << " " << arm << " " << charge  
        << " " << dw23_bin_center << " " << wness_bin_center 
        << " " << value
        << " " << uncertainty;
    return ss.str();
  }
};

int main(int argc, char** argv) {
  if( argc < 2) {
    std::cerr << std::endl << "\tUsage is: " << argv[0] << " <text_file_with_values>" << std::endl;
  } else {
    std::string input_file = argv[1];
    std::cout << std::endl << "Processing file: " << input_file << std::endl;
    std::ifstream in_file(input_file.c_str());
    std::string line = "";
    std::vector<gpr_fit_values> data;
    if(in_file) {
      std::cout << "Checking Data string: " << std::endl;
      while(getline(in_file,line)) {
        if(line[0] == '#') continue;
        gpr_fit_values data_temp;
        std::stringstream ss;
        ss.str(line.c_str());
        ss >> data_temp.index >> data_temp.arm >> data_temp.charge
            >> data_temp.dw23_bin_center >> data_temp.wness_bin_center
            >> data_temp.value >> data_temp.uncertainty;
        std::cout << "\t" << data_temp.GetDataString() << std::endl;
        data.push_back(data_temp);
      }
      std::cout << " Lines successfully extracted: " << data.size() << std::endl;
      std::cout << " Done." << std::endl;
    } else {
      std::cerr << std::endl << "Badly formatted input file string, or problem opening file." << std::endl;
    }
  }
	return 0;
}
