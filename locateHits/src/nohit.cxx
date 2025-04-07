#include "nohit.hxx"

#include "processEntry.cxx"
#include "arguments.cxx"

#include <iostream>
#include <fstream>
#include <string>

/*! \brief fnoa makes a list of the best non-overlapping alignments using the
	output from Foldalign. Foldalign must be run with option -plot_score to
	generate input data for fnoa.

*/
int main(int argc, char* argv[]) {

	// Setup the posible options
	arguments arg(argc, argv); // parse the options

	arg.addString("version", version);

	// Handle -version option
	if ( arg.boolOpt("-version") ) {
		std::cout << program_name << " version: " << version << std::endl;
		std::cout << copyright << std::endl;
		return 0;
	}

	// Handle -help and -h options
	if (arg.boolOpt("-help") || arg.boolOpt("-h")) {
		std::cout << program_name << " version: " << version << std::endl;
		std::cout << copyright << std::endl;
		std::cout << "Usage:\n" << program_name << " [<options>] <file_1> [<file_2>] ...\n";
		std::cout << program_name << " finds the best scoring non-overlapping hits in a FOLDALIGN output file.\n";
		std::cout << "The FOLDALIGN output file must have been made using FOLDALIGN option -plot_score\n";
		std::cout << "The options are:\n";
		arg.printOptions();
		return 0;
	}

	if (bad >= arg.stOpt("-score_cut_off")) {
		std::cerr << program_name << " error!" << std::endl;
		std::cerr << "The score_cut_off is to small compared to the internal ''bad'' score." << std::endl;
		std::cerr << "Please either use a score_cut_off larger than: " << bad << " or lower the bad score in" << std::endl;
		std::cerr << "the file nohit.hxx and recompile" << std::endl;
		return 1;
	}		

	std::string first_line;
	const bool border = arg.boolOpt("-lambda_border");
	try {
		if (arg.numberOfArguments() == 0) {
			while (getline(std::cin, first_line)) {
				processEntry(std::cin, arg.stOpt("-score_cut_off"), border, arg.ltOpt("-max_number_of_hits"), arg.boolOpt("-estimate_parameters"));
			}
		}
		else {
			for(int i = 0; i < arg.numberOfArguments(); i++) {
				std::ifstream input_file;
				input_file.open(arg.argument(i).c_str());
				while (getline(input_file, first_line)) {
					processEntry(input_file, arg.stOpt("-score_cut_off"), border, arg.ltOpt("-max_number_of_hits"), arg.boolOpt("-estimate_parameters"));
				}
				input_file.close();
			}
		}
	}
	catch ( ... ) {
		std::cerr << program_name << ": Unknown error!" << std::endl;
		return 1;
	}
	return 0;
}
