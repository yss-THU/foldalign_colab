#ifndef CONVERTSTRING
#define CONVERTSTRING

#include <string>
#include <sstream>
#include "exception.cxx"

namespace convertString {


	template< typename type >
	inline static type toNumeric(std::string number) {
		type result;
		
		std::stringstream ss(number);
		
		return ss >> result ? result : throw exception(std::string("Not a number: "+number), true);
	}

	inline static std::string toString(std::string textString) {
		return textString;
	}
	
	inline static bool toBool(std::string logic) {
		if (logic.compare("true") == 0) {
			return true;
		}
		return false;
	}


}


#endif
