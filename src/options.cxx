#ifndef OPTIONS
#define OPTIONS

#include<string>
#include<iostream>
#include<iomanip>
#include<stdlib.h>

#include "foldalign.hxx"
#include "helper.cxx"
#include "exception.cxx"

template<class type>
class options {
public:
	options(std::string optionType, const int& numberOfOptions, std::string* opt, 
			std::string* val, std::string* typ, std::string* des,
			type (*convertString2type)(std::string));
	
	options(const options& opt) {copy(opt);}
	options& operator=(const options& opt);
		
	inline type get(const std::string& name) const;
	inline bool exists(const std::string& name) const;
	inline void set(const std::string& name, type value);
	inline void add(std::string name, type value, std::string description = "");
	
	type getDefault(const std::string& name) const;
	std::string getDescription(const std::string& name) const;

	int getSize() const {return size;}

	~options() {clean();}

private:
	options(); // No default constructor
	
	int size;
	
	std::string* names;
	type* values;
	type* defaultValues;
	std::string* descriptions;
	
	void clean();
	void copy(const options& opt);
};

template<class type>
options< type >::options(std::string optionType, const int& numberOfOptions, std::string* opt, 
			std::string* val, std::string* typ, std::string* des,
			type (*convertString2type)(std::string)) {

	size = 0;
	for(int c=0; c < numberOfOptions; c++) {
		if (typ[c].compare(optionType) == 0) {
			size++;
		}
	}

	if (size == 0) {
		return;
	}
	
	names = new std::string[size];
	values = new type[size];
	defaultValues = new type[size];
	descriptions = new std::string[size];
	
	int pos = 0;
	for(int c=0; c < numberOfOptions; c++) {
		if (typ[c].compare(optionType) == 0) {
			names[pos] = opt[c];
			descriptions[pos] = des[c];
			values[pos] = (*convertString2type)(val[c]);
			defaultValues[pos] = values[pos];
			pos++;
		}
	}
}

template<class type>
options<type>& options<type>::operator=(const options& opt) {
	if (this != &opt) {
		clean();		
		copy(opt);
	}
	return *this;
}

template<class type>
void options<type>::copy(const options& opt) {
	size = opt.size;
	names = new std::string[size];
	values = new type[size];
	defaultValues = new type[size];
	descriptions = new std::string[size];
	
	helper::copyArray(names, opt.names, size);
	helper::copyArray(values, opt.values, size);
	helper::copyArray(defaultValues, opt.defaultValues, size);
	helper::copyArray(descriptions, opt.descriptions, size);
}

template< class type >
void options<type>::clean() {
    if (size > 0) {
	delete[] names;
	delete[] values;
	delete[] defaultValues;
	delete[] descriptions;
    }
}

template<class type>
type options< type >::get(const std::string& name) const {

	for(int c=0; c < size; c++) {
		if (name.compare(names[c]) == 0) {
			return values[c];
		}
	}
	std::string error = "Program error! Options.get. Option " + name + " does not exists.";
	throw exception(error, true);
}

template<class type>
void options< type >::set(const std::string& name, type value) {

	for(int c=0; c < size; c++) {
		if (name.compare(names[c]) == 0) {
			values[c] = value;
			return;
		}
	}
	std::string error = "Program error! Options.set. Option " + name + " does not exists.";
	throw exception(error, true);
}

template<class type>
bool options< type >::exists(const std::string& name) const {
	for(int c=0; c < size; c++) {
		if (name.compare(names[c]) == 0) {
			return true;
		}
	}
	return false;
}

template<class type>
void options< type >::add(std::string name, type value, std::string description) {

	for(int c=0; c < size; c++) {
		if (name.compare(names[c]) == 0) {
			std::string error = "Program error! Options.add. Option " + name + " already exists.";
			throw exception(error, true);
		}
	}

	const int newSize = size +1;
	if (size == 0) {
		names = new std::string[newSize];
		values = new type[newSize];
		defaultValues = new type[newSize];
		descriptions = new std::string[newSize];
	
	}
	else {
		helper::expandArray(names, size, newSize);
		helper::expandArray(values, size, newSize);
		helper::expandArray(defaultValues, size, newSize);
		helper::expandArray(descriptions, size, newSize);
	}

	names[size] = name;
	values[size] = value;
	defaultValues[size] = value;
	descriptions[size] = description;

	size++;

	return;
}

template<class type>
type options< type >::getDefault(const std::string& name) const {

	for(int c=0; c < size; c++) {
		if (name.compare(names[c]) == 0) {
			return defaultValues[c];
		}
	}
	std::string error = "Program error! Options.getDefault. Option " + name + " does not exists.";
	throw exception(error, true);
}

template<class type>
std::string options< type >::getDescription(const std::string& name) const {

	for(int c=0; c < size; c++) {
		if (name.compare(names[c]) == 0) {
			return descriptions[c];
		}
	}
	std::string error = "Program error! Options.getDescription. Option " + name + " does not exists.";
	throw exception(error, true);
}
#endif
