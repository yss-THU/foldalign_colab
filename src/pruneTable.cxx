#ifndef PRUNETABLE
#define PRUNETABLE

#include "foldalign.hxx"
#include "readfile.cxx"
#include "helper.cxx"
#include "stack_ssl.cxx"

class pruneTable {
public:
	pruneTable(scoreType pruneStart, scoreType linearPrune, lengthType size, bool NoPrune)
	  : pruneZero(pruneStart), pruneCoefficient(linearPrune), tableSize(size+1), noPrune(NoPrune) {

		table = new scoreType[tableSize];
		buildTable();

	};

	pruneTable& operator=(const pruneTable& p) {
		if (this != &p) {
			clean();
			copy(p);
		}
		return *this;
	};

	pruneTable(const pruneTable& p) {
		copy(p);
	}

	~pruneTable() {
		clean();
	};

	scoreType get(const lengthType length) const {
		return table[length];
	};

	lengthType getLength() const {
		return tableSize;
	};

	void setPruneCoefficient(scoreType value) {
		pruneCoefficient = value;
		buildTable();
	}

	void buildTableFromFile(readfile*& file);
	inline void reSizeTable(const lengthType size);

private:
	scoreType pruneZero;
	scoreType pruneCoefficient;
	lengthType tableSize;
	bool noPrune;

	scoreType* table;

	// Used to store table values read from a file
	struct prunes {
		scoreType score;
		lengthType index;
	};

	void clean() {
		delete[] table;
	};

	void copy(const pruneTable& p) {

		pruneZero = p.pruneZero;
		pruneCoefficient = p.pruneCoefficient;
		tableSize = p.getLength();

		noPrune = p.noPrune;

		table = new scoreType[tableSize];
		for(lengthType i = 0; i < tableSize; i++) {
			table[i] = p.get(i);
		}
	};

	inline void buildTable() {

		if (! noPrune) {
		  table[0] = pruneZero;
		  for(lengthType i = 1; i< tableSize; i++) {
			  table[i] = pruneZero + pruneCoefficient*(i);
		  }
		}
		else {
		  for(lengthType i = 0; i < tableSize; i++) {
		    table[i] = big_neg;
		  }
		}
	}

	inline scoreType pruneScore(const lengthType len_new, const lengthType len_old, const scoreType zero);
	inline void discontiniusPrune(const scoreType startValue, const lengthType start);
	inline void fillPrune(const lengthType start, const scoreType score, const lengthType stop, stack_ssl<prunes, lengthType>& lineStack);
	inline void assignFromStack(scoreType*& store, stack_ssl<prunes, lengthType>& lineStack);
};

inline void pruneTable::buildTableFromFile(readfile*& file) {

	std::string line;
	table[0] = big_neg;

	stack_ssl<prunes, lengthType> lineStack = stack_ssl<prunes, lengthType>();

	scoreType last_score = big_neg;
	lengthType last_index = 0;

	while (file->get_line_failEmpty(line)) {

		int prev = 0;
		lengthType index = helper::getValue(prev, line);
		scoreType score = helper::getValue(prev, line);

		if (index < last_index) {
			std::cerr << "Warning: Length " << index << " in the Pruning: table has a smaller length than the previous entry: " << last_index  << ". This table must be sorted by length value. The value is ignored." << std::endl;
			continue;
		}

		if (index > last_index +1) {
			fillPrune(last_index, last_score, index, lineStack);
		}

		prunes prune = {score, index};
		lineStack.push(prune);

		last_score = score;
		last_index = index;

	}

	last_index++;
	if (last_index > tableSize) {
		helper::expandArray(table, tableSize, last_index);
	}

	assignFromStack(table, lineStack);

	tableSize = last_index;
	pruneZero = table[0];
}

inline void pruneTable::assignFromStack(scoreType*& store, stack_ssl<prunes, lengthType>& lineStack) {

	const lengthType size = lineStack.size();
	for(lengthType i = 0; i < size; i++) {
		prunes prune = lineStack.pop();
		store[prune.index] = prune.score;
	}
	store[0] = store[1];
}

inline void pruneTable::reSizeTable(const lengthType size) {

	// The old table is stored temporarily
	scoreType* old_table = table;
	lengthType old_size = tableSize;

	// The new table is made
	tableSize = size +1;
	table = new scoreType[tableSize];

	// Copy the old table to the new.
	lengthType min_size = old_size < tableSize ? old_size : tableSize;
	helper::copyArray(table, old_table, min_size);

	// If there is no pruning (ie the prune score is big_neg for the last entry)
	// keep it that way
	if (table[tableSize-1] == big_neg) {
		pruneCoefficient = 0;
	}

	// Calculate the pruning scores for the rest of the table
	for(lengthType i = old_size; i < tableSize; i++) {
		table[i] = pruneScore(i, old_size-1, table[old_size-1]);
	}

	delete[] old_table;

}

inline scoreType pruneTable::pruneScore(const lengthType len_new, const lengthType len_old, const scoreType zero) {
	return pruneCoefficient*(len_new - len_old) + zero;
}

inline void pruneTable::discontiniusPrune(const scoreType startValue,
										   const lengthType startLength) {
	for(lengthType i = startLength; i < tableSize; i++) {
		table[i] = pruneScore(i, startLength, startValue);
	}
}

inline void pruneTable::fillPrune(const lengthType start, const scoreType score,
		const lengthType stop, stack_ssl<prunes, lengthType>& lineStack) {

	for(lengthType x = start+1; x < stop; x++) {

		prunes prune = {score, x};
		lineStack.push(prune);

	}
}

#endif
