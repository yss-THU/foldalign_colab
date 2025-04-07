###############################################################################
#                                                                             #
#   Copyright 2004 - 2005 Jakob Hull Havgaard, hull@bioinf.kvl.dk             #
#                                                                             #
#   This file is part of Foldalign                                            #
#                                                                             #
#   Foldalign is free software; you can redistribute it and/or modify         #
#   it under the terms of the GNU General Public License as published by      #
#   the Free Software Foundation; either version 2 of the License, or         #
#   (at your option) any later version.                                       #
#                                                                             #
#   Foldalign is distributed in the hope that it will be useful,              #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
#   GNU General Public License for more details.                              #
#                                                                             #
#   You should have received a copy of the GNU General Public License         #
#   along with Foldalign; if not, write to the Free Software                  #
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA #
#                                                                             #
###############################################################################

all: bin/foldalign bin/locateHits
#options = -pthread -Wall -pedantic -ansi -DDEBUG
#options = -pthread -Wall -pedantic -ansi -g
#options = -pthread -Wall -pedantic -ansi -pg
options = -pthread -Wall -pedantic -ansi -O3
#options = -pthread -Wall -pedantic -ansi -O3 -static

sup_options = -Wno-long-long


#foldalignobjects = \
#                src/arguments.o \
#		src/backtrack.o \
#		src/cell.o \
#		src/constraints.o \
#		src/convertString.o \
#		src/exception.o \
#		src/fold.o \
#		src/foldThreadHandler.o \
#		src/foldK.o \
#		src/globalPruneTable.o \
#		src/helper.o \
#		src/jl.o \
#		src/local.o \
#		src/lockList.o \
#		src/longCell.o \
#		src/longCellPtr.o \
#		src/longCellState.o \
#		src/longTermMemory.o \
#		src/mblcell.o \
#		src/mbllist.o \
#		src/mem_stack.o \
#		src/multistacks.o \
#		src/options.o \
#		src/output.interface.o \
#		src/output.o \
#		src/output.col.o \
#		src/output.stk.o \
#		src/output.summary.o \
#		src/prune.o \
#		src/pruneTable.o \
#		src/readConstraints.o \
#		src/readfile.o \
#		src/results.o \
#		src/runKargs.o \
#		src/scorematrix.o \
#		src/seqs.o \
#		src/sequenceArg.o \
#		src/sequence.o \
#		src/shortTermMemory.o \
#		src/stack.o \
#		src/stack_ssl.o \
#		src/stmSubMatrix.o \
#		src/stmWiWkSubMatrix.o \
#		src/tables.o \
#		src/thread.o \
#		src/threadMaster.o \
#		src/two_link.o

locateHitsobjects = \
		locateHits/src/alignmentCells.o\
		locateHits/src/helper.o \
		locateHits/src/arguments.o \
		locateHits/src/datalist.o \
		locateHits/src/nonoverlap.o \
		locateHits/src/processEntry.o \
		locateHits/src/pValue.o

testobjects = \
		test/src/test.o\
		test/src/testTest.o\
		test/src/testArguments.o \
		test/src/testCell.o\
		test/src/testCombineSimilarityEnergy.o\
		test/src/testConstraints.o\
		test/src/testException.o\
		test/src/testGlobalPrune.o \
		test/src/testHelper.o\
		test/src/testLockList.o\
		test/src/testLongCell.o\
		test/src/testLongTermMemory.o\
		test/src/testMatrix2d.o\
		test/src/testMatrix4d.o\
		test/src/testOptions.o \
		test/src/testPrune.o \
		test/src/testPruneTable.o\
		test/src/testReadConstraints.o\
		test/src/testReadfile.o\
		test/src/testScorematrix.o\
		test/src/testTwo_link.o

cc = g++

bin/foldalign: src/foldalign.o src/revnumber.hxx src/empty src/foldalign.hxx src/arguments.o src/backtrack.o src/cell.o src/constraints.o src/exception.o src/local.o src/longCell.o src/helper.o src/output.col.o src/output.stk.o src/output.summary.o src/scorematrix.o src/seqs.o src/sequence.o  #Makefile
	$(cc) $(options) -o bin/foldalign src/foldalign.o src/arguments.o src/backtrack.o src/cell.o src/constraints.o src/exception.o src/local.o src/longCell.o src/helper.o src/output.col.o src/output.stk.o src/output.summary.o src/scorematrix.o src/seqs.o src/sequence.o
#	$(cc) $(options) src/foldalign.cxx $(foldalignobjects) -o bin/foldalign
#	test/bin/runTest

bin/locateHits: locateHits/src/nohit.cxx $(locateHitsobjects) locateHits/src/nohit.hxx src/revnumber.hxx
	g++ $(options) $(sup_options) locateHits/src/nohit.cxx $(locateHitsobjects) -o bin/locateHits

#src/revnumber.hxx: .git/index src/empty_sub
#	bin/revnumber > src/revnumber.hxx



test: bin/test.foldalign test/lastTest
	bin/test.foldalign

bin/test.foldalign: test/src/foldalignTest.cxx $(testobjects)
	g++ $(options) test/src/foldalignTest.cxx $(testObjects) -o bin/test.foldalign

test/lastTest: bin/foldalign
	test/bin/runTest

#testobjects = \
#	testsrc/test_foldThreadHandler.o

#testsrc/test_foldalign: testsrc/test_foldalign.cxx $(testobjects)
#	$(cc) $(options) testsrc/test_foldalign.cxx -o testsrc/test_foldalign


#testsrc/test_foldThreadHandler.o: testsrc/test_foldThreadHandler.cxx src/foldThreadHandler.o
#	$(cc) $(options) -c testsrc/test_foldThreadHandler.cxx -o testsrc/test_foldThreadHandler.o

#src/foldThreadHandler.o: src/foldThreadHandler.cxx
#	$(cc) $(options) -c src/foldThreadHandler.cxx -o src/foldThreadHandler.o

re:
	touch src/empty
	make

resub:
	touch src/empty_sub
	make

src/empty:
	touch src/empty


src/empty_sub:
	touch src/empty_sub

clean:
	rm -rf src/*.o
	rm -rf locateHits/src/*.o
	rm -rf bin/foldalign bin/locateHits
	touch src/empty_sub

touch:
	touch src/*.cxx

############################### Foldalign ########################

src/foldalign.o: src/foldalign.cxx src/local.o src/backtrack.o src/sequence.o src/scorematrix.o src/seqs.o src/arguments.o src/helper.o src/exception.o src/output.stk.o src/output.summary.o src/cell.o src/longCell.o src/constraints.o
	$(cc) $(options) -c -o src/foldalign.o src/foldalign.cxx
src/arguments.o: src/arguments.cxx src/foldalign.hxx src/helper.o src/exception.o src/options.o src/convertString.o
	$(cc) $(options) -c -o src/arguments.o src/arguments.cxx

src/backtrack.o: src/backtrack.cxx src/helper.o src/output.o src/fold.o src/multistacks.o src/arguments.o src/results.o src/cell.o src/longCell.o src/longCellState.o src/longCellPtr.o src/stack_ssl.o src/mbllist.o src/mblcell.o src/constraints.o src/tables.o src/combineSimilarityEnergy.o
	$(cc) $(options) -c -o src/backtrack.o src/backtrack.cxx

src/cell.o: src/cell.cxx src/mbllist.o src/foldalign.hxx
	$(cc) $(options) -c -o src/cell.o src/cell.cxx

src/combineSimilarityEnergy.o: src/combineSimilarityEnergy.cxx src/foldalign.hxx
	$(cc) $(options) -c -o src/combineSimilarityEnergy.o src/combineSimilarityEnergy.cxx

src/constraints.o: src/constraints.cxx src/foldalign.hxx src/readConstraints.o src/jl.o src/helper.o
	$(cc) $(options) -c -o src/constraints.o src/constraints.cxx

src/convertString.o: src/convertString.cxx src/exception.o
	$(cc) $(options) -c -o src/convertString.o src/convertString.cxx

src/exception.o: src/exception.cxx
	$(cc) $(options) -c -o src/exception.o src/exception.cxx

src/fold.o: src/fold.cxx src/foldalign.hxx src/foldK.o src/sequence.o src/output.o src/helper.o src/arguments.o src/results.o src/longTermMemory.o src/shortTermMemory.o src/exception.o src/stack_ssl.o src/mbllist.o src/longCell.o src/constraints.o src/jl.o src/cell.o src/threadMaster.o
	$(cc) $(options) -c -o src/fold.o src/fold.cxx

src/foldK.o: src/foldK.cxx src/foldalign.hxx src/sequence.o src/output.o src/helper.o src/arguments.o src/results.o src/longTermMemory.o src/shortTermMemory.o src/exception.o src/stack_ssl.o src/mbllist.o src/longCell.o src/constraints.o src/jl.o src/cell.o src/thread.o src/runKargs.o src/foldThreadHandler.o src/stmSubMatrix.o src/tables.o src/combineSimilarityEnergy.o
	$(cc) $(options) -c -o src/foldK.o src/foldK.cxx

src/foldThreadHandler.o: src/foldThreadHandler.cxx src/foldalign.hxx
	$(cc) $(options) -c -o src/foldThreadHandler.o src/foldThreadHandler.cxx

src/globalPruneTable.o: src/globalPruneTable.cxx src/foldalign.hxx src/exception.o
	$(cc) $(options) -c -o src/globalPruneTable.o src/globalPruneTable.cxx

src/helper.o: src/helper.cxx src/foldalign.hxx src/exception.o src/revnumber.hxx
	$(cc) $(options) -c -o src/helper.o src/helper.cxx

src/i.o: src/i.cxx
	$(cc) $(options) -c -o src/i.o src/i.cxx

src/jl.o: src/jl.cxx src/foldalign.hxx
	$(cc) $(options) -c -o src/jl.o src/jl.cxx

src/local.o: src/local.cxx  src/backtrack.o src/fold.o src/cell.o src/longCell.o src/exception.o src/constraints.o
	$(cc) $(options) -c -o src/local.o src/local.cxx

src/lockList.o: src/lockList.cxx src/foldalign.hxx src/exception.o src/helper.o
	$(cc) $(options) -c -o src/lockList.o src/lockList.cxx

src/longCell.o: src/longCell.cxx  src/foldalign.hxx src/mbllist.o
	$(cc) $(options) -c -o src/longCell.o src/longCell.cxx

src/longCellPtr.o: src/longCellPtr.cxx src/foldalign.hxx src/mbllist.o
	$(cc) $(options) -c -o src/longCellPtr.o src/longCellPtr.cxx

src/longCellState.o: src/longCellState.cxx src/foldalign.hxx src/mbllist.o
	$(cc) $(options) -c -o src/longCellState.o src/longCellState.cxx

src/longTermMemory.o: src/longTermMemory.cxx src/longCell.o src/helper.o src/two_link.o src/foldalign.hxx src/exception.o src/lockList.o
	$(cc) $(options) -c -o src/longTermMemory.o src/longTermMemory.cxx

src/matrix2d.o: src/matrix2d.cxx src/foldalign.hxx
	$(cc) $(options) -c -o src/matrix2d.o src/matrix2d.cxx

src/matrix4d.o: src/matrix4d.cxx src/foldalign.hxx
	$(cc) $(options) -c -o src/matrix4d.o src/matrix4d.cxx

src/mblcell.o: src/mblcell.cxx src/cell.o src/mbllist.o src/foldalign.hxx
	$(cc) $(options) -c -o src/mblcell.o src/mblcell.cxx

src/mbllist.o: src/mbllist.cxx src/foldalign.hxx
	$(cc) $(options) -c -o src/mbllist.o src/mbllist.cxx

src/mem_stack.o: src/mem_stack.cxx src/exception.o
	$(cc) $(options) -c -o src/mem_stack.o src/mem_stack.cxx

src/multistacks.o: src/multistacks.cxx src/stack.o
	$(cc) $(options) -c -o src/multistacks.o src/multistacks.cxx

src/options.o: src/options.cxx src/foldalign.hxx src/helper.o src/exception.o src/convertString.o
	$(cc) $(options) -c -o src/options.o src/options.cxx

src/output.col.o: src/output.col.cxx src/arguments.o src/sequence.o src/helper.o src/stack.o src/foldalign.hxx src/output.hxx src/output.interface.o
	$(cc) $(options) -c -o src/output.col.o src/output.col.cxx

src/output.o: src/output.cxx src/arguments.o src/sequence.o src/helper.o src/stack.o src/foldalign.hxx src/output.hxx src/output.interface.o src/output.stk.o src/output.col.o src/output.summary.o
	$(cc) $(options) -c -o src/output.o src/output.cxx

src/output.interface.o: src/output.interface.cxx src/arguments.o src/sequence.o src/helper.o src/stack.o src/foldalign.hxx src/output.hxx
	$(cc) $(options) -c -o src/output.interface.o src/output.interface.cxx

src/output.stk.o: src/output.stk.cxx src/arguments.o src/sequence.o src/helper.o src/stack.o src/foldalign.hxx src/output.hxx src/output.interface.o
	$(cc) $(options) -c -o src/output.stk.o src/output.stk.cxx

src/output.summary.o: src/output.summary.cxx src/arguments.o src/sequence.o src/helper.o src/stack.o src/foldalign.hxx src/output.hxx src/output.interface.o
	$(cc) $(options) -c -o src/output.summary.o src/output.summary.cxx

src/prune.o: src/prune.cxx src/pruneTable.o src/globalPruneTable.o src/foldalign.hxx
	$(cc) $(options) -c -o src/prune.o src/prune.cxx

src/pruneTable.o: src/pruneTable.cxx src/foldalign.hxx src/readfile.o src/helper.o src/stack_ssl.o
	$(cc) $(options) -c -o src/pruneTable.o src/pruneTable.cxx

src/read4dMatrix.o: src/read4dMatrix.cxx src/matrix4d.o src/readfile.o src/scorematrix.o src/helper.o
	$(cc) $(options) -c -o src/read4dMatrix.o src/read4dMatrix.cxx

src/readConstraints.o: src/readConstraints.cxx src/foldalign.hxx src/readfile.o src/helper.o src/combineSimilarityEnergy.o
	$(cc) $(options) -c -o src/readConstraints.o src/readConstraints.cxx

src/readfile.o: src/readfile.cxx src/exception.o
	$(cc) $(options) -c -o src/readfile.o src/readfile.cxx

src/results.o: src/results.cxx src/foldalign.hxx src/mbllist.o
	$(cc) $(options) -c -o src/results.o src/results.cxx

src/runKargs.o: src/runKargs.cxx src/constraints.o src/foldalign.hxx
	$(cc) $(options) -c -o src/runKargs.o src/runKargs.cxx

src/scorematrix.o: src/scorematrix.cxx src/foldalign.hxx src/readfile.o src/helper.o src/stack_ssl.o src/exception.o src/matrix4d.o src/matrix2d.o src/prune.o
	$(cc) $(options) -c -o src/scorematrix.o src/scorematrix.cxx

src/seqs.o: src/seqs.cxx src/foldalign.hxx src/sequence.o src/sequenceArg.o src/helper.o src/scorematrix.o src/arguments.o src/readfile.o src/exception.o
	$(cc) $(options) -c -o src/seqs.o src/seqs.cxx

src/sequence.o: src/sequence.cxx src/foldalign.hxx src/scorematrix.o
	$(cc) $(options) -c -o src/sequence.o src/sequence.cxx

src/sequenceArg.o: src/sequenceArg.cxx src/sequence.o src/arguments.o src/foldalign.hxx
	$(cc) $(options) -c -o src/sequenceArg.o src/sequenceArg.cxx

src/shortTermMemory.o: src/shortTermMemory.cxx src/arguments.o src/helper.o src/scorematrix.o src/exception.o src/foldalign.hxx src/stmWiWkSubMatrix.o
	$(cc) $(options) -c -o src/shortTermMemory.o src/shortTermMemory.cxx

src/stack.o: src/stack.cxx src/exception.o
	$(cc) $(options) -c -o src/stack.o src/stack.cxx

src/stack_ssl.o: src/stack_ssl.cxx
	$(cc) $(options) -c -o src/stack_ssl.o src/stack_ssl.cxx

src/stmSubMatrix.o: src/stmSubMatrix.cxx src/arguments.o src/helper.o src/scorematrix.o src/exception.o src/foldalign.hxx src/shortTermMemory.o src/stmWiWkSubMatrix.o
	$(cc) $(options) -c -o src/stmSubMatrix.o src/stmSubMatrix.cxx

src/stmWiWkSubMatrix.o: src/stmWiWkSubMatrix.cxx src/arguments.o src/helper.o src/scorematrix.o src/exception.o src/foldalign.hxx
	$(cc) $(options) -c -o src/stmWiWkSubMatrix.o src/stmWiWkSubMatrix.cxx

src/tables.o: src/tables.cxx src/foldalign.hxx src/helper.o
	$(cc) $(options) -c -o src/tables.o src/tables.cxx

src/test.readConstraints.o: src/test.readConstraints.cxx src/readConstraints.o
	$(cc) $(options) -c -o src/test.readConstraints.o src/test.readConstraints.cxx

src/thread.o: src/thread.cxx
	$(cc) $(options) -c -o src/thread.o src/thread.cxx

src/threadMaster.o: src/threadMaster.cxx src/foldalign.hxx src/thread.o src/runKargs.o src/foldThreadHandler.o src/foldK.o src/sequence.o src/output.o src/helper.o src/arguments.o src/results.o src/longTermMemory.o src/shortTermMemory.o src/exception.o src/stack_ssl.o src/mbllist.o src/longCell.o src/constraints.o src/jl.o src/cell.o
	$(cc) $(options) -c -o src/threadMaster.o src/threadMaster.cxx

src/two_link.o: src/two_link.cxx src/longCell.o src/foldalign.hxx src/exception.o
	$(cc) $(options) -c -o src/two_link.o src/two_link.cxx

###################### LocateHits ###########################

locateHits/src/alignmentCells.o: locateHits/src/alignmentCells.cxx locateHits/src/nohit.hxx
	$(cc) $(options) -o locateHits/src/alignmentCells.o -c locateHits/src/alignmentCells.cxx

locateHits/src/arguments.o: locateHits/src/arguments.cxx locateHits/src/nohit.hxx locateHits/src/helper.o
	$(cc) $(options) -o locateHits/src/arguments.o -c locateHits/src/arguments.cxx

locateHits/src/datalist.o: locateHits/src/datalist.cxx locateHits/src/nohit.hxx
	$(cc) $(options) -o locateHits/src/datalist.o -c locateHits/src/datalist.cxx

locateHits/src/helper.o: locateHits/src/helper.cxx
	$(cc) $(options) -o locateHits/src/helper.o -c locateHits/src/helper.cxx

locateHits/src/nonoverlap.o: locateHits/src/nonoverlap.cxx locateHits/src/nohit.hxx locateHits/src/alignmentCells.o locateHits/src/datalist.o
	$(cc) $(options) -o locateHits/src/nonoverlap.o -c locateHits/src/nonoverlap.cxx

locateHits/src/pValue.o: locateHits/src/pValue.cxx locateHits/src/nohit.hxx
	$(cc) $(options) $(sup_options) -o locateHits/src/pValue.o -c locateHits/src/pValue.cxx

locateHits/src/processEntry.o: locateHits/src/processEntry.cxx locateHits/src/nohit.hxx locateHits/src/datalist.o locateHits/src/alignmentCells.o locateHits/src/nonoverlap.o locateHits/src/pValue.o
	$(cc) $(options) $(sup_options) -o locateHits/src/processEntry.o -c locateHits/src/processEntry.cxx

locateHits/src/test.datalist.o: locateHits/src/test.datalist.cxx locateHits/src/nohit.hxx locateHits/src/datalist.o locateHits/src/alignmentCells.o locateHits/src/nonoverlap.o
	$(cc) $(options) -o locateHits/src/test.datalist.o -c locateHits/src/test.datalist.cxx

#################### TESTING ###############################

test/src/test.o: test/src/test.cxx
	$(cc) $(options) -o test/src/test.o -c test/src/test.cxx

test/src/testArguments.o: test/src/testArguments.cxx test/src/test.o src/arguments.o src/foldalign.hxx src/options.o src/convertString.o
	$(cc) $(options) -o test/src/testArguments.o -c test/src/testArguments.cxx

test/src/testCell.o: test/src/testCell.cxx test/src/test.o src/cell.o src/foldalign.hxx src/mbllist.o
	$(cc) $(options) -o test/src/testCell.o -c test/src/testCell.cxx

test/src/testCombineSimilarityEnergy.o: test/src/testCombineSimilarityEnergy.cxx src/combineSimilarityEnergy.o src/foldalign.hxx test/src/test.o
	$(cc) $(options) -o test/src/testCombineSimilarityEnergy.o -c test/src/testCombineSimilarityEnergy.cxx

test/src/testConstraints.o: test/src/testConstraints.cxx src/constraints.o test/src/test.o
	$(cc) $(options) -o test/src/testConstraints.o -c test/src/testConstraints.cxx

test/src/testException.o: test/src/testException.cxx test/src/test.o
	$(cc) $(options) -o test/src/testException.o -c test/src/testException.cxx

test/src/testGlobalPrune.o: test/src/testGlobalPrune.cxx test/src/test.o
	$(cc) $(options) -o test/src/testGlobalPrune.o -c test/src/testGlobalPrune.cxx

test/src/testHelper.o: test/src/testHelper.cxx test/src/test.o src/helper.o
	$(cc) $(options) -o test/src/testHelper.o -c test/src/testHelper.cxx

test/src/testLockList.o: test/src/testLockList.cxx test/src/test.o src/lockList.o src/foldalign.hxx src/exception.o
	$(cc) $(options) -o test/src/testLockList.o -c test/src/testLockList.cxx

test/src/testLongCell.o: test/src/testLongCell.cxx test/src/test.o src/longCell.o src/foldalign.hxx src/mbllist.o
	$(cc) $(options) -o test/src/testLongCell.o -c test/src/testLongCell.cxx

test/src/testLongTermMemory.o: test/src/testLongTermMemory.cxx test/src/test.o src/longCell.o src/foldalign.hxx src/helper.o src/two_link.o src/exception.o src/lockList.o src/longTermMemory.o
	$(cc) $(options) -o test/src/testLongTermMemory.o -c test/src/testLongTermMemory.cxx

test/src/testMatrix2d.o: test/src/testMatrix2d.cxx test/src/test.o src/matrix2d.o
	$(cc) $(options) -o test/src/testMatrix2d.o -c test/src/testMatrix2d.cxx

test/src/testMatrix4d.o: test/src/testMatrix4d.cxx test/src/test.o src/matrix4d.o
	$(cc) $(options) -o test/src/testMatrix4d.o -c test/src/testMatrix4d.cxx

test/src/testOptions.o: test/src/testOptions.cxx test/src/test.o src/options.o src/foldalign.hxx src/convertString.o
	$(cc) $(options) -o test/src/testOptions.o -c test/src/testOptions.cxx

test/src/testPrune.o: test/src/testPrune.cxx src/prune.o test/src/test.o
	$(cc) $(options) -o test/src/testPrune.o -c test/src/testPrune.cxx

test/src/testPruneTable.o: test/src/testPruneTable.cxx test/src/test.o src/pruneTable.o
	$(cc) $(options) -o test/src/testPruneTable.o -c test/src/testPruneTable.cxx

test/src/testReadConstraints.o: test/src/testReadConstraints.cxx test/src/test.o src/readConstraints.o src/foldalign.hxx src/jl.o
	$(cc) $(options) -o test/src/testReadConstraints.o -c test/src/testReadConstraints.cxx

test/src/testReadfile.o: test/src/testReadfile.cxx test/src/test.o src/readfile.o src/exception.o
	$(cc) $(options) -o test/src/testReadfile.o -c test/src/testReadfile.cxx

test/src/testScorematrix.o: test/src/testScorematrix.cxx test/src/test.o src/scorematrix.o
	$(cc) $(options) -o test/src/testScorematrix.o -c test/src/testScorematrix.cxx

test/src/testTwo_link.o: test/src/testTwo_link.cxx test/src/test.o src/two_link.o src/longCell.o
	$(cc) $(options) -o test/src/testTwo_link.o -c test/src/testTwo_link.cxx

test/src/testTest.o: test/src/testTest.cxx test/src/test.o src/exception.o
	$(cc) $(options) -o test/src/testTest.o -c test/src/testTest.cxx
