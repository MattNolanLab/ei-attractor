# Makefile for the ARC project

PACK_FILES = ARCtest.sh \
			 args.cc	\
			 args.h		\
			 CellularGA.cc	\
			 CellularGA.h	\
			 Chromosome.cc	\
			 ChromosomeFactory.h	\
			 Chromosome.h	\
			 debug.h		\
			 error.cc		\
			 error.h		\
			 global.h		\
			 Makefile		\
			 mpi			\
			 OneMaxChromosome.cc	\
			 OneMaxChromosome.h		\
			 README			\
			 README.ARC		\
			 SimpleBitArray.h	\
			 testOneMax.cc

BIN = testOneMax *.o
DEPFILE = Makefile.dep

BASE = args.o error.o CellularGA.o Chromosome.o OneMaxChromosome.o \
	   NeuralNetwork.o 
#NNChromosome.o NNOneCharChromosome.o RecNeuralNet.o NNTrainData.o 

USE_MPI ?= 1
DEBUG ?= -DNDEBUG_ADDITIONAL
CXXADD ?= 

ifeq ($(USE_MPI), 1)
	CXX = mpicxx
	CXXFLAGS = -O3 -g -W -Wall $(DEBUG) -DARC_MPI_ENABLED
else
	CXX = g++
	CXXFLAGS = -O3 -g -W -Wall $(DEBUG)
endif

CXXFLAGS += $(CXXADD)


LDFLAGS = -lm

.PHONY:	clean deps

all: testNN

run: testOneMax
	./testOneMax

mpirun: testOneMax
	mpirun -np 2 ./testOneMax

testOneMax: $(BASE) testOneMax.o
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o testOneMax $(BASE) testOneMax.o

testNN: $(BASE) testNN.o
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o testNN $(BASE) testNN.o

kill:
	killall testOneMax

clean:
	rm -f $(BIN)

deps:
	g++ -MM *.cc > $(DEPFILE)

pack:
	tar cvzf xsolan00.tar.gz $(PACK_FILES)


CellularGA.o: CellularGA.cc CellularGA.h global.h debug.h Chromosome.h \
  ChromosomeFactory.h NNTrainData.h NeuralNetwork.h error.h \
  OneMaxChromosome.h SimpleBitArray.h NNChromosome.h
Chromosome.o: Chromosome.cc Chromosome.h debug.h
NNChromosome.o: NNChromosome.cc NNChromosome.h Chromosome.h \
  NeuralNetwork.h debug.h error.h NNTrainData.h global.h
NNOneCharChromosome.o: NNOneCharChromosome.cc NNOneCharChromosome.h \
  Chromosome.h NNChromosome.h NeuralNetwork.h debug.h error.h \
  NNTrainData.h global.h
NNTrainData.o: NNTrainData.cc NNTrainData.h NeuralNetwork.h debug.h \
  error.h
NeuralNetwork.o: NeuralNetwork.cc debug.h NeuralNetwork.h error.h \
  global.h
OneMaxChromosome.o: OneMaxChromosome.cc OneMaxChromosome.h \
  SimpleBitArray.h debug.h Chromosome.h
RecNeuralNet.o: RecNeuralNet.cc RecNeuralNet.h error.h debug.h
args.o: args.cc args.h CellularGA.h global.h debug.h Chromosome.h \
  ChromosomeFactory.h NNTrainData.h NeuralNetwork.h error.h \
  OneMaxChromosome.h SimpleBitArray.h NNChromosome.h
error.o: error.cc error.h
test.o: test.cc
testNN.o: testNN.cc args.h CellularGA.h global.h debug.h Chromosome.h \
  ChromosomeFactory.h NNTrainData.h NeuralNetwork.h error.h \
  OneMaxChromosome.h SimpleBitArray.h NNChromosome.h \
  NNOneCharChromosome.h
testOneMax.o: testOneMax.cc args.h CellularGA.h global.h debug.h \
  Chromosome.h ChromosomeFactory.h NNTrainData.h NeuralNetwork.h error.h \
  OneMaxChromosome.h SimpleBitArray.h NNChromosome.h
