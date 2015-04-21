all: experimenter.cpp
	g++ -g -std=c++11 -O0 experimenter.cpp -o experimenter

clean:
	rm -rf *.o *.aux *.log *.out experimenter
