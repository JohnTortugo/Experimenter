#include <sys/time.h>
#include <sys/signal.h>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
using namespace std;

/// Just to shorten things
typedef pair<string, string> info;

/// Just to shorten things again
#define VALIDATE(_param_, key)	if (_param_.first != key) { \
									fprintf(stderr, "[%s:%d] It was not possible to parse the configuration file.\n", __FUNCTION__, __LINE__); \
									exit(1); \
								}


/// Return the current time in "seconds.microseconds";
double get_wall_time(){
    struct timeval time;

    if (gettimeofday(&time,NULL)) return 0;

	return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

class ResultTable {
private:
	/// tells whether the matrix is initialized
	bool initialized_ = false;

	/// number of rows, cols and depth of the measurement matrix/numvalues.
	int numRows_ = -1, numCols_ = -1, depth_ = -1;

	/// tri-dimensional matrix for storing results
	double* measurements_ = nullptr;

	/// this will store the actual number of measurements for each cell of measurements_
	int* num_values_ = nullptr;

	/// map from row/col "name" to index.
	map<string, int>* rowNames_;
	map<string, int>* colNames_;

	void checkInitialized() {
		if (!initialized_) {
			fprintf(stderr, "[%s:%d] Using an untilized results matrix!!!\n", __FUNCTION__, __LINE__);
			exit(1);
		}
	}

	void setMeasure(int i, int j, int k, double val) {
		this->measurements_[i*numRows_*depth_ + j*depth_ + k] = val;
	}

	double getMeasure(int i, int j, int k) {
		return this->measurements_[i*numRows_*depth_ + j*depth_ + k];
	}

	void addMeasure(int i, int j, double val) {
		int k = this->num_values_[i*numCols_ + j];
		this->measurements_[i*numRows_*depth_ + j*depth_ + k] = val;
		this->num_values_[i*numCols_ + j]++;
	}

	void setNumValues(int i, int j, int val) {
		this->num_values_[i*numCols_ + j] = val;
	}

public:
	ResultTable() : rowNames_(nullptr), colNames_(nullptr)
	{ }

	~ResultTable() {
		delete[] this->measurements_;
		delete[] this->num_values_;
	}

	void initialize(map<string, int>& rowNames, map<string, int>& colNames, int depth) {
		this->rowNames_ = new map<string, int>(rowNames);
		this->colNames_ = new map<string, int>(colNames);
		this->numCols_ = colNames_->size();
		this->numRows_ = rowNames_->size();
		this->depth_ = depth;

		this->measurements_ = new double[this->numRows_ * this->numCols_ * depth];
		this->num_values_ 	= new int[this->numRows_ * this->numCols_];

		for (int i=0; i<rowNames_->size(); i++) {
			for (int j=0; j<colNames_->size(); j++) {
				for (int k=0; k<depth; k++) {
					setMeasure(i, j, k, 0);
				}

				setNumValues(i, j, 0);
			}
		}

		this->initialized_ 	= true;
	}

	void addMeasure(string row, string col, double val) {
		if (rowNames_->find(row) == rowNames_->end() || colNames_->find(col) == colNames_->end()) {
			fprintf(stderr, "[%s:%d] Trying to add a measure in a row/col not existent! <%s, %s>.\n", __FUNCTION__, __LINE__, row.c_str(), col.c_str());
			exit(1);
		}

		int i = (*rowNames_)[row];
		int j = (*colNames_)[col];

		addMeasure(i, j, val);
	}

	double getAverage(string row, string col) {
		int i = (*rowNames_)[row];
		int j = (*colNames_)[col];
		int l = num_values_[i*numCols_ + j];

		double sum = 0;
		for (int ind=0; ind<l; ind++) {
			sum += getMeasure(i, j, ind);
		}

		return sum / l;
	}

};

/// Represent one command to be executed
class Command {
private:
	string impl_;
	string xlabel_;
	string cmd_;

public:
	void impl(info p) {
		VALIDATE(p, "IMPL_NAME");
		this->impl_ = p.second;
	}

	void xlabel(info p) {
		VALIDATE(p, "XLABEL");
		this->xlabel_ = p.second;
	}

	void cmd(info p) {
		VALIDATE(p, "CMD");
		this->cmd_ = p.second;
	}

	string impl() { return this->impl_; }

	string xlabel() { return this->xlabel_; }

	string cmd() { return this->cmd_; }
};

/// Represent each of the experiments to be executed in the project
class Experiment {
private:
	string name_;
	string hypothesis_;
	string description_;
	int num_samples_;

	string plot_type_;
	string plot_title_;
	string plot_xlabel_;
	string plot_ylabel_;

	ResultTable results_;

	vector<Command> cmds_;

	map<string, int> implementations;
	map<string, int> xlabels;

public:
	void name(info p) {
		VALIDATE(p, "EXP_NAME");
		this->name_ = p.second;
	}

	void hypothesis(info p) {
		VALIDATE(p, "HYPOTHESIS");
		this->hypothesis_ = p.second;
	}

	void description(info p) {
		VALIDATE(p, "SHORT_DESCRIPTION");
		this->description_ = p.second;
	}

	void num_samples(info p) {
		VALIDATE(p, "EXECUTIONS");
		this->num_samples_ = std::stoi(p.second);
	}

	void plot_type(info p) {
		VALIDATE(p, "PLOT_TYPE");
		this->plot_type_ = p.second;
	}

	void plot_title(info p) {
		VALIDATE(p, "PLOT_TITLE");
		this->plot_title_ = p.second;
	}

	void plot_xlabel(info p) {
		VALIDATE(p, "PLOT_XLABEL");
		this->plot_xlabel_ = p.second;
	}

	void plot_ylabel(info p) {
		VALIDATE(p, "PLOT_YLABEL");
		this->plot_ylabel_ = p.second;
	}

	string name() { return this->name_ ; }

	string hypothesis() { return this->hypothesis_ ; }

	string description() { return this->description_ ; }

	int num_samples() { return this->num_samples_ ; }

	string plot_type() { return this->plot_type_ ; }

	string plot_title() { return this->plot_title_ ; }

	string plot_xlabel() { return this->plot_xlabel_; }

	string plot_ylabel() { return this->plot_ylabel_; }

	vector<Command>& cmds() { return this->cmds_; }

	void addCommand(Command& cmd) {
		cmds_.push_back(cmd);
	}

	void addResult(string impl, string xlabel, double val) {
		this->results_.addMeasure(impl, xlabel, val);
	}

	void initResults() {
		int numImpls = 0, numXlbls = 0;

		for (auto& c: cmds_) {
			if (implementations.find(c.impl()) == implementations.end()) {
				implementations[c.impl()] = numImpls++;
			}

			if (xlabels.find(c.xlabel()) == xlabels.end()) {
				xlabels[c.xlabel()] = numXlbls++;
			}
		}

		this->results_.initialize(implementations, xlabels, this->num_samples_);
	}

	void execute() {
		cout << endl << "\t\tExecuting experiment: " << name_ << endl;

		for (auto& cmd : cmds_) {
			cout << "\t\t\tExecuting command: " << cmd.cmd() << endl;

			for (int sampleNum=0; sampleNum<num_samples_; sampleNum++) {
				cout << "\t\t\t\tExecuting sample " << sampleNum << endl;
				double startTime 	= get_wall_time();

				int retCode 		= system( cmd.cmd().c_str() );

				double endTime 		= get_wall_time();
				double wallTime 	= endTime - startTime;

				if (retCode == -1) {
					fprintf(stderr, "[%s:%d] Some error occurred while executing command: %s\n", __FUNCTION__, __LINE__, cmd.cmd().c_str());
					exit(1);
				}

				if ( WIFSIGNALED(retCode) && (WTERMSIG(retCode) == SIGINT || WTERMSIG(retCode) == SIGQUIT) ) {
					fprintf(stderr, "[%s:%d] Signal for stopping captured. Aborting!\n", __FUNCTION__, __LINE__); \
					exit(1);
				}

				addResult(cmd.impl(), cmd.xlabel(), wallTime);
			}
		}
	}

	void printResultsToStdout() {
		cout << endl << "\tResults of experiment: " << name_ << endl;

		printf("%13s ", " ");

		for (auto& lbl : xlabels) {
			printf("%13s ", lbl.first.c_str());
		}

		printf("\n");

		for (auto& impl : implementations) {
			printf("%13s ", impl.first.c_str());

			for (auto& lbl : xlabels) {
				printf("%13.5lf ", this->results_.getAverage(impl.first, lbl.first));
			}

			printf("\n");
		}

	}
};

/// Represent the whole project
class Project {
private:
	string title_;
	string abstract_;
	string experimenter_;
	string experimenter_email_;

	string intro_tex_path_;

	vector<Experiment> exps_;
public:
	void title(info p) {
		VALIDATE(p, "TITLE");
		this->title_ = p.second;
	}

	void abstract(info p) {
		VALIDATE(p, "ABSTRACT");
		this->abstract_ = p.second;
	}

	void experimenter(info p) {
		VALIDATE(p, "EXPERIMENTER");
		this->experimenter_ = p.second;
	}

	void experimenter_email(info p) {
		VALIDATE(p, "EXPERIMENTER_EMAIL");
		this->experimenter_email_ = p.second;
	}

	void intro_tex_path(info p) {
		VALIDATE(p, "INTRODUCTION");
		this->intro_tex_path_ = p.second;
	}

	string title() { return this->title_ ; }

	string experimenter() { return this->experimenter_ ; }

	string experimenter_email() { return this->experimenter_email_ ; }

	string intro_tex_path() { return this->intro_tex_path_ ; }

	vector<Experiment>& experiments() {
		return this->exps_;
	}

	void addExperiment(Experiment& expr) {
		exps_.push_back(expr);
	}

	/// Execute all commands in all experiments in the project
	void execute() {
		cout << "Starting execution of project " << title_ << endl;

		for (auto& exp : exps_) {
			exp.initResults();

			exp.execute();
		}
	}

	void printResultsToStdout() {
		cout << endl << endl << "Results of the project \"" << title_ << "\"" << endl;

		for (auto& exp : exps_) {
			exp.printResultsToStdout();
		}
	}
};

info nextInputPair() {
	const int maxSize = 1000;
	char line[maxSize], key[maxSize], value[maxSize];

	/// Read the next line which is not a blank line
	while (fgets(line, maxSize, stdin) != NULL) {
		/// skip blank lines
		if (strlen(line) == 1) continue;

		/// parse the KEY=VALUE
		sscanf(line, "%[^=]=%[^\n]\n", key, value);

		/// finishes the loop and returns the pair
		return make_pair(string(key), string(value));
	}

	fprintf(stderr, "It was not possible to parse the configuration file.");
	exit(1);
}

bool parseInput(Project& proj) {
	/// Parse project header information
	proj.title(nextInputPair());
	proj.abstract(nextInputPair());
	proj.experimenter(nextInputPair());
	proj.experimenter_email(nextInputPair());
	proj.intro_tex_path(nextInputPair());

	/// Lets start parsing the experiments
	info numExpsInfo = nextInputPair();
	VALIDATE(numExpsInfo, "NUM_EXPS");

	/// for each experiment
	for (int i=0; i < std::stoi(numExpsInfo.second); i++ ) {
		Experiment exp;

		/// Experiments metadata
		exp.name(nextInputPair());
		exp.hypothesis(nextInputPair());
		exp.description(nextInputPair());
		exp.num_samples(nextInputPair());

		/// Plot configuration
		exp.plot_type(nextInputPair());
		exp.plot_title(nextInputPair());
		exp.plot_xlabel(nextInputPair());
		exp.plot_ylabel(nextInputPair());

		/// Lets start parsing the commands
		info numCmdsInfo = nextInputPair();
		VALIDATE(numCmdsInfo, "NUM_CMDS");

		/// for each command
		for (int j=0; j < std::stoi(numCmdsInfo.second); j++ ) {
			Command cmd;

			/// Command data
			cmd.impl(nextInputPair());
			cmd.xlabel(nextInputPair());
			cmd.cmd(nextInputPair());

			exp.addCommand(cmd);
		}

		proj.addExperiment(exp);
	}

	return true;
}

/// Execute the experiments and produce a report
int main(int argc, char* argv[]) {
	Project proj;

	if (argc < 2) {
		fprintf(stderr, "You should specify an input configuration file.\n");
		exit(1);
	}

	/// Open and parse the configuration file.
	/// Informations are stored in the "proj" variable.
	if (!freopen(argv[1], "r", stdin) || !parseInput(proj)) {
		fprintf(stderr, "It was not possible to parse the configuration file.\n");
		exit(1);
	}
	fclose(stdin);


	/// Execute the commands, measure their wall time and store in result table
	proj.execute();

	/// Just dump the results table to stdout for debugging
	proj.printResultsToStdout();

	/// Emit latex source code
	proj.printToLatex();

	/// compile latex to pdf
	proj.latexToPdf();

	return 0;
}
