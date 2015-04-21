#include <sys/time.h>
#include <sys/signal.h>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <ctime>
#include <locale>
using namespace std;

/// Just to shorten things
typedef pair<string, string> info;

/// Just to shorten things again
#define VALIDATE(_param_, key)	if (_param_.first != key) { \
									fprintf(stderr, "[%s:%d] It was impossible to parse the configuration file.\n", __FUNCTION__, __LINE__); \
									exit(1); \
								}


/// Return the current time in "seconds.microseconds";
double get_wall_time(){
    struct timeval time;

    if (gettimeofday(&time,NULL)) return 0;

	return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

void ReplaceAll(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where ***'to'*** is a substring of 'from'
    }
}

/// yeah, might be inefficient buttt... itemsToSort is a very small "list"
string getOrderedByValueStr(map<string, int> itemsToSort) {
	vector<string> sortedItems(itemsToSort.size());
	string result;

	for (auto& item : itemsToSort)
		sortedItems[item.second] = item.first;

	for (auto& item : sortedItems)
		result += item + ",";

	return result;
}

string todayDate() {
	std::locale::global(std::locale("pt_BR.utf8"));
	std::time_t t = std::time(NULL);
	char mbstr[100];

	if (std::strftime(mbstr, sizeof(mbstr), "%Y_%m_%d_%H_%M_%S", std::localtime(&t))) {
		return string(mbstr);
	}
	else {
		fprintf(stderr, "[%s:%d] It was impossible to obtain the current date.\n", __FUNCTION__, __LINE__);
		exit(1);
	}
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

	double getSampleStandardDeviation(string row, string col) {
		double average = getAverage(row, col);
		int i = (*rowNames_)[row];
		int j = (*colNames_)[col];
		int l = num_values_[i*numCols_ + j];

		double sum = 0;
		for (int ind=0; ind<l; ind++) {
			double item = getMeasure(i, j, ind);
			double diff = item - average;
			sum += (diff * diff);
		}

		return sqrt( sum / (l-1) );
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
	int id_;
	string name_;
	string hypothesis_;
	string description_;
	int num_samples_;

	string plot_type_;
	string plot_title_;
	string plot_xlabel_;
	string plot_ylabel_;
	string plot_caption_;

	ResultTable results_;

	vector<Command> cmds_;

	map<string, int> implementations;
	map<string, int> xlabels;

	const string EXP_ID_KEY				{ "@@@EXP-ID@@@" };
	const string PLOT_TITLE_KEY			{ "@@@PLOT-TITLE@@@" };
	const string PLOT_XLABEL_KEY		{ "@@@PLOT-XLABEL@@@" };
	const string PLOT_YLABEL_KEY		{ "@@@PLOT-YLABEL@@@" };
	const string PLOT_CAPTION_KEY		{ "@@@PLOT-CAPTION@@@" };
	const string PLOT_NUM_SAMPLES_KEY	{ "@@@NUM-SAMPLES@@@" };
	const string PLOT_XLABELS_LIST_KEY	{ "@@@PLOT-X-LABELS-LIST@@@" };
	const string PLOT_IMPLS_LIST_KEY	{ "@@@PLOT-IMPLS-LIST@@@" };
	const string PLOT_EXP_RES_TABLE_KEY	{ "@@@EXP-RES-TABLE@@@" };
	const string BEG_PLOT_LINE_KEY		{ "@@@BEG-PLOT-LINE@@@" };
	const string END_PLOT_LINE_KEY		{ "@@@END-PLOT-LINE@@@" };
	const string PLOT_LINES_KEY			{ "@@@PLOT-LINES@@@" };
	const string PLOT_ADD_PLOT_LINE_KEY	{ "@@@LINE-COLOR@@@" };
	const string PLOT_ADD_PLOT_FILL_KEY	{ "@@@FILL-COLOR@@@" };
	const string PLOT_ADD_PLOT_COUNTER_KEY	{ "@@@COUNTER@@@" };
	const string PLOT_ADD_PLOT_MARKER_KEY	{ "@@@MARKER@@@" };


public:
	void id(int pid) {
		this->id_ = pid;
	}

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

	void plot_caption(info p) {
		VALIDATE(p, "PLOT_CAPTION");
		this->plot_caption_ = p.second;
	}

	int id() { return this->id_; }

	string name() { return this->name_ ; }

	string hypothesis() { return this->hypothesis_ ; }

	string description() { return this->description_ ; }

	int num_samples() { return this->num_samples_ ; }

	string plot_type() { return this->plot_type_ ; }

	string plot_title() { return this->plot_title_ ; }

	string plot_xlabel() { return this->plot_xlabel_; }

	string plot_ylabel() { return this->plot_ylabel_; }

	string plot_caption() { return this->plot_caption_; }

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
				cout << "\t\t\t\tExecuting sample " << sampleNum << ": ";
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

				cout << " took " << wallTime << " seconds. " << endl;

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

	string getTexResultTable() {
		stringstream results;

		results << "x \t";

		for (int i=0; i<implementations.size(); i++)
			results << "c" << i << " \t c" << i << "_err \t";
		results << endl;

		for (auto& lbl : xlabels) {
			results << lbl.first << " \t";

			for (auto& impl : implementations)
				results << this->results_.getAverage(impl.first, lbl.first) << " \t " <<
						   this->results_.getSampleStandardDeviation(impl.first, lbl.first) << " \t";

			results << endl;
		}

		return results.str();
	}

	/// These "+1"/"-1" all around are for taking care of the "%" latex comment character
	string findAndPathAddPlotLine(string& plotText) {
		auto secBegin 		= plotText.find(BEG_PLOT_LINE_KEY);
		auto secEnd 		= plotText.find(END_PLOT_LINE_KEY);

		string plotLineTmpl = plotText.substr(secBegin+BEG_PLOT_LINE_KEY.size(), secEnd - secBegin-BEG_PLOT_LINE_KEY.size() - 1);

		plotText.replace(secBegin - 1, secEnd + END_PLOT_LINE_KEY.size() + 1 - secBegin, PLOT_LINES_KEY);

		return plotLineTmpl;
	}

	string genPlotTex(string texTemplate) {
		const int NUM_DESIGNS = 18;

		string markers[NUM_DESIGNS] = {	"o", "star", "oplus", "otimes", "square", "triangle", "diamond", "halfdiamond", "halfcircle", "pentagon",
										"*", "star*", "oplus*", "otimes*", "square*", "triangle*", "diamond*", "pentagon*"};

		string colors[NUM_DESIGNS] = {	"red", "blue", "yellow", "green", "violet", "orange", "pink", "cyan", "olive", "magenta",
										"gray", "brown", "purple", "lightgray", "teal", "lime", "black", "darkgray"};

		string plotText = texTemplate;

		ReplaceAll(plotText, EXP_ID_KEY, std::to_string(this->id()));
		ReplaceAll(plotText, PLOT_TITLE_KEY, this->plot_title());
		ReplaceAll(plotText, PLOT_YLABEL_KEY, this->plot_ylabel());
		ReplaceAll(plotText, PLOT_XLABEL_KEY, this->plot_xlabel());
		ReplaceAll(plotText, PLOT_CAPTION_KEY, this->plot_caption());
		ReplaceAll(plotText, PLOT_NUM_SAMPLES_KEY, std::to_string(this->num_samples()));

		string implList 	= getOrderedByValueStr( this->implementations );
		string xlabelsList 	= getOrderedByValueStr( this->xlabels );

		ReplaceAll(plotText, PLOT_XLABELS_LIST_KEY, xlabelsList);
		ReplaceAll(plotText, PLOT_IMPLS_LIST_KEY, implList);

		/// Adds the result table to the .tex
		string texResTable	= getTexResultTable();

		ReplaceAll(plotText, PLOT_EXP_RES_TABLE_KEY, texResTable);

		/// Now we are going to add the "\addplot" lines. First we identify
		/// the position and later we add them.
		/// The function below will find the BEG/END of a line plot. Extract
		/// the addplot template and return it and, replace the beg\end region
		/// with a new marker @@@PLOT_LINES@@@
		string addPlotTemplate = findAndPathAddPlotLine(plotText);
		string plotLines;

		for (int i=0, disId=0; i<implementations.size(); i++, disId++) {
			string newLine = addPlotTemplate;
			disId = disId % NUM_DESIGNS;

			ReplaceAll(newLine, PLOT_ADD_PLOT_LINE_KEY, colors[disId]);
			ReplaceAll(newLine, PLOT_ADD_PLOT_FILL_KEY, colors[disId] + "!60");
			ReplaceAll(newLine, PLOT_ADD_PLOT_COUNTER_KEY, std::to_string(i));
			ReplaceAll(newLine, PLOT_ADD_PLOT_MARKER_KEY, markers[disId]);

			plotLines += newLine;
		}

		/// Add all plot lines at once in the plot.
		ReplaceAll(plotText, PLOT_LINES_KEY, plotLines);

		/// returns the text of the new plot
		return plotText;
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
	const string TITLE_KEY 				{ "@@@TITLE@@@"};
	const string ABSTRACT_KEY 			{ "@@@ABSTRACT@@@"};
	const string EXPERIMENTER_KEY 		{ "@@@EXPERIMENTER@@@"};
	const string EXPERIMENTER_EMAIL_KEY { "@@@EXPERIMENTER-EMAIL@@@"};
	const string INTRODUCTION_KEY		{ "@@@INTRODUCTION@@@" };
	const string ARCH_BRAND_KEY			{ "@@@ARCH-BRAND-NAME@@@" };
	const string OS_UNAME_KEY			{ "@@@OS-UNAME@@@" };
	const string RES_SUBSECS_KEY		{ "@@@RES-SUBSECS@@@" };

	const string EXP_NAME_KEY			{ "@@@EXP-NAME@@@" };
	const string EXP_HYPOTHESIS_KEY		{ "@@@EXP-HYPOTHESIS@@@" };
	const string EXP_SHORT_DESCRIPTION_KEY	{ "@@@EXP-SHORT-DESCRIPTION@@@" };
	const string EXP_PLOT_CODE_KEY		{ "@@@EXP-PLOT-CODE@@@" };


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

	string abstract() { return this->abstract_ ; }

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

	string readFile(string path) {
		ifstream inpFile (path);
		string line;
		stringstream ss;

		if (inpFile.is_open()) {
			while ( getline(inpFile, line) )
				ss << line << "\n";

			inpFile.close();

			return ss.str();
		}
		else {
			fprintf(stderr, "[%s:%d] It was impossible to read the file %s.\n", __FUNCTION__, __LINE__, path.c_str());
			exit(1);
		}
	}

	void writeToFile(string path, string content) {
		ofstream outFile (path);

		if (outFile.is_open()) {
			outFile << content << endl;
			outFile.close();
		}
		else {
			fprintf(stderr, "[%s:%d] It was impossible to read the file %s.\n", __FUNCTION__, __LINE__, path.c_str());
			exit(1);
		}
	}

	string commandOutput(string command) {
		FILE* fp = popen(command.c_str(), "r");
		stringstream ss;
		char buff[4096];

		if (!fp) {
			fprintf(stderr, "[%s:%d] It was impossible to execute the command or capture its output: %s.\n", __FUNCTION__, __LINE__, command.c_str());
			exit(1);
		}

		while(fgets(buff, sizeof(buff), fp)!=NULL)
			ss << buff;

		pclose(fp);

		string aux = ss.str();

		/// Escape for latex special characters
		ReplaceAll(aux, string("_"), string("\\_"));
		ReplaceAll(aux, string("#"), string("\\#"));

		return aux;
	}

	void printToLatex() {
		/// Read templates
		string mainFile_TemplateText 		= readFile("templates/template1/template.tex");
		string resultSection_TemplateText 	= readFile("templates/template1/ressec_template.tex");
		string linePlot_TemplateTex 		= readFile("templates/template1/lineplot_template.tex");
		string barPlot_TemplateTex 			= readFile("templates/template1/barplot_template.tex");
		string reportText					= mainFile_TemplateText;

		/// Produce the architecture representation image
		int rsys = system("lstopo reports/architecture.pdf");

		/// Update information in the main template
		ReplaceAll(reportText, TITLE_KEY, this->title());
		ReplaceAll(reportText, ABSTRACT_KEY, this->abstract());
		ReplaceAll(reportText, EXPERIMENTER_KEY, this->experimenter());
		ReplaceAll(reportText, EXPERIMENTER_EMAIL_KEY, this->experimenter_email());

		/// If the user specified an "intro.tex" file then it is included
		if (this->intro_tex_path() != "") {
			string introText = readFile(this->intro_tex_path());
			ReplaceAll(reportText, INTRODUCTION_KEY, introText);
		}

		/// Append architecture name, frequency, code; append kernel version, machine name, etc.
		ReplaceAll(reportText, OS_UNAME_KEY, commandOutput("uname -a"));
		ReplaceAll(reportText, ARCH_BRAND_KEY, commandOutput("cat /proc/cpuinfo | grep \"model name\" | cut -d':' -f 2 | head -n 1"));

		/// Output a results subsection for each experiment
		string subSecsText;
		for (auto& exp : experiments()) {
			/// fresh new template for the subsection
			string section 		= resultSection_TemplateText;
			string plotTemplate = (exp.plot_type() == "BARS" ? barPlot_TemplateTex : linePlot_TemplateTex);
			string plotText 	= exp.genPlotTex(plotTemplate);

			ReplaceAll(section, EXP_NAME_KEY, 				exp.name());
			ReplaceAll(section, EXP_HYPOTHESIS_KEY, 		exp.hypothesis());
			ReplaceAll(section, EXP_SHORT_DESCRIPTION_KEY, 	exp.description());
			ReplaceAll(section, EXP_PLOT_CODE_KEY, 			plotText);

			/// concatenate the subsection result template filled
			subSecsText += section;
		}

		/// Write all the .tex for the results section in report.tex
		ReplaceAll(reportText, RES_SUBSECS_KEY, subSecsText);

		/// Write the .tex and compile the project
		string today = todayDate();
		writeToFile("reports/report_" + today + ".tex", reportText);

		/// Produce the ".pdf" report and remove temporary files.
		commandOutput("cd reports && pdflatex report_" + today + ".tex");
		commandOutput("cd reports && rm -rf *.aux *.log *.out *.csv");
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

	fprintf(stderr, "It was impossible to parse the configuration file.");
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
		exp.id(i);
		exp.name(nextInputPair());
		exp.hypothesis(nextInputPair());
		exp.description(nextInputPair());
		exp.num_samples(nextInputPair());

		/// Plot configuration
		exp.plot_type(nextInputPair());
		exp.plot_title(nextInputPair());
		exp.plot_xlabel(nextInputPair());
		exp.plot_ylabel(nextInputPair());
		exp.plot_caption(nextInputPair());

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
		fprintf(stderr, "It was impossible to parse the configuration file.\n");
		exit(1);
	}
	fclose(stdin);


	/// Execute the commands, measure their wall time and store in result table
	proj.execute();

	/// Just dump the results table to stdout for debugging
	proj.printResultsToStdout();

	/// Emit latex source code
	proj.printToLatex();

	return 0;
}
