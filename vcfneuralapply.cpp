#include "vcflib/Variant.h"
#include "doublefann.h"
#include "convert.h"
#include <fstream>
#include <iostream>
#include <getopt.h>
#include <vector>
#include "split/join.h"
#include <cmath>

using namespace std;
using namespace vcf;


double mean(const vector<double>& data) {
    double total = 0;
    for (vector<double>::const_iterator i = data.begin(); i != data.end(); ++i) {
        total += *i;
    }
    return total/data.size();
}

double median(vector <double>& data) {
    double median;
    size_t size = data.size();
    // ascending order
    sort(data.begin(), data.end());
    // get middle value
    if (size % 2 == 0) {
        median = (data[size/2-1] + data[size/2]) / 2;
    } else {
        median = data[size/2];
    }
    return median;
}

double variance(const vector <double>& data, const double mean) {
    double total = 0;
    for (vector <double>::const_iterator i = data.begin(); i != data.end(); ++i) {
        total += (*i - mean)*(*i - mean);
    }
    return total / (data.size());
}

double standard_deviation(const vector <double>& data, const double mean) {
    return sqrt(variance(data, mean));
}

struct Stats {
    double mean;
    double stdev;
    Stats(void) : mean(0), stdev(1) { }
};

// ann = fann_create_standard(num_layers, num_input, num_neurons_hidden, num_output);
bool load_ann_metadata(string& ann_metadata_file,
                       vector<string>& fields,
                       map<string, Stats>& stats) {
    ifstream in(ann_metadata_file.c_str());
    if (!in.is_open()) {
        return false;
    }
    string linebuf;
    while (getline(in, linebuf)) {
        // format is: field_name, mean, stdev
        vector<string> m = split(linebuf, "\t ");
        fields.push_back(m[0]);
        Stats& s = stats[m[0]];
        convert(m[1], s.mean);
        convert(m[2], s.stdev);
    }
    in.close();
    return true;
}

bool save_ann_metadata(string& ann_metadata_file,
                       vector<string>& fields,
                       map<string, Stats>& stats) {
    ofstream out(ann_metadata_file.c_str());
    if (!out.is_open()) {
        return false;
    }
    for (vector<string>::iterator f = fields.begin(); f != fields.end(); ++f) {
        Stats& s = stats[*f];
        out << *f << "\t" << s.mean << "\t" << s.stdev << endl;
    }
    out.close();
    return true;
}

void normalize_inputs(vector<double>& record, vector<string>& fields, map<string, Stats>& stats) {
    vector<double>::iterator r = record.begin();
    for (vector<string>::iterator f = fields.begin(); f != fields.end(); ++f, ++r) {
        Stats& s = stats[*f];
        *r = (*r - s.mean) / s.stdev;
    }
}

void read_fields(Variant& var, int ai, vector<string>& fields, vector<double>& record) {
    double td;
    vector<string>::iterator j = fields.begin();
    for (; j != fields.end(); ++j) {
        if (*j == "QUAL") { // special handling...
            td = var.quality;
        } else {
            if (var.info.find(*j) == var.info.end()) {
                td = 0;
            } else {
                if (var.vcf->infoCounts[*j] == 1) { // for non Allele-variant fields
                    convert(var.info[*j][0], td);
                } else {
                    convert(var.info[*j][ai], td);
                }
            }
        }
        record.push_back(td);
    }
}


void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options]" << endl
         << endl
         << "Applies an artificial neural network on a set of fields defined by an associated" << endl
         << "metadata file (ANN_FILE.fields), writes output on stdout, tagged by the ANN result." << endl
         << endl
         << "options:" << endl
         << "    -h, --help              this dialog." << endl
         << "    -v, --vcf-file FILE     specifies the vcf file (or BGZipped, vcf.gz) to use for application" << endl
         << "                            if '-' specified, stdin is used (default)" << endl
         << "    -a, --ann-file FILE     use the ANN in this file (required).  any number may be specified," << endl
         << "                            as an ensemble.  the results are averaged at runtime." << endl
         << "    -o, --output-tag TAG    output the results of the neural network execution as TAG=..." << endl
         << "                            writes results to QUAL by default." << endl
         << "    -i, --info STRING       add descriptive information to the VCF header output TAG definition" << endl
         //<< "    -r, --region          specify a region on which to target the analysis, requires a BGZF" << endl
         //<< "                          compressed file which has been indexed with tabix.  any number of" << endl
         //<< "                          regions may be specified." << endl
         << endl;
}

#define PHRED_MAX 50000.0

long double double2phred(long double prob) {
    if (prob == 1)
        return PHRED_MAX;  // guards against "-0"
    long double p = -10 * (long double) log10(prob);
    if (p < 0 || p > PHRED_MAX) // int overflow guard
        return PHRED_MAX;
    else
        return p;
}

int main(int argc, char** argv)
{   

    unsigned int num_output = 1;
    //double* din = malloc(num_input * num_data * sizeof(double));
    //double* dout = malloc(num_output * num_data * sizeof(double));

    string variantFileName = "-";
    string information;

    bool useQUAL = false;
    bool writeQual = true;

    vector<string> annFiles;

    string outputTag;

    int c;

    if (argc == 1) {
        printSummary(argv);
        exit(1);
    }

    while (true) {
        static struct option long_options[] =
        {  
            /* These options set a flag. */
            //{"verbose", no_argument,       &verbose_flag, 1},
            {"help", no_argument, 0, 'h'},
            {"vcf-file", required_argument, 0, 'v'},
            {"ann-file", required_argument, 0, 'a'},
            {"output-tag", required_argument, 0, 'o'},
            {"info", required_argument, 0, 'i'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hv:a:o:i:",
                         long_options, &option_index);

        if (c == -1)
            break;

        string field;

        switch (c) {

	    case 'v':
                variantFileName = optarg;
                break;

            case 'a':
                annFiles.push_back(optarg);
                break;

            case 'o':
                outputTag = optarg;
                writeQual = false;
                break;

            case 'i':
                information = optarg;
                break;

            case 'h':
                printSummary(argv);
                exit(0);
                break;

            default:
                break;
        }
    }

    if (annFiles.empty()) {
        cerr << "please supply a filename for the neural network, --ann-file" << endl;
        exit(1);
    }

    // now populate the data with VCF input
    VariantCallFile variantFile;
    if (variantFileName == "-") {
        variantFile.open(std::cin);
    } else {
        variantFile.open(variantFileName);
    }

    if (!variantFile.is_open()) {
        return 1;
    } 

    // check that all the fields files are the same
    vector<struct fann*> anns;
    vector<vector<string> > fields;
    vector<map<string, Stats> > stats;

    int i = 0;
    for (vector<string>::iterator annFile = annFiles.begin(); annFile != annFiles.end(); ++annFile, ++i) {
        vector<string> theseFields;
        map<string, Stats> theseStats;
        string annMetadataFile = *annFile + ".meta";
        if (!load_ann_metadata(annMetadataFile, theseFields, theseStats)) {
            cerr << "could not open " << annMetadataFile << " to read metadata" << endl;
            exit(1);
        }

        struct fann *ann = fann_create_from_file(annFile->c_str());
        anns.push_back(ann);
        fields.push_back(theseFields);
        stats.push_back(theseStats);
    }

    fann_type *calc_out;

    Variant var(variantFile);

    if (!outputTag.empty()) {
        string infostr;
        if (!information.empty()) {
            infostr = ".  " + information;
            variantFile.addHeaderLine("##INFO=<ID=" + outputTag + ",Number=A,Type=Float,Description=\"Probability given model described by " + join(annFiles, " ") + infostr + "\">");
        }
    }

    cout << variantFile.header << endl;

    while (variantFile.getNextVariant(var)) {
        if (var.info.find(outputTag) != var.info.end()) {
            var.info[outputTag].clear();
        }
        long double newqual = 0;
        // for alt, get the index,
        for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
            string& alt = *a;
            int altindex = var.getAltAlleleIndex(alt);
            long double fresult = 0; // floating point result, sum of 

            int j = 0;
            for (vector<struct fann*>::iterator ann = anns.begin(); ann != anns.end(); ++ann, ++j) {
                map<string, Stats>& cstats = stats.at(j);
                vector<string>& cfields = fields.at(j);
                vector<double> record;
                read_fields(var, altindex, cfields, record);
                normalize_inputs(record, cfields, cstats);
                vector<fann_type> input; //[fields.size()];
                for (vector<double>::iterator d = record.begin(); d != record.end(); ++d) {
                    input.push_back(*d);
                }
                calc_out = fann_run(*ann, &input[0]);
                fresult += calc_out[0];
            }

            // take mean, and convert to phred
            long double result = double2phred(1 - (fresult / anns.size()));

            if (!outputTag.empty()) {
                var.info[outputTag].push_back(convert(result));
            }
            if (writeQual) {
                if (newqual == 0) {
                    newqual = result;
                } else {
                    newqual = min(newqual, result);
                }
            }
        }
        if (writeQual) {
            var.quality = newqual;
        }

        cout << var << endl;

    }

    for (vector<struct fann*>::iterator ann = anns.begin(); ann != anns.end(); ++ann) {
        fann_destroy(*ann);
    }

    return 0;

}
