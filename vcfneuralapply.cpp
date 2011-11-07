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


void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options]" << endl
         << endl
         << "Applies an artificial neural network on a set of fields defined by an associated" << endl
         << "metadata file (ANN_FILE.fields), writes output on stdout, tagged by the ANN result." << endl
         << endl
         << "options:" << endl
         << "    -h, --help              this dialog." << endl
         << "    -v, --vcf-file FILE     specifies the vcf file (or BGZipped, vcf.gz) to use for training" << endl
         << "                            if '-' specified, stdin is used (default)" << endl
         << "    -a, --ann-file FILE     save the ANN to this file (required).  metadata, specifically" << endl
         << "                            the VCF INFO fields which are used, is saved as FILE.fields." << endl
         << "                            Any number may be specified, the results are averaged at runtime." << endl
         << "    -o, --output-tag TAG    output the results of the neural network execution as TAG=..." << endl
         << "                            writes results to QUAL by default." << endl
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
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hv:a:o:",
                         long_options, &option_index);

        if (c == -1)
            break;

        string field;

        switch (c)
        {

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
    vector<string> fields;

    for (vector<string>::iterator annFile = annFiles.begin(); annFile != annFiles.end(); ++annFile) {
        ifstream fieldsFile;
        fieldsFile.open(string(*annFile + ".fields").c_str());
        string line;
        getline(fieldsFile, line);
        vector<string> theseFields = split(line, "\t");
        if (theseFields.empty()) {
            cerr << "no fields provided, cannot execute neural net" << endl;
            exit(1);
        } else {
            if (theseFields.front() == "QUAL") { // qual is always first according to vcfneuratrain
                useQUAL = true;
                theseFields.erase(theseFields.begin(), theseFields.begin() + 1); // erase the QUAL field, simplifies use
            }
        }
        fieldsFile.close();

        if (fields.empty()) {
            fields = theseFields;
        } else {
            if (fields != theseFields) {
                cerr << "differing field-sets provided" << endl
                     << join(fields, " ") << endl
                     << "... versus ..." << endl
                     << join(theseFields, " ") << endl;
                exit(1);
            }
        }

        struct fann *ann = fann_create_from_file(annFile->c_str());
        anns.push_back(ann);
    }

    fann_type *calc_out;

    Variant var(variantFile);

    if (!outputTag.empty()) {
        variantFile.addHeaderLine("##INFO=<ID=" + outputTag + ",Number=A,Type=Float,Description=\"Probability given model described by " + join(annFiles, " ") + "\">");
        cout << variantFile.header << endl;
    } else {
        cout << variantFile.header; // XXX BUG there shouldn't have to be two ways to write the header!!!
    }

    while (variantFile.getNextVariant(var)) {
        if (var.info.find(outputTag) != var.info.end()) {
            var.info[outputTag].clear();
        }
        double sumqual = 0;
        // for alt, get the index,
        for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
            string& alt = *a;
            int altindex = var.getAltAlleleIndex(alt);

            vector<fann_type> input; //[fields.size()];
            if (useQUAL) {
                input.push_back(var.quality); // QUAL
            }
            double val; // placeholder for conversions
            for (vector<string>::iterator f = fields.begin(); f != fields.end(); ++f) {
                if (var.info[*f].size() == var.alt.size()) {
                    convert(var.info[*f].at(altindex), val); input.push_back(val);
                } else {
                    convert(var.info[*f].front(), val); input.push_back(val);
                }
            }

            long double fresult = 0; // floating point result
            for (vector<struct fann*>::iterator ann = anns.begin(); ann != anns.end(); ++ann) {
                calc_out = fann_run(*ann, &input[0]);
                // rescale to [0,1] and add to the result sum
                fresult += 1 - (1 + calc_out[0]) / 2;
            }

            // convert to phred
            //double result = double2phred(1 - (1 + calc_out[0]) / 2);
            long double result = double2phred(fresult / anns.size());

            if (!outputTag.empty()) {
                var.info[outputTag].push_back(convert(result));
            }
            if (writeQual) {
                sumqual += result;
            }
        }
        if (writeQual) {
            var.quality = sumqual / var.alt.size();
        }

        cout << var << endl;

    }

    for (vector<struct fann*>::iterator ann = anns.begin(); ann != anns.end(); ++ann) {
        fann_destroy(*ann);
    }

    return 0;
}
