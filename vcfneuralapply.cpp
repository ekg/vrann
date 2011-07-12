#include "vcflib/Variant.h"
#include "fann.h"
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
         << "                            the VCF INFO fields which are used, is saved as FILE.fields" << endl
         << "    -o, --output-tag TAG    output the results of the neural network execution as TAG=..." << endl
         << "                            writes results to QUAL by default." << endl
         //<< "    -r, --region          specify a region on which to target the analysis, requires a BGZF" << endl
         //<< "                          compressed file which has been indexed with tabix.  any number of" << endl
         //<< "                          regions may be specified." << endl
         << endl;
}

#define PHRED_MAX 50000.0

long double float2phred(long double prob) {
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
    //float* din = malloc(num_input * num_data * sizeof(float));
    //float* dout = malloc(num_output * num_data * sizeof(float));

    string variantFileName = "-";

    bool useQUAL = false;
    bool writeQual = true;

    vector<string> fields;

    string annFile;

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
                annFile = optarg;
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

    if (annFile.empty()) {
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

    ifstream fieldsFile;
    fieldsFile.open(string(annFile + ".fields").c_str());
    string line;
    getline(fieldsFile, line);
    fields = split(line, "\t");
    if (fields.empty()) {
        cerr << "no fields provided, cannot execute neural net" << endl;
        exit(1);
    } else {
        if (fields.front() == "QUAL") { // qual is always first
            useQUAL = true;
            fields.erase(fields.begin(), fields.begin() + 1); // erase the QUAL field, simplifies use
        }
    }
    fieldsFile.close();

    struct fann *ann = fann_create_from_file(annFile.c_str());

    fann_type *calc_out;

    Variant var(variantFile);
    cout << variantFile.header;

    while (variantFile.getNextVariant(var)) {
        //
        // get the status of the validation (out of band)
        // --- pass if omni mono has variant
        // --- fail otherwise
        // get the parameters we need
        // QUAL
        vector<fann_type> input; //[fields.size()];
        if (useQUAL) {
            input.push_back(var.quality); // QUAL
        }
        float val; // placeholder for conversions
        for (vector<string>::iterator f = fields.begin(); f != fields.end(); ++f) {
            convert(var.info[*f].front(), val); input.push_back(val); // XXX broken for multi-allelic
        }
        calc_out = fann_run(ann, &input[0]);

        // rescale to [0,1] and convert to phred
        float result = float2phred(1 - (1 + calc_out[0]) / 2);

        if (writeQual) {
            var.quality = result;
        } else {
            var.info[outputTag].clear();
            var.info[outputTag].push_back(convert(result));
        }
        cout << var << endl;

    }

    fann_destroy(ann);

    return 0;
}
