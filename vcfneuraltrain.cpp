#include "vcflib/Variant.h"
#include "doublefann.h"
#include "convert.h"
#include <fstream>
#include <iostream>
#include <getopt.h>
#include "split/join.h"

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

struct fann_train_data *read_from_array(vector<double>& din, vector<double>& dout, unsigned int num_data, unsigned int num_input, unsigned int num_output) {
    unsigned int i, j;
    fann_type *data_input, *data_output;
    struct fann_train_data *data =
        (struct fann_train_data *) malloc(sizeof(struct fann_train_data));
    if(data == NULL) {
        fann_error(NULL, FANN_E_CANT_ALLOCATE_MEM);
        return NULL;
    }

    fann_init_error_data((struct fann_error *) data);

    data->num_data = num_data;
    data->num_input = num_input;
    data->num_output = num_output;
    data->input = (double **) calloc(num_data, sizeof(double *));
    if(data->input == NULL) {
        fann_error(NULL, FANN_E_CANT_ALLOCATE_MEM);
        fann_destroy_train(data);
        return NULL;
    }

    data->output = (double **) calloc(num_data, sizeof(double *));
    if(data->output == NULL) {
        fann_error(NULL, FANN_E_CANT_ALLOCATE_MEM);
        fann_destroy_train(data);
        return NULL;
    }

    data_input = (double *) calloc(num_input * num_data, sizeof(double));
    if(data_input == NULL) {
        fann_error(NULL, FANN_E_CANT_ALLOCATE_MEM);
        fann_destroy_train(data);
        return NULL;
    }

    data_output = (double *) calloc(num_output * num_data, sizeof(double));
    if(data_output == NULL) {
        fann_error(NULL, FANN_E_CANT_ALLOCATE_MEM);
        fann_destroy_train(data);
        return NULL;
    }

    for(i = 0; i != num_data; i++) {
        data->input[i] = data_input;
        data_input += num_input;
   
        for(j = 0; j != num_input; j++) {
            data->input[i][j] = din[i*num_input+j];
        }
   
   
        data->output[i] = data_output;
        data_output += num_output;
   
        for(j = 0; j != num_output; j++) {
            data->output[i][j] = dout[i*num_output+j];
        }
    }
    return data;
}

void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options]" << endl
         << endl
         << "Trains an artificial neural network on a given set of fields, saves the network" << endl
         << "for later use in quality score recalibration or variant subsetting." << endl
         << endl
         << "options:" << endl
         << "    -h, --help              this dialog." << endl
         << "    -v, --vcf-file FILE     specifies the vcf file (or BGZipped, vcf.gz) to use for training" << endl
         << "                            if '-' specified, stdin is used (default)" << endl
         << "    -u, --unsupervised      use Self-Organizing Maps to learn from input data" << endl
         << "    -f, --field FIELD       use this field as a training target, any number may be specified" << endl
         << "    -P, --pass-tag TAG      this VCF tag indicates that the record passes the filter which" << endl
         << "                            the neural net will approximate." << endl
         << "    -V, --pass-tag-value V  the tag value which specifies if the given allele passes or fails" << endl
        //<< "    -F, --fail-tag TAG      this VCF tag indicates that the record fails the filter which" << endl
        //<< "                            the neural net will approximate (default: missing PASS tag indicates" << endl
        //<< "                            failure)." << endl
         << "    -a, --ann-file FILE     save the ANN to this file (required).  metadata, specifically" << endl
         << "                            the VCF INFO fields which are used, is saved as FILE.fields" << endl
         << "    -n, --normalize-input   normalize the input data to [-1,1]" << endl
         << "    -l, --layers N          use this many layers in the ANN, default 3" << endl
         << "    -H, --hidden-neurons N  use this many hidden neurons, default 100" << endl
         << "    -e, --target-error N    train until this error is reached, default 0.001" << endl
         << "    -m, --max-epochs N      train no more than this many epochs, default 500000" << endl
         << "    -r, --report-interval N report progress even N epochs, default 10" << endl
         << "    -C, --cascade           use cascade training algorithm to build neural network during training" << endl
         //<< "    -r, --region          specify a region on which to target the analysis, requires a BGZF" << endl
         //<< "                          compressed file which has been indexed with tabix.  any number of" << endl
         //<< "                          regions may be specified." << endl
         << endl;
}


int main(int argc, char** argv)
{   

    unsigned int num_output = 1;
    //double* din = malloc(num_input * num_data * sizeof(double));
    //double* dout = malloc(num_output * num_data * sizeof(double));
    vector<double> din, dout;
    map<string, vector<double> > dins; // to normalize into din

    string variantFileName = "-";

    string passTag;
    string passValue;

    vector<string> fields;

    unsigned int num_layers = 3;
    unsigned int num_neurons_hidden = 100;
    double desired_error = (double) 0.001;
    unsigned int max_epochs = 500000;
    unsigned int epochs_between_reports = 10;

    string annFile;

    bool normalizeInput = false;
    map<string, Stats> stats;
    bool cascadeTrain = false;
    bool unsupervised = false;

    int width = 100;
    int height = 100;
    int num_dimensions = 2;

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
                {"field",  required_argument, 0, 'f'},
                {"pass-tag",  required_argument, 0, 'P'},
                //{"fail-tag",  required_argument, 0, 'F'},
                {"pass-tag-value", required_argument, 0, 'V'},
                {"ann-file", required_argument, 0, 'a'},
                {"normalize-input", no_argument, 0, 'n'},
                {"layers", required_argument, 0, 'l'},
                {"hidden-neurons", required_argument, 0, 'H'},
                {"target-error", required_argument, 0, 'e'},
                {"max-epochs", required_argument, 0, 'm'},
                {"report-interval", required_argument, 0, 'r'},
                {"cascade", no_argument, 0, 'C'},
                {"unsupervised", no_argument, 0, 'u'},
                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hnCv:f:P:a:l:H:e:m:r:V:u:",
                         long_options, &option_index);

        if (c == -1)
            break;

        string field;

        switch (c)
        {

        case 'v':
            variantFileName = optarg;
            break;

        case 'f':
            field = string(optarg);
            fields.push_back(field);
            break;

        case 'C':
            cascadeTrain = true;
            break;

        case 'u':
            unsupervised = true;
            break;

        case 'P':
            passTag = optarg;
            break;

        case 'V':
            passValue = optarg;
            break;

        case 'a':
            annFile = optarg;
            break;

        case 'n':
            normalizeInput = true;
            break;

        case 'l':
            if (!convert(optarg, num_layers)) {
                cerr << "could not parse --layers" << endl;
                exit(1);
            }
            break;

        case 'H':
            if (!convert(optarg, num_neurons_hidden)) {
                cerr << "could not parse --hidden-neurons" << endl;
                exit(1);
            }
            break;

        case 'e':
            if (!convert(optarg, desired_error)) {
                cerr << "could not parse --target-error" << endl;
                exit(1);
            }
            break;

        case 'm':
            if (!convert(optarg, max_epochs)) {
                cerr << "could not parse --max-epochs" << endl;
                exit(1);
            }
            break;

        case 'r':
            if (!convert(optarg, epochs_between_reports)) {
                cerr << "could not parser --report-interval" << endl;
                exit(1);
            }
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
        cerr << "please supply an output filename for the neural network, --ann-file" << endl;
        exit(1);
    }

    if (!unsupervised && passTag.empty()) {
        cerr << "please specify a tag that indicates that the record passes, --pass-tag" << endl;
        exit(1);
    }

    // TODO allow operation which assumes that all not-tagged records are passing or failing?
    //if (failTag.empty()) {
    //    cerr << "please specify a tag that indicates that the record fails, --fail-tag" << endl;
    //    exit(1);
    //}

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

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        if (var.info.find(passTag) == var.info.end()) {
            continue;
        }
        bool missingData = false;
        for (vector<string>::iterator f = fields.begin(); f != fields.end(); ++f) {
            if (var.info.find(*f) == var.info.end()) {
                missingData = true;
            }
        }
        if (missingData) {
            continue;
        }

        //cout << "processing " << var.sequenceName << ":" << var.position << endl;
        //
        // get the status of the validation (out of band)
        // --- pass if tagged
        // --- fail otherwise
        for (int a = 0; a < var.alt.size(); ++a) {
            if (var.info.find(passTag) != var.info.end()) {
                if (var.info[passTag].size() == var.alt.size()) {
                    if (var.info[passTag].at(a) == passValue) {
                        dout.push_back(1); // passing
                    }
                } else if (var.info[passTag].front() == passValue) {
                        dout.push_back(1); // passing
                } else {
                    dout.push_back(0); // failing
                }
            } else {
                continue; // doesn't pass or fail, so we can't use it for training
            }
            double val; // placeholder for conversions
            for (vector<string>::iterator f = fields.begin(); f != fields.end(); ++f) {
                if (*f == "QUAL") {
                    dins["QUAL"].push_back(var.quality); // QUAL
                } else {
                    if (variantFile.infoCounts[*f] == 1) { // for non Allele-variant fields
                        convert(var.info[*f].at(0), val);
                    } else {
                        convert(var.info[*f].at(a), val);
                    }
                    dins[*f].push_back(val);
                }
            }
        }
    }

    unsigned int num_input = dins.size();

    // normalize inputs
    if (normalizeInput) {
        // get normalization vector
        // goal is normalization at 0, sd=1
        int i = 0;
        for (vector<string>::iterator f = fields.begin(); f != fields.end(); ++f, ++i) {
            vector<double>& fv = dins[*f];
            Stats& s = stats[*f];
            // get normalization constants
            s.mean = mean(fv);
            s.stdev = standard_deviation(fv, s.mean);
            // normalize
            for (vector<double>::iterator d = fv.begin(); d != fv.end(); ++d) {
                *d = (*d - s.mean) / s.stdev;
            }
        }
    }

    // now drop the input into the format which fann wants

    vector<vector<double>*> dinsv;
    for (map<string, vector<double> >::iterator i = dins.begin(); i != dins.end(); ++i) {
        dinsv.push_back(&i->second);
    }

    for (int i = 0; i < dinsv.front()->size(); ++i) {
        for (vector<vector<double>*>::iterator v = dinsv.begin(); v != dinsv.end(); ++v) {
            din.push_back((*v)->at(i));
        }
    }

    struct fann *ann;
    if (!unsupervised) {
        ann = fann_create_standard(num_layers, num_input, num_neurons_hidden, num_output);
        fann_set_activation_function_hidden(ann, FANN_SIGMOID_SYMMETRIC);
        fann_set_activation_function_output(ann, FANN_SIGMOID_SYMMETRIC);
    }
 /* else {
    ann = fann_create_som(width, height, num_dimensions);
    ann->som_params->som_topology = FANN_SOM_TOPOLOGY_HEXAGONAL;
    ann->som_params->som_neighborhood = FANN_SOM_NEIGHBORHOOD_DISTANCE;
    ann->som_params->som_learning_decay = FANN_SOM_LEARNING_DECAY_LINEAR;
    ann->som_params->som_learning_rate = 0.01;
    }*/

    unsigned int num_data = din.size() / num_input;
    fann_train_data *data = read_from_array(din, dout, num_data, num_input, num_output);

    if (cascadeTrain) {
        fann_cascadetrain_on_data(ann, data, max_epochs, epochs_between_reports, desired_error);
    } else {
        fann_train_on_data(ann, data, max_epochs, epochs_between_reports, desired_error);
    }

    cerr << "saving neural net to " << annFile << endl;
    fann_save(ann, annFile.c_str());

    fann_destroy(ann);
    fann_destroy_train(data);

    // save the fields we trained against, for use when applying neural net to other data
    string annDescriptionFileName = annFile + ".meta";
    cerr << "writing neural net structure information to " << annDescriptionFileName << endl;
    save_ann_metadata(annDescriptionFileName, fields, stats);

    return 0;
}
