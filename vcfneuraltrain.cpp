#include "vcflib/Variant.h"
#include "fann.h"
#include "convert.h"
#include <fstream>
#include <iostream>
#include <getopt.h>
#include "split/join.h"

using namespace std;
using namespace vcf;


struct fann_train_data *read_from_array(vector<float>& din, vector<float>& dout, unsigned int num_data, unsigned int num_input, unsigned int num_output) {
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
  data->input = (float **) calloc(num_data, sizeof(float *));
  if(data->input == NULL) {
    fann_error(NULL, FANN_E_CANT_ALLOCATE_MEM);
    fann_destroy_train(data);
    return NULL;
  }

  data->output = (float **) calloc(num_data, sizeof(float *));
  if(data->output == NULL) {
    fann_error(NULL, FANN_E_CANT_ALLOCATE_MEM);
    fann_destroy_train(data);
    return NULL;
  }

  data_input = (float *) calloc(num_input * num_data, sizeof(float));
  if(data_input == NULL) {
    fann_error(NULL, FANN_E_CANT_ALLOCATE_MEM);
    fann_destroy_train(data);
    return NULL;
  }

  data_output = (float *) calloc(num_output * num_data, sizeof(float));
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
         << "    -f, --field FIELD       use this field as a training target, any number may be specified" << endl
         << "    -P, --pass-tag TAG      this VCF tag indicates that the record passes the filter which" << endl
         << "                            the neural net will approximate." << endl
         << "    -F, --fail-tag TAG      this VCF tag indicates that the record fails the filter which" << endl
         << "                            the neural net will approximate." << endl
         << "    -a, --ann-file FILE     save the ANN to this file (required).  metadata, specifically" << endl
         << "                            the VCF INFO fields which are used, is saved as FILE.fields" << endl
         << "    -n, --normalize-input   normalize the input data to [-1,1]" << endl
         << "    -l, --layers N          use this many layers in the ANN, default 3" << endl
         << "    -H, --hidden-neurons N  use this many hidden neurons, default 100" << endl
         << "    -e, --target-error N    train until this error is reached, default 0.001" << endl
         << "    -m, --max-epochs N      train no more than this many epochs, default 500000" << endl
         << "    -r, --report-interval N report progress even N epochs, default 10" << endl
         //<< "    -r, --region          specify a region on which to target the analysis, requires a BGZF" << endl
         //<< "                          compressed file which has been indexed with tabix.  any number of" << endl
         //<< "                          regions may be specified." << endl
         << endl;
}


int main(int argc, char** argv)
{   

    unsigned int num_output = 1;
    //float* din = malloc(num_input * num_data * sizeof(float));
    //float* dout = malloc(num_output * num_data * sizeof(float));
    vector<float> din, dout;
    map<string, vector<float> > dins; // to normalize into din

    string variantFileName = "-";

    string passTag;
    string failTag;

    bool useQUAL = false;
    vector<string> fields;

    unsigned int num_layers = 3;
    unsigned int num_neurons_hidden = 100;
    float desired_error = (float) 0.001;
    unsigned int max_epochs = 500000;
    unsigned int epochs_between_reports = 10;

    string annFile;

    bool normalizeInput = false;

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
            {"fail-tag",  required_argument, 0, 'F'},
            {"ann-file", required_argument, 0, 'a'},
            {"normalize-input", no_argument, 0, 'n'},
            {"layers", required_argument, 0, 'l'},
            {"hidden-neurons", required_argument, 0, 'H'},
            {"target-error", required_argument, 0, 'e'},
            {"max-epochs", required_argument, 0, 'm'},
            {"report-interval", required_argument, 0, 'r'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hnv:f:P:F:a:l:H:e:m:r:",
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
                if (field == "QUAL") {
                    useQUAL= true; // hack to support VCF distinction between QUAL and INFO fields
                } else {
                    fields.push_back(field);
                }
                break;

            case 'P':
                passTag = optarg;
                break;

            case 'F':
                failTag = optarg;
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

    if (passTag.empty()) {
        cerr << "please specify a tag that indicates that the record passes, --pass-tag" << endl;
        exit(1);
    }

    // TODO allow operation which assumes that all not-tagged records are passing or failing?
    if (failTag.empty()) {
        cerr << "please specify a tag that indicates that the record fails, --fail-tag" << endl;
        exit(1);
    }

    /*
    while (optind < argc) {
        string field = argv[optind++];
        if (field == "QUAL") {
            useQUAL= true; // hack to support VCF distinction between QUAL and INFO fields
        } else {
            fields.push_back(field);
        }
    }
    */

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
        //
        // get the status of the validation (out of band)
        // --- pass if omni mono has variant
        // --- fail otherwise
        if (var.infoFlags.find(passTag) != var.infoFlags.end()) {
            dout.push_back(1); // passing
        } else if (var.infoFlags.find(failTag) != var.infoFlags.end()) {
            dout.push_back(-1); // not passing
        } else {
            continue; // doesn't pass or fail
        }
        // get the parameters we need
        // QUAL
        if (useQUAL) {
            dins["QUAL"].push_back(var.quality); // QUAL
        }
        float val; // placeholder for conversions
        for (vector<string>::iterator f = fields.begin(); f != fields.end(); ++f) {
            convert(var.info[*f].front(), val); dins[*f].push_back(val);
        }
    }

    unsigned int num_input = dins.size();

    // normalize the inputs to -1, 1
    if (normalizeInput) {
        for (map<string, vector<float> >::iterator i = dins.begin(); i != dins.end(); ++i) {
            vector<float>& in = i->second;
            // get the min and max
            float min = in.front();
            float max = in.front();
            for (vector<float>::iterator j = in.begin(); j != in.end(); ++j) {
                float m = *j;
                if (m < min) {
                    min = m;
                }
                if (m > max) {
                    max = m;
                }
            }
            float scaling = ( 1 - -1 ) / ( max - min );
            // normalize
            for (vector<float>::iterator j = in.begin(); j != in.end(); ++j) {
                float& m = *j;
                m -= min;
                m *= scaling;
            }
        }
    }
    // now drop the input into the format which fann wants

    vector<vector<float>*> dinsv;
    for (map<string, vector<float> >::iterator i = dins.begin(); i != dins.end(); ++i) {
        dinsv.push_back(&i->second);
    }

    for (int i = 0; i < dinsv.front()->size(); ++i) {
        for (vector<vector<float>*>::iterator v = dinsv.begin(); v != dinsv.end(); ++v) {
            din.push_back((*v)->at(i));
        }
    }


    struct fann *ann = fann_create_standard(num_layers, num_input, num_neurons_hidden, num_output);

    fann_set_activation_function_hidden(ann, FANN_SIGMOID_SYMMETRIC);
    fann_set_activation_function_output(ann, FANN_SIGMOID_SYMMETRIC);

    unsigned int num_data = din.size() / num_input;
    //if (din.size() % num_input != 0) cerr << "what???" << endl;
    fann_train_data *data = read_from_array(din, dout, num_data, num_input, num_output);

    //fann_train_on_file(ann, "test.data", max_epochs, epochs_between_reports, desired_error);
    fann_train_on_data(ann, data, max_epochs, epochs_between_reports, desired_error);

    cerr << "saving neural net to " << annFile << endl;
    fann_save(ann, annFile.c_str());

    fann_destroy(ann);

    string annDescriptionFileName = annFile + ".fields";
    cerr << "writing neural net structure information to " << annDescriptionFileName << endl;
    // save the fields we trained against, for use when applying neural net to other data
    ofstream annDescriptionFile;
    annDescriptionFile.open(annDescriptionFileName.c_str());
    if (useQUAL) {
        annDescriptionFile << "QUAL";
        if (!fields.empty()) annDescriptionFile << "\t";
    }
    if (!fields.empty()) {
        annDescriptionFile << join(fields, "\t");
    }
    annDescriptionFile << endl;
    annDescriptionFile.close();

    // TODO free other data, dins, etc.

    return 0;
}
