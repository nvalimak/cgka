/** 
 * TODO Proper builder with FASTQ and FASTA support.
 * TODO Code clean up. Mixture of C and C++ follows...
 */
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <cstring>
#include <cassert>
#include <getopt.h>
#include "bcr-demo.h"
#include "CGkArray.h"

// Include from library RLCSA
#include "bits/deltavector.h"

using namespace std;

#define DEFAULT_SAMPLERATE 16
#define DEFAULT_BLOCKSIZE 16

/**
 * Flags set based on command line parameters
 */
bool verbose = false;
unsigned gk = 0; // K for Gk arrays

void revstr(char *t, ulong n)
{
    char c;
    for (ulong i = 0; i < n / 2; ++i) {
        c = t[i];
        t[i] = t[n - i - 1];
        t[n - i - 1] = c;
    }
}

void write(long length, CSA::DeltaEncoder & de, string const &filename)
{
    ofstream ofs(filename.c_str());
    if (!ofs.good())
    { cerr << "error: unable to write " << filename << endl; abort(); }
    CSA::DeltaVector dv(de, length+1);
    dv.writeTo(ofs);

/* Debug    CSA::DeltaVector::Iterator iter(dv);
    for (long i = 0; i < length+1; ++i)
        cerr << (iter.isSet(i) ? '1' : '0');
    cerr << endl;*/
    ofs.close();
}

void print_usage(char const *name)
{
    cerr << "usage: " << name << " [options] <input> [output]" << endl
         << "Check README or `" << name << " --help' for more information." << endl;
}

void print_help(char const *name)
{
    cerr << "usage: " << name << " [options] <input> [output]" << endl << endl
         << "<input> is the input filename. "
         << "The input must be in plain-text format (i.e. sequences separated by '\\n')." << endl
         << "If no output filename is given, the index is stored as <input>.cgka" << endl << endl
         << "Options:" << endl
         << " -k <int>, --gk <int>          k-mer length (mandatory option)." << endl
         << " -s <int>, --sample-rate <int> Sampling rate for the index, a smaller number " << endl
         << "                               yields a bigger index but can decrease search " << endl
         << "                               time (default: " << DEFAULT_SAMPLERATE << ")." << endl
         << " -h, --help                    Display command line options." << endl
         << " -v, --verbose                 Print progress information." << endl;
}

int atoi_min(char const *value, int min, char const *parameter, char const *name)
{
    int i = atoi(value);
    if (i < min)
    {
        cerr << name << ": argument of " << parameter << " must be greater than or equal to " << min << endl
             << "Check README or `" << name << " --help' for more information." << endl;
        std::exit(1);
    }
    return i;
}

int main(int argc, char **argv) 
{
    /**
     * Parse command line parameters
     */
    if (argc == 1)
    {
        print_usage(argv[0]);
        return 1;        
    }
    unsigned samplerate = 0;
    static struct option long_options[] =
        {
            {"gk",          required_argument, 0, 'k'},
            {"sample-rate", required_argument, 0, 's'},
            {"help",        no_argument,       0, 'h'},
            {"verbose",     no_argument,       0, 'v'},
            {0, 0, 0, 0}
        };
    int option_index = 0;
    int c;
    while ((c = getopt_long(argc, argv, "cR:s:hvk:",
                            long_options, &option_index)) != -1) 
    {
        switch(c) 
        {
        case 'k':
            gk = atoi_min(optarg, 3, "-k", argv[0]); 
            break;
        case 's':
            samplerate = atoi_min(optarg, 1, "-s, --sample-rate", argv[0]); 
            break;
        case 'h':
            print_help(argv[0]);
            return 0;
        case 'v':
            verbose = true; break;
        case '?':
            print_usage(argv[0]);
            return 1;
        default:
            print_usage(argv[0]);
            std::abort();
        }
    }
    
    if (samplerate <= 3)
        cerr << "Warning: small samplerates (-s, --sample-rate) may yield infeasible index sizes" << endl;

    if (argc - optind < 1)
    {
        cerr << argv[0] << ": no input filename given!" << endl;
        print_usage(argv[0]);
        return 1;
    }
        
    if (argc - optind > 2)
        cerr << "Warning: too many filenames given! Ignoring all but first two." << endl;

    if (gk == 0)
    {
        cerr << "error: you must provide the k-mer length (parameter -k <int>)" << endl;
        return 1;
    }

    string inputfile = string(argv[optind++]);
    string outputfile = "";
    if (optind != argc)
        outputfile = string(argv[optind++]);
    else
        outputfile = inputfile; // suffix will be added by TextCollection::save().
    
    istream *fp;
    if (inputfile == "-")
        fp = &std::cin;
    else
        fp = new std::ifstream(inputfile.c_str());

    if (!fp->good())
    {
        cerr << argv[0] << ": unable to read input file " << inputfile << endl;
        exit(1); 
    }

    cerr << std::fixed;
    cerr.precision(2);
    time_t wctime = time(NULL);

    // estimate the total input sequence length
    fp->seekg(0, ios::end);
    long length = fp->tellg();
    if (length == -1)
    {
        cerr << "error: unable to estimate input file size" << endl;
        return 1;
    }
    fp->seekg(0);
    
    /**
     * Build forward/rotation index
     */
    if (verbose)
        cerr << "Building the forward index:" << endl;

    unsigned numberOfTexts = 0, maxTextLength = 0, curLength = 0;
    uchar *s = new uchar[length];
    CSA::DeltaEncoder * de = new CSA::DeltaEncoder(DEFAULT_BLOCKSIZE); // Collects text start positions
    de->setBit(0);
    fp->read((char *)s, length);
    delete fp;
    for (long i = 0; i < length; ++i)
        if (s[i] == '\n') 
        {            
            s[i] = 0;
            if (maxTextLength < curLength)
                maxTextLength = curLength;            
            curLength = 0;
            de->setBit(i+1);
            ++numberOfTexts;
        }
        else
            ++curLength;
    
    write(length, *de, outputfile + ".cgka_map");
    delete de; de = 0;

    uchar *B = bcr_lite(0, 0, length, s);
    delete [] s;
/*    for (long i = 0; i < length; ++i)
        putchar(B[i]? B[i] : '$');
        putchar('\n');*/

    CGkArray *cgka = new CGkArray(B, length, samplerate, numberOfTexts, maxTextLength, gk, verbose);
    // B was already free()'d;

    cgka->save(outputfile);

    delete cgka;
    if (verbose) 
        std::cerr << "Save complete. "
                  << "(total wall-clock time " << std::difftime(time(NULL), wctime) << " s, " 
                  << std::difftime(time(NULL), wctime) / 3600 << " hours)" << endl;
    return 0;
}
 
