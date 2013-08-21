#include <sstream>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <string>
using std::string;
#include <ctime>
#include <cstring>
#include <getopt.h>

#include "CGkArray.h"

void print_usage(char const *name)
{
    cerr << "usage: " << name << " [options] <index>" << endl
         << "Sample program to test out CGkArrays. Check README for more information." << endl;
}

// FIXME Clean up. Use for debugging only.
bool equalVectors(CGkArray::position_vector const & vector1, CGkArray::position_vector const &vector2)
{
    if ( vector1.size() != vector2.size() )
        return false;

    return std::equal ( vector1.begin(), vector1.end(), vector2.begin() );
}

// FIXME Clean up. Use for debugging only.
void revstr(uchar *t, ulong n)
{
    uchar c;
    for (ulong i = 0; i < n / 2; ++i) {
        c = t[i];
        t[i] = t[n - i - 1];
        t[n - i - 1] = c;
    }
}

// FIXME Clean up. Use for debugging only.
void complstr(uchar *t, ulong n)
{
    for (ulong i = 0; i < n; ++i) {
        switch (t[i])
        {
        case 'A': t[i] = 'T'; break;
        case 'T': t[i] = 'A'; break;
        case 'C': t[i] = 'G'; break;
        case 'G': t[i] = 'C'; break;
        default:  t[i] = 'N'; break;
        }
    }
}

int main(int argc, char **argv) 
{
    /**
     * Parse command line parameters
     */
    if (argc <= 1)
    {
        print_usage(argv[0]);
        return 0;
    }

    bool verbose = false; 
    bool debug = false;
    unsigned nqueries = 0;
#ifdef PARALLEL_SUPPORT
    unsigned parallel = 1; /* Disabled for this simple example */
#endif
    static struct option long_options[] =
        {
            {"nqueries",  required_argument, 0, 'q'},
            {"debug",     no_argument,       0, 'D'},
            {"verbose",   no_argument,       0, 'v'},
            {0, 0, 0, 0}
        };
    int option_index = 0;
    int c;
    while ((c = getopt_long(argc, argv, "q:Dv", long_options, &option_index)) != -1) 
    {
        switch(c) 
        {
        case 'q':
            nqueries = atoi(optarg); break;
        case 'v':
            verbose = true; break;
        case 'D':
            debug = true; break;
        case '?':
        case 'h':
            print_usage(argv[0]);
            return 1;
        default:
            print_usage(argv[0]);
            std::abort ();
        }
    }

    // Parse filenames
    if (argc - optind != 1)
    {
        cerr << argv[0] << ": index filename is required" << endl;
        print_usage(argv[0]);
        return 1;
    }
    string indexfile = string(argv[optind++]);

    if (nqueries < 1) 
    {
        cerr << argv[0] << ": parameter -q,--nqueries <int> is mandatory, where <int> is greater than 0." << endl;
        print_usage(argv[0]);
        return 1;
    }   
    
    /**
     * Initialize shared data structures
     */
    if (verbose) cerr << "Loading index " << indexfile << endl;
    CGkArray *tc = new CGkArray(indexfile);

    // Sanity checks
    if (!tc) {
        cerr << "readaligner: could not read index file " << indexfile << endl;
        return 1;
    }

    /**
     * Initialize shared input reader(s)
     */
    /* TODO */

    /**
     * Initialize shared output writer
     */
    /* TODO */

    // Shared counters
    unsigned total_found = 0;
    ulong total_occs = 0;
    time_t wctime = time(NULL);

#ifdef PARALLEL_SUPPORT
    if (parallel != 0)
    {
        if (verbose && parallel == 1) 
            cerr << "Using only one core (default setting for this sample script)" << endl;
        if (verbose && parallel != 1)
            cerr << "Using " << parallel << " cores." << endl;
        omp_set_num_threads(parallel);
    }
    else
        if (verbose) cerr << "Using all " << omp_get_max_threads() << " cores available." << endl;
/* #pragma omp parallel */
#endif

    /**
     * Example on how to compute a so called read-coverage profile
     *
     * Input: Read number (here the first read, n:o 0).
     * Output: k-mer coverage at each position of the read (in backwards order)
     */
    if (debug)
    {
        unsigned readno = 0;
        ulong readpos = tc->getLength(readno) - tc->getGkSize() - 1;
        ulong tmp = tc->initMoveLeft(readno); // Initialize <internal value> for the given read n:o
        CGkArray::sa_range sar = tc->moveLeft(tmp); // Get the SA range of the last k-mer
        while (sar.first <= sar.second) // End of the traversal is signaled by an empty SA range
        { 
            cout << "At read " << readno << " position " << readpos << " coverage is " << tc->countReads(sar) << endl;
            
            // Take one step left (using <internal value>)            
            sar = tc->moveLeft(tmp);
            --readpos;
        }
        assert(readpos == 0); // Arrived at the first position of the read
    }

    /**
     * Read-coverage profile for the reverse complemented strand
     *
     * Input: Read number (here the first read, n:o 0).
     * Output: k-mer coverage for the reverse complemented read
     */
    if (debug)
    {
        // Read number is used here only to retrieve the read sequence
        unsigned readno = 0;
        // Retrieve the read sequence from the index, user must delete [] the buffer
        uchar *read = tc->getRead(readno);
        unsigned l = tc->getLength(readno) - 1;
        ulong readpos = l - tc->getGkSize() + 1;

        // Rev. compl. the given read
        revstr(read, l);
        complstr(read, l);
        
        // Initialize <internal pointer> for an arbitrary pattern
        CGkArray::internal_pointer tmp = tc->initMoveLeft(read); 
        CGkArray::sa_range sar = tc->moveLeft(tmp, read); // Get the SA range of the last k-mer
        while (readpos)
        { 
            --readpos;
            cout << "At reverse complemented read " << readno << " position " << readpos 
                 << " coverage is " << tc->countReads(sar) << endl;
            // Take one step left (using <internal pointer>)            
            sar = tc->moveLeft(tmp, read);
        }
        assert(readpos == 0); // Arrived at the first position of the read
        delete [] read;
    }


    /**
     * Test random positions
     */
    srand(543262346);
    wctime = time(NULL);
    cerr << "Testing " << nqueries << " random positions for Q1..." << endl;
    for (unsigned i = 0; i < nqueries; ++i)
    {
        ulong pos = rand() % tc->getLength();
        if (!tc->isValidTextPos(pos))
            pos -= tc->getGkSize(); // Make it a valid position (at least k nucleotides from the end of the read)

        CGkArray::position_vector occs = tc->reportReads(pos);

        if (debug)
        {
            uchar const *suffix = tc->getSuffix(tc->inverseSA(pos), tc->getGkSize());
            for (CGkArray::position_vector::iterator it = occs.begin(); it != occs.end(); ++it)
                cout << "Pattern " << suffix << " = " << it->first << "," << it->second << endl;
            
            // Search with the uchar *
            CGkArray::sa_range sar = tc->kmerToSARange(suffix);
            CGkArray::position_vector occs2 = tc->reportReads(sar); // query with SA range
            if (!equalVectors(occs,occs2))
            { cerr << "Q1 assert failed: vectors were not equal at i = " << i << endl; abort(); }
            delete [] suffix;
        }
        if (!occs.empty())
            ++total_found;
        total_occs += occs.size();
    }
    cerr << "Number of reported alignments: " << total_occs << endl
         << "Number of reads found: " << total_found << endl;
    cerr << "Wall-clock time: " << std::difftime(time(NULL), wctime) << " seconds (" 
         << std::difftime(time(NULL), wctime) / 3600 << " hours)" << endl;


    srand(543262346);    
    cerr << "Testing " << nqueries << " random positions for Q2..." << endl;
    total_found = 0;
    wctime = time(NULL);
    total_occs = 0;
    for (unsigned i = 0; i < nqueries; ++i)
    {
        ulong pos = rand() % tc->getLength();
        if (!tc->isValidTextPos(pos))
            pos -= tc->getGkSize(); // Make it a valid position (at least k nucleotides from the end of the read)

        ulong occs = tc->countReads(pos);
        if (debug)
        {
            uchar const *suffix = tc->getSuffix(tc->inverseSA(pos), tc->getGkSize());
            cout << "Pattern " << suffix << " = " << occs << endl;

            // Search with the uchar *
            CGkArray::sa_range sar = tc->kmerToSARange(suffix);
            unsigned occs2 = tc->countReads(sar); // query with SA range
            if (occs != occs2)
            { cerr << "Q2 assert failed: counts were not equal at i = " << i << endl; abort(); }
            delete [] suffix;
        }
        if (occs)
            ++total_found;
        total_occs += occs;
    }
    cerr << "Number of reported alignments: " << total_occs << endl
         << "Number of reads found: " << total_found << endl;
    cerr << "Wall-clock time: " << std::difftime(time(NULL), wctime) << " seconds (" 
         << std::difftime(time(NULL), wctime) / 3600 << " hours)" << endl;


    srand(543262346);    
    cerr << "Testing " << nqueries << " random positions for Q3..." << endl;
    total_found = 0;
    wctime = time(NULL);
    total_occs = 0;
    for (unsigned i = 0; i < nqueries; ++i)
    {
        ulong pos = rand() % tc->getLength();
        if (!tc->isValidTextPos(pos))
            pos -= tc->getGkSize(); // Make it a valid position (at least k nucleotides from the end of the read)

        CGkArray::position_vector occs = tc->reportOccs(pos);
        if (debug)
        {
            uchar const *suffix = tc->getSuffix(tc->inverseSA(pos), tc->getGkSize());
            for (CGkArray::position_vector::iterator it = occs.begin(); it != occs.end(); ++it)
                cout << "Pattern " << suffix << " = " << it->first << "," << it->second << endl;

            // Search with the uchar *
            CGkArray::sa_range sar = tc->kmerToSARange(suffix);
            CGkArray::position_vector occs2 = tc->reportOccs(sar); // query with SA range
            if (!equalVectors(occs,occs2))
            { cerr << "Q3 assert failed: vectors were not equal at i = " << i << endl; abort(); }
            delete [] suffix;
        }
        if (!occs.empty())
            ++total_found;
        total_occs += occs.size();
    }
    cerr << "Number of reported alignments: " << total_occs << endl
         << "Number of reads found: " << total_found << endl;
    cerr << "Wall-clock time: " << std::difftime(time(NULL), wctime) << " seconds (" 
         << std::difftime(time(NULL), wctime) / 3600 << " hours)" << endl;


    srand(543262346);    
    cerr << "Testing " << nqueries << " random positions for Q4..." << endl;
    total_found = 0;
    wctime = time(NULL);
    total_occs = 0;
    for (unsigned i = 0; i < nqueries; ++i)
    {
        ulong pos = rand() % tc->getLength();
        if (!tc->isValidTextPos(pos))
            pos -= tc->getGkSize(); // Make it a valid position (at least k nucleotides from the end of the read)

        ulong occs = tc->countOccs(pos);
        if (debug)
        {
            uchar const * suffix = tc->getSuffix(tc->inverseSA(pos), tc->getGkSize());
            cout << "Testing " << suffix << " = " << occs << endl;

            // Search with the uchar *
            CGkArray::sa_range sar = tc->kmerToSARange(suffix);
            unsigned occs2 = tc->countOccs(sar); // query with SA range
            if (occs != occs2)
            { cerr << "Q4 assert failed: counts were not equal at i = " << i << endl; abort(); }
            delete [] suffix;
        }
        if (occs)
            ++total_found;
        total_occs += occs;
    }
    cerr << "Number of reported alignments: " << total_occs << endl
         << "Number of reads found: " << total_found << endl;
    cerr << "Wall-clock time: " << std::difftime(time(NULL), wctime) << " seconds (" 
         << std::difftime(time(NULL), wctime) / 3600 << " hours)" << endl;

    delete tc;
}
