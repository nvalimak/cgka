/******************************************************************************
 *                                                                            *
 *   This program is free software; you can redistribute it and/or modify     *
 *   it under the terms of the GNU Lesser General Public License as published *
 *   by the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                      *
 *                                                                            *
 *   This program is distributed in the hope that it will be useful,          *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *   GNU Lesser General Public License for more details.                      *
 *                                                                            *
 *   You should have received a copy of the GNU Lesser General Public License *
 *   along with this program; if not, write to the                            *
 *   Free Software Foundation, Inc.,                                          *
 *   51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.            *
 *****************************************************************************/
#include "CGkArray.h"

#include <iostream>
#include <map>
#include <unordered_map>
#include <utility>
#include <deque>
#include <stdexcept>
#include <cassert>
#include <cstring> // For strlen()
using std::vector;
using std::pair;
using std::make_pair;
using std::map;

// Alphabet is fixed for getSuffix()
const char CGkArray::ALPHABET_SHIFTED[] = {1, 'A'+1, 'C'+1, 'G'+1, 'N'+1, 'T'+1};

// Save file version info
const uchar CGkArray::versionFlag = 17;

/** sets bit p in e */
#define bitset32(e,p) ((e)[(p)/32] |= (1<<((p)%32)))
/** reads bit p from e */
#define bitget32(e,p) ((((e)[(p)/32] >> ((p)%32))) & 1)


/**
 * Returns a copy of an indexed read.
 *
 * Input: Read number (0-based)
 * Output: Copy of the read. Caller must delete [] the buffer.
 */ 
uchar * CGkArray::getRead(unsigned readno) const
{
    // Position of the trailing '\0' byte in SA is at...
    ulong i = readno;

    unsigned l = getLength(readno);
    uchar *text = new uchar[l]; // Length includes '\0' byte
    text[--l] = 0;
    ulong alphabetrank_i_tmp = 0;
    uchar c  = alphabetrank->access(i, alphabetrank_i_tmp);
    text[--l] = c;
    while (l > 0)
    {
        i = C[c]+alphabetrank_i_tmp-1;
        c = alphabetrank->access(i, alphabetrank_i_tmp);
        text[--l] = c;  
    }
    return text;
}


/**
 * Given suffix i and substring length l, return T[SA[i] ... SA[i]+l].
 *
 * Caller must delete [] the buffer.
 */ 
uchar * CGkArray::getSuffix(ulong dest, unsigned l) const
{
    ++dest;
    uchar *text = new uchar[l+1];
    text[l] = 0;
    for (unsigned i = 0; i < l; ++i)
    {
        unsigned int which = 0;
        int c = 0;
        // ALPHABET_SHIFTED contains chars shifted +1
        for (const char *ch = ALPHABET_SHIFTED; ch < ALPHABET_SHIFTED+sizeof(ALPHABET_SHIFTED); ++ch) {
            c = (char)*ch;
            if (C[c] >= dest) {
                which = dest - C[c-1];
                c--;
                if (c == 0) 
                {
                    text[i] = 0;
                    return text;
                }
                break;
            }
        }
        text[i] = (uchar)c;
        dest = alphabetrank->select(c, which) + 1;
    }
    return text;
}

/**
 * Find the suffix array range for f
 *
 * Input: k-mer
 * Output: Suffix array range
 */
CGkArray::sa_range CGkArray::kmerToSARange(uchar const *kmer) const
{
    ulong smin = 0;
    ulong smax = n-1;
    unsigned i = gk;
    while (i > 0)
    {
        smin = LF(kmer[i-1], smin-1);
        smax = LF(kmer[i-1], smax)-1;
        if (smin > smax)
            return make_pair(1,0); // Not found
        --i;
    }
    return make_pair(smin, smax);
}

/**
 * Move left operation gives an efficient way to step over each k-mer in a read.
 *
 * Each call to this method returns an SA range of the preceding read position,
 * that is, the k-mers in the read are traversed in backwards order.
 * This method allows to compute, for example, so called read-coverage profiles.
 *
 * Input: <internal value>, use initMoveLeft() to initialize.
 * Output: Suffix array range of the preceeding read position
 */
CGkArray::sa_range CGkArray::moveLeft(ulong &i) const
{
    ulong alphabetrank_i_tmp = 0;
    uchar c  = alphabetrank->access(i, alphabetrank_i_tmp);
    if (c == '\0')
        return make_pair(1,0); // Cannot move left.
    
    i = C[c]+alphabetrank_i_tmp-1;
    return make_pair(Blcp->prev(i), Blcp->next(i+1)-1);
}

/**
 * Move left operation for arbitrary patterns
 *
 * Allows an efficient way to step over each k-mer in any given string.
 *
 * Input: <internal pointer>, use initMoveLeft() to initialize, and
 *        String that was used to initialize the internal pointer.
 * Output: Suffix array range of the preceeding position. 
 *         Range can be empty if the k-mer is not found in the index.
 *         After the last step, subsequent calls return an empty SA range.
 */
CGkArray::sa_range CGkArray::moveLeft(CGkArray::internal_pointer &intp, uchar const *pattern) const
{
    sa_range &sar = intp.first;
    unsigned &pos = intp.second;
    if (pos == 0)
        return make_pair(1,0); // End of pattern

    --pos;
    if (sar.first > sar.second)
    { 
        // Previous k-mer was not found in the index:
        // restarting the search from the new pos
        sar = kmerToSARange(pattern + pos);
    }
    else
    {
        // Previous SA range was valid:
        // update with the next symbol
        sar.first  = LF(pattern[pos], sar.first-1);
        sar.second = LF(pattern[pos], sar.second)-1;
        if (sar.first > sar.second)
            return make_pair(1,0); // The k-mer at position pos was not found
        // Truncate the search to k symbols
        sar.first  = Blcp->prev(sar.first);
        sar.second = Blcp->next(sar.second + 1) - 1;
    }
    return sar;
}

/**
 * Initialize move left
 *
 * Returns an internal value that corresponds to the end
 * of the given read.
 *
 * Input: Read number
 * Output: <internal value>
 */
ulong CGkArray::initMoveLeft(unsigned j) const
{
    ulong i = j; // Position of the '\0' terminator of read j
    
    // Move left over k-1 symbols
    unsigned k = gk - 1;
    ulong alphabetrank_i_tmp = 0;
    uchar c  = alphabetrank->access(i, alphabetrank_i_tmp);
    while (k--) 
    {
        i = C[c]+alphabetrank_i_tmp-1;
        c = alphabetrank->access(i, alphabetrank_i_tmp);
    }
    return i;
}

/**
 * Initialize move left for arbitrary pattern
 *
 * Returns <internal pointer> that corresponds to the end
 * of the given read.
 *
 * Input: Arbitrary string, assuming '\0'-terminated
 * Output: <internal pointer> that points to the last k-mer of the given pattern.
 */
CGkArray::internal_pointer CGkArray::initMoveLeft(uchar const *pattern) const
{
    unsigned pos = strlen((char const *)pattern) - gk + 1;
    sa_range sar = kmerToSARange(pattern + pos);
    return make_pair(sar, pos);
}

CGkArray::position_vector CGkArray::reportOccs(sa_range const &range) const
{
    ulong sp = range.first;
    ulong ep = range.second;
    position_vector result;
    result.reserve(ep-sp+1);
        
    ulong tmp_rank_c = 0; // Cache rank value of c.
    for (; sp <= ep; ++sp)
    {
        ulong i = sp;
        ulong dist = 0;
        uchar c = alphabetrank->access(i, tmp_rank_c);
        while (c != '\0' && !sampled->access(i))
        {
            i = C[c]+tmp_rank_c-1; //alphabetrank->rank(c,i)-1;
            c = alphabetrank->access(i, tmp_rank_c);
            ++ dist;
        }
        if (c == '\0')
        {
            // Rank among the end-markers in BWT
            unsigned docId = Doc->access(tmp_rank_c-1);
            result.push_back(make_pair(docId, dist)); 
        }
        else
            result.push_back(textPosToReadPos((*suffixes)[sampled->rank1(i)-1] + dist));
    }
    return result;
}

/**
 * Constructor inits an empty dynamic FM-index.
 * Samplerate defaults to TEXTCOLLECTION_DEFAULT_SAMPLERATE.
 */
CGkArray::CGkArray(uchar * bwt, ulong length, unsigned samplerate_, unsigned numberOfTexts_, 
                   ulong maxTextLength_, unsigned gk_, bool verbose)
    : n(length), samplerate(samplerate_), alphabetrank(0), sampled(0), Blast(0), Blcp(0), gk(gk_),
      suffixes(0), positions(0), textStartPos(0), numberOfTexts(numberOfTexts_), maxTextLength(maxTextLength_), 
      Doc(0)
{
    if (gk < 3)
    {
        cerr << "CGkArray::CGkArray(): error: gk < 3" << endl;
        abort();
    }

    makewavelet(bwt); // Deletes bwt!
    bwt = 0;

    Blcp = buildBlcp();

    // Make sampling tables and B_last (requires B_lcp)
    assert(Blcp != 0);
    maketables(verbose);
}

/**
 * Save index to a file handle
 *
 * Throws a std::runtime_error exception on i/o error.
 * First byte that is saved represents the version number of the save file.
 */
void CGkArray::save(std::string const & filename) const
{
    std::string name = filename + ".cgka";
    std::FILE *file = std::fopen(name.c_str(), "wb");

    // Saving version info:
    if (std::fwrite(&versionFlag, 1, 1, file) != 1)
        throw std::runtime_error("CGkArray::save(): file write error (version flag).");

    if (std::fwrite(&(this->n), sizeof(ulong), 1, file) != 1)
        throw std::runtime_error("CGkArray::save(): file write error (n).");
    if (std::fwrite(&(this->samplerate), sizeof(unsigned), 1, file) != 1)
        throw std::runtime_error("CGkArray::save(): file write error (samplerate).");
    if (std::fwrite(&(this->gk), sizeof(unsigned), 1, file) != 1)
        throw std::runtime_error("CGkArray::save(): file write error (gk).");

    if (std::fwrite(this->C, sizeof(unsigned), 256, file) != 256)
        throw std::runtime_error("CGkArray::save(): file write error (C table).");

    if (std::fwrite(&(this->bwtEndPos), sizeof(ulong), 1, file) != 1)
        throw std::runtime_error("CGkArray::save(): file write error (bwt end position).");
    
    HuffWT::save(alphabetrank, file);
    sampled->save(file);
    Blast->save(file);
    Blcp->save(file);

    suffixes->Save(file);
    positions->Save(file);
    
    if (std::fwrite(&(this->numberOfTexts), sizeof(unsigned), 1, file) != 1)
        throw std::runtime_error("CGkArray::save(): file write error (numberOfTexts).");
    if (std::fwrite(&(this->maxTextLength), sizeof(ulong), 1, file) != 1)
        throw std::runtime_error("CGkArray::save(): file write error (maxTextLength).");

    Doc->save(file);

    fflush(file);
    std::fclose(file);
}



const char ALPHABET_DNA[] = {'A', 'C', 'G', 'T', 'N'};
pair<ulong, ulong> intervalList[256];

/*
 * Based on Algorithm 2 proposed in 
 * Timo Beller, Simon Gog, Enno Ohlebusch, and Thomas Schnattinger:
 * Computing the Longest Common Prefix Array Based on the Burrows-Wheeler Transform
 * SPIRE 2011, LNCS 7024, pp. 197–208, 2011.
 */
void CGkArray::traverseBWT(uint *lcp)
{
    deque<pair<pair<ulong, ulong>, uchar> > intervals;
    intervals.push_back(make_pair(make_pair(0,n-1), 0));

    while (!intervals.empty())
    {
        ulong s = intervals.front().first.first;
        ulong e = intervals.front().first.second;
        uchar l = intervals.front().second;
        intervals.pop_front();
        alphabetrank->getIntervals(s, e);
        for (const char *c = ALPHABET_DNA; c < ALPHABET_DNA + sizeof(ALPHABET_DNA); ++c)
        {
            ulong nmin = intervalList[(int)*c].first;
            ulong nmax = intervalList[(int)*c].second;
            if (nmin > nmax || nmax == ~0lu)
                continue;
            intervalList[(int)*c] = make_pair(1,0);
            if (!bitget32(lcp, nmax+1)) 
            {
                bitset32(lcp, nmax+1);
                if ((unsigned)l+1 < gk)
                    intervals.push_back(make_pair(make_pair(nmin, nmax), l+1));
            }
        }
    }
}
void CGkArray::traverseBWT(uint *lcp, ulong s, ulong e, unsigned l)
{
    deque<pair<pair<ulong, ulong>, uchar> > intervals;
    intervals.push_back(make_pair(make_pair(s,e), l));

    while (!intervals.empty())
    {
        ulong s = intervals.front().first.first;
        ulong e = intervals.front().first.second;
        uchar l = intervals.front().second;
        intervals.pop_front();
        alphabetrank->getIntervals(s,e);
        for (const char *c = ALPHABET_DNA; c < ALPHABET_DNA + sizeof(ALPHABET_DNA); ++c)
        {
            ulong nmin = intervalList[(int)*c].first;
            ulong nmax = intervalList[(int)*c].second;
            if (nmin > nmax || nmax == ~0lu)
                continue;
            intervalList[(int)*c] = make_pair(1,0);
            for (ulong i = nmin; i <= nmax+1; ++i)
                bitset32(lcp, i);
            if (l < gk)
                intervals.push_back(make_pair(make_pair(nmin, nmax), l+1));
        }
    }
}

static_bitsequence_brw32 * CGkArray::buildBlcp()
{
    uint *lcp = new uint[(n+1)/32+1];
    for (ulong i = 0; i < (n+1)/32+1; ++i)
        lcp[i] = 0;
    bitset32(lcp, 0);
    bitset32(lcp, n);
    // Traverse all except suffixes with '\0' in their gk-length prefix
    alphabetrank->setC(C);
    alphabetrank->setList(intervalList);
    for (unsigned i = 0; i < 256; ++i)
        intervalList[i] = make_pair(1,0);
    traverseBWT(lcp);

    // Traverse suffixes with '\0' in their gk-length prefix
    ulong nmin = 0;
    ulong nmax = LF(0, n-1)-1;
    for (unsigned i = 0; i < 256; ++i)
        intervalList[i] = make_pair(1,0);
    traverseBWT(lcp, nmin, nmax, 2);
    for (; nmin <= nmax; ++nmin)
        bitset32(lcp, nmin);

    static_bitsequence_brw32 *b = new static_bitsequence_brw32(lcp, n+1, 16);
    delete [] lcp; 
    return b;
}

/**
 * Load index from a file handle
 *
 * Throws a std::runtime_error exception on i/o error.
 * For more info, see CGkArray::save().
 */
CGkArray::CGkArray(std::string const & filename)
    : n(0), samplerate(0), alphabetrank(0), sampled(0), Blast(0), Blcp(0), gk(0), suffixes(0), positions(0),
      textStartPos(0), numberOfTexts(0), maxTextLength(0), Doc(0) 
{
    // Load textStartPos
    {
        std::ifstream ifs(filename + ".cgka_map");
        if (!ifs.good())
        { std::cerr << "error: unable to read input file " << filename << ".cgka_map" << std::endl; std::abort(); }
        textStartPos = new CSA::DeltaVector(ifs);
    }
    

    std::string name = filename + ".cgka";
    std::FILE *file = std::fopen(name.c_str(), "rb");
    if (!file)
    { std::cerr << "error: unable to read input file " << name << std::endl; std::abort(); }

    uchar verFlag = 0;
    if (std::fread(&verFlag, 1, 1, file) != 1)
        throw std::runtime_error("file read error: incorrect version flag! Please reconstruct the index");
    if (verFlag != CGkArray::versionFlag)
        throw std::runtime_error("CGkArray::CGkArray(): invalid save file version.");

    if (std::fread(&(this->n), sizeof(ulong), 1, file) != 1)
        throw std::runtime_error("CGkArray::CGkArray(): file read error (n).");
    if (std::fread(&samplerate, sizeof(unsigned), 1, file) != 1)
        throw std::runtime_error("CGkArray::CGkArray(): file read error (samplerate).");
    if (std::fread(&(this->gk), sizeof(unsigned), 1, file) != 1)
        throw std::runtime_error("CGkArray::CGkArray(): file read error (gk).");

    if (std::fread(this->C, sizeof(unsigned), 256, file) != 256)
        throw std::runtime_error("CGkArray::CGkArray(): file read error (C table).");

    if (std::fread(&(this->bwtEndPos), sizeof(ulong), 1, file) != 1)
        throw std::runtime_error("CGkArray::CGkArray(): file read error (bwt end position).");

    alphabetrank = HuffWT::load(file);
    sampled = (static_bitsequence_brw32 *)static_bitsequence::load(file);
    Blast = (static_bitsequence_brw32 *)static_bitsequence::load(file);
    Blcp = (static_bitsequence_brw32 *)static_bitsequence::load(file);

    suffixes = new BlockArray(file);
    positions = new BlockArray(file);

    if (std::fread(&(this->numberOfTexts), sizeof(unsigned), 1, file) != 1)
        throw std::runtime_error("CGkArray::CGkArray(): file read error (numberOfTexts).");
    if (std::fread(&(this->maxTextLength), sizeof(ulong), 1, file) != 1)
        throw std::runtime_error("CGkArray::CGkArray(): file read error (maxTextLength).");

    Doc = new ArrayDoc(file);

    
/**
   debug
    for (ulong j = 0; j < n; ++j)
    {
        position_result pr;
        getPosition(pr, j);
        cerr << j << " SA " << pr.first << "," << pr.second << " Blcp " << (Blcp->access(j) ? "1" : "0");
        cerr << " Blast " << (Blast->access(j) ? "1" : "0") << "  ";
        cerr << (char *)getSuffix(j, getGkSize());
        cerr << endl;
    }
*/
    /*cerr << "Blast ";
    for (ulong j = 0; j < n; ++j)
        cerr << (Blast->access(j) ? "1" : "0");
    cerr << endl;
    */
    std::fclose(file);
}


CGkArray::~CGkArray() {
    HuffWT::deleteHuffWT(alphabetrank);
    delete sampled;
    delete suffixes;
    delete positions;
    delete textStartPos; 
    delete Doc;
    delete Blast;
    delete Blcp;
}

void CGkArray::makewavelet(uchar *bwt)
{
    ulong i;
    for (i=0;i<256;i++)
        C[i]=0;
    for (i=0;i<n;++i)
        C[(int)bwt[i]]++;
    
    ulong prev=C[0], temp;
    C[0]=0;
    for (i=1;i<256;i++) {          
        temp = C[i];
        C[i]=C[i-1]+prev;
        prev = temp;
    }
    alphabetrank = HuffWT::makeHuffWT(bwt, n);
    // bwt was already deleted.
    bwt = 0;
}

void CGkArray::maketables(bool verbose)
{
    // Calculate BWT end-marker position (of last inserted text)
    {
        ulong i = 0; // This is the end-marker of first text
        ulong alphabetrank_i_tmp = 0;
        uchar c  = alphabetrank->access(i, alphabetrank_i_tmp);
        while (c != '\0') // This loops over the first text (in backward order) 
                          // identify the end-marker of the last text
        {
            i = C[c]+alphabetrank_i_tmp-1;
            c = alphabetrank->access(i, alphabetrank_i_tmp);
        }

        this->bwtEndPos = i;
    }

    // Init data structures to construct Blast
    unsigned *Bl = new unsigned[n/32+1];
    for (ulong i = 0; i < n/32+1; ++i)
        Bl[i] = 0;
    unordered_map<ulong, ulong> trie;


    // Build up array for text starting positions
//    BlockArray* textStartPos = new BlockArray(numberOfTexts, Tools::CeilLog2(this->n));
//    (*textStartPos)[0] = 0; 

    // Mapping from end-markers to doc ID's:
    BlockArray *endmarkerDocId = new BlockArray(numberOfTexts, Tools::CeilLog2(numberOfTexts));

    unsigned sampleLength = (n%samplerate==0) ? n/samplerate : n/samplerate+1;
    positions = new BlockArray(sampleLength, Tools::CeilLog2(this->n));
    uint *sampledpositions = new uint[n/(sizeof(uint)*8)+1];
    for (ulong i = 0; i < n / (sizeof(uint)*8) + 1; i++)
        sampledpositions[i] = 0;
    
    ulong x,p=bwtEndPos;
    // Keeping track of text position of prev. end-marker seen
    ulong posOfSuccEndmarker = n-1;
    unsigned textId = numberOfTexts;
    ulong alphabetrank_i_tmp =0;
    time_t wctime = time(NULL);
    for (ulong i = n; i > 0;) {
        --i;
        if (verbose && i % 100000000 == 0)
            cerr << "Building positions " << 100.0*(n-i)/n << ", Wall-clock time: " << std::difftime(time(NULL), wctime) << " s." << endl; 

        x=(i==n-1)?0:i+1;
        // Now x == SA[p]

        uchar c = alphabetrank->access(p, alphabetrank_i_tmp);
        //cout << (c ? (char)c : (char)'$'); // debug output T^rev
        if (x % samplerate == 0)
        {
            set_field(sampledpositions,1,p,1);
            (*positions)[x/samplerate] = p;
        }

        if (posOfSuccEndmarker - i > gk)
        {
            ulong r = Blcp->rank1(p);
            ulong& value = trie[r]; // Returns reference to new value 0 if the key-value pair does not exists.
            if (!value)
            {
                value = p+1; // +1 since we assume 0 is equal to "not set"
                bitset32(Bl, p);
            } 
        }                

        if (c == '\0')
        {
            trie.clear(); // Flush the trie
            --textId;
            
            // Record the order of end-markers in BWT:
            ulong endmarkerRank = alphabetrank_i_tmp - 1;
            (*endmarkerDocId)[endmarkerRank] = (textId + 1) % numberOfTexts;
            

            // Store text length and text start position:
            if (textId < (unsigned)numberOfTexts - 1)
            {
//                (*textStartPos)[textId + 1] = x;  // x-1 is text position of end-marker.
                posOfSuccEndmarker = i;
            }

            // LF-mapping from '\0' does not work with this (pseudo) BWT.
            p = textId; // Correct LF-mapping to the last char of the previous text.
        }
        else // Now c != '\0', do LF-mapping:
            p = C[c]+alphabetrank_i_tmp-1;
    }
    trie.clear();
    assert(textId == 0);
    
    if (verbose)
        cerr << "Sampling first phase done. Wall-clock time: " << std::difftime(time(NULL), wctime) << " s." << endl; 

    sampled = new static_bitsequence_brw32(sampledpositions, n, 16);
    delete [] sampledpositions;
    assert(sampled->rank1(n-1) == sampleLength);

    Doc = new ArrayDoc(endmarkerDocId);

    // Record text lengths
/*    textLength = new BlockArray(numberOfTexts, Tools::CeilLog2(maxTextLength));
    for (unsigned i = 0; i < numberOfTexts - 1; ++i)
        (*textLength)[i] = (*textStartPos)[i+1] - (*textStartPos)[i] - 1;
    (*textLength)[numberOfTexts-1] = n - (*textStartPos)[numberOfTexts-1] - 1;
    delete textStartPos;
    textStartPos = 0;*/

    // Suffixes store an offset from the text start position
    suffixes = new BlockArray(sampleLength, Tools::CeilLog2(n+1));
    for(ulong i=0; i<sampleLength; ++i) {
        if (verbose && i % 10000000 == 0)
            cerr << "Building suffixes " << 100.0*(i)/sampleLength 
                 << ", Wall-clock time: " << std::difftime(time(NULL), wctime) << " s." << endl; 

        ulong j = sampled->rank1((*positions)[i]);
        if (j==0) j=sampleLength;
        (*suffixes)[j-1] = (i*samplerate==n)?0:i*samplerate;
    }

    if (verbose)
        cerr << "Sampling second phase done. Wall-clock time: " << std::difftime(time(NULL), wctime) << " s." << endl; 

    Blast = new static_bitsequence_brw32(Bl, n, 16);
    delete [] Bl;
    
    if (verbose)
        cerr << "breakdown of size: " << endl
             << "WT: n/a" << endl
             << "suffixes: " << suffixes->size() << endl
             << "positions: " << positions->size() << endl
             << "sampled: " << sampled->size() << endl
             << "B_last: " << Blast->size() << endl
             << "B_lcp: " << Blcp->size() << endl
             << "Doc: " << Doc->size() << endl
             << "textStartPos: see file *.cgka_map" << endl;
}
