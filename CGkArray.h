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

#ifndef _CGkArray_H_
#define _CGkArray_H_

#include "BlockArray.h"
#include "ArrayDoc.h"
#include "HuffWT.h"

// Include from RLCSA
#include "bits/deltavector.h"

// Include from libcds
#include <basics.h> // Defines W == 32
#include <static_bitsequence.h>

// Libcds includes will collide with #define W.
// Re-defining the word size to ulong:
#undef W
#define W (CHAR_BIT*sizeof(unsigned long))
#undef bitset
#undef bitget

#include <string>
#include <vector>
#include <set>

/**
 * Implementation of the Compressed Gk Arrays
 *
 * Use TextCollectionBuilder to construct.
 *
 */
class CGkArray
{
public:
    /**
     * Data types for results
     */
    // Pair of read number (0-based numbering) and read position (0-based numbering)
    typedef std::pair<unsigned, ulong> position_result;
    // Vector of read number and read position pairs
    typedef std::vector<position_result> position_vector;
    // Range of the suffix array
    typedef std::pair<ulong,ulong> sa_range;
    // Internal pointer for move left (on arbitrary patterns)
    typedef std::pair<sa_range,unsigned> internal_pointer;

    /**
     * Convert from text position to a pair of <read number, read position>
     *
     * Both the read number and position follow a 0-based numbering.
     *
     * Input: Text position
     * Output: Read number + read position
     */
    inline position_result textPosToReadPos(ulong i) const
    {
        CSA::DeltaVector::Iterator iter(*textStartPos);
        unsigned read = iter.rank(i)-1;
        return std::make_pair(read, i - iter.select(read));
    }

    /**
     * Convert from a pair of <read number, read position> to text position
     * 
     * Both the read number and position follow a 0-based numbering.
     *
     * Input: Read number + read position
     * Output: Text position
     */
    inline ulong readPosToTextPos(unsigned read, ulong i) const
    {
        CSA::DeltaVector::Iterator iter(*textStartPos);
        return iter.select(read) + i;
    }

    /**
     * Check if given text position is valid (i.e. k-mer does not span over a '\0' byte).
     *
     * Input: Text position, use readPosToTextPos() to convert
     */
    inline bool isValidTextPos(ulong i) const
    {
        CSA::DeltaVector::Iterator iter(*textStartPos);
        unsigned read = iter.rank(i);
        return (iter.select(read) - i > gk ? true : false);
    }

    // Return total length of text (including 0-terminators).
    inline ulong getLength() const
    { return n; }
    // Return the length of the given read (including 0-terminator)
    ulong getLength(unsigned i) const
    {
        CSA::DeltaVector::Iterator iter(*textStartPos);
        return iter.select(i+1) - iter.select(i);
    }
    // Return the k-mer size used in indexing
    unsigned getGkSize() const
    { return gk; }
    // Return the number of indexed reads
    unsigned getNumberOfReads() const
    { return numberOfTexts; }

    /**
     * Find the suffix array range for the given k-mer
     *
     * Input: k-mer
     * Output: Suffix array range
     */
    sa_range kmerToSARange(uchar const *) const;

    /**
     * Find the suffix array range for the k-mer at the given position
     *
     * Input: text position, use readPosToTextPos() to convert
     * Output: Suffix array range
     */
    sa_range textPosToSARange(ulong x) const
    {
        ulong y = inverseSA(x);
        return make_pair(Blcp->prev(y), Blcp->next(y+1)-1);
    }
    
    /**
     * Move left operation gives an efficient way to step over each k-mer in a read.
     *
     * Each call to this method returns a SA range of the preceding read position,
     * that is, the k-mers of the read are traversed in backwards order.
     * This method allows to compute, for example, so called read-coverage profiles.
     *
     * Input: <internal value>, use initMoveLeft() to initialize.
     * Output: Suffix array range of the preceeding read position. 
     *         After the last step, subsequent calls return an empty SA range.
     */
    sa_range moveLeft(ulong &) const;
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
    sa_range moveLeft(internal_pointer &, uchar const *) const;

    /**
     * Initialize move left
     *
     * Returns an internal value that corresponds to the end
     * of the given read.
     *
     * Input: Read number
     * Output: <internal value> that points to the last k-mer of the given read.
     */
    ulong initMoveLeft(unsigned) const;
    /**
     * Initialize move left for arbitrary pattern
     *
     * Returns <internal pointer> that corresponds to the end
     * of the given read.
     *
     * Input: Arbitrary string, assuming '\0'-terminated
     * Output: <internal pointer> that points to the last k-mer of the given pattern.
     */
    internal_pointer initMoveLeft(uchar const *) const;

    /**
     * Q1 In which reads does k-mer occur?
     *
     * Input: text position, use readPosToTextPos() to convert
     * Output: Vector of read numbers and positions (i.e. position of the last occurrence of f in each read)
     */
    inline position_vector reportReads(ulong x) const
    {
        ulong y = inverseSA(x);
        return reportReads(Blcp->prev(y), Blcp->next(y+1)-1);
    }

    /**
     * Q1 In which reads does k-mer occur?
     *
     * Input: range in the suffix array, use kmerToSARange() to find
     * Output: Vector of read numbers and positions (i.e. position of the last occurrence of f in each read)
     */
    inline position_vector reportReads(sa_range const &range) const
    {
        return reportReads(range.first, range.second);
    }

    /** 
     * Q2 In how many reads does k-mer occur?
     *
     * Input: text position, use readPosToTextPos() to convert
     */
    inline unsigned countReads(ulong x) const
    {
        ulong y = inverseSA(x);
        return countReads(Blcp->prev(y), Blcp->next(y+1)-1);
    }

    /** 
     * Q2 In how many reads does k-mer occur?
     *
     * Input: range in the suffix array, use kmerToSARange() to find
     */
    inline unsigned countReads(sa_range const &range) const
    {
        return countReads(range.first, range.second);
    }

    /**
     * Q3 What are the occurrence positions of k-mer?
     *
     * Input: text position, use readPosToTextPos() to convert
     * Output: Vector of read numbers and positions (i.e. all occurrences of f)
     */
    inline position_vector reportOccs(ulong x) const
    {
        ulong y = inverseSA(x);
        return reportOccs(make_pair(Blcp->prev(y), Blcp->next(y+1)-1));
    }

    /**
     * Q3 What are the occurrence positions of k-mer?
     *
     * Input: range of the suffix array, use kmerToSARange() to find
     * Output: Vector of read numbers and positions (i.e. all occurrences of f)
     */
    position_vector reportOccs(sa_range const &) const;

    /**
     * Q4 What is the number of occurrences of k-mer?
     *
     * Input: text position, use readPosToTextPos() to convert
     */
    inline unsigned countOccs(ulong x) const
    {
        ulong y = inverseSA(x);
        return (Blcp->next(y+1)-1) - Blcp->prev(y) + 1;
    }

    /**
     * Q4 What is the number of occurrences of k-mer?
     *
     * Input: range in the suffix array, use kmerToSARange() to find
     */
    inline unsigned countOccs(sa_range const &range) const
    {
        return range.second - range.first + 1;
    }

    /**
     * Returns a copy of an indexed read.
     *
     * Input: Read number (0-based)
     * Output: Copy of the read. Caller must delete [] the buffer.
     */ 
    uchar * getRead(unsigned) const;

    /**
     * Given suffix i and substring length l, return T[SA[i] ... SA[i]+l].
     *
     * Caller must delete [] the buffer.
     */ 
    uchar * getSuffix(ulong dest, unsigned l) const;

    /**
     * Inverse suffix array
     *
     * Input: Text position, use readPosToTextPos() to convert
     */
    ulong inverseSA(ulong i) const
    {
        ulong skip = samplerate - i % samplerate;
        ulong j;
        if (i / samplerate + 1 >= n / samplerate)
        {
            j = bwtEndPos;
            skip = n - i;
        }
        else
            j = (*positions)[i/samplerate+1];                                                                              
        
        ulong tmp_rank_c = 0; // Cache rank value of c.
        while (skip > 0)
        {
            int c = alphabetrank->access(j, tmp_rank_c); 
            if (c == '\0')
            {
                unsigned endmarkerRank = tmp_rank_c-1; 
                j = Doc->access(endmarkerRank); // LF-mapping for '\0'
                if (j==0)
                    j = numberOfTexts-1;
                else
                    --j;
            }
            else
                j = C[c]+tmp_rank_c-1; 
            skip --;
        }
        return j;
    }
 

    // For given suffix i, return corresponding read number and read position.
    position_result getPosition(ulong i) const
    {
        ulong tmp_rank_c = 0; // Cache rank value of c.
        ulong dist = 0;
        uchar c = alphabetrank->access(i, tmp_rank_c);
        while (c != '\0' && !sampled->access(i))
        {
            i = C[c]+tmp_rank_c-1; 
            c = alphabetrank->access(i, tmp_rank_c);
            ++ dist;
        }
        if (c == '\0')
        {
            // Look-up from the rank among the end-markers in BWT
            return std::make_pair(Doc->access(tmp_rank_c-1), dist); 
        }

        return textPosToReadPos((*suffixes)[sampled->rank1(i)-1] + dist);
    }

    CGkArray(uchar *, ulong, unsigned, unsigned, ulong, unsigned, bool);
    // Index from/to disk
    CGkArray(std::string const &);
    void save(std::string const &) const;
    ~CGkArray();

private:
    // Return C[c] + rank_c(L, i) for given c and i
    inline ulong LF(uchar c, ulong i) const
    {
        if (C[(int)c+1]-C[(int)c] == 0) // FIXME fix alphabet
            return C[(int)c];
        return C[(int)c] + alphabetrank->rank(c, i);
    } 

    // Helper method for Q1
    position_vector reportReads(ulong sp, ulong ep) const
    {
        unsigned nreads = Blast->rank1(ep)-Blast->rank1(sp-1);
        position_vector pv;
        pv.reserve(nreads+1);
        while (sp < ep)
        {
            sp = Blast->next(sp);
            pv.push_back(getPosition(sp++));
        }
        return pv;
    }

    // Helper method for Q2
    inline unsigned countReads(ulong sp, ulong ep) const
    {
        return Blast->rank1(ep) - Blast->rank1(sp-1);
    }
 
    // Required by getSuffix(), assuming DNA alphabet
    static const char ALPHABET_SHIFTED[];

    static const uchar versionFlag;
    ulong n;
    unsigned samplerate;
    unsigned C[256];
    ulong bwtEndPos;
    HuffWT *alphabetrank;

    static_bitsequence_brw32 * sampled;
    static_bitsequence_brw32 * Blast;
    static_bitsequence_brw32 * Blcp;
    unsigned gk;
    BlockArray * suffixes;
    BlockArray * positions;
    CSA::DeltaVector * textStartPos;

    // Total number of texts in the collection
    unsigned numberOfTexts;
    // Length of the longest text
    ulong maxTextLength;

    // Array of document id's in the order of end-markers in BWT
    ArrayDoc *Doc;

    uchar * BWT(uchar *);
    void makewavelet(uchar *);
    void maketables(bool);
    void traverseBWT(uint *);
    void traverseBWT(uint *, ulong, ulong, unsigned);
    static_bitsequence_brw32 * buildBlcp();

    /**
     * Count end-markers in given interval
     */
    inline unsigned CountEndmarkers(ulong sp, ulong ep) const
    {
        if (sp > ep)
            return 0;

        ulong ranksp = 0;
        if (sp != 0)
            ranksp = alphabetrank->rank(0, sp - 1);
    
        return alphabetrank->rank(0, ep) - ranksp;
    }

}; // class CGkArray

#endif
