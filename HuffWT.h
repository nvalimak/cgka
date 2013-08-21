#ifndef _HUFFWT_H_
#define _HUFFWT_H_


#include "BitRank.h"

#include <cstdio>
#include <stdexcept>
#include <utility>

class HuffWT 
{
public:
    class TCodeEntry 
    {
    public:
        unsigned count;
        unsigned bits;
        unsigned code;
        TCodeEntry() {count=0;bits=0;code=0u;};
        void load(std::FILE *file)
        {
            if (std::fread(&count, sizeof(unsigned), 1, file) != 1)
                throw std::runtime_error("TCodeEntry: file read error (Rs).");
            if (std::fread(&bits, sizeof(unsigned), 1, file) != 1)
                throw std::runtime_error("TCodeEntry: file read error (Rs).");
            if (std::fread(&code, sizeof(unsigned), 1, file) != 1)
                throw std::runtime_error("TCodeEntry: file read error (Rs).");
        }
        void save(std::FILE *file)
        {
            if (std::fwrite(&count, sizeof(unsigned), 1, file) != 1)
                throw std::runtime_error("TCodeEntry: file write error (Rs).");
            if (std::fwrite(&bits, sizeof(unsigned), 1, file) != 1)
                throw std::runtime_error("TCodeEntry: file write error (Rs).");
            if (std::fwrite(&code, sizeof(unsigned), 1, file) != 1)
                throw std::runtime_error("TCodeEntry: file write error (Rs).");
        }
    };

private:
    BitRank *bitrank;
    HuffWT *left;
    HuffWT *right;
    TCodeEntry *codetable;
    uchar ch;
    bool leaf;
    unsigned *C;
    std::pair<ulong, ulong> *intervals;

    HuffWT(uchar *, ulong, TCodeEntry *, unsigned);
    void save(std::FILE *);
    HuffWT(std::FILE *, TCodeEntry *);
public:
    static HuffWT * makeHuffWT(uchar *bwt, ulong n);
    static HuffWT * load(std::FILE *);
    static void save(HuffWT *, std::FILE *);
    static void deleteHuffWT(HuffWT *);
    ~HuffWT(); 

    // C needs to be an array of [0..255]
    void setC(unsigned *C_)
    {
        C = C_;
        if (left) left->setC(C_);
        if (right) right->setC(C_);
    }
    // list needs to be an array of [0..255]
    void setList(std::pair<ulong, ulong> *list)
    {
        intervals = list;
        if (left) left->setList(list);
        if (right) right->setList(list);
    }
    // Updates the list of intervals. Let [i,j] be the interval of substring w.
    // Then the new intervals cover cw intervals for all symbols c.
    //  Beller et al. Computing the Longest Common Prefix Array Based on the {B}urrows-{W}heeler
    //  Transform. Proc. 18th Intl. Symp. String Processing and Information Retrieval, 2011.
    void getIntervals(ulong i, ulong j)
    {        
        if (leaf)
        {
            intervals[ch] = std::make_pair(C[(int)ch] + i, C[(int)ch] + j);
            return;
        }
        ulong a1 = bitrank->rank(i-1);
        ulong b1 = bitrank->rank(j)-1;
        ulong a0 = i-a1;
        ulong b0 = j-b1-1;
        if (b1 >= a1)
            right->getIntervals(a1,b1);
        if (b0 >= a0)
            left->getIntervals(a0,b0);
    }

    inline ulong rank(uchar c, ulong i) const { // returns the number of characters c before and including position i
        HuffWT const *temp=this;
        if (codetable[c].count == 0) return 0;
        unsigned level = 0;
        unsigned code = codetable[c].code;
        while (!temp->leaf) {
            if ((code & (1u<<level)) == 0) {
                i = i-temp->bitrank->rank(i); 
                temp = temp->left; 
            }
            else { 
                i = temp->bitrank->rank(i)-1; 
                temp = temp->right;
            }
            ++level;
        } 
        return i+1;
    };   

    inline ulong select(uchar c, ulong i, unsigned level = 0) const 
    {
        if (leaf)
            return i-1;

//        if (codetable[c].count == 0) return 0;
//        unsigned level = 0;
        unsigned code = codetable[c].code;
            if ((code & (1u<<level)) == 0) {
                i = left->select(c, i, level+1);
                i = bitrank->select0(i+1); 
            }
            else 
            {
                i = right->select(c, i, level+1);
                i = bitrank->select(i+1); 
            }
        return i;
    };   

    inline bool IsCharAtPos(uchar c, ulong i) 
    {
        HuffWT const *temp=this;
        if (codetable[c].count == 0) return false;
        unsigned level = 0;
        unsigned code = codetable[c].code;      
        while (!temp->leaf) {
            if ((code & (1u<<level))==0) {
                if (temp->bitrank->IsBitSet(i)) return false;
                i = i-temp->bitrank->rank(i); 
                    temp = temp->left; 
            }
            else { 
                if (!temp->bitrank->IsBitSet(i)) return false;         
                i = temp->bitrank->rank(i)-1; 
                temp = temp->right;
            }
            ++level;
        } 
        return true;
    }
    inline uchar access(ulong i) const 
    {
        HuffWT const *temp=this;
        while (!temp->leaf) {
            if (temp->bitrank->IsBitSet(i)) {
                i = temp->bitrank->rank(i)-1;
                temp = temp->right;
            }
            else {
                i = i-temp->bitrank->rank(i); 
                temp = temp->left;      
            }         
        }
        return (int)temp->ch;
    }

    inline uchar access(ulong i, ulong &rank) const
    {
        HuffWT const *temp=this;
        while (!temp->leaf) {
            if (temp->bitrank->IsBitSet(i)) {
                i = temp->bitrank->rank(i)-1;
                temp = temp->right;
            }
            else {
                i = i-temp->bitrank->rank(i); 
                temp = temp->left;      
            }         
        }
        rank = i+1;
        return (int)temp->ch;
    }
};
#endif
