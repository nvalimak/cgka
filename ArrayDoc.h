
#ifndef _ARRAYDOC_H_
#define _ARRAYDOC_H_

#include <vector>

class ArrayDoc {
public:
    ArrayDoc(BlockArray *input)
        : data(input)
    {

    }
    ArrayDoc(FILE *fp)
        : data(0)
    {
        data = new BlockArray(fp);
    }

    ~ArrayDoc()
    {
        delete data;
    }

    void save(FILE *fp)
    {
        data->Save(fp);
    }
    
    inline unsigned access(unsigned i)
    {
        return (*data)[i];
    }
    inline std::vector<int> accessAll(unsigned i, unsigned j)
    {
        std::vector<int> res;
        res.reserve(j-i+1);

        for (; i <= j; ++i)
            res.push_back((*data)[i]);

        return res;
    }
    
    std::vector<int> access(unsigned i, unsigned j, unsigned min, unsigned max)
    {
        std::vector<int> res;
        res.reserve(j-i+1);

        for (; i <= j; ++i)
            if ((*data)[i] >= min && (*data)[i] <= max)
                res.push_back((*data)[i]);

        return res;
    }
    

    unsigned count(unsigned i, unsigned j, unsigned min, unsigned max)
    {
        unsigned c = 0;
        for (; i <= j; ++i)
            if ((*data)[i] >= min && (*data)[i] <= max)
                ++c;
        return c;
    }
    
    ulong size() const
    {
        return data->spaceInBits()/8;
    }
    
private:
    BlockArray *data;
};

#endif
