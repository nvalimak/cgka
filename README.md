Compressed Gk Arrays
====

Implementation of Compressed Gk arrays proposed in:

* Niko Välimäki and Eric Rivals: Scalable and Versatile k-mer Indexing 
for High-Throughput Sequencing Data. In Proc. 9th International
Symposium on Bioinformatics Research and Applications (ISBRA 2013),
Springer Verlag, LNBI 7875, pp. 237-248, May 20-22, 2012.


Acknowledgments
----
This is joint work with <a href="http://www.lirmm.fr/~rivals/">Eric Rivals</a> (CRNS/LIRMM, France).
The BWT construction algorithm is based on <a href="https://github.com/lh3/ropebwt">ropebwt/BCR</a> by Heng Li.


TODO
----

* Support for FASTQ and FASTA inputs.
* Use BCR+LCP for the B_lcp construction. (See e.g. Markus J. Bauer, 
    Anthony J. Cox, Giovanna Rosone, Marinella Sciortino: Lightweight LCP 
    Construction for Next-Generation Sequencing Datasets. WABI 2012: 326-337)


Basic usage
----

1) Compile the software by issuing the command `make'.

2) Construct an index for the sequences by `./builder -v -k 8 -s 4 input.txt', 
   where parameter -k determines the k-mer length, and -s determines the 
   sampling rate. See builder.cpp for an example how the index is constructed.
   FASTQ and FASTA inputs are not yet supported.

3) Run an example script with 100 random position queries using
   `./cgkquery -v -q 100 input.txt'.


Brief summary of the CGkArray.h interface
----

For construction, see builder.cpp for an example.

Queries can be issued using either

    1) text position
    2) k-mer  (uchar pointer to a string of length k)

For text positions, you can use methods

    1) textPosToReadPos()
    2) readPosToTextPos()

to convert between text position and a pair of <read number, read position>. 
All numberings are 0-based.

For k-mer queries, use the method kmerToSARange() to recover the corresponding
suffix array range first, and then issue the wanted query. Queries take the SA range
as their only parameter.

Queries from Q1 to Q4 are supported. See the paper for details.
Queries return either a pair <read number, read position> or a vector of said pairs.

cgkarray.cpp contains an example for computing a read-coverage profile
by traversing all the k-mers in a given read.
