//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_GAP_ALIGN_H__
#define PRODDL_GAP_ALIGN_H__

/*  A GLOBAL ALIGNMENT PROGRAM (GAP):

copyright (c) 1992 Xiaoqiu Huang
The distribution of the program is granted provided no charge is made
and the copyright notice is included.
E-mail: huang@cs.mtu.edu

Proper attribution of the author as the source of the software would
be appreciated: "On global sequence alignment" (CABIOS).
Xiaoqiu Huang
Department of Computer Science
Michigan Technological University
Houghton, MI 49931

The GAP program computes a global alignment of two sequences
without penalizing terminal gaps. It delivers the alignment in
linear space, so long sequences can be aligned. 

Users supply scoring parameters. In the simplest form, users just
provide 3 integers: ms, q and r, where ms is the score of a mismatch
and the score of an i-symbol indel is -(q + r * i). Each match
automatically receives score 10. This simple scoring scheme may be
used for DNA sequences. NOTE: all scores are integers.

In general, users can define an alphabet of characters appearing
in the sequences and a matrix that gives the substitution score
for each pair of symbols in the alphabet. The 127 ASCII characters
are eligible. The alphabet and matrix are given in a file, where
the first line lists the characters in the alphabet and the lower
triangle of the matrix comes next. An example file looks as follows:

ARNDC	       
13
-15  19
-10 -22  11
-20 -10 -20  18
-10 -20 -10 -20  12

Here the -22 at position (3,2) is the score of replacing N by R.
This general scoring scheme is useful for protein sequences where the
set of protein characters and Dayhoff matrix are specified in the file.

The GAP program is written in C and runs under Unix systems on
Sun workstations and under DOS systems on PCs.
We think that the program is portable to many machines.

Sequences to be analyzed are stored in separate files.
An input file contains all characters of a sequence, separated by
newline characters, in linear order. No other characters are allowed.
Since upper case and lower case characters are different, use the same
case consistently. A sample sequence file of 4 lines is shown below.

GAATTCTAATCTCCCTCTCAACCCTACAGTCACCCATTTGGTATATTAAA
GATGTGTTGTCTACTGTCTAGTATCCCTCAAGTAGTGTCAGGAATTAGTC
ATTTAAATAGTCTGCAAGCCAGGAGTGGTGGCTCATGTCTGTAATTCCAG
CACTGGAGAGGTAGAAGTG

To find the best alignment of two sequences in files A and B,
use a command of form

gap  A  B  gs  ms  q  r > result

where gap is the name of the object code, gs is the minimum length
of any gap in the short sequence receiving a constant gap penalty,
ms is a negative integer specifying mismatch weight, q and r are
non-negative integers specifying gap-open and gap-extend penalties,
respectively. Output alignment is saved in the file "result".

For using a scoring matrix defined in file S, use a command of form

gap  A  B  gs  S  q  r > result

Note that ms is replaced by the file S.

Acknowledgments
The functions diff2() and display() evolved from those written by Gene Myers.
We made the following modifications: similarity weights (integer), instead of
distance weights (float), are used, terminal gaps are not penalized, and
any gap of length at least gs in the short sequence is given a constant
penalty.
*/

#include <blitz/array.h>

#include <boost/noncopyable.hpp>

namespace PRODDL {
   

  class GapAlign : public boost::noncopyable {

    ///////////////////  Public Interface  /////////////////////

  public:

    enum { maxAlphabetSymbol = 128 };

    typedef blitz::Array<int,2> IntMatrix;

    typedef blitz::Array<int,1> Ints;

    typedef blitz::Array<char,1> SeqData;

    struct Totals {

      int score;
      int n_matches;
      int n_mismatches;
      int alignment_length;

      // derived values

      int percent_matches;
      int n_gaps;

      Totals() {}

      Totals(int _score, int _no_mat, int _no_mis, int _al_len):
	score(_score), n_matches(_no_mat), n_mismatches(_no_mis), 
	alignment_length(_al_len) {

	percent_matches = _al_len > 0 ? (100*_no_mat)/_al_len : 0;
	n_gaps = _al_len - _no_mat - _no_mis;
      }

    };

    GapAlign(int gapConstLen=1000, int gapOpen=1, int gapExtend=10, int mismatchWeight=-10) {
      init(gapConstLen,gapOpen,gapExtend,mismatchWeight);
    } 

    GapAlign(int gapConstLen, int gapOpen, int gapExtend,
	const SeqData& alphabet,const IntMatrix& subsMatr ) {
         init(gapConstLen,gapOpen,gapExtend,alphabet,subsMatr );
    }

    void init(int gapConstLen=1000, int gapOpen=1, int gapExtend=10, int mismatchWeight=-10 );

    void init(int gapConstLen, int gapOpen, int gapExtend, 
	      const SeqData& alphabet,const IntMatrix& subsMatr );


    // Align sequences 'A' and 'B' and write aligned sequences ('A' and 'B'
    // elements with gap symbol inserted) to 'a' and 'b'.
    // The size of 'a' and 'b' must be at least 'alignmentStorage(A,B)' each.
    // The end of data in 'a' and 'b' will be marked by 'endSymbol' symbol.
    // The optional aggregate results (similarity score, number of matches etc)
    // will be available by calling 'totals()' after 'align()' has been called.

    void
    align(const SeqData& A, const SeqData& B,
	  SeqData& a, SeqData& b);


    const Totals&
    totals() const {

      return m_totals;

    }

    // Storage size guaranteed to accomodate each sequence
    // after inclusion of gap symbols and end symbol.

    int alignmentStorage(const SeqData& A, const SeqData& B) {
      return A.size() + B.size() + 1;
    }


    ///////////////// End of Public Interface ///////////////////


  protected:

    void init_common(int gapConstLen, int gapOpen, int gapExtend);

    /* Append "Delete k" op */

    // diff2 work methods

    inline int gap(int k)  { return (k) <= 0 ? 0 : q+r*(k); }	/* k-symbol indel score */

    inline int gap2(int k)  { return (k) <= 0 ? 0 : ((k) <= gaplen ? q+r*(k) : pay); }

    inline 
    void DEL(int k)
    { 
      al_len += k;	
      if (last < 0)	
	last = S(iS-1) -= (k);
      else
	last = S(iS++) = -(k);
    }
	
    /* Append "Insert k" op */

    inline 
    void INS(int k)	
    { 
      al_len += k;	
      if (last > 0)	
	last = S(iS-1) += (k);
      else				
	last = S(iS++) = (k);
    }
    /* Append "Replace" op */

    inline 
    void REP() 	
    { 
      last = S(iS++) = 0; 
      al_len += 1;	
    }

    void
    GLOBAL(const SeqData& A, const SeqData& B);



    /* diff2(A,B,M,N,tb,te,sc,sr,ec,er) returns the score of an optimum conversion
       between A[1..M] and B[1..N] that begins(ends) with a delete if tb(te) is zero
       and appends such a conversion to the current script. If sc = 0, then
       the beginning deletion is not penalized; if sr = 0, the beginning insertion is
       not penalized; if ec = 0, the ending deletion is not charged; if er = 0;
       then the ending insertion is not charged. Any insertion of length at least
       gaplen is given a constant cost */

    int diff2(const char *A, const char *B, int M, int N, 
	      int tb, int te, int sc, int sr, int ec, int er);


    /* Alignment format routine */

    // Take alignmnet operations in 'S' and construct
    // in 'a' and 'b' aligned sequences by including
    // gap symbols 'gapSymbol' (e.g. '-') and end symbol 'endSymbol' (e.g. '\0').
    // The size of 'a' and 'b' must be at least 'alignmentStorage(A,B)' each.
    // The end of data in 'a' and 'b' will be marked by 'endSymbol' symbol.

    void
    format(const SeqData& A, const SeqData& B, const Ints& S,
	   SeqData& a, SeqData& b);



    //////////////////// Data Members  /////////////////////////

  protected:

    // Aggregate results of alignment will be stored here and updated after each
    // call to 'align()'

    Totals m_totals;

    // symbols to use for gaps and end marker in aligned format of sequences

    char gapSymbol, endSymbol;

    //////////// variables to be assigned by initialization routine ////////////

    int match, mismh;		/* max and min substitution weights */
    IntMatrix v;	/* substitution scores */
    int  q, r;       /* gap penalties */
    int  qr;         /* qr = q + r */
    int  gaplen;     /* minimum length for constant-cost insertion */
    int  pay;	/* constant-cost for long insertion */

    ///////////////  diff2() work data  ///////////////////

    Ints CC, DD;			/* saving matrix scores */
    Ints RR, SS;		 	/* saving start-points */
    Ints S;			/* alignment operations created by diff */

    int  zero;				/* int type zero        */

    int iS;				/* Current script append index into S */
    int  last;				/* Last script op appended */

    int no_mat; 				/* number of matches */ 
    int no_mis; 				/* number of mismatches */ 
    int al_len; 				/* length of alignment */

    // Debugging bound checking info for diff2

    const char *begA, *endA, *begB, *endB;

  }; // class GapAlign

} // namespace PRODDL


#endif // PRODDL_GAP_ALIGN_H__
