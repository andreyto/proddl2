//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//

#include "PRODDL/gap_align.hpp"

#include "PRODDL/Common/logger.hpp"

namespace PRODDL {

  void GapAlign::init(int gapConstLen, int gapOpen, int gapExtend, int mismatchWeight ) {
    init_common(gapConstLen,gapOpen,gapExtend);
    int ms = mismatchWeight;
    ATLOG_SWITCH_1(dbg::assertion(dbg::error,DBG_ASSERTION(mismatchWeight < 0)));
    v.resize(maxAlphabetSymbol,maxAlphabetSymbol);
    match = 10;
    mismh = ms;
    /* set match and mismatch weights */
    for (int i = 0; i < v.rows() ; i++ )
      for (int j = 0; j < v.columns() ; j++ ) {
	if (i == j )
	  v(i,j) = 10;
	else
	  v(i,j) = mismh;
      }
  }

  void GapAlign::init(int gapConstLen, int gapOpen, int gapExtend, 
		      const SeqData& alphabet,const IntMatrix& subsMatr ) {
    init_common(gapConstLen,gapOpen,gapExtend);
    v.resize(maxAlphabetSymbol,maxAlphabetSymbol);
    v = 0;
    match = 0;
    mismh = 0;
    for (int i = 0; i < alphabet.size() ; i++ ) {
      int i_s = alphabet(i);
      ATLOG_SWITCH_1(dbg::assertion(dbg::error,DBG_ASSERTION( i_s != gapSymbol && i_s != endSymbol )));
      for (int j = 0; j < alphabet.size() ; j++ ) { 
	int j_s = alphabet(j);
	int val = subsMatr(i,j);
	if ( val > match ) match = val;
	if ( val < mismh ) mismh = val;
	v(i_s,j_s) = val;
      }
    }
  }

  void
  GapAlign::align(const SeqData& A, const SeqData& B,
		  SeqData& a, SeqData& b) {
    if( A.size() >= B.size() ) {
      GLOBAL(A,B);
      format(A,B,S,a,b);
    }
    else {
      GLOBAL(B,A);
      format(B,A,S,b,a);
    }
  }

  void 
  GapAlign::init_common(int gapConstLen, int gapOpen, int gapExtend) 
  {
    zero = 0;
    gaplen = gapConstLen;
    ATLOG_SWITCH_1(dbg::assertion(dbg::error,DBG_ASSERTION(gapConstLen > 0)));
    q = gapOpen;
    ATLOG_SWITCH_1(dbg::assertion(dbg::error,DBG_ASSERTION(gapOpen >= 0)));
    r = gapExtend;
    ATLOG_SWITCH_1(dbg::assertion(dbg::error,DBG_ASSERTION(gapExtend >= 0)));
    pay = q + r * gaplen;
    qr = q + r;
    gapSymbol = '-';
    endSymbol = '\0';
  }      

  void
  GapAlign::GLOBAL(const SeqData& A, const SeqData& B)
  { 

    int M = A.size();
    int N = B.size();

    begA = A.dataFirst();
    endA = begA + M;
    begB = B.dataFirst();
    endB = begB + N;

    /* allocate space for all vectors */

    int j = (N + 1);
    CC.resize(j);
    DD.resize(j);
    RR.resize(j);
    SS.resize(j);
    int i = (M + 1);
    S.resize(i + j);
    iS = 0;
    last = 0;
    al_len = 0;
    no_mat = 0;
    no_mis = 0;
    //diff2() assumes utit-based offsets
    int score = diff2(A.dataFirst()-1,B.dataFirst()-1,A.size(),B.size(),
		      q,q,0,0,0,0);
    m_totals = Totals(score,no_mat,no_mis,al_len);
  }

  int 
  GapAlign::diff2(const char *A, const char *B, int M, int N, 
		  int tb, int te, int sc, int sr, int ec, int er)
  { 

    int   midi, midj, type;	/* Midpoint, type, and cost */
    int midc;
    int  ss,cc;

    int iV;

    int i,j;

    int c, e, d, s;
    int t;
    int  g, temp;

    /* Boundary cases: M <= 1 or N == 0 */

    if (N <= 0)
      { if (M > 0) DEL(M);
	if ( !sc || !ec )
	  return 0;
	else
	  return - gap(M);
      }
    if (M <= 1)
      { if (M <= 0)
	  { INS(N);
	    if ( !sr || !er )
	      return 0;
	    else
	      return - gap2(N);
	  }
	int midc = - (sc * (tb + r) + er * gap2(N) );
	int midj = -1;
	if ( midc < ( c =  - (ec * (te + r) + sr * gap2(N) ) ) )
	  { midc = c;
	    midj = 0;
	  }
	ATLOG_SWITCH_2(dbg::assertion(dbg::error,DBG_ASSERTION( begA <= (A+1) && (A+1) < endA )));
	iV = A[1];
	for (j = 1; j <= N; j++)
	  { 
	    ATLOG_SWITCH_2(dbg::assertion(dbg::error,DBG_ASSERTION( begB <= (B+j) && (B+j) < endB )));
	    c = v(iV,static_cast<int>(B[j])) - ( sr * gap2(j-1) + er * gap2(N-j) );
	    if (c > midc)
	      { midc = c;
		midj = j;
	      }
	  }
	if (midj == -1)
	  { DEL(1); INS(N); }
	else
	  if (midj == 0)
	    { INS(N); DEL(1); }
	  else
	    { if (midj > 1) INS(midj-1);
	      REP();
	      ATLOG_SWITCH_2(dbg::assertion(dbg::error,DBG_ASSERTION( begA <= (A+1) && (A+1) < endA && begB <= (B+midj) && (B+midj) < endB )));
	      if ( A[1] == B[midj] )
		no_mat += 1;
	      else
		no_mis += 1;
	      if (midj < N) INS(N-midj);
	    }
	return midc;
      }

    /* Divide: Find optimum midpoint (midi,midj) of cost midc */

    midi = M/2;			/* Forward phase:                          */
    CC(0) = 0;			/*   Compute C(M/2,k) & D(M/2,k) for all k */
    t = - q * sr;
    if ( N <= gaplen )
      for (j = 1; j <= N; j++)
	{ CC(j) = t = (t-r) * sr;
	  DD(j) = t-q;
	}
    else
      { for (j = 1; j <= gaplen; j++)
	  { CC(j) = t = (t-r) * sr;
	    DD(j) = t-q;
	  }
	for (j = gaplen+1; j <= N; j++)
	  { CC(j) = t = -pay * sr;
	    DD(j) = t - q;
	  }
      }
    if ( !ec ) DD(N) += q;
    t = -tb * sc;
    for (i = 1; i <= midi; i++)
      { s = CC(0);
	CC(0) = c = t = (t-r) * sc;
	e = t-q;
	g = t - pay;
	ATLOG_SWITCH_2(dbg::assertion(dbg::error,DBG_ASSERTION( begA <= (A+i) && (A+i) < endA )));
	iV = A[i];
	for (j = 1; j <= N; j++)
	  { if ((c = c - qr) > (e = e - r)) e = c;
	    if ( j == N && !ec )
	      { if ((c = CC(j) ) > (d = DD(j) )) d = c;}
	    else
	      if ((c = CC(j) - qr) > (d = DD(j) - r)) d = c;
	    ATLOG_SWITCH_2(dbg::assertion(dbg::error,DBG_ASSERTION( begB <= (B+j) && (B+j) < endB )));
	    c = s+v(iV,static_cast<int>(B[j]));
	    if (c < d) c = d;
	    if (c < e) c = e;
	    if ( j - gaplen > 0 )
	      { if ( g < ( temp = CC(j-gaplen-1) - pay ) )
		  g = temp;
		if ( c < g ) c = g;
	      }
	    s = CC(j);
	    CC(j) = c;
	    DD(j) = d;
	  }
      }
    DD(0) = CC(0);

    RR(N) = 0;			/* Reverse phase:                          */
    t = -q * er;			/*   Compute R(M/2,k) & S(M/2,k) for all k */
    if ( N <= gaplen )
      for (j = N-1; j >= 0; j--)
	{ RR(j) = t = (t-r) * er;
	  SS(j) = t-q;
	}
    else
      { temp = N - gaplen;
	for (j = N-1; j >= temp; j--)
	  { RR(j) = t = (t-r) * er;
	    SS(j) = t-q;
	  }
	for (j = temp-1; j >= 0; j--)
	  { RR(j) = t = -pay * er;
	    SS(j) = t - q;
	  }
      }
    if ( !sc ) SS(0) += q;
    t = -te * ec;
    for (i = M-1; i >= midi; i--)
      { s = RR(N);
	RR(N) = c = t = (t-r) * ec;
	g = t - pay;
	e = t-q;
	ATLOG_SWITCH_2(dbg::assertion(dbg::error,DBG_ASSERTION( begA <= (A+i+1) && (A+i+1) < endA )));
	iV = A[i+1];
	for (j = N-1; j >= 0; j--)
	  { if ((c = c - qr) > (e = e - r)) e = c;
	    if ( !j && !sc )
	      { if ((c = RR(j) ) > (d = SS(j) )) d = c;}
	    else
	      if ((c = RR(j) - qr) > (d = SS(j) - r)) d = c;
	    ATLOG_SWITCH_2(dbg::assertion(dbg::error,DBG_ASSERTION( begB <= (B+j+1) && (B+j+1) < endB )));
	    c =  s+v(iV,static_cast<int>(B[j+1]));
	    if (c < d) c = d;
	    if (c < e) c = e;
	    if ( j + gaplen < N )
	      { if ( g < ( temp = RR(j+gaplen+1) - pay ) )
		  g = temp;
		if ( c < g ) c = g;
	      }
	    s = RR(j);
	    RR(j) = c;
	    SS(j) = d;
	  }
      }
    SS(N) = RR(N);

    midc = CC(0)+RR(0);		/* Find optimal midpoint */
    midj = 0;
    type = 1;
    for (j = 0; j <= N; j++)
      if ((c = CC(j) + RR(j)) >= midc)
	if (c > midc || CC(j) != DD(j) && RR(j) == SS(j))
	  { midc = c;
	    midj = j;
	  }
    for (j = N; j >= 0; j--)
      { if ( j == N )
	  d = q * ec;
	else
	  if ( j == 0 )
	    d = q * sc;
	  else
	    d = q;
	if ((c = DD(j) + SS(j) + d) > midc)
	  { midc = c;
	    midj = j;
	    type = 2;
	  }
      }

    /* Conquer: recursively around midpoint */

    cc = midj == N ? ec : 1;
    ss = midj == 0 ? sc : 1;
    if (type == 1)
      { 
	diff2(A,B,midi,midj,tb,q,sc,sr,cc,1);
	diff2(A+midi,B+midj,M-midi,N-midj,q,te,ss,1,ec,er);
      }
    else
      { 
	diff2(A,B,midi-1,midj,tb,zero,sc,sr,cc,1);
	DEL(2);
	diff2(A+midi+1,B+midj,M-midi-1,N-midj,zero,te,ss,1,ec,er);
      }
    return midc;
  }

  void
  GapAlign::format(const SeqData& A, const SeqData& B, const Ints& S,
		   SeqData& a, SeqData& b)
  { 

    int i = 0, j = 0, op = 0;
    // max length of each sequence after alignmnet insertions is always <= maxLine
    int maxLine = alignmentStorage(A,B);

    ATLOG_SWITCH_1(dbg::assertion(dbg::error,DBG_ASSERTION( a.size() >= maxLine && b.size() >= maxLine )));
    ATLOG_SWITCH_1(dbg::assertion(dbg::error,DBG_ASSERTION( &a != &A && &b != &B && &a != &B && &b != &A && &a != &b )));

    int ia = 0, ib = 0;

    int M = A.size();
    int N = B.size();

    iS = 0;

    while (i < M || j < N)
      { if (op == 0 && S(iS) == 0)
	  { op = S(iS++);
	    a(ia++) = A(i++);
	    b(ib++) = B(j++);
	  }
	else
	  { if (op == 0)
	      op = S(iS++);
	    if (op > 0)
	      { a(ia++) = '-';
		b(ib++) = B(j++);
		op--;
	      }
	    else
	      { a(ia++) = A(i++);
		b(ib++) = '-';
		op++;
	      }
	  }
      }
    a(ia) = '\0';
    b(ib) = '\0';
  }


} // namespace PRODDL
