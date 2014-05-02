//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_MPI_TAGS_H__
#define PRODDL_MPI_TAGS_H__


// Central repository for MPI tags

namespace PRODDL {

  namespace Mpi {


    const int TAG_START = 101;
    const int TAG_EXIT =                   TAG_START + 1;
    const int TAG_TINY_VECT =              TAG_START + 2;
    const int TAG_ARRAY_DATA =             TAG_START + 3;
    const int TAG_ARRAY_STRUCT =           TAG_START + 4;
    const int TAG_ROTATION =               TAG_START + 5;
    const int TAG_TRANVALUES =             TAG_START + 6;

    

  }} // namespace PRODDL::Mpi

#endif // PRODDL_MPI_TAGS_H__
