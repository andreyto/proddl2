//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_COMMON_ALGOR_H__
#define PRODDL_COMMON_ALGOR_H__

// This header combines several commonly used headers
// It also injects some names from blitz specific namespaces into PRODDL namespaces

#include "PRODDL/Common/math.hpp"

#include "PRODDL/Blitz/bztinyvec.hpp"

#include "PRODDL/Blitz/bzarray_iter.hpp"

#include "PRODDL/Common/bz_vect_ext.hpp"


namespace PRODDL { namespace Math {

  using namespace blitz_ext;

}}

#endif // PRODDL_COMMON_ALGOR_H__
