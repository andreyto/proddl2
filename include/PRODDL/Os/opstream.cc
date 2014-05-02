//  Copyright (C) 2001  Petter Urkedal (petter.urkedal@matfys.lth.se)

//  This file is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.

//  This file is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.

//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

//  As a special exception, you may use this file as part of a free
//  software library without restriction.  Specifically, if other files
//  instantiate templates or use macros or inline functions from this
//  file, or you compile this file and link it with other files to
//  produce an executable, this file does not by itself cause the
//  resulting executable to be covered by the GNU General Public
//  License.  This exception does not however invalidate any other
//  reasons why the executable file might be covered by the GNU General
//  Public License.

//  $Id: opstream.cc,v 1.2 2001/09/26 17:00:47 petter Exp $


#include <iostream>
#include <stdexcept>
#include "opstream.h"

#include <stdio.h>
#include <unistd.h>

namespace more {
namespace io {

  template<typename CharT, typename Traits>
    basic_opstreambuf<CharT, Traits>::
    basic_opstreambuf(char const* cmd)
	: m_fp(popen(cmd, "w")), m_buf(new CharT[s_nbuf])
    {
	if (m_buf)  setp(m_buf, m_buf + s_nbuf);
    }

  template<typename CharT, typename Traits>
    basic_opstreambuf<CharT, Traits>::
    ~basic_opstreambuf()
    {
	sync();
	delete m_buf;
	if (m_fp)  pclose((FILE*)m_fp);
    }

  template<typename CharT, typename Traits>
    basic_opstreambuf<CharT, Traits>::int_type
    basic_opstreambuf<CharT, Traits>::
    overflow(int_type c)
    {
	if (pbase() && pbase() != pptr()) {
	    if (fwrite(pbase(), sizeof(CharT), pptr() - pbase(), (FILE*)m_fp)
		!= (size_t)(pptr() - pbase()))
		return Traits::eof();
	    setp(pbase(), epptr());
	}
	if (!Traits::eq_int_type(c, Traits::eof())) {
	    if (fputc(c, (FILE*)m_fp) == EOF)
		return Traits::eof();
	    return c;
	}
	else
	    return Traits::not_eof(c);
    }

  template<typename CharT, typename Traits>
    int
    basic_opstreambuf<CharT, Traits>::
    sync()
    {
	overflow();
	fflush((FILE*)m_fp);
	return 0;
    }

  template class basic_opstreambuf< char, std::char_traits<char> >;
  template class basic_opstream< char, std::char_traits<char> >;

}
} // namespace more
