//  Copyright (C) 2001  Petter Urkedal (petter.urkedal@matfys.lth.se)
//
//  This file is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This file is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//  As a special exception, you may use this file as part of a free
//  software library without restriction.  Specifically, if other files
//  instantiate templates or use macros or inline functions from this
//  file, or you compile this file and link it with other files to
//  produce an executable, this file does not by itself cause the
//  resulting executable to be covered by the GNU General Public
//  License.  This exception does not however invalidate any other
//  reasons why the executable file might be covered by the GNU General
//  Public License.
//
//  $Id: ipstream.cc,v 1.2 2001/09/07 23:36:30 petter Exp $


#include <more/io/ipstream.h>
#include <stdexcept>
#include <cwchar>

namespace more {
namespace io {

  template<typename Char, typename Traits>
    basic_ipstreambuf<Char, Traits>::
    basic_ipstreambuf(char const* cmd)
	: m_fp(popen(cmd, "r")) //, m_buf(new char_type[s_nbuf])
    {
	setg(m_buf, m_buf, m_buf);
    }

  template<typename Char, typename Traits>
    basic_ipstreambuf<Char, Traits>::
    ~basic_ipstreambuf()
    {
	if (m_fp) pclose((FILE*)m_fp);
    }

  template<typename Char, typename Traits>
    basic_ipstreambuf<Char, Traits>::int_type
    basic_ipstreambuf<Char, Traits>::underflow()
    {
	if (m_fp == 0)
	    throw std::logic_error("more::io::basic_ipstreambuf: "
				   "File is not open.");
	std::size_t nread = std::fread(eback(), sizeof(char_type), s_nbuf,
				       (FILE*)m_fp);
	if (nread == 0)
	    return traits_type::eof();
	else {
	    setg(eback(), eback(), eback() + nread);
	    return *eback();
	}
    }

  template class basic_ipstreambuf<char, std::char_traits<char> >;
  template class basic_ipstreambuf<wchar_t, std::char_traits<wchar_t> >;

}
}
