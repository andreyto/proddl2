//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
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
//  $Id: ipstream.h,v 1.2 2001/12/15 19:23:49 petter Exp $


#include <ios>
#include <iostream>

namespace more {
namespace io {

  /** \class basic_ipstreambuf ipstream.h more/io/ipstream.h
   ** IO buffer of an input pipe stream. */
  template<typename Char, typename Traits>
    struct basic_ipstreambuf
	: std::basic_streambuf<Char, Traits>
    {
	typedef Char char_type;
	typedef Traits traits_type;
 	typedef typename Traits::int_type int_type;
	typedef typename Traits::pos_type pos_type;
	typedef typename Traits::off_type off_type;

	basic_ipstreambuf(char const*);
	virtual ~basic_ipstreambuf();
	bool is_open() const { return m_fp; }

      protected:
	virtual int_type underflow();

      private:
	static int const s_nbuf = 256;

	void* m_fp;
	char_type m_buf[s_nbuf];
    };

  /** \class basic_ipstream ipstream.h more/io/ipstream.h
   ** An input pipe stream. */
  template< typename Char, typename Traits=std::char_traits<Char> >
    struct basic_ipstream : std::basic_istream<Char, Traits>
    {
	typedef basic_ipstreambuf<Char, Traits> streambuf_type;
	basic_ipstream(char const* cmd)
	    : std::basic_istream<Char, Traits>(m_buf=new streambuf_type(cmd))
	{
	    if (!is_open()) setstate(std::ios_base::failbit);
	}

	~basic_ipstream() { delete m_buf; }

	bool is_open() const { return m_buf->is_open(); }

      private:
	streambuf_type* m_buf;
    };

  typedef basic_ipstream<char> ipstream;

}
}
