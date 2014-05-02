//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
//  Copyright (C) 2001--2002  Petter Urkedal (petter.urkedal@matfys.lth.se)
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
//  $Id: fstream.h,v 1.2 2002/02/03 23:44:50 petter Exp $


#ifndef MORE_IO_FSTREAM_H
#define MORE_IO_FSTREAM_H

#include <ios>
#include <streambuf>

namespace more {
namespace io {

  /** \class basic_ifstream fstream.h more/io/fstream.h
   **
   ** Basically the same as the \c std:: one, but also provides
   ** on-the-fly decompression of \c bzip2 or \c gzip compressed
   ** files. */
  template<typename Char, typename Traits = std::char_traits<Char> >
    struct basic_ifstream : std::basic_istream<Char, Traits>
    {
	typedef Char char_type;
	typedef Traits traits_type;
	typedef typename traits_type::int_type int_type;
	typedef typename traits_type::pos_type pos_type;
	typedef typename traits_type::off_type off_type;
      private:
	typedef std::basic_streambuf<Char, Traits> streambuf_type;

      public:
	basic_ifstream()
	    : std::basic_istream<Char, Traits>(0), m_buf(0) {}

	basic_ifstream(char const* name,
		       std::ios_base::openmode mode = std::ios_base::in)
	    : std::basic_istream<Char, Traits>(0), m_buf(0)
	{
	    open(name, mode);
	}

	~basic_ifstream() { close(); }

	streambuf_type* rdbuf() { return m_buf; }
	bool is_open() { return m_buf; }
	streambuf_type* open(char const*, std::ios_base::openmode);
	streambuf_type* close() { delete m_buf; m_buf = 0; return 0; }

      private:
	streambuf_type* m_buf;
	bool m_is_pipe;
    };

  typedef basic_ifstream<char, std::char_traits<char> > ifstream;

}
}

#endif
