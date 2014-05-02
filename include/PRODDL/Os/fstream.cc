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
//  $Id: fstream.cc,v 1.3 2002/02/03 23:44:50 petter Exp $


#include <more/io/fstream.h>
#include <more/io/ipstream.h>
#include <more/io/filesys.h>
#include <fstream>
#include <cwchar>

namespace more {
namespace io {

  template<typename Char, typename Traits>
    basic_ifstream<Char, Traits>::streambuf_type*
    basic_ifstream<Char, Traits>::
    open(char const* c_name, std::ios_base::openmode mode = std::ios_base::in)
    {
	typedef basic_ipstreambuf<Char, Traits> ipstreambuf_type;
	typedef std::basic_filebuf<Char, Traits> filebuf_type;
	close();

	enum { plain, gz, bz2 } format = plain;
	std::string name = c_name;

	int n = name.length();
	if (n >= 3 && name.substr(n - 3) == ".gz")
	    format = gz;
	else if (n >= 4 && name.substr(n - 4) == ".bz2")
	    format = bz2;
	else if (!file_status(name).exists()) {
	    if (file_status(name + ".gz").exists()) {
		name += ".gz";
		format = gz;
	    }
	    else if (file_status(name + ".bz2").exists()) {
		name += ".bz2";
		format = bz2;
	    }
	    else
		return 0;
	}

	switch (format) {
	  case plain:
	    m_buf = new filebuf_type();
	    static_cast<filebuf_type*>(m_buf)->open(name.c_str(), mode);
	    break;
	  case gz:
	    m_buf = new ipstreambuf_type(("gunzip -c "+name).c_str());
	    break;
	  case bz2:
	    m_buf = new ipstreambuf_type(("bunzip2 -c "+name).c_str());
	    break;
	}
	init(m_buf);
	return m_buf;
    }

  template class basic_ifstream<char, std::char_traits<char> >;
  template class basic_ifstream<wchar_t, std::char_traits<wchar_t> >;

}
}
