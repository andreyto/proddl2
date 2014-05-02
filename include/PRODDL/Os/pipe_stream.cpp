//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_OS_PIPE_STREAM_H__
# error This file must be included from pipe_stream.hpp
#endif

#include <cstdio>

#include "PRODDL/Os/pipe_stream_except.hpp"

namespace PRODDL {

  template<typename Char, typename Traits>
    basic_iopipe_streambuf<Char, Traits>::
    basic_iopipe_streambuf(const char* systemCommand, bool isInput)
	: m_fp( popen(systemCommand, isInput ? "r" : "w"))
    {
    }

  template<typename Char, typename Traits>
    basic_iopipe_streambuf<Char, Traits>::
    ~basic_iopipe_streambuf()
    {
	if (m_fp) {
	  if( pclose((FILE*)m_fp) != 0 )
	    throw pipe_stream_error("Closing pipe stream returned an error");
	}
    }

  template<typename Char, typename Traits>
  basic_ipipe_streambuf<Char, Traits>::
  basic_ipipe_streambuf(const char* systemCommand) :
    basic_iopipe_streambuf<Char, Traits>(systemCommand,true) 
    {
      if (m_buf) setg(m_buf, m_buf, m_buf);
    }

  template<typename Char, typename Traits>
    typename basic_ipipe_streambuf<Char, Traits>::base_class::int_type
    basic_ipipe_streambuf<Char, Traits>::underflow()
    {
	if (m_fp == 0)
	  throw pipe_stream_error("Input pipe stream is not open.");
	std::size_t nread = std::fread(eback(), sizeof(typename base_class::char_type), s_nbuf,
				       (FILE*)m_fp);
	if (nread == 0)
	    return traits_type::eof();
	else {
	    setg(eback(), eback(), eback() + nread);
	    return *eback();
	}
    }


  template<typename Char, typename Traits>
  basic_opipe_streambuf<Char, Traits>::
  basic_opipe_streambuf(char const* systemCommand) :
    basic_iopipe_streambuf<Char, Traits>(systemCommand,false) 
    {
      if (m_buf) setp(m_buf, m_buf + s_nbuf);
    }


  template<typename CharT, typename Traits>
    basic_opipe_streambuf<CharT, Traits>::
    ~basic_opipe_streambuf()
    {
	sync();
    }

  template<typename CharT, typename Traits>
    typename basic_opipe_streambuf<CharT, Traits>::base_class::int_type
    basic_opipe_streambuf<CharT, Traits>::
    overflow(typename base_class::int_type c)
    {
	if (pbase() && pbase() != pptr()) {
	    if (std::fwrite(pbase(), sizeof(CharT), pptr() - pbase(), (FILE*)m_fp)
		!= (size_t)(pptr() - pbase()))
		return Traits::eof();
	    setp(pbase(), epptr());
	}
	if (!Traits::eq_int_type(c, Traits::eof())) {
	    if (std::fputc(c, (FILE*)m_fp) == EOF)
		return Traits::eof();
	    return c;
	}
	else
	    return Traits::not_eof(c);
    }

  template<typename CharT, typename Traits>
    int
    basic_opipe_streambuf<CharT, Traits>::
    sync()
    {
	overflow();
	std::fflush((FILE*)m_fp);
	return 0;
    }


} // namespace PRODDL
 
