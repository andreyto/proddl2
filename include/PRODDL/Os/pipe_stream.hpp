//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_OS_PIPE_STREAM_H__
#define PRODDL_OS_PIPE_STREAM_H__

// Pipe input and output streams with corresponding stream buffers.
// Modelled after those from MORE library.

#include <ios>
#include <iostream>

namespace PRODDL {

  template<typename Char, typename Traits>
  class basic_iopipe_streambuf : public std::basic_streambuf<Char, Traits>
  {
  public:
    typedef Char char_type;
    typedef Traits traits_type;
    typedef typename Traits::int_type int_type;
    typedef typename Traits::pos_type pos_type;
    typedef typename Traits::off_type off_type;

    basic_iopipe_streambuf(const char* systemCommand,bool isInput);
    virtual ~basic_iopipe_streambuf();

    bool is_open() const { return m_fp != 0; }

  protected:
    static int const s_nbuf = 256;

    void* m_fp;

    char_type m_buf[s_nbuf];
  };

  template<typename Char, typename Traits>
  class basic_ipipe_streambuf : public basic_iopipe_streambuf<Char,Traits>
  {
  public:

    typedef basic_iopipe_streambuf<Char,Traits> base_class;
    
    basic_ipipe_streambuf(const char* systemCommand);

    virtual ~basic_ipipe_streambuf() 
    {};

  protected:
    virtual typename base_class::int_type underflow();

  };


  template< typename Char, typename Traits=std::char_traits<Char> >
  class basic_ipipe_stream : public std::basic_istream<Char, Traits>
  {
  public:
    typedef basic_ipipe_streambuf<Char, Traits> streambuf_type;
    basic_ipipe_stream(const char* systemCommand)
      : std::basic_istream<Char, Traits>(m_buf=new streambuf_type(systemCommand))
    {
      if (!is_open()) setstate(std::ios_base::failbit);
    }

    ~basic_ipipe_stream() { delete m_buf; }

    bool is_open() const { return m_buf->is_open(); }

  private:
    streambuf_type* m_buf;

  };

  typedef basic_ipipe_stream<char> ipipe_stream;


  template<typename Char, typename Traits>
  class basic_opipe_streambuf : public basic_iopipe_streambuf<Char,Traits>
  {

  public:

    typedef basic_iopipe_streambuf<Char,Traits> base_class;

    basic_opipe_streambuf(const char* systemCommand);

    virtual ~basic_opipe_streambuf();

  protected:

    virtual int sync();

    virtual typename base_class::int_type overflow(typename base_class::int_type size = traits_type::eof());

  };

  template< typename Char, typename Traits = std::char_traits<Char> >
    class basic_opipe_stream
        : public std::basic_ostream<Char, Traits>
    {
    public:
	typedef basic_opipe_streambuf<Char, Traits> streambuf_type;

	basic_opipe_stream(const char* systemCommand) : 
	  std::basic_ostream<Char, Traits>(m_buf=new streambuf_type(systemCommand))
	{
	    if (!is_open()) setstate(std::ios_base::failbit);
	}

	~basic_opipe_stream() { delete m_buf; }

	bool is_open() const { return m_buf->is_open(); }

      private:
	streambuf_type* m_buf;
    };

  typedef basic_opipe_stream<char> opipe_stream;

} // namespace PRODDL


#include "PRODDL/Os/pipe_stream.cpp"

#endif // PRODDL_OS_PIPE_STREAM_H__
