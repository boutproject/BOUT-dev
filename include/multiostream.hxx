#ifndef __MULTIOSTREAM_H__
#define __MULTIOSTREAM_H__

#include <streambuf>
#include <vector>
#include <algorithm>

/// Template class to split streams
/*!
  from http://accu.org/index.php/journals/260
*/
template<typename char_type, typename traits = std::char_traits<char_type> >
class multioutbuf : public std::basic_streambuf<char_type, traits> {
  private:
  typedef std::vector<std::basic_ostream<char_type, traits>* >
  stream_container;
  typedef typename stream_container::iterator
  iterator;
  stream_container streams_;
  
  public:
  void add(std::basic_ostream<char_type,
	   traits>& str) {
    iterator pos = std::find(streams_.begin(),
                        streams_.end(), &str);

    // Already been added
    if(pos != streams_.end()) {
      return;
    }
    
    streams_.push_back(&str);

  }

  void remove(std::basic_ostream<char_type,
                               traits>& str) {
    iterator pos = std::find(streams_.begin(),
                        streams_.end(), &str);

    if(pos != streams_.end()) {
      streams_.erase(pos);
    }
    
  }
protected:
  virtual std::streamsize xsputn(
                    const char_type* sequence,
                    std::streamsize num) {
    iterator current = streams_.begin();
    iterator end = streams_.end();

    for(; current != end; ++current) {
      (*current)->write(sequence, num);
      (*current)->flush();
    }

    return num;
  }

  virtual int overflow(int c) {
    iterator current = streams_.begin();
    iterator end = streams_.end();

    for(; current != end; ++current) {
      (*current)->put(static_cast<char>(c));
      (*current)->flush();
    }

    return c;
  }
};

template<typename char_type, typename traits>
class multioutbuf_init {
private:
  multioutbuf<char_type, traits> buf_;

public:
  multioutbuf<char_type, traits>* buf() {
    return &buf_;
  }
};

template<typename char_type, typename traits
               = std::char_traits<char_type> >
class multiostream : private
       multioutbuf_init<char_type, traits>, 
      public
       std::basic_ostream<char_type, traits> {
 private:
  typedef ::multioutbuf_init<char_type, traits> multioutbuf_init;

 public:
  multiostream() : multioutbuf_init(), std::basic_ostream<char_type,
 traits>(multioutbuf_init::buf()) {}

  void add(std::basic_ostream<char_type,
                              traits>& str) {
    multioutbuf_init::buf()->add(str);
  }

  void remove(std::basic_ostream<char_type,
                              traits>& str) {
    multioutbuf_init::buf()->remove(str);
  }
};

typedef multiostream<char> cmultiostream;
typedef multiostream<wchar_t> wmultiostream;

#endif // __MULTIOSTREAM_H__

