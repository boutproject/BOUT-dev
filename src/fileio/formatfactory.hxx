
class FormatFactory;

#ifndef __FORMATFACTORY_H__
#define __FORMATFACTORY_H__

#include "bout/dataformat.hxx"
#include <bout/sys/uncopyable.hxx>

class FormatFactory : private Uncopyable {
public:
  /// Return a pointer to the only instance
  static FormatFactory* getInstance();

  std::unique_ptr<DataFormat> createDataFormat(const char *filename = nullptr,
                                               bool parallel = true);

private:
  static FormatFactory* instance; ///< The only instance of this class (Singleton)
  
  int matchString(const char *str, int n, const char **match);
};

#endif // __FORMATFACTORY_H__

