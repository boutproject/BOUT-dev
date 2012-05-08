
class FormatFactory;

#ifndef __FORMATFACTORY_H__
#define __FORMATFACTORY_H__

#include "dataformat.hxx"

#include <uncopyable.hxx>

class FormatFactory : private Uncopyable {
public:
  /// Return a pointer to the only instance
  static FormatFactory* getInstance();
  
  DataFormat* createDataFormat(const char *filename = NULL);
private:
  static FormatFactory* instance; ///< The only instance of this class (Singleton)
  
  int matchString(const char *str, int n, const char **match);
};

// For backwards compatability
DataFormat* data_format(const char *filename = NULL);

#endif // __FORMATFACTORY_H__

