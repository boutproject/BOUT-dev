
class CASHKARPScheme;

#ifndef __CASHKARP_SCHEME_H__
#define __CASHKARP_SCHEME_H__

#include <bout/rkscheme.hxx>
#include <utils.hxx>

class CASHKARPScheme : public RKScheme {
public:
  CASHKARPScheme(Options* options);
};

#endif // __CASHKARP_SCHEME_H__
