
class CASHKARPScheme;

#ifndef BOUT_CASHKARP_SCHEME_H
#define BOUT_CASHKARP_SCHEME_H

#include <bout/rkscheme.hxx>
#include <bout/utils.hxx>

class CASHKARPScheme : public RKScheme {
public:
  CASHKARPScheme(Options* options);
};

namespace {
RegisterRKScheme<CASHKARPScheme> registerrkschemecashkarp(RKSCHEME_CASHKARP);
}

#endif // BOUT_CASHKARP_SCHEME_H
