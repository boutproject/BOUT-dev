/*
 * Extend the Domain class, to allow extra data
 * to be attached to each domain
 */

class QuiltDomain;

#ifndef __QUILTDOMAIN_H__
#define __QUILTDOMAIN_H__

#include "../domain.hxx"

class QuiltDomain : public Domain {
public:
  QuiltDomain(int NX, int NY) : Domain(NX,NY) {}
  ~QuiltDomain();
  
  int proc;
  
  QuiltDomain *parent;
};

#endif // __QUILTDOMAIN_H__
