// From Scott Meyers' "Effective C++, third edition"

#ifndef __UNCOPYABLE_H__
#define __UNCOPYABLE_H__

/// Inherit from this class (private) to prevent copying
class Uncopyable {
protected:
  Uncopyable() {}
  ~Uncopyable() {}
private:
  Uncopyable(const Uncopyable&);
  Uncopyable& operator=(const Uncopyable&);
};

#endif // __UNCOPYABLE_H__
