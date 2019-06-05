// From Scott Meyers' "Effective C++, third edition"

#ifndef __UNCOPYABLE_H__
#define __UNCOPYABLE_H__

/// Inherit from this class (private) to prevent copying
class Uncopyable {
protected:
  Uncopyable() = default;
  ~Uncopyable() = default;
  Uncopyable(const Uncopyable &) = delete;
  Uncopyable &operator=(const Uncopyable &) = delete;
};

#endif // __UNCOPYABLE_H__
