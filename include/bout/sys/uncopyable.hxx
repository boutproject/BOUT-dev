// From Scott Meyers' "Effective C++, third edition"

#ifndef BOUT_UNCOPYABLE_H
#define BOUT_UNCOPYABLE_H

/// Inherit from this class (private) to prevent copying
class Uncopyable {
protected:
  Uncopyable() = default;
  ~Uncopyable() = default;

public:
  Uncopyable(const Uncopyable&) = delete;
  Uncopyable& operator=(const Uncopyable&) = delete;
};

#endif // BOUT_UNCOPYABLE_H
