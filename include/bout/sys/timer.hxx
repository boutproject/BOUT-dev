
class Timer;

#ifndef __TIMER_H__
#define __TIMER_H__

#include <map>
#include <string>

class Timer {
public:
  Timer();
  Timer(const std::string &label);
  ~Timer();

  double getTime();
  double resetTime();
  
  static double getTime(const std::string &label);
  static double resetTime(const std::string &label);
  
  static void cleanup(); ///< Frees memory
  
private:
  /// Structure to contain timing information
  struct timer_info {
    double time;    ///< Total time
    bool running;   ///< Is the timer currently running?
    double started; ///< Start time
  };
  
  static std::map<std::string, timer_info*> info;
  
  static timer_info *getInfo(const std::string &label);
  
  timer_info* timing;
};

#endif // __TIMER_H__
