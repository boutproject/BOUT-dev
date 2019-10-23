
class Timer;

#ifndef __TIMER_H__
#define __TIMER_H__

#include <map>
#include <string>

/*!
 * Timing class for performance benchmarking and diagnosis
 *
 * To record the time spent in a particular function, create a Timer object
 * when you wish to start timing
 * 
 *     void someFunction() {
 *       Timer timer("test"); // Starts timer
 *       
 *     } // Timer stops when goes out of scope
 *
 * Each time this function is called, the total time spent in someFunction
 * will be accumulated. To get the total time spent use getTime()
 * 
 *     Timer::getTime("test"); // Returns time in seconds as double
 *
 * To reset the timer, use resetTime
 * 
 *     Timer::resetTime("test"); // Timer reset to zero, returning time as double
 */
class Timer {
public:
  /*!
   * Create a timer. This constructor is equivalent to Timer("")
   */
  Timer();
  
  /*!
   * Create a timer, continuing from last time if the same label
   * has already been used
   */
  Timer(const std::string &label);
  
  /*!
   * Stop the timer
   */
  ~Timer();
  
  /*!
   * Get the time in seconds for time particular Timer object
   *
   *     Timer timer("test");
   *     // Some calculation
   *     output << timer.getTime();
   *     // timer still counting
   */
  double getTime();
  
  /*!
   * Get the time in seconds, reset timer to zero
   */
  double resetTime();
  
  /*!
   * The total time in seconds 
   */
  static double getTime(const std::string &label);
  
  /*!
   * The total time in seconds, resets the timer to zero
   */
  static double resetTime(const std::string &label);
  
  /*!
   * Clears all timers, freeing memory
   */
  static void cleanup();
  
private:
  /// Structure to contain timing information
  struct timer_info {
    double time;    ///< Total time
    bool running;   ///< Is the timer currently running?
    double started; ///< Start time
    unsigned int counter; ///< Number of Timer objects associated with this
                          ///< timer_info
  };
  
  static std::map<std::string, timer_info*> info;
  
  static timer_info *getInfo(const std::string &label);
  
  timer_info* timing;
};

#endif // __TIMER_H__
