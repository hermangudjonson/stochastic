
/*
 * Timekeeper class taken from branching project for help recording times
 */

#ifndef TIMEKEEPER_H
#define TIMEKEEPER_H

#include <vector>

/*
 * TimeKeeper is a helper class for keeping track of time event logic 
 * in various listener classes
 */
class TimeKeeper {
public:
  /* two methods of construction:
   * -- keeping time by event
   * -- keeping times given 
   */
  TimeKeeper(std::vector<double> &t,double PRC=1e-15) : prc(PRC) {times_set = false;}
  TimeKeeper(std::vector<double> &t,std::vector<double> &ts,double PRC=1e-15) : prc(PRC)
  {
    t = ts;
    times_set = true;
  }
  
   // for use in init, time is starting time
  void init_times(std::vector<double> &t,double time)
  {
    // if event recording, need to add starting time
    if (!times_set) {
      t.push_back(time);
    }
    
    // increment tindex until we reach first recording time
    while (t[tindex] < time && !AlmostEqual(t[tindex],time,prc)) {++tindex;}
  }

  // for new event, bool returned indicates if we need a new record item
  // add new time if needed as well
  bool new_entry(std::vector<double> &t,double time)
  {
    if (!times_set) {
      // update by event
      if (!AlmostEqual(time,t.back(),prc)){
	t.push_back(time); // add new time entry
	return true;
      }
    }
    return false;
  }

  // step tindex if record time is less than or equal to next event time
  // if true is returned, tindex has been incremented and record items should be updated
  // allowing tindex to go past time vector
  bool step_time(std::vector<double> &t,double time)
  {
    if (tindex < int(t.size())) {
      if (t[tindex] < time && !AlmostEqual(t[tindex],time,prc)) {
	++tindex;
	return true;
      }
    }
    return false;
  }

  // condition for valid tindex
  bool in_range(std::vector<double> &t) {return tindex < int(t.size());}
  
  // condition for recording items from new event
  bool record(std::vector<double> &t,double time)
  {
    if (tindex < int(t.size()))
      return (t[tindex] > time || AlmostEqual(t[tindex],time,prc));
    else
      return false;
  }
  int tindex = 0;
private:
  double prc;
  bool times_set;
  bool AlmostEqual(double a,double b,double EPS=1.0e-15){return std::abs(a - b) < EPS;}
};

#endif
