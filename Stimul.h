#ifndef __STIMUL__
#define __STIMUL__

#include "Cell.h"

#include <vector>

namespace ryr
{

class Stimul
{

public:

  Stimul(): start_(0), end_(10000000) {}
  std::vector<unsigned> cell_id_;
  double start_, end_, interval_, current_, duration_;
  std::vector<cell::MainCell*> cells_;

  void update_voltage(double time, double dt) {
    if(time > start_) {
      time -= int(time / interval_) * interval_;

      if(time < duration_) {
        for(std::vector<cell::MainCell*>::iterator it = cells_.begin();
            it != cells_.end(); ++it) {
          if(*it) {
            (*it)->v() += current_ * dt;
          }
        }
      }
    }
  }

};

}

#endif
