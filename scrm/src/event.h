#ifndef scrm_src_event
#define scrm_src_event

#include "forest.h"
class Forest;

struct Event {
  const double start_height;
  const double end_height; 
  const std::vector<Node*> contemporaries;
};

class EventIterator {
 public:
  EventIterator();
  EventIterator(Forest *forest, const double &start_height);
  ~EventIterator();

  struct Event next();

 private:
  Forest* forest_;
  std::vector<Node*> contemporaries_;
  double current_height_;
  std::vector<Node*>::iterator twig_iterator_;
};

#endif
