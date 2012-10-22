#ifndef scrm_src_event
#define scrm_src_event

#include "forest.h"
class Forest;

class Event {
 public:
  Event();
  Event(Forest* forest, 
        double start_height, 
        double end_height, 
        std::vector<Node*> contemporaries);
  ~Event() { };

  double start_height() { return this->start_height_; };
  double end_height() { return this->end_height_; };
  std::vector<Node*> contemporaries() { return this->contemporaries_; };

  Node* getRandomContemporary();
  void removeFromContemporaries(Node* node);


 private:
  Forest* forest_;
  double start_height_;
  double end_height_; 
  std::vector<Node*> contemporaries_;
};

class EventIterator {
 public:
  EventIterator();
  EventIterator(Forest *forest, const double &start_height);
  ~EventIterator();

  Event next();
  
#ifdef UNITTEST
  friend class TestEvent;
#endif

 private:
  Forest* forest_;
  std::vector<Node*> contemporaries_;
  double start_height_;
  std::vector<Node*>::iterator twig_iterator_;
};

#endif
