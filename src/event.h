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

  double start_height() const { return this->start_height_; };
  double end_height()   const { return this->end_height_; };
  double length()       const { return (end_height() - start_height()); };

  std::vector<Node*> contemporaries() const { return this->contemporaries_; };

  Node* getRandomContemporary() const;
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
  EventIterator(Forest *forest, Node const* start_node);
  ~EventIterator() {};

  void next();
  bool good() const;

  Event operator*() const { return current_event_; }
  Event operator++() { next(); return current_event_; }
  Event operator++(int) { Event tmp = current_event_;
                          next();
                          return tmp; }
  
#ifdef UNITTEST
  friend class TestEvent;
#endif

 private:
  Forest* forest_;
  Event current_event_;
  std::vector<Node*> contemporaries_;
  NodeIterator twig_iterator_;
  bool good_;
};

#endif
