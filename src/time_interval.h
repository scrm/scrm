#ifndef scrm_src_time_interval
#define scrm_src_time_interval

#include "forest.h"

class Forest;

class TimeInterval {
 public:
  friend class TimeIntervalIterator;

#ifdef UNITTEST
  friend class TestTimeInterval;
#endif

  TimeInterval();
  TimeInterval(Forest* forest, 
        double start_height, 
        double end_height, 
        std::vector<Node*> *contemporaries);
  ~TimeInterval() { };

  double start_height() const { return this->start_height_; };
  double end_height()   const { return this->end_height_; };
  double length()       const { return (end_height() - start_height()); };

  size_t numberOfContemporaries() const { return contemporaries_->size(); }
  Node* getRandomContemporary() const;

 private:
  bool checkContemporaries() const;
  const std::vector<Node*> contemporaries() const { return *contemporaries_; };
  
  Forest* forest_;
  double start_height_;
  double end_height_; 
  std::vector<Node*> *contemporaries_;
};



class TimeIntervalIterator {
 public:
  TimeIntervalIterator();
  TimeIntervalIterator(Forest *forest, Node const* start_node);
  ~TimeIntervalIterator() {};

  void next();
  bool good() const;

  TimeInterval operator*() const { return current_event_; }
  TimeInterval operator++() { next(); return current_event_; }
  TimeInterval operator++(int) { TimeInterval tmp = current_event_;
                          next();
                          return tmp; }

  // Splits the current interval in two parts by adding a node inside the interval;
  // Only affects the event after the next "next()" which than represents the
  // second part of the interval.
  void splitCurrentInterval(Node* splitting_node, Node* del_node = NULL) {
    this->inside_node_ = splitting_node;
    if (del_node != NULL) removeFromContemporaries(del_node);
    //current_event_.end_height_ = splitting_node->height();
  };

  void recalculateInterval();
  void removeFromContemporaries(Node* node);
  void addToContemporaries(Node* node);

#ifdef UNITTEST
  friend class TestTimeInterval;
#endif

 private:
  Forest* forest_;
  TimeInterval current_event_;
  std::vector<Node*> contemporaries_;
  NodeIterator node_iterator_;
  bool good_;

  Node* inside_node_;

};

#endif
