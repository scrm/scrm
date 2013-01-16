#ifndef scrm_node
#define scrm_node

#include <cstddef>
#include <cfloat>
#include <stdexcept>
#include <iostream>
#include <cassert>


class Node
{
 public:                       

#ifdef UNITTEST
  friend class TestNode;
#endif

  Node();
  Node(double height);
  Node(double height, bool active);
  Node(double height, bool active, size_t last_update);
  Node(double height, bool active, size_t last_update, size_t samples_below);
  Node(double height, bool active, size_t last_update, size_t samples_below, double length_below);
  ~Node();

  //Getters & Setters
  double height() const { return this->height_; }
  void set_height(const double &height) { this->height_ = height; }

  double parent_height() const;
  double height_above() const { return this->parent_height() - this->height(); }

  bool active() const { return this->active_; }
  void set_active(bool active) { this->active_ = active; }
  void activate() { this->set_active(true); }
  void deactivate() { this->set_active(false); }

  Node *parent() const; 
  void set_parent(Node *parent) { this->parent_ = parent; }; 

  Node *higher_child() const { return this->higher_child_; }
  void set_higher_child(Node *higher_child) { this->higher_child_ = higher_child; }

  Node *lower_child() const { return this->lower_child_; }
  void set_lower_child(Node *lower_child) { this->lower_child_ = lower_child; }

  size_t last_update() const { if ( active() || is_fake() ) return 0; 
    return(last_update_); }
  void set_last_update(size_t position) { this->last_update_ = position; }

  size_t samples_below() const { return samples_below_; }
  void set_samples_below(size_t samples) { samples_below_ = samples; }

  double length_below() const { return length_below_; }
  void set_length_below(double length) { length_below_ = length; }

  void change_child(Node* from, Node* to);
  int  numberOfChildren() const;

  bool is_fake() const; 
  bool is_root() const; 
  bool is_ultimate_root() const;
  bool in_sample() const;

 private:
  void init(double heigh=-1, bool active=true, size_t last_update = 0, size_t samples_below=0, double length_below=0);

  double height_;        // The total height of the node
  bool   active_;        // Indicates if the branch above is active,
  // i.e. on the local tree
  size_t last_update_;   // The sequence position on which the branch above the node
  // was last checked for recombination events. Ignored for
  // active nodes, which are always up to date, and for fake
  // nodes, which don't recombine.
  size_t samples_below_; // the number of sampled nodes in the subtree below this node
  double length_below_;  // the total length of active branches in the subtree below this node

  //The tree structure
  Node *parent_;
  Node *higher_child_;
  Node *lower_child_;
};

//Inlining this functions slighly increases preformance
inline bool Node::is_fake() const {
  return ( this->height() == FLT_MAX ); 
}

inline bool Node::is_ultimate_root() const {
  return ( this->is_fake() && this->parent_ == NULL );
}


#endif
