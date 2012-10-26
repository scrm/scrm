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
   ~Node();
                         
   //Getters & Setters
   double height() const { return this->height_; }
   void set_height(const double &height) { this->height_ = height; }

   double parent_height() const;
   double height_above() { return this->parent_height() - this->height(); }

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

   void change_child(Node* from, Node* to);
  
   bool uncoalesed() const { return this->uncoalesed_; }
   void set_uncoalesed(const bool &uncoalesed) { this->uncoalesed_ = uncoalesed; }

   bool is_fake() const; 
   bool is_root() const; 
   bool is_ultimate_root() const;

  private:
   void init(double heigh=-1, bool active=true);

   double height_;     //The total height of the node
   bool   active_;     //Indicates if the branch above is active, 
                       //i.e. on the local tree
   bool   uncoalesed_; //The nodes marks an unimplemented coalesence event

   //The tree structure
   Node *parent_;
   Node *higher_child_;
   Node *lower_child_;
};

#endif
