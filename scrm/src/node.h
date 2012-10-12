#ifndef scrm_node
#define scrm_node

#include <cstddef>
#include <cfloat>

class Node
{
  public:                       
   Node();
   Node(double height);
   Node(double height, bool active);
   ~Node();
                         
   //Getters & Setters
   double height() { return this->height_; }
   void set_height(const double &height) { this->height_ = height; }

   double parent_height();
   double height_above() { return this->parent_height() - this->height(); }

   bool active() { return this->active_; }
   void set_active(bool active) { this->active_ = active; }
   void activate() { this->set_active(true); }
   void deactivate() { this->set_active(false); }

   Node *parent() { return this->parent_; }
   void set_parent(Node *parent) { this->parent_ = parent; }; 

   Node *higher_child() { return this->higher_child_; }
   void set_higher_child(Node *higher_child) { this->higher_child_ = higher_child; }

   Node *lower_child() { return this->lower_child_; }
   void set_lower_child(Node *lower_child) { this->lower_child_ = lower_child; }

  private:
   void init(double heigh=-1, bool active=false);

   double height_;
   //double height_above_;
   bool active_;

   Node *parent_;
   Node *higher_child_;
   Node *lower_child_;
};

#endif
