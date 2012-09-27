#ifndef scrm_node
#define scrm_node

#include <cstddef>
using namespace std;

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

   double height_above() { return this->height_above_; }
   void set_height_above(const double &height) { this->height_above_ = height; }
  

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
   double height_above_;
   bool active_;

   Node *parent_;
   Node *higher_child_;
   Node *lower_child_;
};

#endif
