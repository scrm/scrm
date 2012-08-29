#ifndef scrm_node
#define scrm_node

class Node
{
  public:                       
   Node();
   Node(int height);
   ~Node();                     
                         
   //Getters & Setters
   int height() { return this->height_; }
   void set_height(int height) { this->height_ = height; }

   Node *parent() { return this->parent_; }
   void set_parent(Node *parent) { this->parent_ = parent; }

   Node *higher_child() { return this->higher_child_; }
   void set_higher_child(Node *higher_child) { this->higher_child_ = higher_child; }

   Node *lower_child() { return this->lower_child_; }
   void set_lower_child(Node *lower_child) { this->lower_child_ = lower_child; }

  private:
   int height_;  
   Node *parent_;
   Node *higher_child_;
   Node *lower_child_;
};

#endif
