#ifndef scrm_node
#define scrm_node

class Node
{
  public:                       
   Node();                      
   ~Node();                     
                         
   //Getters & Setters
   int get_height() { return this->height; }
   void set_height(int height) { this->height = height; }

   Node *get_parent() { return this->parent; }
   void set_parent(Node *parent) { this->parent = parent; }

   Node *get_higher_child() { return this->higher_child; }
   void set_higher_child(Node *higher_child) { this->higher_child = higher_child; }

   Node *get_lower_child() { return this->lower_child; }
   void set_lower_child(Node *lower_child) { this->lower_child = lower_child; }

  private:
   int height;  
   Node *parent;
   Node *higher_child;
   Node *lower_child;
};

#endif
