#include "param.h"
#include "forest.h"
//#include <vector> this is included in forest.h
#include <valarray>



class seg_data{
 public:
  seg_data();
  ~seg_data(){};

  seg_data(Forest * forest, double max_length, int max_num_mut);
  vector <double> positions;
  vector < valarray<int> > haplotypes;	
};


class SegDataContainer {
 public:
  SegDataContainer();
  SegDataContainer(Param const* user_para, Forest const* forest);
  ~SegDataContainer();

  friend std::ostream& operator<<(std::ostream& stream, const SegDataContainer &sdc);

  bool seg_bool() const { return user_para_->seg_bool(); };
  int numseg;

  double total_seq_length() const { return forest_->model().loci_length(); };
  int nsam() const { return forest_->model().sample_size(); };

  void append_new_seg_data(Forest *forest);
  void print_to_file();

 private:
  vector <seg_data*> seg_datas_;

  Param const* user_para_;
  Forest const* forest_; 

  size_t remaining_max_num_mut;
};


std::valarray <int> find_haplotypes(Node *node, int nsam);
void traversal(Node *node, std::valarray <int> &haplotype);
