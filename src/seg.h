#include "param.h"
#include "forest.h"
//#include <vector> this is included in forest.h
#include <valarray>



class seg_data{
	public:
		seg_data();
		~seg_data(){};
		seg_data(Forest * forest,double max_length, int max_num_mut);
		vector <double> positions;
		vector < valarray<int> > haplotypes;	
};


class seg_data_container{
	public:
		bool seg_bool;
		string FILENAME;
		int remaining_max_num_mut;
		int numseg;
		double total_seq_length;	
		int nsam;
		seg_data_container();
		seg_data_container(scrm::param user_para);
		~seg_data_container();
		void append_new_seg_data(Forest *forest);
		void print_to_file();
	private:
		vector <seg_data*> seg_datas;
};


std::valarray <int> find_haplotypes(Node *node, int nsam);
void traversal(Node *node, std::valarray <int> &haplotype);
