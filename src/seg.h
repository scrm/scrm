//seg.hpp
#include "param.h"
#include "forest.h"
class seg_data{
	public:
		seg_data();
		~seg_data(){};
		seg_data(Forest * forest,double max_length, int max_num_mut);
		vector<double> positions;
	
};


class seg_data_container{
	public:
		bool seg_bool;
		string FILENAME;
		int remaining_max_num_mut;
		double total_seq_length;	
		seg_data_container();
		seg_data_container(scrm::param user_para);
		~seg_data_container();
		void append_new_seg_data(Forest *forest);
		void print_to_file();
	private:
		vector <seg_data*> seg_datas;
};

//double unifRand();
//int poisson_rand_var(double lambda);
