#include "forest.h"

void Forest::find_descndnt(){
	int nsam = this->writable_model()->sample_size();
	for(size_t i = 0; i < this->getNodes()->size(); ++i) {
		//std::valarray<int> tip_descndnt(nsam);
		std::vector<int> tip_descndnt (nsam,0); 
		this->nodes()->get_copy(i)->descndnt = tip_descndnt;
	}
	for(size_t i = 0; i < nsam; ++i) {
		Node * parent = this->nodes()->get_copy(i);//->parent();
		unsigned int label = parent->label();	
		parent->descndnt[label-1]=1;
		while ( !parent->is_root()){
			parent = parent->parent();
			parent->descndnt[label-1]=1;
		//while ( !parent->is_root() && parent->local() ){
			//if (parent->descndnt.size() == nsam){
				//parent->descndnt = parent->descndnt + tip_descndnt;
			//}
			//else{
				//parent->descndnt = tip_descndnt;
			//}
			//parent = parent->parent();
		}
	}
	
}






void Forest::exp_mut_num(int total_mut){
	vector <double> cum_bl;
	double total_brchlen=0;
	for(size_t i = 0; i < this->getNodes()->size(); ++i) {
		if (this->getNodes()->get(i)->local()){
			total_brchlen = total_brchlen + this->getNodes()->get(i)->height_above();
			cum_bl.push_back(total_brchlen);
		}
	}

	for (int mut_i=0;mut_i<total_mut;mut_i++){
		unsigned int brch_index=0;
		double u=random_generator_->sample()*total_brchlen;
		while (u>cum_bl[brch_index]){
			 brch_index =brch_index + 1;	
		}
		this->nodes()->get_copy(brch_index)->set_mut_num(this->nodes()->get_copy(brch_index)->mut_num()+1);
	}
}

void Forest::seg_data(string treefile, int in_total_mut){
	int total_mut;
	double total_bl= this->local_root()->length_below();

	if (in_total_mut!=0){
		total_mut = in_total_mut;
		//give the probability ...
		}
	else{
		total_mut=random_generator_->samplePoisson(total_bl * this->writable_model()->mutation_rate());
	}
	
	
	this->find_descndnt();
	this->exp_mut_num(total_mut);
	std::ofstream tree_file;
	tree_file.open (treefile.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
	tree_file<<"segsites: "<<total_mut<<"\n";
	int nsam = this->writable_model()->sample_size();

	for (unsigned int tip_i=0;tip_i<nsam;tip_i++){
		//tree_file<<mt_tree.tip_name[tip_i]<<" ";
		//cout<<mt_tree.tip_name[tip_i]<<" ";
		for (unsigned int node_i=0;node_i<this->getNodes()->size();node_i++){
			if (this->getNodes()->get(node_i)->mut_num()>0){
				for (int num_repeat=0;num_repeat<this->getNodes()->get(node_i)->mut_num();num_repeat++ ){				
					//site_data_file<<mt_tree.descndnt2[node_i][tip_i] <<" ";
					tree_file<<this->getNodes()->get(node_i)->descndnt[tip_i] ;
				}
			}
		}
		tree_file<<"\n";
	}

	tree_file.close();	
}
			
