#include "seg.h"
#include "tree_point.h"


seg_data_container::seg_data_container(scrm::param user_para){
	this->seg_bool = user_para.seg_bool;
	this->FILENAME = user_para.treefile; 
	this->total_seq_length = user_para.nsites;
	this->remaining_max_num_mut = user_para.total_mut;
	this->nsam = user_para.nsam;
	this->numseg = 0 ;
}


void seg_data_container::append_new_seg_data(Forest *forest){
	if (this->seg_bool){
		seg_data* seg_data_ptr=new seg_data(forest, total_seq_length,remaining_max_num_mut);
		seg_datas.push_back(seg_data_ptr);
		remaining_max_num_mut = remaining_max_num_mut - seg_datas.back()->positions.size();
		numseg = numseg + seg_datas.back()->positions.size();
	}
}


void seg_data_container::print_to_file(){
	if (this->seg_bool){
	  std::ofstream tree_file;
	  tree_file.open (FILENAME.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
	  tree_file << "segsites: "<< numseg <<"\n";
	  if (numseg>0){
		  tree_file << "Positions: ";
		  for (size_t i=0; i< seg_datas.size();i++){
			  for (size_t j=0; j< seg_datas[i]->positions.size();j++){
				  tree_file << seg_datas[i]->positions[j] / total_seq_length << " ";
			  }
		  }
		  tree_file << "\n";
	
		  for (size_t k=0; k< nsam; k++){
			  for (size_t i=0; i< seg_datas.size();i++){
				for (size_t j=0; j< seg_datas[i]->haplotypes.size();j++){
				  tree_file << seg_datas[i]->haplotypes[j][k];
				}
			  }
			  tree_file <<"\n";
		  }
		  tree_file <<"\n";
	  }
	  tree_file.close();	
	}
}


seg_data_container::~seg_data_container(){
  for (size_t i=0; i< seg_datas.size();i++){
	  delete seg_datas[i];
  }

}

seg_data::seg_data(
Forest * forest,
double max_length,
int max_num_mut /*! if max_num_mut < 0, then it is infinity */
){
	double position_at = forest->current_base();
	int remaining_num_mut = max_num_mut;
	position_at = position_at + forest->random_generator()->sampleExpo(forest->local_tree_length() * forest->writable_model() ->mutation_rate() );
	double max_base = min(forest->next_base(),max_length);
	//forest->printTree_cout();
	while ( (remaining_num_mut > 0 || max_num_mut < 0 ) && position_at < max_base){
		TreePoint mut_point = forest->samplePoint();
		haplotypes.push_back(find_haplotypes(mut_point.base_node(), forest->writable_model()->sample_size()));
		positions.push_back(position_at);
		position_at = position_at + forest->random_generator()->sampleExpo(forest->local_tree_length() * forest->writable_model() ->mutation_rate() );
		//cout<<mut_point.base_node()<<" ";	
		remaining_num_mut--;
	}	//cout<<endl;
}

std::valarray <int> find_haplotypes(Node *node, int nsam){
	std::valarray <int> haplotype(nsam);
	traversal(node, haplotype);
	return haplotype;
}


void traversal(Node *node, std::valarray <int> &haplotype){
	if (node->first_child() == NULL && ((node->label())>0)){
		haplotype[node->label()-1]=1;
	}
	else if (node->first_child()->local() && node->second_child()->local()){
		Node *left = tracking_local_node(node->first_child());
		traversal(left, haplotype);
		Node *right = tracking_local_node(node->second_child());
		traversal(right, haplotype);
	}
	else if (!node->first_child()->local() ){
		traversal(node->second_child(), haplotype);
	}
	else{
		traversal(node->first_child(), haplotype);
	}
	
}




//void Forest::find_descndnt(){
	//int nsam = this->writable_model()->sample_size();
	//for(size_t i = 0; i < this->getNodes()->size(); ++i) {
		////std::valarray<int> tip_descndnt(nsam);
		//std::vector<int> tip_descndnt (nsam,0); 
		//this->nodes()->get_copy(i)->descndnt = tip_descndnt;
	//}
	//for(size_t i = 0; i < nsam; ++i) {
		//Node * parent = this->nodes()->get_copy(i);//->parent();
		//unsigned int label = parent->label();	
		//parent->descndnt[label-1]=1;
		//while ( !parent->is_root()){
			//parent = parent->parent();
			//parent->descndnt[label-1]=1;
		////while ( !parent->is_root() && parent->local() ){
			////if (parent->descndnt.size() == nsam){
				////parent->descndnt = parent->descndnt + tip_descndnt;
			////}
			////else{
				////parent->descndnt = tip_descndnt;
			////}
			////parent = parent->parent();
		//}
	//}
	
//}






//void Forest::exp_mut_num(int total_mut){
	//vector <double> cum_bl;
	//double total_brchlen=0;
	//for(size_t i = 0; i < this->getNodes()->size(); ++i) {
		//if (this->getNodes()->get(i)->local()){
			//total_brchlen = total_brchlen + this->getNodes()->get(i)->height_above();
			//cum_bl.push_back(total_brchlen);
		//}
	//}

	//for (int mut_i=0;mut_i<total_mut;mut_i++){
		//unsigned int brch_index=0;
		//double u=random_generator_->sample()*total_brchlen;
		//while (u>cum_bl[brch_index]){
			 //brch_index =brch_index + 1;	
		//}
		//this->nodes()->get_copy(brch_index)->set_mut_num(this->nodes()->get_copy(brch_index)->mut_num()+1);
	//}
//}

//void Forest::seg_data(string treefile, int in_total_mut){
	//int total_mut;
	//double total_bl= this->local_root()->length_below();

	//if (in_total_mut!=0){
		//total_mut = in_total_mut;
		////give the probability ...
		//}
	//else{
		//total_mut=random_generator_->samplePoisson(total_bl * this->writable_model()->mutation_rate());
	//}
	
	
	//this->find_descndnt();
	//this->exp_mut_num(total_mut);
	//std::ofstream tree_file;
	//tree_file.open (treefile.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
	//tree_file<<"segsites: "<<total_mut<<"\n";
	//int nsam = this->writable_model()->sample_size();

	//for (unsigned int tip_i=0;tip_i<nsam;tip_i++){
		////tree_file<<mt_tree.tip_name[tip_i]<<" ";
		////cout<<mt_tree.tip_name[tip_i]<<" ";
		//for (unsigned int node_i=0;node_i<this->getNodes()->size();node_i++){
			//if (this->getNodes()->get(node_i)->mut_num()>0){
				//for (int num_repeat=0;num_repeat<this->getNodes()->get(node_i)->mut_num();num_repeat++ ){				
					////site_data_file<<mt_tree.descndnt2[node_i][tip_i] <<" ";
					//tree_file<<this->getNodes()->get(node_i)->descndnt[tip_i] ;
				//}
			//}
		//}
		//tree_file<<"\n";
	//}

	//tree_file.close();	
//}

