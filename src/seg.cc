#include "seg.h"
#include "tree_point.h"


seg_data_container::seg_data_container(Param user_para, Forest const* forest){
	this->seg_bool = user_para.seg_bool;
	this->FILENAME = user_para.treefile; 
	this->total_seq_length = forest->model().loci_length();
	this->remaining_max_num_mut = user_para.total_mut;
	this->nsam = forest->model().sample_size(); 
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


/*
std::ostream& Forest::generateSegData(std::ostream &output, int total_mut) {
	double total_bl = this->local_root()->length_below();

	if (total_mut == 0) {
		total_mut = random_generator_->samplePoisson(total_bl * this->model().mutation_rate());
	}
	
	this->find_descndnt();
	this->exp_mut_num(total_mut);
	output << "segsites: " << total_mut << std::endl;
	int nsam = this->model().sample_size();

	for (size_t tip_i=0; tip_i<nsam; tip_i++){
		for (size_t node_i=0; node_i < this->getNodes()->size(); ++node_i){
			if (this->getNodes()->get(node_i)->mut_num() > 0){
				for (int num_repeat=0; num_repeat < this->getNodes()->get(node_i)->mut_num(); num_repeat++){				
					output << this->getNodes()->get(node_i)->descndnt[tip_i] ;
				}
			}
		}
		output << std::endl;
	}
  return output;
}
*/
			

