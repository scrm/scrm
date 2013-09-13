/*
 * scrm is an implementation of the Sequential-Coalescent-with-Recombination Model.
 * 
 * Copyright (C) 2013 Paul R. Staab, Sha (Joe) Zhu and Gerton Lunter
 * 
 * This file is part of scrm.
 * 
 * scrm is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "seg.h"
#include "tree_point.h"


SegDataContainer::SegDataContainer(Param const* user_para, Forest const* forest){
  this->user_para_ = user_para;
  this->forest_ = forest;
  this->remaining_max_num_mut_ = forest->model().mutation_exact_number();
  this->numseg = 0;
}


void SegDataContainer::append_new_seg_data(Forest *forest){
  if (this->seg_bool()){
    SegDataBlock* seg_data_ptr=new SegDataBlock(forest, remaining_max_num_mut_);
    seg_datas_.push_back(seg_data_ptr);
    remaining_max_num_mut_ = remaining_max_num_mut_ - seg_datas_.back()->positions.size();
    numseg = numseg + seg_datas_.back()->positions.size();
  }
}


std::ostream& operator<< (std::ostream& stream, const SegDataContainer &sdc) {
  if ( sdc.seg_bool() ) {
    stream << "segsites: "<< sdc.numseg << std::endl;
    if (sdc.numseg>0){
      stream << "positions: ";
      for (size_t i=0; i< sdc.seg_datas_.size();i++){
        for (size_t j=0; j< sdc.seg_datas_[i]->positions.size();j++){
          stream << sdc.seg_datas_[i]->positions[j] / sdc.total_seq_length() << " ";
        }
      }
      stream << std::endl;

      for (size_t k=0; k< sdc.nsam(); k++){
        for (size_t i=0; i< sdc.seg_datas_.size();i++){
          for (size_t j=0; j< sdc.seg_datas_[i]->haplotypes.size();j++){
            stream << sdc.seg_datas_[i]->haplotypes[j][k];
          }
        }
        stream <<"\n";
      }
      //stream <<"\n";
    }
  }
  return stream;
}


SegDataContainer::~SegDataContainer(){
  for (size_t i=0; i< seg_datas_.size();i++){
    delete seg_datas_[i];
  }
}

SegDataBlock::SegDataBlock(
    Forest * forest /*! Current ARG */,
    //double seq_length /*! Sequence length */,
    int max_num_mut /*! Maximal number of mutation on the local genealogy. If max_num_mut < 0, then it is infinity */
    ){
  double position_at = forest->current_base();
  int remaining_num_mut = max_num_mut;
  position_at += forest->random_generator()->sampleExpo(forest->local_tree_length() * 
                 forest->writable_model()->mutation_rate() );
  double max_base = min(forest->next_base(),double(forest->model().loci_length()));
  //forest->printTree_cout();
  while ( (remaining_num_mut > 0 || max_num_mut < 0 ) && position_at < max_base){
    TreePoint mut_point = forest->samplePoint();
    haplotypes.push_back(find_haplotypes(mut_point.base_node(), forest->model().sample_size()));
    positions.push_back(position_at);
    position_at += forest->random_generator()->sampleExpo(forest->local_tree_length() * 
                 forest->writable_model()->mutation_rate() );
    //cout<<mut_point.base_node()<<" ";	
    remaining_num_mut--;
  }	//cout<<endl;
}

std::valarray <int> find_haplotypes(Node *node, int nsam){
  std::valarray <int> haplotype(nsam);
  traversal(node, haplotype);
  return haplotype;
}


void traversal(Node *node, std::valarray<int> &haplotype){
  //std::cout << "start " << node << std::endl;
  if (node->in_sample()){
    haplotype[node->label()-1]=1;
  }
  else if (node->second_child() == NULL) {
    traversal(node->first_child(), haplotype);
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
  else if (!node->second_child()->local()) {
    traversal(node->first_child(), haplotype);
  }
  else {
    assert( 0 );
  }
  //std::cout << "end" << std::endl;
}
