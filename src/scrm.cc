/*
 * scrm is an implementation of the Sequential-Coalescent-with-Recombination Model.
 * 
 * Copyright (C) 2013, 2014 Paul R. Staab, Sha (Joe) Zhu and Gerton Lunter
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

#include <iostream>
#include <ctime>

#include "param.h"
#include "forest.h"
#include "random/random_generator.h"
#include "random/mersenne_twister.h"

#include <string>
void Node::extract_bl_and_label ( std::string::iterator in_it ){
    std::string::iterator bl_start = in_it;
    for (; (*(bl_start-1)) != ':'; --bl_start ){ }
    double bl = strtod( &(*bl_start), NULL );
    this->set_bl ( bl );

    std::string::iterator label_start = (bl_start-2);
    while ( (*(label_start)) != ',' &&  (*(label_start)) != '(' &&  (*(label_start)) != ')' ){
        --label_start;
    }
    
    this->set_label ( ( (*(label_start)) == ')' ? 0 : strtol ( &(*(label_start+1)), NULL , 10)) );
    //if ( this->in_sample() ) this->set_height ( 0 );
        //if ( this->in_sample() ) this->set_height ( this->bl()   );
}

Node* Forest::readNewickNode( std::string &in_str, std::string::iterator &it, size_t parenthesis_balance, Node* const parent ){
    Node * node = new Node( 0, 0.0 );
    node->set_parent ( parent );

    this->nodes()->push_back( node );
    //std::string::iterator it = current_it;
    
    std::cout << "Node " << node <<" starts from [ " << std::endl;//\""<< *it <<"\" " ;
    for (  ; it != in_str.end(); ++it){
        std::cout << "\""<<(*it) <<"\"" ;
        if      ( (*it) == '(' )  { // Start of a internal node, extract a new node
            parenthesis_balance++;
            Node* child_1 = this->readNewickNode ( in_str, it = (it+1),parenthesis_balance, node );
            node->set_first_child ( child_1 );
            if ( node->first_child() != NULL){ node->set_height ( node->first_child()->height() + 40000*node->first_child()->bl()  ) ;}
        }
        else if ( (*(it+1)) == ',' ){ //
            node->extract_bl_and_label(it);            
            std::cout << " "<< parent<<" has first child node "<<node<<" with branch length " << node->bl() <<", and with the label "<< node->label()<<" height "<< node->height()<<" ]  Node " << node << " closed return " << *it  << std::endl;
            return node;
            //return it;
        }
        else if ( (*(it)) == ',' ){ //             
            Node* child_2 = this->readNewickNode ( in_str, it=(it + 1), parenthesis_balance, node );
            ////node->set_second_child ( this->nodes()->first() );
            //std::cout <<" Node " << node << " closed" << std::endl;
            //if ( parent != NULL ) 
            node->set_second_child ( child_2 );
            
            //return it;
        }
        else if ( (*(it+1)) == ')'  ){
            // Before return, extract the branch length for the second node
            node->extract_bl_and_label(it);
            std::cout << " " << parent << " has second child node " << node << " with branch length " << node->bl() <<", and with the label "<< node->label()<<" height "<< node->height()<<" ] Node " << node << " closed return "<< *it << std::endl;
            //return (it);
            return node;
        }
        //else if ( (*(it)) == ')' ){
            //parenthesis_balance--;    
            //if ( parenthesis_balance == 0 ){
                
                ////return it;    
            //}
        //}
        //else if ( (*(it)) == ')'  && parenthesis_balance == 0 ){
            
            //return (it);
             ////Before return, extract the branch length for the second node
            ////node->extract_bl_and_label(it);
            ////std::cout << " first child node finishes at "<< (*it) << ", with branch length " << node->bl() <<", and with the label "<< node->label();
            ////if ( parent != NULL ) parent->set_second_child ( node );
            ////return (it);
        //}
		else if ( (*(it)) == ';' ) {
            std::cout <<" Node " << node << " closed" << std::endl;            
            //return (it-1);
            return node;
        }
        else {
                continue;
        }
        
    }
    
}

void Forest::readNewick( std::string &in_str ){
  this->set_current_base(0.0);
  this->segment_count_ = 1;

std::string::iterator it = in_str.begin();
  (void)this->readNewickNode( in_str, it );
  this->set_local_root( this->nodes()->first() );
  this->set_primary_root(this->nodes()->first() );

    std::cout<< std::endl<<"there are "<< this->nodes()->size()<<std::endl;
    
    (void)this->nodes()->sorted();
assert(this->printNodes());
    //assert(this->printTree());

    this->sampleNextBase();
    this->calcSegmentSumStats();
}


#ifndef UNITTEST
int main(int argc, char *argv[]){
  try {
      Param user_para(argc, argv);
    Model model;
    user_para.parse(model);
    std::string tre_str ( "((1:1,2:1):3,3:4);" );
    //std::string tre_str ( "((11:1.1,22:22.22):6,(33:3.33,44:4.44):5.55);" );
std::cout << tre_str << std::endl;
    //std::string tre_str ( "(6:6,(3:3,4:4):5);" );

    MersenneTwister rg = MersenneTwister(1);

    Forest forest = Forest(&model, &rg);

    std::ostream *output = &std::cout;

    forest.readNewick ( tre_str );
    
    forest.printLocusSumStats(*output);
    rg.clearFastFunc();
    return EXIT_SUCCESS;
    //assert(forest.printTree());
    //// Organize output
    //std::ostream *output = &std::cout;

    //// Parse command line arguments
    //Param user_para(argc, argv);
    //Model model;
    //user_para.parse(model);

    //// Print help if user asked for it
    //if (user_para.help()) {
      //user_para.printHelp(*output); 
      //return EXIT_SUCCESS;
    //}
    //if (user_para.version()) {


      //*output << "scrm " << VERSION << std::endl;
      //return EXIT_SUCCESS;
    //}

    //MersenneTwister rg = MersenneTwister(user_para.random_seed());
    //*output << user_para << std::endl;
    //*output << rg.seed() << std::endl;

    //// Create the forest
    //Forest forest = Forest(&model, &rg);

    //// Loop over the independent loci/chromosomes
    //for (size_t rep_i=0; rep_i < model.loci_number(); ++rep_i) {

      //// Mark the start of a new independent sample
      //*output << std::endl << "//" << std::endl;

      //// Now set up the ARG, and sample the initial tree
      //forest.buildInitialTree();

      //while (forest.next_base() < model.loci_length()) { 
        //// Sample next genealogy
        //forest.sampleNextGenealogy();
      //}
      
      //forest.printLocusSumStats(*output);
      //forest.clear();
    //}

    // Clean-up and exit
    rg.clearFastFunc();
    return EXIT_SUCCESS;
  }
  catch (const std::exception &e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    std::cerr << "Try 'scrm --help' for more information." << std::endl;
    return EXIT_FAILURE;
  }
}
#endif
