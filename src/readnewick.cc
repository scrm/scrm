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

void Forest::readInitialTree( std::string in_str, std::string::iterator current_it ){
    size_t num_b = 0;
	for ( size_t i = 0; i < in_str.size(); i++){
		if      (in_str[i] == '(') num_b++;
		else if (in_str[i] == ')') num_b--;
        else continue;
	}
	if ( num_b != 0 ) throw std::invalid_argument(in_str + "Parenthesis not balanced!" );
    //ReadNewick* basic_info = new ReadNewick (tree_str);
    //this->build_inintial_nodes ( basic_info );
    //this->connect_graph () ;
    //delete basic_info;
}



void Forest::build_inintial_nodes ( ReadNewick * basic_info ){
    for ( size_t i = 0 ; i < Tree_info->brchlens.size(); i++ ){
        bool is_tip = ( Tree_info->node_labels[i] == Tree_info->node_contents[i]);
        Node * node = new Node ( Tree_info->node_labels[i],
                                 Tree_info->node_contents[i],
                                 strtod( Tree_info->brchlens[i].c_str(), NULL),
                                 is_tip );
        this->nodes_.add ( node );
    }    
}

void Forest::connect_graph (){
    
}

ReadNewick::ReadNewick ( std::string &in_str ) {
    size_t found_bl = this->in_str.find(':');
    for ( size_t i_str_len = 1; i_str_len < this->in_str.size(); ){
        if ( this->in_str[i_str_len]=='e' && ( this->in_str[i_str_len+1]=='-' || this->in_str[i_str_len+1]=='+' ) ){
            i_str_len++;
        }
        else if ( start_of_tax_name( this->in_str, i_str_len) ){
            size_t str_start_index = i_str_len;
            string label = extract_label( this->in_str, i_str_len );
            //this->node_labels.push_back(label);
            
            string node_content = ( this->in_str[str_start_index-1]==')' ) ? extract_One_node_content ( this->in_str, str_start_index-1 )
                                                                            : label ;
            //node_contents.push_back(node_content);

            i_str_len += label.size();

            string brchlen;
            if ( found_bl != string::npos ){
                size_t found=min(min( this->in_str.find(",",i_str_len+1), this->in_str.find(")",i_str_len+1)), this->in_str.size());
                brchlen = this->in_str.substr(i_str_len+1,found-i_str_len-1);
            }
            found_bl = this->in_str.find(":", found_bl+1);
            //brchlens.push_back(brchlen);
            NewickNode new_node ();
            newick_nodes_.push_back ( new_node );
        }
        else {
            i_str_len++;
        }
    }
}

string GraphReader::extract_label( string &in_str, size_t i ){
    string label( in_str.substr ( i , end_of_label_or_bl ( in_str, i ) + 1-i ) );
    return label;
}


string GraphReader::extract_One_node_content( string &in_str, size_t back_parenthesis_index ){
    size_t front_parenthesis_index = Parenthesis_balance_index_backwards( in_str, back_parenthesis_index );
    return in_str.substr ( front_parenthesis_index, back_parenthesis_index - front_parenthesis_index + 1);
}


void ReadNewick::connect_graph(){
    for ( auto it = nodes_.iterator(); it.good(); ++it){
        dout << " node " << (*it) << ": " << (*it)->label <<"\t"<<(*it)->node_content<<endl;
        if ( (*it)->node_content[0] != '(' ) continue;
        
        char child_node1[(*it)->node_content.length()];
        for ( size_t i_content_len = 1; i_content_len < (*it)->node_content.length(); ){
            if ((*it)->node_content[i_content_len]=='(' ||  start_of_tax_name((*it)->node_content,i_content_len) ){    
                size_t j_content_len = ((*it)->node_content[i_content_len] == '(') ? Parenthesis_balance_index_forwards( (*it)->node_content, i_content_len ) + 1:
                                                                                               i_content_len;
                int child1_node_content_i = 0;
                for ( ; j_content_len < (*it)->node_content.length(); j_content_len++){
                    child_node1[child1_node_content_i] = (*it)->node_content[j_content_len];
                    char stop = (*it)->node_content[j_content_len+1];
                    if ( stop == ',' || stop == ')' || stop == ':'){
                        child_node1[child1_node_content_i+1]='\0';
                        break;
                    }
                    child1_node_content_i++;
                }
                string child_node1_str = child_node1;        
                i_content_len = j_content_len + 2;
                for ( size_t j = 0; j < nodes_.size(); j++){
                    if (child_node1_str == this->nodes_.at(j)->label) {
                        (*it)->add_child( this->nodes_.at(j) );
                        //dout << "node " << &this->nodes_[i] << " has child "<< &this->nodes_[j]<<endl;
                    }
                }
            }
            else { i_content_len++;}
        }    
    }    
}


/*! \brief Identify if its the start of the taxon name in a newick string, should be replaced by using (isalpha() || isdigit())  */
bool ReadNewick::start_of_tax_name( string &in_str, size_t i ){
	if      (  in_str[i-1] == '('  &&   in_str[i] != '(' ) return true;
    else if (  in_str[i-1] == ','  &&   in_str[i] != '(' ) return true; 
    else if ( (in_str[i-1] == ')') && ( in_str[i] != ')' || in_str[i]!=':' || in_str[i]!=',' || in_str[i]!=';' ) ) return true;
    else return false;
}


size_t ReadNewick::Parenthesis_balance_index_backwards( string &in_str, size_t i ){
	size_t j = i;
	int num_b = 0;
	for ( ; j > 0 ; j-- ){
		if      ( in_str[j] == '(' ) num_b--;
		else if ( in_str[j] == ')' ) num_b++;
		else continue;
		if ( num_b == 0 ) break;
	}
	return j;
}


size_t ReadNewick::Parenthesis_balance_index_forwards( string &in_str, size_t i ){
	size_t j = i;
	int num_b = 0;
	for ( ; j < in_str.size(); j++ ){
		if      ( in_str[j] == '(' ) num_b++;		
		else if ( in_str[j] == ')' ) num_b--;
        else continue;
		if ( num_b == 0 ) break;
	}
	return j;
}

size_t ReadNewick::end_of_label_or_bl( string &in_str, size_t i ){
    for ( size_t j = i; j < in_str.size(); j++){
        if      ( in_str[j+1] == ',' )    return j;
        else if ( in_str[j+1] == ')' )    return j;
        else if ( in_str[j+1] == ':' )    return j;
        else if ( in_str[j+1] == ';' )    return j;
        else continue;
    }
}


/*! \brief Checking Parenthesis of a (extended) Newick string */
void ReadNewick::check_Parenthesis( string &in_str ){
	int num_b = 0;
	for ( size_t i = 0; i < in_str.size(); i++){
		if      (in_str[i] == '(') num_b++;
		else if (in_str[i] == ')') num_b--;
        else continue;
	}
	if ( num_b != 0 ) throw std::invalid_argument(in_str + "Parenthesis not balanced!" );
}
