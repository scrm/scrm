void ReadNewick::read ( std::string tre_str) 
    this->check_Parenthesis(tre_str);
    this->check_labeled( tre_str );
	// check & sign, this should be illigal for hybrid-Lambda, 
	

		size_t found_bl = net_str.find(':');
		for (size_t i_str_len=1;i_str_len<net_str.size();){
			if (net_str[i_str_len]=='e' && (net_str[i_str_len+1]=='-' || net_str[i_str_len+1]=='+')){
				i_str_len++;
			}
			else{
				if ( start_of_tax_name(net_str,i_str_len) ){
					size_t str_start_index = i_str_len;
					string label = extract_label(net_str,i_str_len);
					labels.push_back(label);

					string node_content;
					if ( net_str[str_start_index-1]==')' ){
						size_t rev_dummy_i = Parenthesis_balance_index_backwards( net_str, str_start_index-1 );
						size_t substr_len = str_start_index-rev_dummy_i;
						node_content = net_str.substr(rev_dummy_i, substr_len );
					}
					else {
						node_content=label;
					}
					i_str_len += label.size();

					node_contents.push_back(node_content);
					string brchlen;
					if ( found_bl != string::npos ){
						size_t found=min(min(net_str.find(",",i_str_len+1),net_str.find(")",i_str_len+1)),net_str.size());
						brchlen = net_str.substr(i_str_len+1,found-i_str_len-1);					
					}
					found_bl = net_str.find(":", found_bl+1);
					brchlens.push_back(brchlen);
				}
				else {
					i_str_len++;
				}
			}
		}
			
		//int label_counter = brchlens.size();
		for ( size_t new_i_label=0 ; new_i_label < brchlens.size(); new_i_label++ ){
			Node empty_node;
			NodeContainer.push_back(empty_node);
			NodeContainer[new_i_label].label = labels[new_i_label];
			NodeContainer[new_i_label].node_content = node_contents[new_i_label];
            NodeContainer[new_i_label].set_brchlen1( strtod(brchlens[new_i_label].c_str(), NULL) );
		}

		for ( size_t i = 1; i < NodeContainer.size()-1; i++ ){
			size_t j;
			for ( j = i+1; j < NodeContainer.size()-1; j++ ){
				if ( NodeContainer[j].label==NodeContainer[i].label ){
					if ( NodeContainer[j].node_content[0] == '(' ){
						NodeContainer[i].node_content = NodeContainer[j].node_content;
					}
					NodeContainer[i].set_brchlen2 ( NodeContainer[j].brchlen1() );
					break;
				}
			}
			if ( NodeContainer[j].label == NodeContainer[i].label ) NodeContainer.erase(NodeContainer.begin()+j);			
		}
		
    this->extract_tax_and_tip_names();

    this->connect_graph();
    this->NodeContainer.back().find_tip();
    this->NodeContainer.back().find_hybrid_descndnt();
    this->NodeContainer.back().CalculateRank();
    this->max_rank = NodeContainer.back().rank();	
    this->enumerate_internal_branch( this->NodeContainer.back() );
    this->init_descendant();
    this->init_node_clade();
    //this->rewrite_descendant();
    this->check_isNet();
    this->check_isUltrametric();
	//dout<<"Net constructed"<<endl;
}


void ReadNewick::connect_graph(){
    for ( size_t i = 0; i < NodeContainer.size(); i++ ){
        if ( NodeContainer[i].node_content[0] != '(' ) continue;
        
        char child_node1[NodeContainer[i].node_content.length()];
        for ( size_t i_content_len = 1; i_content_len < NodeContainer[i].node_content.length(); ){
            if (NodeContainer[i].node_content[i_content_len]=='(' ||  start_of_tax_name(NodeContainer[i].node_content,i_content_len) ){	
                size_t j_content_len = (NodeContainer[i].node_content[i_content_len] == '(') ? Parenthesis_balance_index_forwards( NodeContainer[i].node_content, i_content_len ) + 1:
                                                                                               i_content_len;
                int child1_node_content_i = 0;
                for ( ; j_content_len < NodeContainer[i].node_content.length(); j_content_len++){
                    child_node1[child1_node_content_i] = NodeContainer[i].node_content[j_content_len];
                    char stop = NodeContainer[i].node_content[j_content_len+1];
                    if ( stop == ',' || stop == ')' || stop == ':'){
                        child_node1[child1_node_content_i+1]='\0';
                        break;
                    }
                    child1_node_content_i++;
                }
                string child_node1_str = child_node1;		
                i_content_len = j_content_len + 2;
                for ( size_t j = 0; j < NodeContainer.size(); j++){
                    if (child_node1_str == NodeContainer[j].label) NodeContainer[i].add_child( &NodeContainer[j] );
                }
            }
            else { i_content_len++;}
        }	
    }    
}




/*! \brief Identify if its the start of the taxon name in a newick string, should be replaced by using (isalpha() || isdigit())  */
bool ReadNewick::start_of_tax_name( string in_str, size_t i ){
	//bool start_bool = false;
	//if ( (in_str[i]!='(' && in_str[i-1]=='(') || (in_str[i-1]==',' && in_str[i]!='(') || ( (in_str[i-1]==')') && ( in_str[i]!=')' || in_str[i]!=':' || in_str[i]!=',' || in_str[i]!=';' ) ) ) {
		//start_bool=true;	
	//}	
	//return 	start_bool;
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

