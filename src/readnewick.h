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

#include "forest.h"

class NewickNode{
    NewickNode();
    ~NewickNode();
    string label;
    string node_content;
    double bl;
};

class ReadNewick {
    ReadNewick( std::string );
    ~ReadNewick(){}

    vector < NewickNode > newick_nodes_;
    void build_forest ( Forest * forest );    
    
    //void read ( std::string tre_str);
    void connect_graph(){
    bool start_of_tax_name( string in_str, size_t i );
    size_t Parenthesis_balance_index_backwards( string &in_str, size_t i );
    size_t Parenthesis_balance_index_forwards( string &in_str, size_t i );
    void check_Parenthesis( string &in_str );
}
