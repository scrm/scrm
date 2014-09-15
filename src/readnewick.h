#include "forest.h"
#include <vector>
class ReadNewick {
    ReadNewick( std::string );
    ~ReadNewick(){}
    
    
    vector<string> labels;
    vector<string> node_contents;
    vector<string> brchlens;
    void build_forest ( Forest * forest );
    
    
    void read ( std::string tre_str);
    void connect_graph(){
bool start_of_tax_name( string in_str, size_t i );
size_t Parenthesis_balance_index_backwards( string &in_str, size_t i );
size_t Parenthesis_balance_index_forwards( string &in_str, size_t i );
void check_Parenthesis( string &in_str );
}
