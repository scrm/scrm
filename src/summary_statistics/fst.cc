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

#include "fst.h"

void Fst::calculate(const Forest &forest) {
    std::cout << " ok, compute fst "<< std::endl;
    std::cout << " sample_size_.size() "<< sample_size_.size() <<std::endl;
    // at least one mutation in the tree
    // computing within differences
    double Hw = 0, Hb = 0;
    // first add to Hw for population A
    size_t pop_shift = 0;
    for ( size_t pop = 0; pop < sample_size_.size(); pop++){
        size_t sample_size_tmp = sample_size_[pop];
        double Hw_tmp = 0;
        double Npairs = 0;
        for( size_t sample_i = 0 ; sample_i < (sample_size_tmp-1); sample_i++ ){
            for( size_t sample_j = (sample_i+1); sample_j < sample_size_tmp ; sample_j++ ){
                Npairs++;
                for ( size_t base = 0; base < seg_sites_->haplotypes_.size(); base++ ){
                    Hw_tmp += abs( seg_sites_->haplotypes_[base][sample_i+pop_shift] - seg_sites_->haplotypes_[base][sample_j+pop_shift] );
                }
            }
        }
        //cout << "Differnce in population " <<pop<<" "<< Hw_tmp<<endl;
        Hw += Hw_tmp ;
        pop_shift += sample_size_tmp;
    }

    
    double Npairs = 0;
    size_t pop_i_shift = 0;
    size_t pop_j_shift = sample_size_[0];
    for ( size_t pop_i = 0; pop_i < (sample_size_.size() - 1); pop_i++ ){
        
        size_t sample_size_pop_i = sample_size_[pop_i];
        
        for( size_t sample_i = 0 ; sample_i < sample_size_pop_i; sample_i++ ){
    
            for ( size_t pop_j = (pop_i+1); pop_j < sample_size_.size() ; pop_j++ ){
                size_t sample_size_pop_j = sample_size_[pop_j];
                for( size_t sample_j = 0 ; sample_j < sample_size_pop_j; sample_j++ ){
                    Npairs++;
                    for ( size_t base = 0; base < seg_sites_->haplotypes_.size(); base++ ){
                        Hb += abs( seg_sites_->haplotypes_[base][sample_i+pop_i_shift] - seg_sites_->haplotypes_[base][sample_j+pop_j_shift] ) ;
                    }
                }
        
            }
        }
        pop_i_shift += sample_size_pop_i;
        pop_j_shift += sample_size_[pop_i + 1];
    }
    
    //cout << "Hw = " << Hw <<endl;
    //cout << "Hb = " << Hb <<endl;
    //cout << "fst = " << 1.0 - (Hw/Hb) <<endl;
    fst_ = ( 1 -  (Hw/Hb) ); 
}

void Fst::printLocusOutput(std::ostream &output) const {
    std::cout << " ok, print fst "<< std::endl;
  output << "Fst: " << fst_ << std::endl;
}
