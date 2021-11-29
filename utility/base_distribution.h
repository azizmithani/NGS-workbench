/* 
 * File:   base_distribution.h
 * Author: Aziz Mithani
 *
 * Created on 08 May 2011, 15:12
 */

#ifndef BASE_DISTRIBUTION_H
#define	BASE_DISTRIBUTION_H

#include <map>
#include <string>
#include "dna.h"
#include "position_base_distribution.h"

using namespace std;

class base_distribution {
public:
    base_distribution(dna::SEQUENCE_LENGTH_LIST* sequence_lengths, bool allocate_space);
    virtual ~base_distribution();

    void initialise(bool allocate_space);
    void write(string out_file);
    position_base_distribution** get_base_distribution_list(string sequence_name);
    void set_base_distribution_list(string sequence_name, position_base_distribution** seq_base_distribution_ptr);
    static dna::SEQUENCE_LENGTH_LIST* read_sequence_lengths(const string& base_distribution_file);
    static base_distribution* read(const string& base_distribution_file);
    static base_distribution* read(const string& base_distribution_file, dna::SEQUENCE_LENGTH_LIST* sequence_lengths);

    void set_sequence_lengths(dna::SEQUENCE_LENGTH_LIST* sequence_lengths);
    dna::SEQUENCE_LENGTH_LIST* get_sequence_lengths() const;
 
    
private:
    
    dna::SEQUENCE_LENGTH_LIST* sequence_lengths;
    map<string, position_base_distribution::SEQ_BASE_DISTRIBUTION_LIST*> base_distribution_list;
};

#endif	/* BASE_DISTRIBUTION_H */

