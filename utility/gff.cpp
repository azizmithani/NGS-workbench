/* 
 * File:   GFFList.cpp
 * Author: aziz
 * 
 * Created on May 10, 2011, 11:39 AM
 */

#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include "dna.h"
#include "gff.h"
#include "gff_entry.h"
#include "utility.h"


const string gff::FEATURE_EXON = "exon";
const string gff::FEATURE_GENE = "gene";


gff::gff() {
}

gff::gff(const string& gff_file) {

    // read the gff file
    read_gff_file(gff_file);
}

gff::gff(const gff& old_gff) {
    GFF_LIST::const_iterator the_gff_list = old_gff.gff_list.begin();
    for (; the_gff_list != old_gff.gff_list.end(); the_gff_list++) {
        string sequence_name = the_gff_list->first;

        // allocate the space for GFF entries for this sequence
        SEQ_GFF_LIST* new_seq_gff_list_ptr = new vector<gff_entry*>;

        SEQ_GFF_LIST* seq_gff_list_ptr = the_gff_list->second;
        SEQ_GFF_LIST::iterator the_seq_gff_list = seq_gff_list_ptr->begin();
        for (; the_seq_gff_list != seq_gff_list_ptr->end(); the_seq_gff_list++) {
            // call the copy constructor
            gff_entry* gff_entry_ptr = new gff_entry(**the_seq_gff_list);

            new_seq_gff_list_ptr->push_back(gff_entry_ptr);

        }
        this->gff_list[sequence_name] = new_seq_gff_list_ptr;
    }
}

gff::~gff() {
    GFF_LIST::iterator the_gff_list = this->gff_list.begin();
    for (; the_gff_list != this->gff_list.end(); the_gff_list++) {

        // get the GFF List for this sequence
        SEQ_GFF_LIST* seq_gff_list_ptr = the_gff_list->second;
        // and the iterator
        SEQ_GFF_LIST::iterator the_seq_gff_list = seq_gff_list_ptr->begin();
        // delete GFF entries for all genes
        for (; the_seq_gff_list != seq_gff_list_ptr->end(); the_seq_gff_list++) {
            // delete the GFF entry
            delete *the_seq_gff_list;
        }

        delete seq_gff_list_ptr;
    }

}

gff::SEQ_GFF_LIST* gff::get_gff_list(const string& sequence_name) {
    if (sequence_name.empty()) {
        return NULL;
    } else {
        return this->gff_list[sequence_name];
    }
}

void gff::set_gff_list(const string& sequence_name, SEQ_GFF_LIST* seq_gff_list_ptr) {
    this->gff_list[sequence_name] = seq_gff_list_ptr;
}

void gff::read_gff_file(const string& gff_file) {
    // open the input stream
    ifstream ifs_gff(gff_file.c_str());
    string str_gff = "";
    string previous_sequence_name = "";
    SEQ_GFF_LIST* seq_gff_list_ptr;
    while (ifs_gff.good()) {
        // read the next line from the GFF file
        getline(ifs_gff, str_gff);

        if (str_gff.empty()) { // ignore empty lines
            continue;
        }
        if (str_gff[0] == '#') { // ignore header
            continue;
        }

        vector<string> entries;
        utility::str_split(str_gff, entries, "\t");

        gff_entry* gff_entry_ptr = new gff_entry(entries, str_gff);

        if (gff_entry_ptr->get_sequence_name().compare(previous_sequence_name) != 0) { // new reference is starting
            if (!previous_sequence_name.empty()) { // this is not the first reference
                this->gff_list[previous_sequence_name] = seq_gff_list_ptr;
            }

            // initialise the GFF list for the new sequence
            seq_gff_list_ptr = new vector<gff_entry*>;
            // save the sequence name
            previous_sequence_name = gff_entry_ptr->get_sequence_name();
        }

        seq_gff_list_ptr->push_back(gff_entry_ptr);

    } // end while ifs_gff is good
    this->gff_list[previous_sequence_name] = seq_gff_list_ptr;

}

gff_entry* gff::get_gff_entry(SEQ_GFF_LIST* seq_gff_list_ptr, int position) {

    SEQ_GFF_LIST::iterator the_seq_gff_list = seq_gff_list_ptr->begin();
    for (; the_seq_gff_list != seq_gff_list_ptr->end(); the_seq_gff_list++) {
        gff_entry* gff_entry_ptr = *the_seq_gff_list;
        if (position >= gff_entry_ptr->get_start() - OFFSET && position <= gff_entry_ptr->get_end() + OFFSET) {
            return *the_seq_gff_list;
        }
    }

    return NULL;
}

void gff::write(const string& out_file) {
    write(out_file, false);
}

void gff::write(const string& out_file, bool only_valid_entries) {

    ofstream ofs_gff(out_file.c_str());
    ofs_gff << "##gff-version 3" << endl;

    GFF_LIST::iterator the_gff_list = this->gff_list.begin();
    for (; the_gff_list != this->gff_list.end(); the_gff_list++) {
        string sequence_name = the_gff_list->first;
        SEQ_GFF_LIST* seq_gff_list_ptr = the_gff_list->second;
        SEQ_GFF_LIST::iterator the_seq_gff_list = seq_gff_list_ptr->begin();
        for (; the_seq_gff_list != seq_gff_list_ptr->end(); the_seq_gff_list++) {
            gff_entry* gff_entry_ptr = *the_seq_gff_list;

            if (only_valid_entries && !gff_entry_ptr->valid()) {
                continue;
            }
            if (!gff_entry_ptr->modified()) {
                ofs_gff << gff_entry_ptr->get_gff_entry() << endl;
            } else {
                vector<string> entries;
                utility::str_split(gff_entry_ptr->get_gff_entry(), entries, "\t");

                ofs_gff << sequence_name << "\t" << entries[1] << "\t" << entries[2] << "\t" << gff_entry_ptr->get_start() << "\t" << gff_entry_ptr->get_end() << "\t" << entries[5] << "\t" << entries[6] << "\t" << entries[7] << "\t" << entries[8] << endl;
            }

        }
    }


    ofs_gff.close();
}

void gff::reset() {

    GFF_LIST::iterator the_gff_list = this->gff_list.begin();
    for (; the_gff_list != this->gff_list.end(); the_gff_list++) {
        string sequence_name = the_gff_list->first;
        SEQ_GFF_LIST* seq_gff_list_ptr = the_gff_list->second;
        SEQ_GFF_LIST::iterator the_seq_gff_list = seq_gff_list_ptr->begin();
        for (; the_seq_gff_list != seq_gff_list_ptr->end(); the_seq_gff_list++) {
            gff_entry* gff_entry_ptr = *the_seq_gff_list;

            gff_entry_ptr->set_valid(true);
        }
    }

}

string gff::create_gff_entry(const string& fasta_name, int sequence_start, int sequence_end, const string& sequence_name, const string& fasta_tag) {
    return fasta_name + "\t.\tgene\t" + utility::to_string(sequence_start) + "\t" + utility::to_string(sequence_end) + "\t.\t+\t.\tID=" + sequence_name + ";fasta_tag=" + fasta_tag;
}

void gff::filter_genes_with_problematic_snps_in_the_table(const string& snp_file, const string& out_file) {

    // open the file streams
    ifstream ifs_snp(snp_file.c_str());

    // go through the list of snp file to save snps each reference (chromosome)
    string str_snp = "";
    string sequence_name = "";
    string current_sequence_name = ""; // reference (chromosome) being currently read
    gff::SEQ_GFF_LIST* seq_gff_list_ptr = NULL;
    while (ifs_snp.good()) {

        // read the snp 
        getline(ifs_snp, str_snp);

        if (str_snp.empty()) { // ignore empty lines
            continue;
        }

        vector<string> entries;
        utility::str_split(str_snp, entries, "\t");

        bool is_valid = true;       
        if (entries.size() < 4) { 
            is_valid = false;
        } else {
            for (int i = 3; i < entries.size(); i++) {
                if (entries[i].empty()) {
                    is_valid = false;
                    break;
                } else if (!dna::is_nucleotide(entries[i])) {
                    is_valid = false;
                    break;                    
                }
            }
        }
        
        if (!is_valid) {
            // invalid entry

            // get the sequence name
            sequence_name = entries[0];

            if (sequence_name.compare(current_sequence_name) != 0) { // new reference
                seq_gff_list_ptr = this->gff_list[sequence_name];
                current_sequence_name = sequence_name;
            }

            int snp_position = atoi(entries[1].c_str());
            gff_entry* gff_entry_ptr = get_gff_entry(seq_gff_list_ptr, snp_position);

            if (gff_entry_ptr == NULL) { // gene not found, ignore this SNP position
                continue;
            }
            
            gff_entry_ptr->set_valid(false);
        }

    } // end while


    this->write(out_file, true);

    // clean up
    seq_gff_list_ptr = NULL;

    // close the stream
    ifs_snp.close();

}

void gff::filter_genes_with_problematic_snps( const string& snp_file, const string& out_file) {

    // open the file streams
    ifstream ifs_snp(snp_file.c_str());

    // go through the list of snp file to save snps each reference (chromosome)
    string str_snp = "";
    string sequence_name = "";
    string current_sequence_name = ""; // reference (chromosome) being currently read
    gff::SEQ_GFF_LIST* seq_gff_list_ptr = NULL;
    while (ifs_snp.good()) {

        // read the snp 
        getline(ifs_snp, str_snp);

        if (str_snp.empty()) { // ignore empty lines
            continue;
        }

        vector<string> entries;
        utility::str_split(str_snp, entries, "\t");

        if (entries.size() < 4 || entries[3].empty() || (!dna::is_nucleotide(entries[3]) && !dna::is_extended_nucleotide(entries[3]))) {
            // invalid entry           

            // get the sequence name
            sequence_name = entries[0];

            if (sequence_name.compare(current_sequence_name) != 0) { // new reference
                seq_gff_list_ptr = this->gff_list[sequence_name];
                current_sequence_name = sequence_name;
            }

            int snp_position = atoi(entries[1].c_str());
            gff_entry* gff_entry_ptr = get_gff_entry(seq_gff_list_ptr, snp_position);

            gff_entry_ptr->set_valid(false);
        }

    } // end while


    this->write(out_file, true);

    // clean up
    seq_gff_list_ptr = NULL;

    // close the stream
    ifs_snp.close();

}

void gff::read_gene_coordinate_file(const string& coordinate_file) {
    // open the input stream
    ifstream ifs_coordinate(coordinate_file.c_str());
    string str_coordinate = "";
    string previous_sequence_name = "";
    SEQ_GFF_LIST* seq_gff_list_ptr;
    while (ifs_coordinate.good()) {
        // read the next line from the GFF file
        getline(ifs_coordinate, str_coordinate);

        if (str_coordinate.empty()) { // ignore empty lines
            continue;
        }

        vector<string> entries;
        utility::str_split(str_coordinate, entries, "\t");

        gff_entry* gff_entry_ptr = new gff_entry(entries[1], atoi(entries[2].c_str()), atoi(entries[3].c_str()), entries[0]);

        if (gff_entry_ptr->get_sequence_name().compare(previous_sequence_name) != 0) { // new reference is starting
            if (!previous_sequence_name.empty()) { // this is not the first reference
                this->gff_list[previous_sequence_name] = seq_gff_list_ptr;
            }

            // initialise the GFF list for the new sequence
            seq_gff_list_ptr = new vector<gff_entry*>;
            // save the sequence name
            previous_sequence_name = gff_entry_ptr->get_sequence_name();
        }

        seq_gff_list_ptr->push_back(gff_entry_ptr);

    } // end while ifs_gff is good
    this->gff_list[previous_sequence_name] = seq_gff_list_ptr;

}

bool gff::filter_features(set<string>& features_to_keep) {
    
    if (features_to_keep.empty()) {
        return false;
    }
    
    bool gff_changed = false;
    GFF_LIST::iterator the_gff_list = this->gff_list.begin();
    for (; the_gff_list != this->gff_list.end(); the_gff_list++) {
        set<SEQ_GFF_LIST::iterator> gff_entries_to_be_deleted;
        SEQ_GFF_LIST* seq_gff_list_ptr = the_gff_list->second;
        SEQ_GFF_LIST::iterator the_seq_gff_list = seq_gff_list_ptr->begin();
        for (; the_seq_gff_list != seq_gff_list_ptr->end(); the_seq_gff_list++) {
            gff_entry* gff_entry_ptr = *the_seq_gff_list;

            if (features_to_keep.find(gff_entry_ptr->get_feature()) == features_to_keep.end()) { // this feature is not to be kept
                    // delete the GFF entry
                delete *the_seq_gff_list;
                gff_entries_to_be_deleted.insert(the_seq_gff_list);
                gff_changed = true;
            }
        }
        
        for (set<SEQ_GFF_LIST::iterator>::reverse_iterator it = gff_entries_to_be_deleted.rbegin(); it != gff_entries_to_be_deleted.rend(); ++it) {
                seq_gff_list_ptr->erase(*it);
        }
        
    }
    
    return gff_changed;
}

void gff::get_sequence_names(set<string>& sequence_names) {
    GFF_LIST::iterator the_seq_gff_list = gff_list.begin();
    for(; the_seq_gff_list != gff_list.end(); the_seq_gff_list++) {
        sequence_names.insert(the_seq_gff_list->first);
    }
}

void gff::create_index_array(map<string, int*>&  index_array) {

    //map<string, int*> index_array;// = new map<string, int*>();
    GFF_LIST::iterator the_gff_list = this->gff_list.begin();
    for (; the_gff_list != this->gff_list.end(); the_gff_list++) {
        string sequence_name = the_gff_list->first;
        SEQ_GFF_LIST* seq_gff_list_ptr = the_gff_list->second;
        
        gff_entry* gff_entry_ptr = *(seq_gff_list_ptr->rbegin());
        int* seq_index_array = new int[gff_entry_ptr->get_end()+1];        
        int gene_idx = 1;
        SEQ_GFF_LIST::iterator the_seq_gff_list = seq_gff_list_ptr->begin();
        for (; the_seq_gff_list != seq_gff_list_ptr->end(); the_seq_gff_list++) {
            gff_entry* gff_entry_ptr = *the_seq_gff_list;
            for (int idx = gff_entry_ptr->get_start() - OFFSET; idx <= gff_entry_ptr->get_end() + OFFSET; idx++) {
                seq_index_array[idx]= gene_idx;
            }
            gene_idx++;
        }
        index_array[sequence_name] = seq_index_array;
    }
    
}

void gff::delete_index_array(map<string, int*>&  index_array) {

    //map<string, int*> index_array;// = new map<string, int*>();
    map<string, int*>::iterator index_array_ptr = index_array.begin();
    for (; index_array_ptr != index_array.end(); index_array_ptr++) {
        delete index_array_ptr->second;
    }

    
}
