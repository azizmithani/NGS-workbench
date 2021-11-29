/* 
 * File:   snp.cpp
 * Author: aziz
 * 
 * Created on May 24, 2011, 10:57 AM
 */

#include <string>
#include <fstream>
#include <iostream>
#include "snp.h"
#include "gff.h"
#include "gff_entry.h"
#include "utility.h"

snp::snp() {
}

snp::snp(int position, char base) {
    this->position = position;
    this->base = base;
}

snp::snp(const snp& orig) {
    this->position = orig.position;
    this->base = orig.base;
}

snp::~snp() {
}

void snp::setBase(char base) {
    this->base = base;
}

char snp::get_base() const {
    return base;
}

void snp::set_position(int position) {
    this->position = position;
}

int snp::get_position() const {
    return position;
}

bool snp::operator==(const snp& the_snp) const {
    return (this->position == the_snp.position & this->base == the_snp.base);
}

bool snp::operator<(const snp& the_snp) const {
    if (this->position < the_snp.position) {
        return true;
    } else if (this->position == the_snp.position) {
        if (this->base < the_snp.base) {
            return true;
        }
    }
    return false;
}

bool snp::operator>(const snp& the_snp) const {
    if (this->position > the_snp.position) {
        return true;
    } else if (this->position == the_snp.position) {
        if (this->base > the_snp.base) {
            return true;
        }
    }
    return false;
}

void snp::extract_snps(const string& snp_file, const string& coordinate_file, const string& out_file, bool ignore_Ns_in_the_reference) {

    // create gff list from the coordinate file
    gff* gff_ptr = new gff();
    gff_ptr->read_gene_coordinate_file(coordinate_file);

    // open the file streams
    ifstream ifs_snp(snp_file.c_str());
    ofstream ofs_out(out_file.c_str());

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

        // get the sequence name
        sequence_name = entries[0];

        if (ignore_Ns_in_the_reference && entries[2].compare("N") == 0) { // ignore positions with "N" in the reference, if suggested so
            continue;
        } else if (entries[2].size() > 1 || (!dna::is_nucleotide(entries[2]) && !dna::is_extended_nucleotide(entries[2]))) { // ignore headers
            continue;
        }

        if (sequence_name.compare(current_sequence_name) != 0) { // new reference
            seq_gff_list_ptr = gff_ptr->get_gff_list(sequence_name);
            current_sequence_name = sequence_name;
        }

        int snp_position = atoi(entries[1].c_str());
        gff_entry* gff_entry_ptr = gff::get_gff_entry(seq_gff_list_ptr, snp_position);

        if (gff_entry_ptr == NULL) { // GFF entry not found, ignore the SNP
            continue;
        }

        ofs_out << str_snp << endl;

    } // end while


    // clean up
    seq_gff_list_ptr = NULL;

    // close the stream
    ifs_snp.close();

}

void snp::check_polyploid_snp_list(dna::BASE_LIST* polyploid_snp_list, base_distribution* polyploid_base_distribution_list) {

    dna::BASE_LIST::iterator the_snp_list = polyploid_snp_list->begin();
    for (; the_snp_list != polyploid_snp_list->end(); the_snp_list++) {
        string sequence_name = the_snp_list->first;
        dna::SEQ_BASE_LIST* seq_snp_list_ptr = the_snp_list->second;

        // go through all SNPs
        if (seq_snp_list_ptr == NULL) {
            continue;
        }

        position_base_distribution** seq_base_distribution_list = polyploid_base_distribution_list->get_base_distribution_list(sequence_name);
        int sequence_length = (*(polyploid_base_distribution_list->get_sequence_lengths()))[sequence_name];

        for (int position = 0; position < sequence_length; position++) {
            char snp_base = seq_snp_list_ptr[position];

            if (snp_base == dna::BASE_BLANK) { // no snp
                continue;
            }

            // get the coverage for this position
            int coverage = seq_base_distribution_list[position]->get_coverage();

            set<char> new_bases;

            if (check_polyploid_base(seq_base_distribution_list[position]->get_coverage_A(), coverage)) {
                new_bases.insert('A');
            }
            if (check_polyploid_base(seq_base_distribution_list[position]->get_coverage_T(), coverage)) {
                new_bases.insert('T');
            }
            if (check_polyploid_base(seq_base_distribution_list[position]->get_coverage_G(), coverage)) {
                new_bases.insert('G');
            }
            if (check_polyploid_base(seq_base_distribution_list[position]->get_coverage_C(), coverage)) {
                new_bases.insert('C');
            }

            if (!new_bases.empty()) {
                char new_base = (*dna::get_reverse_nucleotide_map())[utility::to_string(new_bases)];
                if (snp_base != new_base) {
                    seq_snp_list_ptr[position] = new_base;
                }
            }
        } // end for each position

    } // end for each sequence

}

void snp::identify_missing_polyploid_snps(dna::BASE_LIST* polyploid_snp_list, base_distribution* polyploid_base_distribution_list, dna::BASE_LIST* control_snp_list, dna::BASE_LIST* polyploid_reference_base_list, dna::BASE_LIST* control_reference_base_list) {

    dna::BASE_LIST::iterator the_control_snp_list = control_snp_list->begin();
    for (; the_control_snp_list != control_snp_list->end(); the_control_snp_list++) {
        string sequence_name = the_control_snp_list->first;
        dna::SEQ_BASE_LIST* seq_control_snp_list_ptr = the_control_snp_list->second;

        if (seq_control_snp_list_ptr == NULL) {
            continue;
        }

        // also get polyploid snp list for this sequence
        dna::SEQ_BASE_LIST* seq_polyploid_snp_list_ptr = (*polyploid_snp_list)[sequence_name];
        // and the reference base lists
        dna::SEQ_BASE_LIST* seq_polyploid_reference_base_list_ptr = (*polyploid_reference_base_list)[sequence_name];
        dna::SEQ_BASE_LIST* seq_control_reference_base_list_ptr = (*control_reference_base_list)[sequence_name];


        position_base_distribution** seq_base_distribution_list = polyploid_base_distribution_list->get_base_distribution_list(sequence_name);
        int sequence_length = (*(polyploid_base_distribution_list->get_sequence_lengths()))[sequence_name];

        // go through all SNPs
        for (int position = 0; position < sequence_length; position++) {

            if (seq_control_snp_list_ptr[position] == dna::BASE_BLANK) { // no snp
                continue;
            }

            // get the base in polyploid
            char snp_base = seq_polyploid_snp_list_ptr[position];
            if (snp_base != dna::BASE_BLANK) { //  SNP has been called at this position, move to the next position
                continue;
            }

            // set the reference base
            seq_polyploid_reference_base_list_ptr[position] = seq_control_reference_base_list_ptr[position];

            // get the coverage for this position
            int coverage = seq_base_distribution_list[position]->get_coverage();

            if (coverage == 0) {
                seq_polyploid_snp_list_ptr[position] = dna::BASE_ZERO_COVERAGE;
                continue;
            } else if (coverage < POLYPLOID_COVERAGE_THRESHOLD) {
                seq_polyploid_snp_list_ptr[position] = dna::BASE_LOW_COVERAGE;
                continue;
            }

            set<char> new_bases;

            if (check_polyploid_base(seq_base_distribution_list[position]->get_coverage_A(), coverage)) {
                new_bases.insert('A');
            }
            if (check_polyploid_base(seq_base_distribution_list[position]->get_coverage_T(), coverage)) {
                new_bases.insert('T');
            }
            if (check_polyploid_base(seq_base_distribution_list[position]->get_coverage_G(), coverage)) {
                new_bases.insert('G');
            }
            if (check_polyploid_base(seq_base_distribution_list[position]->get_coverage_C(), coverage)) {
                new_bases.insert('C');
            }

            if (!new_bases.empty()) {
                char new_base = (*dna::get_reverse_nucleotide_map())[utility::to_string(new_bases)];
                //                if (new_base != seq_polyploid_reference_base_list_ptr[position]) {
                seq_polyploid_snp_list_ptr[position] = new_base;
                //                } else {
                //                    cout << position  + 1 << endl;
                //                }
            }
        } // end for each position

    } // end for each sequence

}

void snp::identify_missing_polyploid_snps(dna::BASE_LIST* polyploid_snp_list, base_distribution* polyploid_base_distribution_list, dna::BASE_LIST* polyploid_reference_base_list, dna::BASE_LIST* control_reference_base_list) {

    dna::BASE_LIST::iterator the_snp_list = polyploid_snp_list->begin();
    for (; the_snp_list != polyploid_snp_list->end(); the_snp_list++) {
        string sequence_name = the_snp_list->first;
        dna::SEQ_BASE_LIST* seq_snp_list_ptr = the_snp_list->second;

        if (seq_snp_list_ptr == NULL) {
            continue;
        }

        // get the reference base lists for this sequence
        dna::SEQ_BASE_LIST* seq_polyploid_reference_base_list_ptr = (*polyploid_reference_base_list)[sequence_name];
        dna::SEQ_BASE_LIST* seq_control_reference_base_list_ptr = (*control_reference_base_list)[sequence_name];

        // and the base distribution
        position_base_distribution** seq_base_distribution_list = polyploid_base_distribution_list->get_base_distribution_list(sequence_name);
        int sequence_length = (*(polyploid_base_distribution_list->get_sequence_lengths()))[sequence_name];

        // go through all positions
        for (int position = 0; position < sequence_length; position++) {

            if (seq_snp_list_ptr[position] != dna::BASE_BLANK) { //SNP has been called at this position, move to the next position
                continue;
            }

            // get the coverage for this position
            int coverage = seq_base_distribution_list[position]->get_coverage();

            if (coverage == 0) {
                continue;
            } else if (coverage < POLYPLOID_COVERAGE_STRICT_THRESHOLD) {
                continue;
            } else {
                // check coverage on 5 base-pairs on either sides
                bool low_coverage = false;
                for (int p = position - 5; p <= position + 5; p++) {
                    if (p < 0) {
                        continue;
                    } else if (p == position) { // we have already checked the coverage at the position
                        continue;
                    } else if (p >= sequence_length) {
                        break;
                    }
                    if (seq_base_distribution_list[p]->get_coverage() < POLYPLOID_BASE_COVERAGE_THRESHOLD) {
                        low_coverage = true;
                        break;
                    }
                } // end for

                if (low_coverage) {
                    continue;
                }
            }

            set<char> new_bases;

            if (check_polyploid_base_strict(seq_base_distribution_list[position]->get_coverage_A(), coverage)) {
                new_bases.insert('A');
            }
            if (check_polyploid_base_strict(seq_base_distribution_list[position]->get_coverage_T(), coverage)) {
                new_bases.insert('T');
            }
            if (check_polyploid_base_strict(seq_base_distribution_list[position]->get_coverage_G(), coverage)) {
                new_bases.insert('G');
            }
            if (check_polyploid_base_strict(seq_base_distribution_list[position]->get_coverage_C(), coverage)) {
                new_bases.insert('C');
            }

            if (!new_bases.empty()) {
                char new_base = (*dna::get_reverse_nucleotide_map())[utility::to_string(new_bases)];
                if (new_base != seq_control_reference_base_list_ptr[position]) { // its a SNP
                    // set the SNP base
                    seq_snp_list_ptr[position] = new_base;
                    // set the reference base
                    seq_polyploid_reference_base_list_ptr[position] = seq_control_reference_base_list_ptr[position];
                    //                } else {
                    //                    cout << position  + 1 << endl;
                }
            }
        } // end for each position

    } // end for each sequence

}

void snp::check_diploid_snp_list(dna::BASE_LIST* diploid_snp_list, base_distribution* diploid_base_distribution_list, dna::BASE_LIST* polyploid_snp_list, dna::BASE_LIST* diploid_reference_base_list, dna::BASE_LIST* polyploid_reference_base_list) {

    dna::BASE_LIST::iterator the_polyploid_snp_list = polyploid_snp_list->begin();
    for (; the_polyploid_snp_list != polyploid_snp_list->end(); the_polyploid_snp_list++) {
        // get sequence name
        string sequence_name = the_polyploid_snp_list->first;
        // get snp list for this sequence
        dna::SEQ_BASE_LIST* seq_polyploid_snp_list_ptr = the_polyploid_snp_list->second;

        if (seq_polyploid_snp_list_ptr == NULL) {
            continue;
        }

        // also get diploid snp list for this sequence
        dna::SEQ_BASE_LIST* seq_diploid_snp_list_ptr = (*diploid_snp_list)[sequence_name];
        // and the reference base lists for this sequence
        dna::SEQ_BASE_LIST* seq_diploid_reference_base_list_ptr = (*diploid_reference_base_list)[sequence_name];
        dna::SEQ_BASE_LIST* seq_polyploid_reference_base_list_ptr = (polyploid_reference_base_list == NULL ? NULL : (*polyploid_reference_base_list)[sequence_name]);
        // and the base distribution list
        position_base_distribution** seq_base_distribution_list = diploid_base_distribution_list->get_base_distribution_list(sequence_name);
        int sequence_length = (*(diploid_base_distribution_list->get_sequence_lengths()))[sequence_name];

        // go through all SNPs
        for (int position = 0; position < sequence_length; position++) {

            if (seq_polyploid_snp_list_ptr[position] == dna::BASE_BLANK) { // no snp
                continue;
            }

            // get the base in diploid
            char snp_base = seq_diploid_snp_list_ptr[position];
            if (snp_base != dna::BASE_BLANK) { //  SNP has been called at this position
                // check if its an ambiguous base
                string str_snp_base(1, snp_base);
                if (dna::is_extended_nucleotide(str_snp_base)) {
                    seq_diploid_snp_list_ptr[position] = dna::BASE_AMBIGUOUS;
                }
                //move to the next position
                continue;
            }

            // we are now checking the position which is not called a SNP in the diploid.
            // check if multiple bases are present at this position (possibly, heterozygous)

            if (seq_polyploid_reference_base_list_ptr != NULL) {
                // set the reference base
                seq_diploid_reference_base_list_ptr[position] = seq_polyploid_reference_base_list_ptr[position];
            }

            // get the coverage for this position
            int coverage = seq_base_distribution_list[position]->get_coverage();

            if (coverage == 0) {
                seq_diploid_snp_list_ptr[position] = dna::BASE_ZERO_COVERAGE;
                continue;
            } else if (coverage < DIPLOID_COVERAGE_THRESHOLD) {
                seq_diploid_snp_list_ptr[position] = dna::BASE_LOW_COVERAGE;
                continue;
            }

            set<char> new_bases;

            if (check_diploid_base(seq_base_distribution_list[position]->get_coverage_A(), coverage)) {
                new_bases.insert('A');
            }
            if (check_diploid_base(seq_base_distribution_list[position]->get_coverage_T(), coverage)) {
                new_bases.insert('T');
            }
            if (check_diploid_base(seq_base_distribution_list[position]->get_coverage_G(), coverage)) {
                new_bases.insert('G');
            }
            if (check_diploid_base(seq_base_distribution_list[position]->get_coverage_C(), coverage)) {
                new_bases.insert('C');
            }

            if (new_bases.empty()) {
                // ideally, this should not happen
                seq_diploid_snp_list_ptr[position] = seq_diploid_reference_base_list_ptr[position];
            } else if (new_bases.size() == 1) {
                seq_diploid_snp_list_ptr[position] = *new_bases.begin();
            } else if (new_bases.size() > 1) {
                seq_diploid_snp_list_ptr[position] = dna::BASE_AMBIGUOUS;
            }

        } // end for each position

    } // end for each sequence


}

void snp::check_diploid_snp_list(dna::BASE_LIST* diploid_snp_list, base_distribution* diploid_base_distribution_list, dna::BASE_LIST* diploid_reference_base_list) {

    dna::BASE_LIST::iterator the_diploid_reference_base_list = diploid_reference_base_list->begin();
    for (; the_diploid_reference_base_list != diploid_reference_base_list->end(); the_diploid_reference_base_list++) {
        // get sequence name
        string sequence_name = the_diploid_reference_base_list->first;
        // get reference base list for this sequence
        dna::SEQ_BASE_LIST* seq_diploid_reference_base_list_ptr = the_diploid_reference_base_list->second;

        if (seq_diploid_reference_base_list_ptr == NULL) {
            continue;
        }

        // also get diploid snp list for this sequence
        dna::SEQ_BASE_LIST* seq_diploid_snp_list_ptr = (*diploid_snp_list)[sequence_name];
        // and the base distribution list
        position_base_distribution** seq_base_distribution_list = diploid_base_distribution_list->get_base_distribution_list(sequence_name);
        int sequence_length = (*(diploid_base_distribution_list->get_sequence_lengths()))[sequence_name];

        // go through all SNPs
        for (int position = 0; position < sequence_length; position++) {

            if (seq_diploid_reference_base_list_ptr[position] == dna::BASE_BLANK) { // no snp
                continue;
            }

            // get the base in diploid
            char snp_base = seq_diploid_snp_list_ptr[position];
            if (snp_base != dna::BASE_BLANK) { //  SNP has been called at this position
                // check if its an ambiguous base
                string str_snp_base(1, snp_base);
                if (dna::is_extended_nucleotide(str_snp_base)) {
                    seq_diploid_snp_list_ptr[position] = dna::BASE_AMBIGUOUS;
                }
                //move to the next position
                continue;
            }

            // we are now checking the position which is not called a SNP in the diploid.

            // get the coverage for this position
            int coverage = seq_base_distribution_list[position]->get_coverage();

            if (coverage == 0) {
                seq_diploid_snp_list_ptr[position] = dna::BASE_ZERO_COVERAGE;
                continue;
            } else if (coverage < DIPLOID_COVERAGE_THRESHOLD) {
                seq_diploid_snp_list_ptr[position] = dna::BASE_LOW_COVERAGE;
                continue;
            }

            set<char> new_bases;

            if (check_diploid_base(seq_base_distribution_list[position]->get_coverage_A(), coverage)) {
                new_bases.insert('A');
            }
            if (check_diploid_base(seq_base_distribution_list[position]->get_coverage_T(), coverage)) {
                new_bases.insert('T');
            }
            if (check_diploid_base(seq_base_distribution_list[position]->get_coverage_G(), coverage)) {
                new_bases.insert('G');
            }
            if (check_diploid_base(seq_base_distribution_list[position]->get_coverage_C(), coverage)) {
                new_bases.insert('C');
            }

            if (new_bases.empty()) {
                // ideally, this should not happen
                seq_diploid_snp_list_ptr[position] = seq_diploid_reference_base_list_ptr[position];
            } else if (new_bases.size() == 1) {
                seq_diploid_snp_list_ptr[position] = *new_bases.begin();
            } else if (new_bases.size() > 1) {
                seq_diploid_snp_list_ptr[position] = dna::BASE_AMBIGUOUS;
            }

        } // end for each position

    } // end for each sequence


}

bool snp::check_polyploid_base(int base_coverage, int total_coverage) {
    if (base_coverage >= POLYPLOID_BASE_COVERAGE_THRESHOLD && (double) base_coverage / (double) total_coverage >= POLYPLOID_BASE_PROPORTION_THRESHOLD) {
        return true;
    } else {
        return false;
    }

}

bool snp::check_polyploid_base_strict(int base_coverage, int total_coverage) {
    if (base_coverage >= POLYPLOID_BASE_COVERAGE_STRICT_THRESHOLD && (double) base_coverage / (double) total_coverage >= POLYPLOID_BASE_PROPORTION_THRESHOLD) {
        return true;
    } else {
        return false;
    }

}

bool snp::check_diploid_base(int base_coverage, int total_coverage) {
    if (base_coverage >= DIPLOID_BASE_COVERAGE_THRESHOLD && (double) base_coverage / (double) total_coverage >= DIPLOID_BASE_PROPORTION_THRESHOLD) {
        return true;
    } else {
        return false;
    }

}

string snp::to_string() {
    return utility::to_string(position) + ":" + utility::to_string(base);
}

void snp::snps_per_gene(const string& snp_file, const string& gff_or_coordinate_file, const string& out_file, bool ignore_Ns_in_the_reference, bool is_gff_file, bool is_snp_table, bool count_snp_positions) {

    // create gff list from the coordinate file
    gff* gff_ptr;
    if (is_gff_file) {
        gff_ptr = new gff(gff_or_coordinate_file);
    } else {
        gff_ptr = new gff();
        gff_ptr->read_gene_coordinate_file(gff_or_coordinate_file);
    }

    // count the number of lines
    int n_lines = 1;
    if (is_snp_table && !count_snp_positions) {
        string str_snp;
        ifstream ifs_snp_lines(snp_file.c_str());
        while (ifs_snp_lines.good()) {

            // read the snp 
            getline(ifs_snp_lines, str_snp);

            if (str_snp.empty()) { // ignore empty lines
                continue;
            }

            vector<string> entries;
            utility::str_split(str_snp, entries, "\t");
            n_lines = entries.size() - 3;
            break;
        }
        ifs_snp_lines.close();
    }

    // open the file streams
    ifstream ifs_snp(snp_file.c_str());
    ofstream ofs_out(out_file.c_str());

    // variable to hold snp counts
    map <string, int*> snp_counts;
    // initialise the memory
    set<string> sequence_names;
    gff_ptr->get_sequence_names(sequence_names);
    set<string>::iterator the_sequence_name = sequence_names.begin();
    for (; the_sequence_name != sequence_names.end(); the_sequence_name++) {
        string sequence_name = *the_sequence_name;

        // get the GFF List for this sequence
        gff::SEQ_GFF_LIST* seq_gff_list_ptr = gff_ptr->get_gff_list(sequence_name);
        // and the iterator
        gff::SEQ_GFF_LIST::iterator the_seq_gff_list = seq_gff_list_ptr->begin();
        // delete GFF entries for all genes
        for (; the_seq_gff_list != seq_gff_list_ptr->end(); the_seq_gff_list++) {
            gff_entry* gff_entry_ptr = *the_seq_gff_list;

            string gff_id = gff_entry_ptr->get_id();
            snp_counts[gff_id] = new int[n_lines];
            for (int c = 0; c < n_lines; c++) {
                (snp_counts[gff_id])[c] = 0;
            }
        }
    }


    // go through the list of snp file to save snps each reference (chromosome)
    string str_snp = "";
    string sequence_name = "";
    string current_sequence_name = ""; // reference (chromosome) being currently read
    gff::SEQ_GFF_LIST* seq_gff_list_ptr = NULL;
    map <int, string> line_names;
    while (ifs_snp.good()) {

        // read the snp 
        getline(ifs_snp, str_snp);

        if (str_snp.empty()) { // ignore empty lines
            continue;
        }

        vector<string> entries;
        utility::str_split(str_snp, entries, "\t");

        // get the sequence name
        sequence_name = entries[0];

        if (ignore_Ns_in_the_reference && entries[2].compare("N") == 0) { // ignore positions with "N" in the reference, if suggested so
            continue;
        } else if (entries[2].size() > 1 || (!dna::is_nucleotide(entries[2]) && !dna::is_extended_nucleotide(entries[2]))) { // ignore headers
            for (int col = 3; col < entries.size(); col++) {
                line_names[col] = entries[col];
            }

            continue;
        }

        if (sequence_name.compare(current_sequence_name) != 0) { // new reference
            seq_gff_list_ptr = gff_ptr->get_gff_list(sequence_name);
            current_sequence_name = sequence_name;
        }

        int snp_position = atoi(entries[1].c_str());
        gff_entry* gff_entry_ptr = gff::get_gff_entry(seq_gff_list_ptr, snp_position);

        if (gff_entry_ptr == NULL) { // GFF entry not found, ignore the SNP
            continue;
        }

        string reference_base = entries[2];
        // get the snp counts corresponding to this gene 
        map <string, int*>::iterator the_snp_count = snp_counts.find(gff_entry_ptr->get_id());
        if (count_snp_positions) {
            *(the_snp_count->second) += 1;
        } else {
            for (int col = 3; col < entries.size(); col++) {
                if (entries[col].empty() || entries[col].compare(reference_base) == 0) { // no SNP
                    continue;
                } else {
                    *(the_snp_count->second + col - 3) += 1;
                }

            }
        }

    } // end while


    // clean up
    seq_gff_list_ptr = NULL;

    // close the stream
    ifs_snp.close();

    // write the output
    ofs_out << gff_entry::header_short();
    if (count_snp_positions || !is_snp_table) {
        ofs_out << "\tSNPs per Gene" << endl;
    } else {
        for (int col = 0; col < n_lines; col++) {
            string header = line_names[col + 3];
            if (header.empty()) {
                header = utility::to_string(col + 1);
            }
            ofs_out << "\t" << header;
        }
        ofs_out << endl;
    }
    //    set<string> sequence_names;
    //    gff_ptr->get_sequence_names(sequence_names);
    //    set<string>::iterator the_sequence_name = sequence_names.begin();
    the_sequence_name = sequence_names.begin();
    for (; the_sequence_name != sequence_names.end(); the_sequence_name++) {
        string sequence_name = *the_sequence_name;

        // get the GFF List for this sequence
        gff::SEQ_GFF_LIST* seq_gff_list_ptr = gff_ptr->get_gff_list(sequence_name);
        // and the iterator
        gff::SEQ_GFF_LIST::iterator the_seq_gff_list = seq_gff_list_ptr->begin();
        // delete GFF entries for all genes
        for (; the_seq_gff_list != seq_gff_list_ptr->end(); the_seq_gff_list++) {
            gff_entry* gff_entry_ptr = *the_seq_gff_list;

            string gff_id = gff_entry_ptr->get_id();
            map <string, int*>::iterator the_snp_count = snp_counts.find(gff_id);
            ofs_out << gff_entry_ptr->to_string_short();
            for (int line = 0; line < n_lines; line++) {
                ofs_out << "\t" << *(the_snp_count->second + line);
            }
            ofs_out << endl;

        }
    }

    ofs_out.close();

}

void snp::get_consensus_base(const string& snp_file, const string& out_file) {

    // open the file streams
    ifstream ifs_snp(snp_file.c_str());
    ofstream ofs_out(out_file.c_str());

    // get the nucleotide map
    dna::NUCLEOTIDE_MAP* nucleotide_map_ptr = dna::get_nucleotide_map();
    dna::REVERSE_NUCLEOTIDE_MAP* reverse_nucleotide_map_ptr = dna::get_reverse_nucleotide_map();

    // go through the snp file 
    string str_snp = "";
    while (ifs_snp.good()) {

        // read the snp 
        getline(ifs_snp, str_snp);

        if (str_snp.empty()) { // ignore empty lines
            continue;
        }

        vector<string> entries;
        utility::str_split(str_snp, entries, "\t");

        string consensus_base;
        // accept the alternate base if "N" in the reference, or reference base is not part of the consensus
        if (entries[2].compare("N") == 0 || entries[8].substr(0, 3).compare("1/1") == 0) {
            consensus_base = entries[3];
        } else if (entries[8].substr(0, 3).compare("0/1") == 0) { // reference base is to included in the consensus
            // get the bases corresponding to the alternate base (alternate base could be an ambiguous character)
            set<char> bases = (*nucleotide_map_ptr)[(entries[3])[0]];
            // insert the reference base
            bases.insert((entries[2])[0]);

            // get the consensus base
            consensus_base = (*reverse_nucleotide_map_ptr)[utility::to_string(bases)];
        } else {
            continue;
        }

        ofs_out << entries[0] << "\t" << entries[1] << "\t" << entries[2] << "\t" << consensus_base << endl;

    } // end while


    // close the stream
    ifs_snp.close();
    ofs_out.close();
}

void snp::get_consensus_base_GATK(const string& snp_file, const string& out_file) {

    // open the file streams
    ifstream ifs_snp(snp_file.c_str());
    ofstream ofs_out(out_file.c_str());

    // get the nucleotide map
    dna::NUCLEOTIDE_MAP* nucleotide_map_ptr = dna::get_nucleotide_map();
    dna::REVERSE_NUCLEOTIDE_MAP* reverse_nucleotide_map_ptr = dna::get_reverse_nucleotide_map();

    // go through the snp file 
    string str_snp = "";
    while (ifs_snp.good()) {

        // read the snp 
        getline(ifs_snp, str_snp);

        if (str_snp.empty()) { // ignore empty lines
            continue;
        } else if (str_snp[0] == '#') {
            continue;
        }

        vector<string> entries;
        utility::str_split(str_snp, entries, "\t");
        
        set<char> bases;
        string consensus_base;
        // accept the alternate base if "N" in the reference, or reference base is not part of the consensus
        int pos = entries[9].find_first_of(':');
        if (entries[9].substr(0, pos).find_first_of('0') != string::npos) {
            bases.insert((entries[3])[0]);
        }
        
        vector<string> alt_bases;
        utility::str_split(entries[4], alt_bases, ",");
        vector<string>::iterator it_alt_bases = alt_bases.begin();
        for (; it_alt_bases != alt_bases.end(); it_alt_bases++){
            bases.insert((*it_alt_bases)[0]);
        }

        consensus_base = (*reverse_nucleotide_map_ptr)[utility::to_string(bases)];
        
        ofs_out << entries[0] << "\t" << entries[1] << "\t" << entries[3] << "\t" << consensus_base << endl;

    } // end while


    // close the stream
    ifs_snp.close();
    ofs_out.close();
}

void snp::compare_snps(const string& snp_table_file, const string& polyploid_snp_file) {//, const string& out_file) {

    // read SNP list
    dna::SNP_LIST polyploid_snp_list;
    dna::read_snp_file(polyploid_snp_file, polyploid_snp_list);


    // open the file stream
    ifstream ifs_snp(snp_table_file.c_str());

    // get the nucleotide map
    dna::NUCLEOTIDE_MAP* nucleotide_map_ptr = dna::get_nucleotide_map();
    dna::REVERSE_NUCLEOTIDE_MAP* reverse_nucleotide_map_ptr = dna::get_reverse_nucleotide_map();

    // variable to store the header (e.g., name of the diploid)
    map<int, string> headers;
    // count of shared snps
    map<int, int> shared_snps_count;
    int positions_ignored = 0;
    int total_positions = 0;

    // go through the snp file 
    string str_snp = "";
    string sequence_name = "";
    string current_sequence_name = ""; // reference (chromosome) being currently read
    dna::SEQ_SNP_LIST* seq_polyploid_snp_list;
    while (ifs_snp.good()) {

        // read the snp 
        getline(ifs_snp, str_snp);

        if (str_snp.empty()) { // ignore empty lines
            continue;
        }

        vector<string> entries;
        utility::str_split(str_snp, entries, "\t");

        if (entries[2].size() > 1) { // header

            for (int col = 3; col < entries.size(); col++) {
                headers[col] = entries[col];
            }
            continue;
        }

        // increment the SNP positions count
        total_positions++;

        // get the sequence name
        sequence_name = entries[0];

        if (sequence_name.compare(current_sequence_name) != 0) { // new reference
            seq_polyploid_snp_list = polyploid_snp_list[sequence_name];
            current_sequence_name = sequence_name;
        }

        // get the position
        int position = atoi(entries[1].c_str());
        // get the polyploid base at this position
        string polyploid_base = (*seq_polyploid_snp_list)[position];

        if (polyploid_base.empty() || (!dna::is_nucleotide(polyploid_base) && !dna::is_extended_nucleotide(polyploid_base))) {
            // snp not called in polyploid, or problematic position (e.g. low coverage), ignore
            positions_ignored++;
            //cout << position  << "\t" << polyploid_base << endl;
            continue;
        }

        // get the set of all bases for the polyploid base
        set<char> bases = (*nucleotide_map_ptr)[polyploid_base[0]];
        for (int col = 3; col < entries.size(); col++) {
            string base = entries[col];
            if (bases.find(base[0]) != bases.end()) {
                shared_snps_count[col]++;
            }
        }

    } // end while

    // close the stream
    ifs_snp.close();


    cout << "Accession\tShared SNPs" << endl;
    for (int col = 0; col < shared_snps_count.size(); col++) {
        string header = headers[col + 3];
        if (header.empty()) {
            header = utility::to_string(col + 1);
        }
        cout << header << "\t" << shared_snps_count[col + 3] << endl;
    }

    cout << "Total positions: " << total_positions << endl;
    cout << "Positions ignored: " << positions_ignored << endl;
    //    ofstream ofs_out(out_file.c_str());
    //    ofs_out.close();
}
