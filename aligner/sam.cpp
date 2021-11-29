/* 
 * File:   sam.cpp
 * Author: Aziz Mithani
 * 
 * Created on May 7, 2011, 4:55 PM
 */

#include <fstream>
#include <iostream>
#include <vector>
#include <set>
#include "dna.h"
#include "sam.h"
#include "utility.h"
#include "position_base_distribution.h"

//sam::sam() {
//    initialise();
//}
//

sam::sam(const string& sam_file) {
    initialise();
    this->sam_file = sam_file;
    read_sequence_lengths();
}

sam::sam(const string& sam_file, const string& tmp_dir, bool paired_data, const string& unique_mapping_tag, const string& repeat_mapping_tag) {
    initialise();
    this->sam_file = sam_file;
    this->b_paired_data = paired_data;
    this->unique_mapping_tag = unique_mapping_tag;
    this->repeat_mapping_tag = repeat_mapping_tag;
    this->tmp_dir = tmp_dir;
    read_sequence_lengths();
}

sam::~sam() {

    // delete the base distribution list
    delete this->the_base_distribution;
    // delete sequence lengths
    delete this->sequence_lengths;
    //this->sequence_lengths.clear();
}

void sam::initialise() {
    this->sam_file = "";
    this->unique_mapping_tag = "";
    this->repeat_mapping_tag = "";
    this->tmp_dir = "";
    this->b_paired_data = true;
    this->b_ignore_non_uniquely_mapped_reads = true;
    this->b_ignore_unpaired_reads = false;
    this->mapping_quality_threshold = 20;
    this->the_base_distribution = NULL;
    this->sequence_lengths = new dna::SEQUENCE_LENGTH_LIST;

    this->path_sam_tools = ""; ///Volumes/SAN001/src/samtools-0.1.12/";
    this->path_sam_tools_pl = ""; ///Volumes/SAN001/src/samtools-0.1.12/misc/";
}

void sam::set_tmp_dir(string tmp_dir) {
    this->tmp_dir = tmp_dir;
}

string sam::get_tmp_dir() const {
    return tmp_dir;
}

void sam::set_base_distribution(base_distribution* the_base_distribution) {
    this->the_base_distribution = the_base_distribution;
}

base_distribution* sam::get_base_distribution() const {
    return the_base_distribution;
}

void sam::set_sequence_lengths(dna::SEQUENCE_LENGTH_LIST* sequence_lengths) {
    this->sequence_lengths = sequence_lengths;
}

dna::SEQUENCE_LENGTH_LIST* sam::get_sequence_lengths() const {
    return sequence_lengths;
}

void sam::set_mapping_quality_threshold(int mapping_quality_threshold) {
    this->mapping_quality_threshold = mapping_quality_threshold;
}

int sam::get_mapping_quality_threshold() const {
    return mapping_quality_threshold;
}

void sam::set_ignore_unpaired_reads(bool b_ignore_unpaired_reads) {
    this->b_ignore_unpaired_reads = b_ignore_unpaired_reads;
}

bool sam::ignore_unpaired_reads() const {
    return b_ignore_unpaired_reads;
}

void sam::set_unique_mapping_tag(string unique_mapping_tag) {
    this->unique_mapping_tag = unique_mapping_tag;
}

string sam::get_unique_mapping_tag() const {
    return unique_mapping_tag;
}

void sam::set_ignore_non_uniquely_mapped_reads(bool b_ignore_non_uniquely_mapped_reads) {
    this->b_ignore_non_uniquely_mapped_reads = b_ignore_non_uniquely_mapped_reads;
}

bool sam::ignore_non_uniquely_mapped_reads() const {
    return b_ignore_non_uniquely_mapped_reads;
}

void sam::set_paired_data(bool b_paired_data) {
    this->b_paired_data = b_paired_data;
}

bool sam::paired_data() const {
    return b_paired_data;
}

void sam::set_sam_file(string sam_file) {
    this->sam_file = sam_file;
}

string sam::get_sam_file() const {
    return sam_file;
}

void sam::set_read_count(int read_count) {
    this->read_count = read_count;
}

int sam::get_read_count() const {
    return read_count;
}

void sam::read_sequence_lengths() {
    if (sam_file.empty()) {
        return;
    } else if (!utility::file_exists(sam_file)) {
        return;
    }

    if (sequence_lengths == NULL) {
        sequence_lengths = new dna::SEQUENCE_LENGTH_LIST;
    }
    sequence_lengths->clear();
    // open the input file stream
    ifstream ifs_sam(sam_file.c_str());

    string str_sam = "";
    while (ifs_sam.good()) {
        getline(ifs_sam, str_sam);

        if (str_sam.empty()) { // ignore empty lines
            continue;
        } else if (str_sam[0] == '@') { // header
            if (str_sam.substr(1, 2).compare("SQ") == 0) {
                // get sequence lengths from the header
                vector<string> entries;
                utility::str_split(str_sam, entries, "\t");

                string sequence_name = entries[1].substr(3, entries[1].length());
                int sequence_length = atoi(entries[2].substr(3, entries[2].length()).c_str());

                (*sequence_lengths)[sequence_name] = sequence_length;
            }
        } else {
            break;
        }

    } // end while

    ifs_sam.close();
}

int sam::extract_aligned_and_unaligned_reads(const string& raw_sam_file, const string& out_aligned_reads_file, const string& out_unaligned_reads_file) {
    /*
     * Returns the number of aligned reads
     */

    // open the input file stream
    ifstream ifs_sam(raw_sam_file.c_str());
    // open the output file streams
    ofstream ofs_aligned(out_aligned_reads_file.c_str());
    ofstream ofs_unaligned(out_unaligned_reads_file.c_str());

    string str_sam = "";
    int aligned_reads = 0;
    while (ifs_sam.good()) {
        getline(ifs_sam, str_sam);

        if (str_sam.empty()) { // ignore empty lines
            continue;
        } else if (str_sam[0] == '@') { // header
            // write to both aligned and unaligned files
            ofs_unaligned << str_sam << endl;
            ofs_aligned << str_sam << endl;
            continue;
        }


        string delimiter = "\t";
        string::size_type last_pos = 0;
        string::size_type pos = 0;
        for (int i = 0; i < 3; i++) {
            // Skip delimiters.  Note the "not_of"
            last_pos = str_sam.find_first_not_of(delimiter, pos);
            // Find next "non-delimiter".
            pos = str_sam.find_first_of(delimiter, last_pos);
        }
        string sequence_name = str_sam.substr(last_pos, pos - last_pos);

        if (sequence_name[0] == '*') { // unaligned reads
            ofs_unaligned << str_sam << endl;
        } else { // aligned reads
            ofs_aligned << str_sam << endl;
            aligned_reads++;
        }
    } // end while (ifs_sam.good())

    ifs_sam.close();
    ofs_aligned.close();
    ofs_unaligned.close();

    return aligned_reads;
}

int sam::filter(const string& out_filtered_sam_file) {

    return filter(out_filtered_sam_file, -1, "");
}

int sam::filter(const string& out_filtered_sam_file, int max_insert_size, const string& gff_file) {
    /*
     * Removes unwanted reads
     */

    // open the input file stream
    ifstream ifs_sam(sam_file.c_str());
    // open the output file stream
    ofstream ofs_filtered(out_filtered_sam_file.c_str());

    // read gff file, if specified
    gff* gff_ptr;
    if (!gff_file.empty()) {
        gff_ptr = new gff(gff_file);
    }

    string read_1 = "";
    string read_2 = "";
    int filtered_reads_count = 0;
    string sequence_name_1 = "";
    string current_sequence_name = ""; // reference (chromosome) being currently read
    gff::SEQ_GFF_LIST* seq_gff_list_ptr = NULL;
    while (ifs_sam.good()) {
        getline(ifs_sam, read_1);

        if (read_1.empty()) { // ignore empty lines
            continue;
        } else if (read_1[0] == '@') { // header
            // write to the output 
            ofs_filtered << read_1 << endl;
            continue;
        }

        // break read into components
        vector<string> entries_1;
        utility::str_split(read_1, entries_1, "\t");

        // ignore supplementary alignment (this will be actually a supplementary alignment of previous read pair - 2nd read)
        if (is_supplementary_alignment(entries_1)) {
            continue;
        }

        // check if this read should be printed
        bool read_1_mapped = is_mapped(entries_1); // read mapped?
        bool is_valid_read = true;
        bool is_unique_1 = true;
        if (read_1_mapped) { // only check if the read is mapped
            is_valid_read = is_high_quality(entries_1); // is mapping quality high?
            //            is_unique_1 = is_unique(entries_1);
            //            if (require_both_reads_unique) {
            //                is_valid_read = is_valid_read && is_unique_1;
            //            }
        }

        bool read_2_mapped = false;
        if (this->b_paired_data) {
            // get the 2nd read in the pair
            getline(ifs_sam, read_2);

            // break read into components
            vector<string> entries_2;
            utility::str_split(read_2, entries_2, "\t");

            // ignore supplementary alignment (this will be actually a supplementary alignment of the current read pair - 1st read)
            while (is_supplementary_alignment(entries_2)) {
                // get the 2nd read in the pair
                getline(ifs_sam, read_2);

                // break read into components
                entries_2.clear();
                utility::str_split(read_2, entries_2, "\t");
            }

            if (!is_valid_read) { // first read in the pair is not valid, no need to check the second read
                continue;
            }

            // check if this read should be printed
            read_2_mapped = is_mapped(entries_2); // read mapped?

            if (read_2_mapped) {// only check if the read is mapped
                // check if this read should be printed
                is_valid_read = is_high_quality(entries_2);
                //                bool is_unique_2 = is_unique(entries_2);
                //                if (require_both_reads_unique) {
                //                    is_valid_read = is_valid_read && is_unique_2;
                //                } else {
                //                    // one may be unique, other may be repeat
                //                    bool is_repeat_1 = true;
                //                    bool is_repeat_2 = true;
                //                    if (!is_unique_1) {
                //                        is_repeat_1 = is_repeat(entries_1);
                //                    }
                //                    if (!is_unique_2) {
                //                        is_repeat_2 = is_repeat(entries_2);
                //                    }
                //                    is_valid_read = is_valid_read && ((is_unique_1 && is_unique_2) || (is_unique_1 && is_repeat_2) ||
                //                            (is_repeat_1 && is_unique_2));
                //                }
            }

            if (!is_valid_read) { // this read in the pair is not valid, continue
                continue;
            }

            // check for correct orientation. ignore F1F2 and R1R2. Don't check if reads mapped on different contigs/genes/chromosomes is required
            if (max_insert_size == -1) {
                is_valid_read = is_correct_orientation(entries_1, entries_2);
                if (!is_valid_read) { // this read pair is not in the correct orientation
                    continue;
                }
            }

            if (max_insert_size == 0) { // insert size to be ignored when filtering
                is_valid_read = true;
            } else if (max_insert_size > 0 && (abs(atoi(entries_1[8].c_str())) > max_insert_size)) {
                is_valid_read = false;
            } else if (gff_ptr != NULL) {
                // get the sequence name
                sequence_name_1 = entries_1[2];
                string sequence_name_2 = entries_2[2];

                if (sequence_name_1.compare(sequence_name_2) != 0) { // different references, don't print this read pair
                    is_valid_read = false;
                } else {
                    if (sequence_name_1.compare(current_sequence_name) != 0) { // new reference
                        seq_gff_list_ptr = gff_ptr->get_gff_list(sequence_name_1);
                        current_sequence_name = sequence_name_1;
                    }

                    if (seq_gff_list_ptr == NULL) { // reference not found in the GFF file, ignore 
                        continue;
                    }

                    int read_start_1 = atoi(entries_1[3].c_str());
                    int read_start_2 = atoi(entries_2[3].c_str());
                    gff_entry* gff_entry_1_ptr = gff::get_gff_entry(seq_gff_list_ptr, read_start_1);
                    gff_entry* gff_entry_2_ptr = gff::get_gff_entry(seq_gff_list_ptr, read_start_2);

                    if (gff_entry_1_ptr == NULL || gff_entry_2_ptr == NULL) { // GFF entry not found, don't print the read pair
                        is_valid_read = (max_insert_size == -1 ? false : true);
                    } else {
                        if (gff_entry_1_ptr == gff_entry_2_ptr) { // check if reads are mapped on the same gene or different genes
                            is_valid_read = (max_insert_size == -1 ? true : false); 
                        } else {
                            is_valid_read = (max_insert_size == -1 ? false : true);
                        }
                    }

                } // if sequence name 1 != sequence name 2
            }

        } else {
            b_ignore_unpaired_reads = false;
        }

        if (b_ignore_unpaired_reads) {
            if (is_valid_read && read_1_mapped && read_2_mapped) {
                ofs_filtered << read_1 << endl;
                ofs_filtered << read_2 << endl;
                // increment the filtered reads count
                filtered_reads_count += 2;
            }
        } else {
            if (is_valid_read && read_1_mapped) {
                ofs_filtered << read_1 << endl;
                // increment the filtered reads count
                filtered_reads_count++;
            }
            if (is_valid_read && read_2_mapped) {
                ofs_filtered << read_2 << endl;
                // increment the filtered reads count
                filtered_reads_count++;
            }
        }
    } // end while (ifs_sam.good())

    ofs_filtered.close();
    //    if (this->b_paired_data) {
    //        return filtered_reads_count / 2;
    //    } else {
    //        return filtered_reads_count;
    //    }

    return filtered_reads_count;

}

// uses an array index for gff instead of searching the gff entry pointers. speeds up things considerably

int sam::filter_new(const string& out_filtered_sam_file, int max_insert_size, const string& gff_file) {
    /*
     * Removes unwanted reads
     */

    // open the input file stream
    ifstream ifs_sam(sam_file.c_str());
    // open the output file stream
    ofstream ofs_filtered(out_filtered_sam_file.c_str());

    // read gff file, if specified
    gff* gff_ptr;
    map<string, int*> index_array;
    if (!gff_file.empty()) {
        gff_ptr = new gff(gff_file);
        // create index for each gene
        gff_ptr->create_index_array(index_array);
    }


    string read_1 = "";
    string read_2 = "";
    int filtered_reads_count = 0;
    string sequence_name_1 = "";
    string current_sequence_name = ""; // reference (chromosome) being currently read
    //gff::SEQ_GFF_LIST* seq_gff_list_ptr = NULL;
    int* seq_index_array = NULL;
    while (ifs_sam.good()) {
        getline(ifs_sam, read_1);

        if (read_1.empty()) { // ignore empty lines
            continue;
        } else if (read_1[0] == '@') { // header
            // write to the output 
            ofs_filtered << read_1 << endl;
            continue;
        }

        // break read into components
        vector<string> entries_1;
        utility::str_split(read_1, entries_1, "\t");

        // ignore supplementary alignment (this will be actually a supplementary alignment of previous read pair - 2nd read)
        if (is_supplementary_alignment(entries_1)) {
            continue;
        }

        // check if this read should be printed
        bool read_1_mapped = is_mapped(entries_1); // read mapped?
        bool is_valid_read = true;
        //bool is_unique_1 = true;
        if (read_1_mapped) { // only check if the read is mapped
            is_valid_read = is_high_quality(entries_1); // is mapping quality high?
        }

        bool read_2_mapped = false;
        if (this->b_paired_data) {
            // get the 2nd read in the pair
            getline(ifs_sam, read_2);

            // break read into components
            vector<string> entries_2;
            utility::str_split(read_2, entries_2, "\t");

            // ignore supplementary alignment (this will be actually a supplementary alignment of the current read pair - 1st read)
            while (is_supplementary_alignment(entries_2)) {
                // get the 2nd read in the pair
                getline(ifs_sam, read_2);

                // break read into components
                entries_2.clear();
                utility::str_split(read_2, entries_2, "\t");
            }

            if (!is_valid_read) { // first read in the pair is not valid, no need to check the second read
                continue;
            }

            // check if this read should be printed
            read_2_mapped = is_mapped(entries_2); // read mapped?

            if (read_2_mapped) {// only check if the read is mapped
                // check if this read should be printed
                is_valid_read = is_high_quality(entries_2);
            }

            if (!is_valid_read) { // this read in the pair is not valid, continue
                continue;
            }

            // check for correct orientation. ignore F1F2 and R1R2. Don't check if reads mapped on different contigs/genes/chromosomes is required
            if (max_insert_size == -1) {
                is_valid_read = is_correct_orientation(entries_1, entries_2);
                if (!is_valid_read) { // this read pair is not in the correct orientation
                    continue;
                }
            }

            if (max_insert_size == 0) { // insert size to be ignored when filtering
                is_valid_read = true;
            } else if (max_insert_size > 0 && (abs(atoi(entries_1[8].c_str())) > max_insert_size)) {
                is_valid_read = false;
            } else if (gff_ptr != NULL) {
                // get the sequence name
                sequence_name_1 = entries_1[2];
                string sequence_name_2 = entries_2[2];

                if (sequence_name_1.compare(sequence_name_2) != 0) { // different references, don't print this read pair
                    is_valid_read = (max_insert_size == -1 ? false : true);
                } else {
                    if (sequence_name_1.compare(current_sequence_name) != 0) { // new reference
                        seq_index_array = index_array[sequence_name_1];
                        current_sequence_name = sequence_name_1;
                    }

                    if (seq_index_array == NULL) { // reference not found in the GFF file, ignore 
                        continue;
                    }

                    int read_start_1 = atoi(entries_1[3].c_str());
                    int read_start_2 = atoi(entries_2[3].c_str());

                    //                    if (entries_1[0].compare("K00114:181:H3LG5BBXX:8:1101:1010:39501") == 0) {
                    //                        cout << seq_index_array[read_start_1] << "\t" << seq_index_array[read_start_2] << endl;
                    //                    }

                    if (seq_index_array[read_start_1] == seq_index_array[read_start_2]) {
                        is_valid_read = (max_insert_size == -1 ? true : false);
                    } else {
                        is_valid_read = (max_insert_size == -1 ? false : true);
                    }
                } // if sequence name 1 != sequence name 2
            }

        } else {
            b_ignore_unpaired_reads = false;
        }

        //        if (entries_1[0].compare("K00114:181:H3LG5BBXX:8:1101:1010:39501") == 0) {
        //            cout << is_valid_read << endl;
        //        }


        if (b_ignore_unpaired_reads) {
            if (is_valid_read && read_1_mapped && read_2_mapped) {
                ofs_filtered << read_1 << endl;
                ofs_filtered << read_2 << endl;
                // increment the filtered reads count
                filtered_reads_count += 2;
            }
        } else {
            if (is_valid_read && read_1_mapped) {
                ofs_filtered << read_1 << endl;
                // increment the filtered reads count
                filtered_reads_count++;
            }
            if (is_valid_read && read_2_mapped) {
                ofs_filtered << read_2 << endl;
                // increment the filtered reads count
                filtered_reads_count++;
            }
        }


    } // end while (ifs_sam.good())

    ofs_filtered.close();
    //    if (this->b_paired_data) {
    //        return filtered_reads_count / 2;
    //    } else {
    //        return filtered_reads_count;
    //    }

    if (!gff_file.empty()) {
        seq_index_array = NULL;
        gff_ptr->delete_index_array(index_array);
        delete gff_ptr;
    }
    return filtered_reads_count;

}

/*
 int sam::filter(const string& out_filtered_sam_file, int max_insert_size, const string& gff_file) {

    filter(out_filtered_sam_file, max_insert_size, gff_file, true);
 }

 int sam::filter(const string& out_filtered_sam_file, int max_insert_size, const string& gff_file, bool require_both_reads_unique) {
    //
    // Removes unwanted reads
    //

    // open the input file stream
    ifstream ifs_sam(sam_file.c_str());
    // open the output file stream
    ofstream ofs_filtered(out_filtered_sam_file.c_str());

    // read gff file, if specified
    gff* gff_ptr;
    if (!gff_file.empty()) {
        gff_ptr = new gff(gff_file);
    }

    string read_1 = "";
    string read_2 = "";
    int filtered_reads_count = 0;
    string sequence_name_1 = "";
    string current_sequence_name = ""; // reference (chromosome) being currently read
    gff::SEQ_GFF_LIST* seq_gff_list_ptr = NULL;
    while (ifs_sam.good()) {
        getline(ifs_sam, read_1);

        if (read_1.empty()) { // ignore empty lines
            continue;
        } else if (read_1[0] == '@') { // header
            // write to the output 
            ofs_filtered << read_1 << endl;
            continue;
        }

        // break read into components
        vector<string> entries_1;
        utility::str_split(read_1, entries_1, "\t");

        // check if this read should be printed
        bool read_1_mapped = is_mapped(entries_1); // read mapped?
        bool is_valid_read = true;
        bool is_unique_1 = true;
        if (read_1_mapped) { // only check if the read is mapped
            is_valid_read = is_high_quality(entries_1); // is mapping quality high?
            is_unique_1 = is_unique(entries_1);
            if (require_both_reads_unique) {
                is_valid_read = is_valid_read && is_unique_1;
            }
        }

        bool read_2_mapped = false;
        if (this->b_paired_data) {
            // get the 2nd read in the pair
            getline(ifs_sam, read_2);

            if (!is_valid_read) { // first read in the pair is not valid, no need to check the second read
                continue;
            }
            // break read into components
            vector<string> entries_2;
            utility::str_split(read_2, entries_2, "\t");

            // check if this read should be printed
            read_2_mapped = is_mapped(entries_2); // read mapped?

            if (read_2_mapped) {// only check if the read is mapped
                // check if this read should be printed
                is_valid_read = is_high_quality(entries_2);
                bool is_unique_2 = is_unique(entries_2);
                if (require_both_reads_unique) {
                    is_valid_read = is_valid_read && is_unique_2;
                } else {
                    // one may be unique, other may be repeat
                    bool is_repeat_1 = true;
                    bool is_repeat_2 = true;
                    if (!is_unique_1) {
                        is_repeat_1 = is_repeat(entries_1);
                    }
                    if (!is_unique_2) {
                        is_repeat_2 = is_repeat(entries_2);
                    }
                    is_valid_read = is_valid_read && ((is_unique_1 && is_unique_2) || (is_unique_1 && is_repeat_2) ||
                            (is_repeat_1 && is_unique_2));
                }
            }

            if (!is_valid_read) { // this read in the pair is not valid, continue
                continue;
            }

            if (max_insert_size == 0) { // insert size to be ignored when filtering
                is_valid_read = true;
            } else if (max_insert_size > 0 && (abs(atoi(entries_1[8].c_str())) > max_insert_size)) {
                is_valid_read = false;
            } else if (gff_ptr != NULL) {
                // get the sequence name
                sequence_name_1 = entries_1[2];
                string sequence_name_2 = entries_2[2];

                if (sequence_name_1.compare(sequence_name_2) != 0) { // different references, don't print this read pair
                    is_valid_read = false;
                } else {
                    if (sequence_name_1.compare(current_sequence_name) != 0) { // new reference
                        seq_gff_list_ptr = gff_ptr->get_gff_list(sequence_name_1);
                        current_sequence_name = sequence_name_1;
                    }

                    if (seq_gff_list_ptr == NULL) { // reference not found in the GFF file, ignore 
                        continue;
                    }

                    int read_start_1 = atoi(entries_1[3].c_str());
                    int read_start_2 = atoi(entries_2[3].c_str());
                    gff_entry* gff_entry_1_ptr = gff::get_gff_entry(seq_gff_list_ptr, read_start_1);
                    gff_entry* gff_entry_2_ptr = gff::get_gff_entry(seq_gff_list_ptr, read_start_2);

                    if (gff_entry_1_ptr == NULL || gff_entry_2_ptr == NULL) { // GFF entry not found, don't print the read pair
                        is_valid_read = false;
                    } else if (gff_entry_1_ptr != gff_entry_2_ptr) { // reads mapped on different genes, don't print the read pair
                        is_valid_read = false;
                    }
                } // if sequence name 1 != sequence name 2
            }

        } else {
            b_ignore_unpaired_reads = false;
        }

        if (b_ignore_unpaired_reads) {
            if (is_valid_read && read_1_mapped && read_2_mapped) {
                ofs_filtered << read_1 << endl;
                ofs_filtered << read_2 << endl;
                // increment the filtered reads count
                filtered_reads_count += 2;
            }
        } else {
            if (is_valid_read && read_1_mapped) {
                ofs_filtered << read_1 << endl;
                // increment the filtered reads count
                filtered_reads_count++;
            }
            if (is_valid_read && read_2_mapped) {
                ofs_filtered << read_2 << endl;
                // increment the filtered reads count
                filtered_reads_count++;
            }
        }
    } // end while (ifs_sam.good())

    ofs_filtered.close();
    //    if (this->b_paired_data) {
    //        return filtered_reads_count / 2;
    //    } else {
    //        return filtered_reads_count;
    //    }

    return filtered_reads_count;

}
 
 */

bool sam::is_high_quality(vector<string>& entries) {
    return (atoi(entries[4].c_str()) < mapping_quality_threshold) ? false : true;
}

bool sam::is_unique(vector<string>& entries) {
    if (b_ignore_non_uniquely_mapped_reads && !unique_mapping_tag.empty() && entries[11].find(unique_mapping_tag) == string::npos) {
        return false;
    } else {
        return true;
    }
}

bool sam::is_repeat(vector<string>& entries) {
    if (!repeat_mapping_tag.empty() && entries[11].find(repeat_mapping_tag) == string::npos) {
        return false;
    } else {
        return true;
    }
}

bool sam::is_mapped(vector<string>& entries) {
    if (atoi(entries[1].c_str()) & FLAG_SEQUENCE_UNMAPPED) {
        return false;
    } else {
        return true;
    }

}

bool sam::is_supplementary_alignment(vector<string>& entries) {
    if (atoi(entries[1].c_str()) & FLAG_SUPPLEMENTARY_ALIGNMENT) {
        return true;
    } else {
        return false;
    }

}

bool sam::is_correct_orientation(vector<string>& entries_1, vector<string>& entries_2) {
    if ((atoi(entries_1[2].c_str()) & FLAG_SEQUENCE_ON_REVERSE_STRAND) && (atoi(entries_2[2].c_str()) & FLAG_SEQUENCE_ON_REVERSE_STRAND)) {//R1R2
        return false;
    } else if ((atoi(entries_1[2].c_str()) & ~FLAG_SEQUENCE_ON_REVERSE_STRAND) && (atoi(entries_2[2].c_str()) & ~FLAG_SEQUENCE_ON_REVERSE_STRAND)) {//F1F2
        return false;
    } else {
        return true;
    }

}

void sam::calculate_base_distribution(const string & out_base_distribution_file) {
    /*
     * Calculate base distributions
     */

    //initialise base distributions list
    the_base_distribution = new base_distribution(this->sequence_lengths, true);

    // open the input file stream
    ifstream ifs_sam(sam_file.c_str());

    string str_sam = "";
    while (ifs_sam.good()) {
        getline(ifs_sam, str_sam);

        if (str_sam.empty()) { // ignore empty lines
            continue;
        } else if (str_sam[0] == '@') { // header
            continue;
        }

        vector<string> entries;
        utility::str_split(str_sam, entries, "\t");

        string sequence_name = entries[2];
        int read_start = atoi(entries[3].c_str());
        string read = entries[9];
        string base_qualities = entries[10];

        // find the occurrence of "D","I" or "S" in 6th column (5th with 0-based indexing)
        int position;
        if ((position = entries[5].find("D")) != string::npos || (position = entries[5].find("I")) != string::npos || (position = entries[5].find("S")) != string::npos) {
            // process this read

            string new_read = "";
            string new_base_qualities = "";
            string cigar = entries[5];
            string length = "";
            int last_non_digit_position = -1;
            int sum = 0;
            int i = 0;
            int start = 0;
            int total_deleted = 0;
            while (i < cigar.length()) {
                for (i = start; i < cigar.length(); i++) {
                    if (isdigit(cigar[i])) {
                        length += cigar[i];
                    } else if (cigar[i] == 'S') {
                        if (i == cigar.length() - 1) {
                            for (int j = read.length() - 1; j >= read.length() - atoi(length.c_str()); j--) {
                                read[j] = '-';
                                base_qualities[j] = '-';
                            }
                        } else {
                            read = read.substr(atoi(length.c_str())); //,read.length());
                            base_qualities = base_qualities.substr(atoi(length.c_str()));
                        }
                        length = "";
                        last_non_digit_position = i;
                    } else if (cigar[i] == 'D') {
                        int deletion_size = atoi(cigar.substr(last_non_digit_position + 1, i - last_non_digit_position).c_str());
                        last_non_digit_position = i;
                        int deletion_start = sum;
                        start = i + 1;
                        length = "";
                        //sum+= DelLength;

                        string gap(deletion_size, '.');
                        new_read += read.substr(0, deletion_start) + gap;
                        read = read.substr(deletion_start, read.length());

                        new_base_qualities += base_qualities.substr(0, deletion_start) + gap;
                        base_qualities = base_qualities.substr(deletion_start, base_qualities.length());

                        total_deleted += deletion_start;
                        sum = 0;
                        break;
                    } else if (cigar[i] == 'I') {
                        int insertion_size = atoi(cigar.substr(last_non_digit_position + 1, i - last_non_digit_position).c_str());
                        last_non_digit_position = i;
                        int insertion_start = sum;
                        start = i + 1;
                        sum = 0;
                        length = "";

                        new_read += read.substr(0, insertion_start);
                        read = read.substr(insertion_start + insertion_size, read.length());

                        new_base_qualities += base_qualities.substr(0, insertion_start);
                        base_qualities = base_qualities.substr(insertion_start + insertion_size, base_qualities.length());

                        total_deleted += insertion_start + insertion_size;

                        break;
                    } else if (cigar[i] == 'M' || cigar[i] == 'N') {
                        sum += atoi(length.c_str());
                        length = "";
                        last_non_digit_position = i;
                    }
                }
            } // end while
            new_read += read.substr(0, read.length());
            new_base_qualities += base_qualities.substr(0, base_qualities.length());

            read = new_read;
            base_qualities = new_base_qualities;
        }

        // get the base distribution for the sequence to which the read is mapped
        position_base_distribution** seq_base_distribution_ptr = the_base_distribution->get_base_distribution_list(sequence_name);
        int sequence_length = (*sequence_lengths)[sequence_name];
        for (int i = 0; i < read.length(); i++) {
            if (read[i] == '-') { // soft-clip .. ignore this position in this read
                continue;
            } else if (read[i] == '.') { // gap due to deletion .. ignore
                continue;
            }

            if (read_start + i > sequence_length)
                break;
            // increment the frequency for the base at this position
            //cout << read_start + i << "\t" << seq_base_distribution_ptr[read_start + i - 1]->to_string() << endl;

            // Ignore this base if the sequencing quality is low
            int base_quality = int(base_qualities[i]) - 33; // base qualities are ascii - 33
            if (base_quality < BASE_QUALITY_THRESHOLD) {
                continue;
            }

            seq_base_distribution_ptr[read_start + i - 1]->increment_coverage(read[i]);
            //the_base_distribution->set_base_distribution_list(sequence_name, seq_base_distribution_ptr);
        } // end for (each character in the read) 

    } // end while (ifs_sam.good())

    ifs_sam.close();

    // write to file, if a file is provided
    if (!out_base_distribution_file.empty()) {
        the_base_distribution->write(out_base_distribution_file);
    }
}

void sam::to_bam(const string& out_bam_file, bool full_run) {

    // file for output by samtools which we don't need 
    string tmp_file = this->tmp_dir + utility::extract_filename(this->sam_file) + "_to_bam.tmp";

    string command;
    // covert the SAM file to BAM file
    string unsorted_bam_file = out_bam_file + ".unsorted";
    command = "samtools view -bS -o " + unsorted_bam_file + " " + this->sam_file + " 2> " + tmp_file;
    system(command.c_str());

    // sort the bam file by coordinates
    command = "samtools sort " + unsorted_bam_file + " " + utility::remove_extension(out_bam_file) + " 2> " + tmp_file;
    system(command.c_str());

    // remove the unsorted BAM file
    utility::remove_file(unsorted_bam_file);

    if (full_run) {
        // index the bam file
        command = "samtools index " + out_bam_file + " 2> " + tmp_file;
        system(command.c_str());

    }
    // convert the sorted BAM file to back to SAM file
    command = "samtools view -h -o " + this->sam_file + " " + out_bam_file + " 2> " + tmp_file;
    system(command.c_str());

    // remove the temporary output file
    utility::remove_file(tmp_file);

}

// Changed on 5-Oct-2011 11:50 to use mpileup command in samtools 0.1.18 (r982:295).
// Previous version renamed to generate_variant_lists_old and kept below

void sam::generate_variant_lists(const string& reference_file, const string& bam_file, bool report_ambiguous_bases,
        int min_coverage_threshold, int max_coverage_threshold, int snp_quality_threshold, int indel_quality_threshold,
        const string& out_pileup_file, const string& out_snp_file, const string& out_insertion_file, const string& out_deletion_file,
        const string& out_accepted_snp_file, const string& out_accepted_insertion_file, const string& out_accepted_deletion_file,
        int* snp_count, int* snp_count_n, int* insertion_count, int* deletion_count, int* accepted_snp_count, int* accepted_insertion_count, int* accepted_deletion_count,
        bool old_indel_format, bool generate_pileup) { // old indel format is used for genome fixing

    // file for output by samtools which we don't need 
    string tmp_file = this->tmp_dir + utility::extract_filename(this->sam_file) + "_generate_variant_lists.tmp";

    string command;
    // create the index file
    command = "samtools faidx " + reference_file + " 2> " + tmp_file;
    system(command.c_str());

    // generate the pileup, if required
    if (generate_pileup) {
        generate_pileup_file(reference_file, bam_file, tmp_file, out_pileup_file, max_coverage_threshold);
    }

    //string pileup_file_no_n = out_pileup_file + ".noN";
    string pileup_file_filtered = out_pileup_file + ".filtered";

    vector< string > snp_list_n;

    // open the input file stream
    ifstream ifs_pileup(out_pileup_file.c_str());
    // get the list of positions with Ns and write the positions without Ns to the output file
    // also remove the indels caused by very few (misaligned) reads
    string str_pileup = "";
    while (ifs_pileup.good()) {
        getline(ifs_pileup, str_pileup);

        if (str_pileup.empty()) { // ignore empty lines
            continue;
        } else if (str_pileup[0] == '#') { // ignore comments/header
            continue;
        }

        vector<string> entries;
        utility::str_split(str_pileup, entries, "\t");
        // Get those positions where the reference has an 'N' and the concensus base does not have an "N", and meets the quality threshold
        if (entries[3].compare("N") == 0 && entries[4].compare("N") != 0 && entries[4].compare(".") != 0 && atof(entries[5].c_str()) >= snp_quality_threshold && entries[7].find("INDEL") == string::npos) {
            // extract the coverage
            int coverage_start = entries[7].find("DP=");
            int coverage_end = entries[7].find(";", coverage_start);
            int coverage = atoi(entries[7].substr(coverage_start + 3, coverage_end - coverage_start - 3).c_str());
            if (coverage >= min_coverage_threshold) {
                if (!report_ambiguous_bases) {
                    if (entries[4].find(",") != string::npos) {
                        continue;
                    } else if (dna::is_extended_nucleotide(entries[4])) {
                        continue;
                    }
                }
                string str_snp = utility::to_string(entries, 0, 1) + "\t" + entries[3] + "\t" + get_consensus_base(entries[4], ",") + "\t" + utility::to_string(entries, 5, 9);
                snp_list_n.push_back(str_snp);

            }
            //            snp_list_n.push_back(utility::to_string(entries, 0, 7));

        }
    } // end while ifsPileup.good()

    ifs_pileup.close();



    // filter the pileup file to remove error-prone variant calls caused by factors not considered in the statistical model of samtools.
    // NOTE: we use "-N 3 -W 1" so that dense SNPs are not filtered out. We use -w 5 to filter out snps within 5 bp of an indel
    command = "vcfutils.pl varFilter -1 0 -4 0 -d " + utility::to_string(min_coverage_threshold) + " -D " + utility::to_string(max_coverage_threshold) + " " + out_pileup_file + " > " + pileup_file_filtered;
    system(command.c_str());

    // open output streams for variations
    ofstream ofs_snp(out_snp_file.c_str());
    ofstream ofs_accepted_snp(out_accepted_snp_file.c_str());
    ofstream ofs_insertion(out_insertion_file.c_str());
    ofstream ofs_deletion(out_deletion_file.c_str());
    ofstream ofs_accepted_insertion(out_accepted_insertion_file.c_str());
    ofstream ofs_accepted_deletion(out_accepted_deletion_file.c_str());

    // initialise
    int local_snp_count = 0;
    int local_snp_count_n = 0;
    int local_insertion_count = 0;
    int local_deletion_count = 0;
    int local_accepted_snp_count = 0;
    int local_accepted_insertion_count = 0;
    int local_accepted_deletion_count = 0;

    // open the input file stream
    ifstream ifs_pileup_filtered(pileup_file_filtered.c_str());
    // get the list of positions with SNPs and Indels
    string str_pileup_filtered = "";
    while (ifs_pileup_filtered.good()) {
        getline(ifs_pileup_filtered, str_pileup_filtered);

        if (str_pileup_filtered.empty()) { // ignore empty lines
            continue;
        } else if (str_pileup_filtered[0] == '#') { // ignore comments/header
            continue;
        }

        vector<string> entries;
        utility::str_split(str_pileup_filtered, entries, "\t");

        if (entries[7].find("INDEL") != string::npos) { // indel
            string str_indel;
            bool is_insertion = false;
            if (old_indel_format) {
                str_indel = utility::to_string(entries, 0, 1);
                string indel_bases;
                string long_seq, short_seq;
                if (entries[3].length() < entries[4].length()) { // insertion
                    long_seq = entries[4];
                    short_seq = entries[3];
                    is_insertion = true;
                } else { // deletion
                    long_seq = entries[3];
                    short_seq = entries[4];
                    is_insertion = false;
                }
                indel_bases = long_seq.substr(1, long_seq.length() - 1);
                string suffix = "";
                if (short_seq.length() > 1) {
                    suffix = short_seq.substr(1);
                }
                indel_bases.erase(indel_bases.length() - suffix.length());
                str_indel += "\t" + indel_bases + "\t" + utility::to_string(entries, 5, entries.size() - 1);
            } else {
                str_indel = utility::to_string(entries, 0, 1) + "\t" + utility::to_string(entries, 3, entries.size() - 1);
                if (entries[3].length() < entries[4].length()) { // insertion
                    is_insertion = true;
                } else { // deletion
                    is_insertion = false;
                }
            }
            if (is_insertion) {
                ofs_insertion << str_indel << endl;
                local_insertion_count++;
            } else {
                ofs_deletion << str_indel << endl;
                local_deletion_count++;
            }
            if (entries[4].compare(".") != 0 && entries[4].find(",") == string::npos && entries[9].substr(0, 3).compare("0/1") != 0 && entries[9].substr(0, 3).compare("1/0") != 0) { // only consider positions where there is no ambiguity
                if (atoi(entries[5].c_str()) >= indel_quality_threshold) {
                    if (entries[3].length() < entries[4].length()) { // insertion
                        ofs_accepted_insertion << str_indel << endl;
                        local_accepted_insertion_count++;
                    } else { // deletion
                        ofs_accepted_deletion << str_indel << endl;
                        local_accepted_deletion_count++;
                    }
                }
            }

        } else if (entries[4].compare(".") != 0 && entries[3].compare("N") != 0) { // different base than the reference and consensus is not N
            string str_output = utility::to_string(entries, 0, 1) + "\t" + entries[3] + "\t" + get_consensus_base(entries[4], ",") + "\t" + utility::to_string(entries, 5, 9);
            // add position to the snp list
            ofs_snp << str_output << endl;
            local_snp_count++;
            if (atoi(entries[5].c_str()) >= snp_quality_threshold) { // check if the SNP quality is meets the threshold
                if (!report_ambiguous_bases && (entries[4].find(",") != string::npos || entries[9].substr(0, 3).compare("0/1") == 0 || entries[9].substr(0, 3).compare("1/0") == 0)) { // ambiguous bases to be reported?
                    // no and this is an ambiguous base -> do nothing
                    continue;
                } else {
                    // this position meets the criteria. add it to the snp list
                    ofs_accepted_snp << str_output << endl;
                    local_accepted_snp_count++;
                }
            }
        }
    } // end while ifsPileupFiltered.good()
    ifs_pileup_filtered.close();

    // add positions with Ns in the reference to the output
    vector<string>::iterator the_snp = snp_list_n.begin();
    for (; the_snp != snp_list_n.end(); the_snp++) {
        ofs_snp << *the_snp << endl;
        ofs_accepted_snp << *the_snp << endl;
    }
    // set and increase the SNP count for 'N'
    local_snp_count_n = snp_list_n.size();
    local_snp_count += snp_list_n.size();
    local_accepted_snp_count += snp_list_n.size();

    // remove intermediate files
    utility::remove_file(pileup_file_filtered);

    // remove the index file created for pileup
    //   utility::remove_file(reference_file + ".fai");

    // remove the temporary output file
    utility::remove_file(tmp_file);
    // close output streams
    ofs_snp.close();
    ofs_insertion.close();
    ofs_deletion.close();
    ofs_accepted_snp.close();
    ofs_accepted_insertion.close();
    ofs_accepted_deletion.close();

    // set the variables if required
    if (snp_count != NULL)
        *snp_count = local_snp_count;
    if (snp_count_n != NULL)
        *snp_count_n = local_snp_count_n;
    if (insertion_count != NULL)
        *insertion_count = local_insertion_count;
    if (deletion_count != NULL)
        *deletion_count = local_deletion_count;
    if (accepted_snp_count != NULL)
        *accepted_snp_count = local_accepted_snp_count;
    if (accepted_insertion_count != NULL)
        *accepted_insertion_count = local_accepted_insertion_count;
    if (accepted_deletion_count != NULL)
        *accepted_deletion_count = local_accepted_deletion_count;

}

// New version (see above) created on 5-Oct-2011 11:50 to use mpileup command in samtools 0.1.18 (r982:295)

void sam::generate_variant_lists_old(const string& reference_file, const string& bam_file, bool report_ambiguous_bases,
        int min_coverage_threshold, int max_coverage_threshold, int snp_quality_threshold, int indel_quality_threshold,
        double indel_consensus_threshold, double error_dependency_coefficient, const string& out_pileup_file,
        const string& out_snp_file, const string& out_insertion_file, const string& out_deletion_file,
        const string& out_accepted_snp_file, const string& out_accepted_insertion_file, const string& out_accepted_deletion_file,
        int* snp_count, int* snp_count_n, int* insertion_count, int* deletion_count, int* accepted_snp_count, int* accepted_insertion_count, int* accepted_deletion_count) {

    // file for output by samtools which we don't need 
    string tmp_file = this->tmp_dir + utility::extract_filename(this->sam_file) + "_generate_variant_lists.tmp";

    string command;
    // create the index file
    command = "samtools faidx " + reference_file + " 2> " + tmp_file;
    system(command.c_str());

    // create the pileup file
    string optional_parameters = "";
    if (error_dependency_coefficient >= 0.0) {
        optional_parameters += " -T " + utility::to_string(error_dependency_coefficient) + " ";
    }
    command = "samtools pileup " + optional_parameters + " -cBf " + reference_file + " " + bam_file + " > " + out_pileup_file;
    system(command.c_str());


    string pileup_file_no_n = out_pileup_file + ".noN";
    string pileup_file_filtered = out_pileup_file + ".filtered";

    vector< string > snp_list_n;

    // open the output file stream
    ofstream ofs_pileup(pileup_file_no_n.c_str());

    // open the input file stream
    ifstream ifs_pileup(out_pileup_file.c_str());
    // get the list of positions with Ns and write the positions without Ns to the output file
    // also remove the indels caused by very few (misaligned) reads
    string str_pileup = "";
    while (ifs_pileup.good()) {
        getline(ifs_pileup, str_pileup);

        if (str_pileup.empty()) { // ignore empty lines
            continue;
        }

        vector<string> entries;
        utility::str_split(str_pileup, entries, "\t");

        if (entries[2].compare("N") == 0 && entries[3].compare("N") != 0) {
            if (atoi(entries[7].c_str()) >= min_coverage_threshold) {
                if (!report_ambiguous_bases && !dna::is_nucleotide(entries[3])) {
                    continue;
                } else {
                    snp_list_n.push_back(utility::to_string(entries, 0, 7));
                }
            }
            //            snp_list_n.push_back(utility::to_string(entries, 0, 7));

        } else if (entries[2].compare("*") && atoi(entries[5].c_str()) == 0) { // this is a false indel
            continue;
        } else {
            ofs_pileup << str_pileup << endl;
        }
    } // end while ifsPileup.good()

    ifs_pileup.close();
    ofs_pileup.close();



    // filter the pileup file to remove error-prone variant calls caused by factors not considered in the statistical model of samtools.
    // NOTE: we use "-N 3 -W 1" so that dense SNPs are not filtered out. We use -w 5 to filter out snps within 5 bp of an indel
    command = "samtools.pl varFilter -d " + utility::to_string(min_coverage_threshold) + " -D " + utility::to_string(max_coverage_threshold) + " -w 5 -N 3 -W 1 " + pileup_file_no_n + " > " + pileup_file_filtered;
    system(command.c_str());


    // open output streams for variations
    ofstream ofs_snp(out_snp_file.c_str());
    ofstream ofs_accepted_snp(out_accepted_snp_file.c_str());
    ofstream ofs_insertion(out_insertion_file.c_str());
    ofstream ofs_deletion(out_deletion_file.c_str());
    ofstream ofs_accepted_insertion(out_accepted_insertion_file.c_str());
    ofstream ofs_accepted_deletion(out_accepted_deletion_file.c_str());

    // initialise
    int local_snp_count = 0;
    int local_snp_count_n = 0;
    int local_insertion_count = 0;
    int local_deletion_count = 0;
    int local_accepted_snp_count = 0;
    int local_accepted_insertion_count = 0;
    int local_accepted_deletion_count = 0;

    // open the input file stream
    ifstream ifs_pileup_filtered(pileup_file_filtered.c_str());
    // get the list of positions with SNPs and Indels
    string str_pileup_filtered = "";
    while (ifs_pileup_filtered.good()) {
        getline(ifs_pileup_filtered, str_pileup_filtered);

        if (str_pileup_filtered.empty()) { // ignore empty lines
            continue;
        }

        vector<string> entries;
        utility::str_split(str_pileup_filtered, entries, "\t");

        if (entries[2].compare("*") == 0) { // indel
            if (entries[8].compare("*") != 0) {
                string str_output = entries[0] + "\t" + entries[1] + "\t" + entries[8].substr(1, entries[8].length()) + "\t" + entries[5] + "\t" + entries[10] + "\t" + entries[7];
                if (entries[8].at(0) == '+') { // insertion
                    ofs_insertion << str_output << endl;
                    local_insertion_count++;
                    if (atoi(entries[5].c_str()) >= indel_quality_threshold && atoi(entries[10].c_str()) >= min_coverage_threshold && atof(entries[10].c_str()) / atof(entries[7].c_str()) >= indel_consensus_threshold) {
                        ofs_accepted_insertion << str_output << endl;
                        local_accepted_insertion_count++;
                    }
                } else { // deletion
                    ofs_deletion << str_output << endl;
                    local_deletion_count++;
                    if (atoi(entries[5].c_str()) >= indel_quality_threshold && atoi(entries[10].c_str()) >= min_coverage_threshold && atof(entries[10].c_str()) / atof(entries[7].c_str()) >= indel_consensus_threshold) {
                        ofs_accepted_deletion << str_output << endl;
                        local_accepted_deletion_count++;
                    }
                }
            } else if (entries[9].compare("*") != 0) {
                string str_output = entries[0] + "\t" + entries[1] + "\t" + entries[9].substr(1, entries[9].length()) + "\t" + entries[5] + "\t" + entries[11] + "\t" + entries[7];
                if (entries[9].at(0) == '+') { // insertion
                    ofs_insertion << str_output << endl;
                    local_insertion_count++;
                    if (atoi(entries[5].c_str()) >= indel_quality_threshold && atoi(entries[11].c_str()) >= min_coverage_threshold && atof(entries[11].c_str()) / atof(entries[7].c_str()) >= indel_consensus_threshold) {
                        ofs_accepted_insertion << str_output << endl;
                        local_accepted_insertion_count++;
                    }
                } else { // deletion
                    ofs_deletion << str_output << endl;
                    local_deletion_count++;
                    if (atoi(entries[5].c_str()) >= indel_quality_threshold && atoi(entries[11].c_str()) >= min_coverage_threshold && atof(entries[11].c_str()) / atof(entries[7].c_str()) >= indel_consensus_threshold) {
                        ofs_accepted_deletion << str_output << endl;
                        local_accepted_deletion_count++;
                    }
                }
            } else {
                // ignore
            }

            //            // check the indel
            //            if (atoi(entries[5].c_str()) >= indel_quality_threshold) {
            //                if (entries[8].compare("*") != 0 && atoi(entries[10].c_str()) >= min_coverage_threshold && atof(entries[10].c_str()) / atof(entries[7].c_str()) > indel_consensus_threshold) {
            //                    string str_output = entries[0] + "\t" + entries[1] + "\t" + entries[8] + "\t" + entries[5] + "\t" + entries[10] + "\t" + entries[7];
            //                    if (entries[8].at(0) == '+') { // insertion
            //                        ofs_accepted_insertion << str_output << endl;
            //                    } else { // deletion
            //                        ofs_accepted_deletion << str_output << endl;
            //                    }
            //                } else if (entries[9].compare("*") != 0 && atoi(entries[11].c_str()) >= min_coverage_threshold && atof(entries[11].c_str()) / atof(entries[7].c_str()) > indel_consensus_threshold) {
            //                    string str_output = entries[0] + "\t" + entries[1] + "\t" + entries[9] + "\t" + entries[5] + "\t" + entries[11] + "\t" + entries[7];
            //                    if (entries[9].at(0) == '+') { // insertion
            //                        ofs_accepted_insertion << str_output << endl;
            //                    } else { // deletion
            //                        ofs_accepted_deletion << str_output << endl;
            //                    }
            //                }
            //            }

        } else if (entries[2].compare(entries[3]) != 0 && entries[3].compare("N") != 0) { // different base than the reference and consensus is not N
            string str_output = utility::to_string(entries, 0, 7);
            // add position to the snp list
            ofs_snp << str_output << endl;
            local_snp_count++;
            if (atoi(entries[5].c_str()) >= snp_quality_threshold) { // check if the SNP quality is meets the threshold
                if (!report_ambiguous_bases && !dna::is_nucleotide(entries[3])) { // ambiguous bases to be reported?
                    // no and this is an ambiguous base -> do nothing
                    continue;
                } else {
                    // this position meets the criteria. add it to the snp list
                    ofs_accepted_snp << str_output << endl;
                    local_accepted_snp_count++;
                }
            }
        }
    } // end while ifsPileupFiltered.good()
    ifs_pileup_filtered.close();

    // add positions with Ns in the reference to the output
    vector<string>::iterator the_snp = snp_list_n.begin();
    for (; the_snp != snp_list_n.end(); the_snp++) {
        ofs_snp << *the_snp << endl;
        ofs_accepted_snp << *the_snp << endl;
    }
    // set and increase the SNP count for 'N'
    local_snp_count_n = snp_list_n.size();
    local_snp_count += snp_list_n.size();
    local_accepted_snp_count += snp_list_n.size();

    // remove intermediate files
    utility::remove_file(pileup_file_filtered);
    utility::remove_file(pileup_file_no_n);

    // remove the index file created for pileup
    //   utility::remove_file(reference_file + ".fai");

    // remove the temporary output file
    utility::remove_file(tmp_file);
    // close output streams
    ofs_snp.close();
    ofs_insertion.close();
    ofs_deletion.close();
    ofs_accepted_snp.close();
    ofs_accepted_insertion.close();
    ofs_accepted_deletion.close();

    // set the variables if required
    if (snp_count != NULL)
        *snp_count = local_snp_count;
    if (snp_count_n != NULL)
        *snp_count_n = local_snp_count_n;
    if (insertion_count != NULL)
        *insertion_count = local_insertion_count;
    if (deletion_count != NULL)
        *deletion_count = local_deletion_count;
    if (accepted_snp_count != NULL)
        *accepted_snp_count = local_accepted_snp_count;
    if (accepted_insertion_count != NULL)
        *accepted_insertion_count = local_accepted_insertion_count;
    if (accepted_deletion_count != NULL)
        *accepted_deletion_count = local_accepted_deletion_count;

}

void sam::generate_pileup_file(const string& reference_file, const string& bam_file, const string& tmp_file, const string& out_pileup_file, int max_coverage_threshold) {
    // generate the mpileup file
    // step 1
    string raw_bcf_file = out_pileup_file + ".raw.bcf";
    string str_max_coverage_threshold = utility::to_string(max_coverage_threshold);
    string command = "samtools mpileup -uBA -d " + str_max_coverage_threshold + " -L " + str_max_coverage_threshold + " -f " + reference_file + " " + bam_file + " 2> " + tmp_file + " | bcftools view -bcg - > " + raw_bcf_file + " 2> " + tmp_file;
    system(command.c_str());
    // step 2
    command = "bcftools view " + raw_bcf_file + " > " + out_pileup_file;
    system(command.c_str());

    // remove intermediate file
    utility::remove_file(raw_bcf_file);

}

int sam::get_read_end(int read_start, string cigar) {

    int sum = 0;
    string length = "";
    for (int i = 0; i < cigar.length(); i++) {
        if (isdigit(cigar[i])) {
            length += cigar[i];
        } else if (cigar[i] == 'S' || cigar[i] == 'H') {
            // ignore
            length = "";
        } else if (cigar[i] == 'D') {
            sum += atol(length.c_str());
            length = "";
        } else if (cigar[i] == 'I') {
            length = "";
        } else if (cigar[i] == 'M' || cigar[i] == 'N') {
            sum += atol(length.c_str());
            length = "";
        }
    } // end for

    return read_start + sum - 1;
}

void sam::extract_read_names(const string& sam_file, set<string>* read_names) {

    extract_read_names(sam_file, NULL, false, read_names);

}

void sam::extract_read_names(const string& sam_file, gff* gff_ptr, bool modified_genes_only, set<string>* read_names) {

    if (read_names == NULL) {
        read_names = new set<string>;
    }

    // open the input file stream
    ifstream ifs_sam(sam_file.c_str());

    string str_sam = "";
    string sequence_name = "";
    string previous_sequence_name = "";
    gff_entry* gff_entry_ptr = NULL;
    gff::SEQ_GFF_LIST* seq_gff_list_ptr = NULL;
    gff::SEQ_GFF_LIST::iterator the_seq_gff_list;
    while (ifs_sam.good()) {
        getline(ifs_sam, str_sam);

        if (str_sam.empty()) { // ignore empty lines
            continue;
        } else if (str_sam[0] == '@') { // ignore header
            continue;
        }

        vector<string> entries;
        utility::str_split(str_sam, entries, "\t");

        if (gff_ptr != NULL) {
            // get the sequence names
            sequence_name = entries[2];
            // read coordinates
            int read_start = atoi(entries[3].c_str());
            int read_end = get_read_end(read_start, entries[5]);

            // get the SNP list and GFF data for this sequence
            bool different_sequence_name;
            if ((different_sequence_name = (sequence_name.compare(previous_sequence_name) != 0)) || read_start > gff_entry_ptr->get_end()) {
                // get the new data for this sequence
                if (different_sequence_name) {
                    // coordinate list for this sequence
                    seq_gff_list_ptr = gff_ptr->get_gff_list(sequence_name);

                    the_seq_gff_list = seq_gff_list_ptr->begin();
                    if (the_seq_gff_list != seq_gff_list_ptr->end()) {
                        gff_entry_ptr = *the_seq_gff_list;
                        the_seq_gff_list++;
                    }

                    previous_sequence_name = sequence_name;
                }
                while (gff_entry_ptr != NULL && read_start > gff_entry_ptr->get_end()) {
                    if (the_seq_gff_list != seq_gff_list_ptr->end()) {
                        gff_entry_ptr = *the_seq_gff_list;
                        the_seq_gff_list++;
                    } else {
                        gff_entry_ptr = NULL;
                    }
                }

                if (gff_entry_ptr == NULL) {
                    // this read lies after the last gene and since the
                    // sam file is sorted all subsequent reads with lie after
                    // the last gene. exit the loop
                    break;
                } else if (read_end < gff_entry_ptr->get_start()) {
                    // this read does not lie within a valid gene
                    continue;
                }
            } else if (read_end < gff_entry_ptr->get_start()) {
                continue;
            }// end if different sequence names or read starts after current gene

            if (modified_genes_only && !gff_entry_ptr->modified()) { // we want only reads lying within the modified genes
                // this one is not a modified gene, ignore it
                continue;
            }
            // save the read name
            read_names->insert(entries[0]);
        } else {
            // save the read name
            read_names->insert(entries[0]);
        }
    } // end while ifs_sam.good

}

int sam::get_read_end_coordinates(int read_start, string cigar) {

    int sum = 0;
    string length = "";
    for (int i = 0; i < cigar.length(); i++) {
        if (isdigit(cigar[i])) {
            length += cigar[i];
        } else if (cigar[i] == 'S' || cigar[i] == 'H') {
            // ignore
            length = "";
        } else if (cigar[i] == 'D') {
            sum += atoi(length.c_str());
            length = "";
        } else if (cigar[i] == 'I') {
            length = "";
        } else if (cigar[i] == 'M' || cigar[i] == 'N') {
            sum += atoi(length.c_str());
            length = "";
        }
    } // end for

    return read_start + sum - 1;
}

void sam::get_optional_fields(vector<string>& entries, map<string, string>& optional_fields) {
    for (int i = OPTIONAL_FIELD_START; i < entries.size(); i++) {
        vector<string> f;
        utility::str_split(entries[i], f, ":");
        optional_fields[f[0]] = f[2];
    }

    //return Fields;
}

bool sam::get_read_snps(const set<snp>& read_region_snp_list, int& read_start, const string& field_md, const string& read_sequence, const string& base_qualities, set< snp >& read_snp_list) {

    bool low_quality_snps = false;
    // returns true if all snp bases are good quality and false if bad quality snps (which are ignored) are present	
    set< snp >::iterator the_read_region_snp_list = read_region_snp_list.begin();
    for (; the_read_region_snp_list != read_region_snp_list.end(); the_read_region_snp_list++) {

        snp the_snp = *the_read_region_snp_list;

        // Get the relative position of this SNP on the read 
        int snp_relative_position = get_relative_position(the_snp.get_position(), read_start, field_md);

        // ignore this base if there's is a deletion in this read at SNP position (SNPRelativePosition == -1)
        if (snp_relative_position == -1) {
            continue;
        }

        // Ignore this base if the sequencing quality is low
        int base_quality = int(base_qualities[snp_relative_position]) - 33; // base qualities are ascii - 33
        if (base_quality < BASE_QUALITY_THRESHOLD) {
            low_quality_snps = true;
            continue;
        }

        // Get the base on the read
        char read_base = read_sequence[snp_relative_position];

        // ignore this base if it's not part of the consensus 
        set<char> bases = (*dna::get_nucleotide_map())[the_snp.get_base()];
        if (bases.find(read_base) == bases.end()) {
            continue;
        }

        // Add the SNP to this read's SNP list
        snp new_snp(the_snp.get_position(), read_base);
        read_snp_list.insert(new_snp);

    }
    return low_quality_snps;

}

int sam::get_relative_position(int position, long read_start, string field_md) {

    bool in_deletion = false;
    string length = "";
    long sum = read_start;
    int relative_position = (int) (position - read_start);
    int gap = 0;
    for (int i = 0; i < field_md.length(); i++) {
        if (isdigit(field_md[i])) {
            if (in_deletion)
                in_deletion = false;
            length += field_md[i];
        } else if (field_md[i] == 'A' || field_md[i] == 'C' || field_md[i] == 'G' || field_md[i] == 'T') {
            if (!in_deletion) {
                sum += atol(length.c_str());
                if (sum >= position)
                    break;
            } else {
                if (sum == position) {
                    return -1;
                }
                gap++;
            }

            sum++; // increment the sum by one to include the current position 
            length = "";
        } else if (field_md[i] == '^') {
            sum += atol(length.c_str());
            if (sum > position)
                break;
            else if (sum == position) {
                return -1;
            }
            length = "";
            in_deletion = true;
        }

    } // end for
    return relative_position - gap;
}

char sam::get_consensus_base(string& str_bases, string separator) {
    vector<string> bases;
    utility::str_split(str_bases, bases, separator);
    set<char> ordered_bases;
    for (vector<string>::iterator it = bases.begin(); it != bases.end(); it++) {
        ordered_bases.insert((*it)[0]);
    }
    return (*dna::get_reverse_nucleotide_map())[utility::to_string(ordered_bases)];

}

void sam::mapping_summary(const string& header_file, const string& pileup_file, const string& gff_file, const string& out_file, int coverage_threshold) {

    map < string, int* > coverage;
    map < string, int > ref_length;

    ifstream ifs_header(header_file.c_str());
    string str_header = "";
    while (ifs_header.good()) {
        getline(ifs_header, str_header);

        if (str_header.empty()) { // ignore empty lines
            continue;
        } else if (str_header.find("@SQ") != string::npos) { // get sequence length from the header
            vector<string> entries;
            utility::str_split(str_header, entries, "\t");

            string reference = entries[1].substr(3, entries[1].length());
            string strLength = entries[2].substr(3, entries[2].length());
            int length = atoi(strLength.c_str());

            // memory initialisation
            int *ref_coverage = new int[length];
            for (int i = 0; i < length; i++) {
                *(ref_coverage + i) = 0;
            }

            ref_length[reference] = length;
            coverage[reference] = ref_coverage;
        } else if (str_header.find("@") == string::npos) {
            break;
        }
    } // end while

    ifs_header.close();

    ifstream ifs_pileup(pileup_file.c_str());
    string str_pileup = "";
    string current_sequence_name = ""; // reference (chromosome) being currently read
    int* ref_coverage;
    while (ifs_pileup.good()) {
        getline(ifs_pileup, str_pileup);

        if (str_pileup.empty()) { // ignore empty lines
            continue;
        } else if (str_pileup[0] == '#') { //ignore header
            continue;
        } else {
            vector<string> entries;
            utility::str_split(str_pileup, entries, "\t");

            string sequence_name = entries[0];
            int pos = atoi(entries[1].c_str());

            if (sequence_name.compare(current_sequence_name) != 0) { // new reference
                ref_coverage = coverage[sequence_name];
                current_sequence_name = sequence_name;
            }

            int coverage_start = entries[7].find("DP=");
            int coverage_end = entries[7].find(";", coverage_start);
            int pos_coverage = atoi(entries[7].substr(coverage_start + 3, coverage_end - coverage_start - 3).c_str());


            *(ref_coverage + pos) = pos_coverage;
        }
    } // end while

    ifs_pileup.close();

    // read gff file
    gff* gff_ptr;
    gff_ptr = new gff(gff_file);
    // open the output stream
    ofstream ofs_out(out_file.c_str());
    ofs_out << gff_entry::header_short() << "\t" << "Mapping Start" << "\t" << "Mapping End" << "\t" << "Covered Positions" << "\t" << "Covered Proportion" << endl;

    map < string, int >::iterator it_ref_length = ref_length.begin();
    for (; it_ref_length != ref_length.end(); it_ref_length++) {
        string sequence_name = it_ref_length->first;
        // get the coverages for this sequence
        ref_coverage = coverage[sequence_name];
        // get the GFF List for this sequence
        gff::SEQ_GFF_LIST* seq_gff_list_ptr = gff_ptr->get_gff_list(sequence_name);
        // get the iterator to the GFF entries in the current sequence
        gff::SEQ_GFF_LIST::iterator the_seq_gff_list = seq_gff_list_ptr->begin();
        for (; the_seq_gff_list != seq_gff_list_ptr->end(); the_seq_gff_list++) {
            gff_entry* gff_entry_ptr = *the_seq_gff_list;

            int gff_start = gff_entry_ptr->get_start() - gff::OFFSET;
            int gff_end = gff_entry_ptr->get_end() + gff::OFFSET;

            int covered_positions = 0;
            int mapping_start = -1;
            int mapping_end = -1;
            for (int i = gff_start; i <= gff_end; i++) {
                int pos_coverage = *(ref_coverage + i);
                if (pos_coverage >= coverage_threshold) {
                    if (mapping_start < 0) { // this is the first non-zero position for this gene
                        mapping_start = i;
                        //mapping_end = i;
                    } else {
                        mapping_end = i;
                    }
                    // increment the number of covered position
                    covered_positions++;
                }
            }

            double covered_proportion = (double) covered_positions / (gff_end - gff_start - gff::OFFSET * 2 + 1);
            ofs_out << gff_entry_ptr->to_string_short() << "\t" << mapping_start << "\t" << mapping_end << "\t" << covered_positions << "\t" << covered_proportion << endl;

        }
    }

    ofs_out.close();
}

void sam::soft_clip_summary(const string& sam_file, const string& out_file, int mapping_quality_threshold, const string& out_sam_file, const string& gff_file) {

    map < string, int* > soft_clip_counts;
    map < string, int > ref_length;

    ifstream ifs_sam(sam_file.c_str());
    ofstream ofs_sam;
    if (!out_sam_file.empty()) {
        ofs_sam.open(out_sam_file.c_str());
    }
    string str_sam = "";
    string current_sequence_name = ""; // reference (chromosome) being currently read
    while (ifs_sam.good()) {
        getline(ifs_sam, str_sam);

        if (str_sam.empty()) { // ignore empty lines
            continue;
        } else if (str_sam[0] == '@') { // header
            if (str_sam.find("@SQ") != string::npos) { // get sequence length from the header
                vector<string> entries;
                utility::str_split(str_sam, entries, "\t");

                string reference = entries[1].substr(3, entries[1].length());
                string strLength = entries[2].substr(3, entries[2].length());
                int length = atoi(strLength.c_str());

                // memory initialisation
                int* ref_soft_clip_counts = new int[length];
                for (int i = 0; i < length; i++) {
                    *(ref_soft_clip_counts + i) = 0;
                }

                ref_length[reference] = length;
                soft_clip_counts[reference] = ref_soft_clip_counts;
            }

            // output the header if output file is specified
            if (!out_sam_file.empty()) {
                ofs_sam << str_sam << endl;
            }

            continue;
        }

        // we will be here when we see an alignment .. split the alignment into individual fields
        vector<string> entries;
        utility::str_split(str_sam, entries, "\t");

        int mapping_qual = atoi(entries[4].c_str());
        if (mapping_qual < mapping_quality_threshold || entries[5].find_first_of('S') == string::npos) { // ignore if low mapping quality or if soft clip is not present
            continue;
        }

        string sequence_name = entries[2];
        int* ref_soft_clip_counts;
        if (sequence_name.compare(current_sequence_name) != 0) { // new reference
            ref_soft_clip_counts = soft_clip_counts[sequence_name];
            current_sequence_name = sequence_name;
        }

        string cigar = entries[5];
        int alignment_start = atoi(entries[3].c_str());
        vector<int> clip_sizes;
        vector<int> clip_positions_read;
        vector<int> clip_positions_genome;
        extract_soft_clips(cigar, alignment_start, clip_sizes, clip_positions_read, clip_positions_genome);

        vector<int>::const_iterator the_genomic_position = clip_positions_genome.begin();
        for (; the_genomic_position != clip_positions_genome.end(); the_genomic_position++) {
            int position = *the_genomic_position;
            *(ref_soft_clip_counts + position) += 1;
        }

        if (!out_sam_file.empty()) {
            ofs_sam << str_sam << endl;
        }

    } // end while

    ifs_sam.close();
    if (!out_sam_file.empty()) {
        ofs_sam.close();
    }

    if (gff_file.empty()) {
        ofstream ofs_out(out_file.c_str());
        map < string, int >::iterator it_ref_length = ref_length.begin();
        for (; it_ref_length != ref_length.end(); it_ref_length++) {
            string sequence_name = it_ref_length->first;
            int length = it_ref_length->second;
            // get the soft clip counts for this sequence
            int* ref_soft_clip_counts = soft_clip_counts[sequence_name];

            for (int i = 0; i < length; i++) {
                int count = *(ref_soft_clip_counts + i);
                if (count > 0) {
                    ofs_out << sequence_name << "\t" << i + 1 << "\t" << count << endl;
                }
            }

        }

        ofs_out.close();
    } else {
        // read gff file
        gff* gff_ptr;
        gff_ptr = new gff(gff_file);
        // open the output stream
        ofstream ofs_out(out_file.c_str());
        ofs_out << gff_entry::header_short() << "\t" << "Left Soft Clip Position" << "\t" << "Read Count (Left)" << "\t" << "Right Soft Clip Position" << "\t" << "Read Count (Right)" << endl;

        map < string, int >::iterator it_ref_length = ref_length.begin();
        for (; it_ref_length != ref_length.end(); it_ref_length++) {
            string sequence_name = it_ref_length->first;
            // get the coverages for this sequence
            int* ref_soft_clip_counts = soft_clip_counts[sequence_name];
            // get the GFF List for this sequence
            gff::SEQ_GFF_LIST* seq_gff_list_ptr = gff_ptr->get_gff_list(sequence_name);
            // get the iterator to the GFF entries in the current sequence
            gff::SEQ_GFF_LIST::iterator the_seq_gff_list = seq_gff_list_ptr->begin();
            for (; the_seq_gff_list != seq_gff_list_ptr->end(); the_seq_gff_list++) {
                gff_entry* gff_entry_ptr = *the_seq_gff_list;

                int gff_start = gff_entry_ptr->get_start() - gff::OFFSET;
                int gff_end = gff_entry_ptr->get_end() + gff::OFFSET;

                int left_sc_read_count = 0;
                int right_sc_read_count = 0;
                int left_sc_position = -1;
                int right_sc_position = -1;
                for (int i = gff_start; i <= gff_end; i++) {
                    int pos_count = *(ref_soft_clip_counts + i);
                    if (pos_count > 0) {
                        if (left_sc_position < 0) { // this is the first non-zero position for this gene
                            left_sc_position = i;
                            left_sc_read_count = pos_count;
                            //mapping_end = i;
                        } else {
                            right_sc_position = i;
                            right_sc_read_count = pos_count;
                        }
                    }
                }

                ofs_out << gff_entry_ptr->to_string_short() << "\t" << left_sc_position << "\t" << left_sc_read_count << "\t" << right_sc_position << "\t" << right_sc_read_count << endl;

            }
        }

        ofs_out.close();
    }
}

void sam::extract_soft_clips(const string& cigar, int alignment_start, vector<int>& clip_sizes, vector<int>& clip_positions_read, vector<int>& clip_positions_genome) {

    int position_genome = 0;
    int position_read = 0;
    string str_length = "";
    bool sc_at_read_start = true; // flag to see if the soft clip is at the beginning of a read or at the end
    for (int i = 0; i < cigar.length(); i++) {
        if (isdigit(cigar[i])) {
            str_length += cigar[i];
        } else if (cigar[i] == 'S') {
            int length = atoi(str_length.c_str());
            clip_sizes.push_back(length);
            clip_positions_read.push_back(position_read);
            if (sc_at_read_start) {
                clip_positions_genome.push_back(alignment_start + position_genome);
            } else {
                clip_positions_genome.push_back(alignment_start - 1 + position_genome);
            }
            position_read += length;
            str_length = "";
        } else if (cigar[i] == 'D') {
            position_genome += atoi(str_length.c_str());
            str_length = "";
        } else if (cigar[i] == 'I') {
            position_read += atoi(str_length.c_str());
            str_length = "";
        } else if (cigar[i] == 'M' || cigar[i] == 'N') {
            int length = atoi(str_length.c_str());
            position_genome += length;
            position_read += length;
            str_length = "";
            sc_at_read_start = false;
        }
    } // end for

}

void sam::identify_soft_clip_at_mapping_boundary(const string& mapping_summary_file, const string& soft_clip_summary_file, const string& out_file, int read_coverage_threshold, double gene_coverage_threshold) {

    //    map < string, int* > coverage;
    //    map < string, int > ref_length;

    ifstream ifs_map(mapping_summary_file.c_str());
    ifstream ifs_soft_clip(soft_clip_summary_file.c_str());
    string str_map = "";
    string str_soft_clip = "";

    // save and ignore the header line
    getline(ifs_map, str_map);
    getline(ifs_soft_clip, str_soft_clip);
    string header_map = str_map;
    string header_soft_clip = str_soft_clip;

    // open the output stream
    ofstream ofs_out(out_file.c_str());
    ofs_out << header_soft_clip << endl;
    while (ifs_map.good()) {
        getline(ifs_map, str_map);
        getline(ifs_soft_clip, str_soft_clip);

        if (str_map.empty()) { // ignore empty lines
            continue;
        }

        vector<string> entries_map;
        utility::str_split(str_map, entries_map, "\t");

        vector<string> entries_soft_clip;
        utility::str_split(str_soft_clip, entries_soft_clip, "\t");

        if (entries_map[1] != entries_soft_clip[1]) {
            cerr << "Incorrect file format";
            exit(1);
        }

        ofs_out << entries_map[0] + "\t" + entries_map[1] + "\t" + entries_map[2];
        if (entries_map[3] == entries_soft_clip[3]) { // left soft clip is at the start of mapping 
            if (atoi(entries_soft_clip[4].c_str()) >= read_coverage_threshold && atof(entries_map[6].c_str()) >= gene_coverage_threshold) { // adequate read coverage and gene coverage
                ofs_out << "\t" << entries_soft_clip[3] << "\t" << entries_soft_clip[4];
            } else {
                ofs_out << "\t" << -1 << "\t" << 0;
            }
        } else {
            ofs_out << "\t" << -1 << "\t" << 0;
        }
        if (entries_map[4] == entries_soft_clip[5]) {// right soft clip is at the end of mapping 
            if (atoi(entries_soft_clip[6].c_str()) >= read_coverage_threshold && atof(entries_map[6].c_str()) >= gene_coverage_threshold) { // adequate read coverage and gene coverage
                ofs_out << "\t" << entries_soft_clip[5] << "\t" << entries_soft_clip[6];
            } else {
                ofs_out << "\t" << -1 << "\t" << 0;
            }
        } else {
            ofs_out << "\t" << -1 << "\t" << 0;
        }
        ofs_out << endl;

    } // end while

    ifs_map.close();
    ifs_soft_clip.close();
    ofs_out.close();
}

void sam::process_soft_clips(const string& sam_file, const string& out_sam_file, const string& soft_clip_file) {

    map < string, set<int>* > left_soft_clip_positions;
    map < string, set<int>* > right_soft_clip_positions;

    set<int>* ref_left_soft_clip_positions; // = new vector<int>;
    set<int>* ref_right_soft_clip_positions; // = new vector<int>;
    string current_sequence_name = ""; // reference (chromosome) being currently read
    string str_soft_clip = "";
    ifstream ifs_soft_clip(soft_clip_file.c_str());
    // ignore the header line
    getline(ifs_soft_clip, str_soft_clip);
    while (ifs_soft_clip.good()) {
        getline(ifs_soft_clip, str_soft_clip);

        if (str_soft_clip.empty()) { // ignore empty lines
            continue;
        }

        vector<string> entries;
        utility::str_split(str_soft_clip, entries, "\t");

        string sequence_name = entries[0];
        if (sequence_name.compare(current_sequence_name) != 0) { // new reference
            // get the left soft clip positions for this reference
            if (left_soft_clip_positions.find(sequence_name) == left_soft_clip_positions.end()) {
                left_soft_clip_positions[sequence_name] = new set<int>();
            }
            ref_left_soft_clip_positions = left_soft_clip_positions[sequence_name];

            // get the right soft clip positions for this reference
            if (right_soft_clip_positions.find(sequence_name) == right_soft_clip_positions.end()) {
                right_soft_clip_positions[sequence_name] = new set<int>();
            }
            ref_right_soft_clip_positions = right_soft_clip_positions[sequence_name];

            // save the sequence name
            current_sequence_name = sequence_name;
        }

        int left_soft_clip_pos = atoi(entries[3].c_str());
        int right_soft_clip_pos = atoi(entries[5].c_str());
        if (left_soft_clip_pos > 0) {
            ref_left_soft_clip_positions->insert(left_soft_clip_pos);
        }
        if (right_soft_clip_pos > 0) {
            ref_right_soft_clip_positions->insert(right_soft_clip_pos);
        }
    } // end while   

    ifstream ifs_sam(sam_file.c_str());
    ofstream ofs_sam(out_sam_file.c_str());
    string str_sam = "";
    current_sequence_name = ""; // reference (chromosome) being currently read
    while (ifs_sam.good()) {
        getline(ifs_sam, str_sam);

        if (str_sam.empty()) { // ignore empty lines
            continue;
        } else if (str_sam[0] == '@') { // header
            ofs_sam << str_sam << endl;
            continue;
        }

        // we will be here when we see an alignment .. split the alignment into individual fields
        vector<string> entries;
        utility::str_split(str_sam, entries, "\t");

        if (entries[5].find_first_of('S') == string::npos) { // soft clip is not present
            // output the alignment as it is
            ofs_sam << str_sam << endl;
            continue;
        }

        // get sequence name
        string sequence_name = entries[2];
        // get alignment start
        int alignment_start = atoi(entries[3].c_str());

        if (sequence_name.compare(current_sequence_name) != 0) { // new reference
            // get the left soft clip positions for this reference
            ref_left_soft_clip_positions = left_soft_clip_positions[sequence_name];
            // get the right soft clip positions for this reference
            ref_right_soft_clip_positions = right_soft_clip_positions[sequence_name];
            // save the sequence name
            current_sequence_name = sequence_name;
        }

        if (ref_left_soft_clip_positions->find(alignment_start) != ref_left_soft_clip_positions->end()) {
            // this alignment contains a soft clip in the beginning 
            process_left_soft_clip(entries);

        } else {
            int alignment_end = get_read_end(alignment_start, entries[5]);
            if (ref_right_soft_clip_positions->find(alignment_end) != ref_right_soft_clip_positions->end()) {
                // this alignment contains a soft clip in the end
                process_right_soft_clip(entries);
            }

        }


        // output the alignment as it is
        ofs_sam << utility::to_string(entries, 0, entries.size() - 1) << endl;

    } // end while

    ifs_sam.close();
    ofs_sam.close();


}

void sam::process_right_soft_clip(vector<string>& entries) {

    string cigar = entries[5];
    int sc_pos = cigar.find_last_of('S');
    int m_pos = cigar.find_last_of('M');

    if (m_pos > sc_pos) {
        // this alignment has a soft clip in the beginning. Ignore
        return;

    }
    // extract the soft clip size
    int sc_size = atoi(cigar.substr(m_pos + 1, sc_pos - m_pos).c_str());

    string str_length = "";
    int m_size_start_pos = 0;
    for (int i = m_pos - 1; i >= 0; i--) {
        if (isdigit(cigar[i])) {
            str_length += cigar[i];
            m_size_start_pos = i;
        } else {
            break;
        }
    }
    int m_size = atoi(cigar.substr(m_size_start_pos, m_pos - m_size_start_pos).c_str());
    // create a new cigar    
    cigar = cigar.substr(0, m_size_start_pos) + utility::to_string(m_size + sc_size) + "M";

    // replace the old cigar with new one
    entries[5] = cigar;
}

void sam::process_left_soft_clip(vector<string>& entries) {

    string cigar = entries[5];
    int sc_pos = cigar.find_first_of('S');
    int m_pos = cigar.find_first_of('M');

    if (sc_pos > m_pos) {
        // this alignment has a soft clip in the end. Ignore
        return;

    }
    // extract the soft clip size
    int sc_size = atoi(cigar.substr(0, sc_pos).c_str());
    // extract the M size
    int m_size = atoi(cigar.substr(sc_pos + 1, m_pos - sc_pos).c_str());

    // create a new cigar    
    cigar = utility::to_string(m_size + sc_size) + "M" + cigar.substr(m_pos + 1);

    // replace the old cigar with new one
    entries[5] = cigar;
    // also update the alignment start position
    entries[3] = utility::to_string(atoi(entries[3].c_str()) - sc_size);
}

// this function is internal -- to be called from aligner 

void sam::process_soft_clips(const string& out_sam_file, const string& out_pileup_file, const string& bam_file, const string& reference_file,
        int min_coverage_threshold, int max_coverage_threshold, double gene_coverage_threshold, const string& accepted_out_sc_summary_file) {

    // file for output by samtools which we don't need 
    string tmp_file = this->tmp_dir + utility::extract_filename(this->sam_file) + "_generate_variant_lists.tmp";

    string command;
    // create the index file
    command = "samtools faidx " + reference_file + " 2> " + tmp_file;
    system(command.c_str());

    // generate the pileup
    generate_pileup_file(reference_file, bam_file, tmp_file, out_pileup_file, max_coverage_threshold);

    string gff_file = reference_file + ".gff";

    // calculate the mapping summary
    string out_mapping_summary_file = this->sam_file + ".map.sum";
    mapping_summary(this->sam_file, out_pileup_file, gff_file, out_mapping_summary_file, min_coverage_threshold);

    // calculate soft clip summary
    string out_sc_summary_file = this->sam_file + ".soft.clip.sum";
    soft_clip_summary(this->sam_file, out_sc_summary_file, this->mapping_quality_threshold, "", gff_file);
    // identify soft clips that are present at mapping boundaries
    identify_soft_clip_at_mapping_boundary(out_mapping_summary_file, out_sc_summary_file, accepted_out_sc_summary_file, min_coverage_threshold, gene_coverage_threshold);

    // process the soft clips
    process_soft_clips(this->sam_file, out_sam_file, accepted_out_sc_summary_file);

    // clean up
    utility::remove_file(out_mapping_summary_file);
    utility::remove_file(out_sc_summary_file);
    //utility::remove_file(accepted_out_sc_summary_file);
}