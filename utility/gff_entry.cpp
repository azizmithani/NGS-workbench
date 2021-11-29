/* 
 * File:   gff_entry.cpp
 * Author: Aziz Mithani
 * 
 * Created on May 10, 2011, 11:29 AM
 */

#include "gff_entry.h"
#include "utility.h"

const string gff_entry::FEATURE_ID_TAG = "ID=";

gff_entry::gff_entry(string sequence_name, int start, int end, string str_gff_entry) {
    this->sequence_name = sequence_name;
    this->start = start;
    this->end = end;
    this->str_gff_entry = str_gff_entry;
    this->gff_entry_modified = false;
    this->gff_entry_valid = true;
    this->feature = "gene";
    this->strand = '.';
}

gff_entry::gff_entry(vector<string>& entries, string str_gff_entry) {
    this->sequence_name = entries[0];
    this->start = atoi(entries[3].c_str());
    this->end = atoi(entries[4].c_str());
    this->str_gff_entry = str_gff_entry;
    this->gff_entry_modified = false;
    this->gff_entry_valid = true;
    this->feature = entries[2];
    this->strand = (entries[6])[0];
}

gff_entry::gff_entry(const gff_entry& old_gff_entry) {
    this->sequence_name = old_gff_entry.sequence_name;
    this->start = old_gff_entry.start;
    this->end = old_gff_entry.end;
    this->str_gff_entry = old_gff_entry.str_gff_entry;
    this->gff_entry_modified = old_gff_entry.gff_entry_modified;
    this->gff_entry_valid = old_gff_entry.gff_entry_valid;
    this->feature = old_gff_entry.feature;
    this->strand = old_gff_entry.strand;
}

gff_entry::~gff_entry() {
}

int gff_entry::get_end() const {
    return end;
}

int gff_entry::get_start() const {
    return start;
}

void gff_entry::set_sequence_name(string sequence_name) {
    this->sequence_name = sequence_name;
}

string gff_entry::get_sequence_name() const {
    return sequence_name;
}

string gff_entry::get_feature() const {
    return feature;
}

char gff_entry::get_strand() const {
    return strand;
}

void gff_entry::set_modified(bool modified) {
    this->gff_entry_modified = modified;
}

bool gff_entry::modified() const {
    return gff_entry_modified;
}

void gff_entry::set_valid(bool gff_entry_valid) {
    this->gff_entry_valid = gff_entry_valid;
}

bool gff_entry::valid() const {
    return gff_entry_valid;
}

string gff_entry::get_gff_entry() const {
    return str_gff_entry;
}

void gff_entry::update_start(int start) {
    update_start(start, true);
}

void gff_entry::update_start(int start, bool set_modified_flag) {
    if (this->start == start) {
        return;
    }
    this->start = start;
    if (set_modified_flag) {
        this->gff_entry_modified = true;
    }
}

void gff_entry::update_end(int end) {
    update_end(end, true);
}

void gff_entry::update_end(int end, bool set_modified_flag) {
    if (this->end == end) {
        return;
    }
    this->end = end;
    if (set_modified_flag) {
        this->gff_entry_modified = true;
    }
}

string gff_entry::get_id() {
    return sequence_name + ":" + feature + ":" + utility::to_string(start) + ":" + utility::to_string(end);
}

string gff_entry::to_string_short() {
    return sequence_name + "\t" + utility::to_string(start) + "\t" + utility::to_string(end);
}

string gff_entry::header_short() {
    // column names should correspond to to_string_short()
    return "Sequence Name\tStart\tEnd";
}

string gff_entry::get_feature_id() {
    int start_pos = str_gff_entry.find(FEATURE_ID_TAG,0);
    if (start_pos == string::npos) {
        return "";
    }
    int end_pos = str_gff_entry.find_first_of(';',start_pos);
    if (end_pos == string::npos) {
        end_pos = str_gff_entry.length();
    }
   
    return str_gff_entry.substr(start_pos + FEATURE_ID_TAG.size(), end_pos - start_pos - FEATURE_ID_TAG.size());
    
}
