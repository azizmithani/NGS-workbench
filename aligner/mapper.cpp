/* 
 * File:   mapper.cpp
 * Author: Aziz Mithani
 *
 * Created on May 7, 2011, 6:30 PM
 */

#include "mapper.h"
#include "utility.h"

//mapper::mapper() {
//}

//mapper::mapper(int max_insert_size){
//    this->id = "";    
//    this->tmp_dir = "";
//    this->max_insert_size = max_insert_size;
//}

//mapper::mapper(string id, string tmp_dir, int max_insert_size){
//    this->id = id;
//    this->max_insert_size = max_insert_size;
//    this->tmp_dir = utility::check_dir(tmp_dir);
//    // Check if the temporary directory already exists
//    if (!utility::dir_exists(this->tmp_dir)) {
//        // if not create it
//        utility::make_dir(this->tmp_dir);
//    }
//    
//}

mapper::mapper(string id, string tmp_dir){
    this->id = id;
    this->tmp_dir = utility::check_dir(tmp_dir);
    // Check if the temporary directory already exists
    if (!utility::dir_exists(this->tmp_dir)) {
        // if not create it
        utility::make_dir(this->tmp_dir);
    }
    
}

mapper::~mapper() {
}

string mapper::get_id() const {
    return this->id;
}

void mapper::set_id(string id) {
    this->id = id;
}

void mapper::set_tmp_dir(string tmp_dir) {
    this->tmp_dir = tmp_dir;
}

string mapper::get_tmp_dir() const {
    return tmp_dir;
}

//void mapper::set_max_insert_size(int max_insert_size) {
//    this->max_insert_size = max_insert_size;
//}
//
//int mapper::get_max_insert_size() const {
//    return this->max_insert_size;
//}
