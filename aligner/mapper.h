/* 
 * File:   mapper.h
 * Author: aziz
 *
 * Created on May 7, 2011, 5:43 PM
 */

#ifndef MAPPER_H
#define	MAPPER_H

#include <string>
#include "sam.h"

using namespace std;

class mapper {
public:
    //mapper();
    //mapper(int max_insert_size);
    //mapper(string id,  string tmp_dir, int max_insert_size);
    mapper(string id,  string tmp_dir);
    virtual ~mapper();
    
    void set_id(string id);    
    string get_id() const;
    void set_tmp_dir(string tempDirectory);
    string get_tmp_dir() const;
//    void set_max_insert_size(int max_insert_size);
//    int get_max_insert_size() const;

    virtual sam* align_reads(const string& reference_file, const string& left_reads_file, const string& right_reads_file, const string& outAlignedSAMFile, const string& outunAlignedSAMFile) = 0;
    virtual sam* align_reads(const string& reference_file, const string& left_reads_file, const string& right_reads_file, const string& outAlignedSAMFile, const string& outunAlignedSAMFile, bool mem) = 0;
protected:
    string id;
//    int max_insert_size;
    string tmp_dir;
    bool b_keep_intermediary_files;

};

#endif	/* MAPPER_H */

