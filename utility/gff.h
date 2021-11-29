/* 
 * File:   gff.h
 * Author: Aziz Mithani
 *
 * Created on May 10, 2011, 11:39 AM
 */

#ifndef GFF_H
#define	GFF_H

#include <map>
#include <set>
#include <string>
#include "gff_entry.h"

using namespace std;

class gff {
public:
    gff();
    gff(const gff& old_gff);
    gff(const string& gff_file);
    virtual ~gff();

    // type definition
    typedef vector<gff_entry*> SEQ_GFF_LIST;
    typedef map<string, SEQ_GFF_LIST* > GFF_LIST;

    // constants
    static const int OFFSET = 20;
    static const string FEATURE_EXON;
    static const string FEATURE_GENE;

    void read_gff_file(const string& gff_file);
    SEQ_GFF_LIST* get_gff_list(const string& sequence_name);
    void set_gff_list(const string& sequence_name, SEQ_GFF_LIST* seq_gff_list_ptr);
    static gff_entry* get_gff_entry(SEQ_GFF_LIST* seq_gff_list_ptr, int position);
    static string create_gff_entry(const string& fasta_name, int sequence_start, int sequence_end, const string& sequence_name, const string& fasta_tag);
    void write(const string& out_file);
    void write(const string& out_file, bool valid_only);
    void reset(); // reset the valid flag to true in all gff entries
    void filter_genes_with_problematic_snps_in_the_table(const string& snp_file, const string& out_file);
    void filter_genes_with_problematic_snps(const string& snp_file, const string& out_file);
    void read_gene_coordinate_file(const string& coordinate_file);
    bool filter_features(set<string>& features_to_keep);
    void get_sequence_names(set<string>& sequence_names);
    void create_index_array(map<string, int*>&  index_array);
    void delete_index_array(map<string, int*>&  index_array);
    
private:
    GFF_LIST gff_list;
};

#endif	/* GFF_H */

