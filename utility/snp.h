/* 
 * File:   snp.h
 * Author: aziz
 *
 * Created on May 24, 2011, 10:57 AM
 */

#ifndef SNP_H
#define	SNP_H

#include <string>
#include "dna.h"
#include "base_distribution.h"

using namespace std;

class snp {
public:
    snp();
    snp(int position, char base);
    snp(const snp& orig);
    virtual ~snp();
    void setBase(char base);
    char get_base() const;
    void set_position(int position);
    int get_position() const;
    bool operator==(const snp& the_snp) const ; // compare the two class objects
    bool operator<(const snp& the_snp) const ; // compare the two class objects
    bool operator>(const snp& the_snp) const ; // compare the two class objects

    static void extract_snps(const string& snp_file, const string& coordinate_file, const string& out_file, bool ignore_Ns_in_the_reference);
    static void check_polyploid_snp_list(dna::BASE_LIST* polyploid_snp_list, base_distribution* polyploid_base_distribution_list);
    static void check_diploid_snp_list(dna::BASE_LIST* diploid_snp_list, base_distribution* diploid_base_distribution_list, dna::BASE_LIST* polyploid_snp_list, dna::BASE_LIST* reference_base_list, dna::BASE_LIST* polyploid_reference_base_list);
    static void check_diploid_snp_list(dna::BASE_LIST* diploid_snp_list, base_distribution* diploid_base_distribution_list, dna::BASE_LIST* diploid_reference_base_list);
    static bool check_polyploid_base(int base_coverage, int total_coverage);
    static bool check_polyploid_base_strict(int base_coverage, int total_coverage);
    static bool check_diploid_base(int base_coverage, int total_coverage);
    static void identify_missing_polyploid_snps(dna::BASE_LIST* polyploid_snp_list, base_distribution* polyploid_base_distribution_list, dna::BASE_LIST* control_snp_list, dna::BASE_LIST* polyploid_reference_base_list, dna::BASE_LIST* control_reference_base_list);
    static void identify_missing_polyploid_snps(dna::BASE_LIST* polyploid_snp_list, base_distribution* polyploid_base_distribution_list, dna::BASE_LIST* polyploid_reference_base_list, dna::BASE_LIST* control_reference_base_list);
    string to_string();
    static void snps_per_gene(const string& snp_file, const string& gff_or_coordinate_file, const string& out_file, bool ignore_Ns_in_the_reference, bool is_gff_file, bool is_snp_table, bool count_snp_positions);
    static void get_consensus_base(const string& snp_file, const string& out_file);
    static void get_consensus_base_GATK(const string& snp_file, const string& out_file);
    static void compare_snps(const string& snp_table_file, const string& polyploid_snp_file);
    
private:
    static const int POLYPLOID_COVERAGE_THRESHOLD = 5;
    static const int POLYPLOID_COVERAGE_STRICT_THRESHOLD = 10;
    static const int DIPLOID_COVERAGE_THRESHOLD = 5;
    static const int POLYPLOID_BASE_COVERAGE_THRESHOLD = 3;
    static const int POLYPLOID_BASE_COVERAGE_STRICT_THRESHOLD = 5;
    static const int DIPLOID_BASE_COVERAGE_THRESHOLD = 5;
    static const double POLYPLOID_BASE_PROPORTION_THRESHOLD = 0.05;
    static const double DIPLOID_BASE_PROPORTION_THRESHOLD = 0.30;
    
    int position;
    char base;

};

#endif	/* SNP_H */

