/* 
 * File:   sam.h
 * Author: Aziz Mithani
 *
 * Created on May 7, 2011, 4:55 PM
 */

#ifndef SAM_H
#define	SAM_H

#include <string>
#include <set>
#include <vector>
#include <map>
#include "base_distribution.h"
#include "dna.h"
#include "snp.h"

using namespace std;

class sam {
public:
    //sam();
    sam(const string& sam_file);
    sam(const string& sam_file, const string& tmp_dir, bool paired_data, const string& unique_mapping_tag, const string& repeat_mapping_tag);
    virtual ~sam();

    // getter/setter
    void initialise();
    void set_tmp_dir(string tmp_dir);
    string get_tmp_dir() const;
    void set_base_distribution(base_distribution* the_base_distribution);
    base_distribution* get_base_distribution() const;
    void set_sequence_lengths(dna::SEQUENCE_LENGTH_LIST* sequence_lengths);
    dna::SEQUENCE_LENGTH_LIST* get_sequence_lengths() const;
    void set_mapping_quality_threshold(int mapping_quality_threshold);
    int get_mapping_quality_threshold() const;
    void set_ignore_unpaired_reads(bool b_ignore_unpaired_reads);
    bool ignore_unpaired_reads() const;
    void set_unique_mapping_tag(string unique_mapping_tag);
    string get_unique_mapping_tag() const;
    void set_repeat_mapping_tag(string unique_mapping_tag);
    string get_repeat_mapping_tag() const;
    void set_ignore_non_uniquely_mapped_reads(bool b_ignore_non_uniquely_mapped_reads);
    bool ignore_non_uniquely_mapped_reads() const;
    void set_paired_data(bool b_paired_data);
    bool paired_data() const;
    void set_sam_file(string sam_file);
    string get_sam_file() const;
    void set_read_count(int read_count);
    int get_read_count() const;

    void read_sequence_lengths();
    static int extract_aligned_and_unaligned_reads(const string& raw_sam_file, const string& out_aligned_reads_file, const string& out_unaligned_reads_file);
    int filter(const string& out_filtered_sam_file);
    int filter(const string& out_filtered_sam_file, int max_insert_size, const string& gff_file);
    int filter_new(const string& out_filtered_sam_file, int max_insert_size, const string& gff_file);
    //    int filter(const string& out_filtered_sam_file, int max_insert_size, const string& gff_file, bool require_both_unique_reads);
    void calculate_base_distribution(const string& out_base_distribution_file);
    void to_bam(const string& out_bam_file, bool full_run);
    void generate_variant_lists(const string& reference_file, const string& bam_file, bool report_ambiguous_bases,
            int min_coverage_threshold, int max_coverage_threshold, int snp_quality_threshold, int indel_quality_threshold,
            const string& out_pileup_file, const string& out_snp_file, const string& out_insertion_file, const string& out_deletion_file,
            const string& out_accepted_snp_file, const string& out_accepted_insertion_file, const string& out_accepted_deletion_file,
            int* snp_count, int* snp_count_n, int* insertion_count, int* deletion_count, int* accepted_snp_count, int* accepted_insertion_count, int* accepted_deletion_count,
            bool old_indel_format, bool generate_pileup);
    void generate_variant_lists_old(const string& reference_file, const string& bam_file, bool report_ambiguous_bases,
            int min_coverage_threshold, int max_coverage_threshold, int snp_quality_threshold, int indel_quality_threshold,
            double indel_consensus_threshold, double error_dependency_coefficient, const string& out_pileup_file,
            const string& out_snp_file, const string& out_insertion_file, const string& out_deletion_file,
            const string& out_accepted_snp_file, const string& out_accepted_insertion_file, const string& out_accepted_deletion_file,
            int* snp_count, int* snp_count_n, int* insertion_count, int* deletion_count, int* accepted_snp_count, int* accepted_insertion_count, int* accepted_deletion_count);
    static void extract_read_names(const string& sam_file, set<string>* read_names);
    static void extract_read_names(const string& sam_file, gff* gff_ptr, bool modified_genes_only, set<string>* read_names);
    static int get_read_end_coordinates(int read_start, string cigar);
    static void get_optional_fields(vector<string>& entries, map<string, string>& optional_fields);
    static bool get_read_snps(const set<snp>& read_region_snp_list, int& read_start, const string& field_md, const string& read_sequence, const string& base_qualities, set< snp >& read_snp_list);
    static int get_relative_position(int position, long read_start, string field_md);
    static char get_consensus_base(string& str_bases, string separator);
    static void mapping_summary(const string& header_file, const string& pileup_file, const string& gff_file, const string& out_file, int coverage_threshold);
    static void soft_clip_summary(const string& sam_file, const string& out_file, int mapping_quality_threshold, const string& out_sam_file, const string& gff_file);
    static void identify_soft_clip_at_mapping_boundary(const string& mapping_summary_file, const string& soft_clip_summary_file, const string& out_file, int read_coverage_threshold, double gene_coverage_threshold);
    static void process_soft_clips(const string& sam_file, const string& out_sam_file, const string& soft_clip_file);
    void process_soft_clips(const string& out_sam_file, const string& out_pileup_file, const string& bam_file, const string& reference_file, int min_coverage_threshold,
            int max_coverage_threshold, double gene_coverage_threshold, const string& accepted_out_sc_summary_file);
    static void generate_pileup_file(const string& reference_file, const string& bam_file, const string& tmp_file, const string& out_pileup_file, int max_coverage_threshold);

    static const int BASE_QUALITY_THRESHOLD = 20;
    static const int FLAG_MATE_UNMAPPED = 8;

private:

    string path_sam_tools;
    string path_sam_tools_pl;
    string tmp_dir;

    int read_count;
    string sam_file;
    bool b_paired_data;
    bool b_ignore_non_uniquely_mapped_reads;
    string unique_mapping_tag;
    string repeat_mapping_tag;
    bool b_ignore_unpaired_reads;
    int mapping_quality_threshold;
    dna::SEQUENCE_LENGTH_LIST* sequence_lengths;
    base_distribution* the_base_distribution;

    bool is_high_quality(vector<string>& entries);
    bool is_unique(vector<string>& entries);
    bool is_repeat(vector<string>& entries);
    bool is_mapped(vector<string>& entries);
    bool is_supplementary_alignment(vector<string>& entries);
    bool is_correct_orientation(vector<string>& entries_1, vector<string>& entries_2);
    static int get_read_end(int read_start, string cigar);
    static void extract_soft_clips(const string& cigar, int alignment_start, vector<int>& clip_sizes, vector<int>& clip_positions_read, vector<int>& clip_positions_genome);
    static void process_right_soft_clip(vector<string>& entries);
    static void process_left_soft_clip(vector<string>& entries);

    // constants
    static const int FLAG_SEQUENCE_UNMAPPED = 4;
    static const int FLAG_SUPPLEMENTARY_ALIGNMENT = 2048;
    static const int BASE_DISTRIBUTION_COLUMNS = 5;
    static const int OPTIONAL_FIELD_START = 11;
    static const int FLAG_SEQUENCE_ON_REVERSE_STRAND = 16;
    
    //static const int FLAG_MATE_ON_REVERSE_STRAND = 32;
    
};

#endif	/* SAM_H */

