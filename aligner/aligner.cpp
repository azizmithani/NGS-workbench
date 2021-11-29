/* 
 * File:   aligner.cpp
 * Author: Aziz Mithani
 * 
 * Created on 08 May 2011, 08:01
 */

#include <iostream>
#include "aligner.h"
#include "mapper.h"
#include "utility.h"

aligner::aligner(string id, string tmp_dir, bool b_keep_intermediary_files) {
    this->id = id;
    this->tmp_dir = utility::check_dir(tmp_dir);
    this->b_keep_intermediary_files = b_keep_intermediary_files;

    // Check if the temporary directory already exists
    if (!utility::dir_exists(this->tmp_dir)) {
        // if not create it
        utility::make_dir(this->tmp_dir);
    }
    // set filenames
    string file_prefix = this->tmp_dir + this->id;
    this->aligned_sam_file = file_prefix + "_aligned.sam";
    this->unaligned_sam_file = file_prefix + "_unaligned.sam";
    this->filtered_sam_file = file_prefix + ".sam";
    this->base_distribution_file = file_prefix + ".base.dist";
    this->bam_file = file_prefix + ".bam";
    this->pileup_file = file_prefix + ".pileup";
    this->snp_file = file_prefix + ".snp";
    this->indel_file = file_prefix + ".indel";
    this->insertion_file = file_prefix + ".ins";
    this->deletion_file = file_prefix + ".del";
    this->accepted_snp_file = file_prefix + ".snp.accepted";
    this->accepted_insertion_file = file_prefix + ".ins.accepted";
    this->accepted_deletion_file = file_prefix + ".del.accepted";
    this->accepted_soft_clip_file = file_prefix +  + ".soft.clip.accepted";

    // defaults
    this->b_report_ambiguous_bases = false;
    this->b_ignore_unpaired_reads = false;
    this->min_coverage_threshold = 5;
    this->max_coverage_threshold = 100000;
    this->snp_quality_threshold = 20;
    this->indel_quality_threshold = 100;
    this->gene_coverage_threshold = 0.60;
}

aligner::~aligner() {
}

string aligner::get_id() const {
    return this->id;
}

void aligner::set_id(string id) {
    this->id = id;
}

void aligner::set_tmp_dir(string tmp_dir) {
    this->tmp_dir = tmp_dir;
}

string aligner::get_tmp_dir() const {
    return tmp_dir;
}

void aligner::set_mapper(mapper* the_mapper) {
    this->the_mapper = the_mapper;

    if (this->the_mapper->get_id().empty()) {
        this->the_mapper->set_id(this->id);
    }
    if (this->the_mapper->get_tmp_dir().empty()) {
        this->the_mapper->set_tmp_dir(this->tmp_dir);
    }
}

mapper* aligner::get_mapper() const {
    return this->the_mapper;
}

void aligner::set_filtered_sam_file(string filtered_sam_file) {
    this->filtered_sam_file = filtered_sam_file;
}

string aligner::get_filtered_sam_file() const {
    return filtered_sam_file;
}

void aligner::set_unaligned_sam_file(string unaligned_sam_file) {
    this->unaligned_sam_file = unaligned_sam_file;
}

string aligner::get_unaligned_sam_file() const {
    return unaligned_sam_file;
}

void aligner::set_aligned_sam_file(string aligned_sam_file) {
    this->aligned_sam_file = aligned_sam_file;
}

string aligner::get_aligned_sam_file() const {
    return aligned_sam_file;
}

void aligner::set_base_distribution_file(string base_distribution_file) {
    this->base_distribution_file = base_distribution_file;
}

string aligner::get_base_distribution_file() const {
    return base_distribution_file;
}

void aligner::set_keep_intermediary_files(bool b_keep_intermediary_files) {
    this->b_keep_intermediary_files = b_keep_intermediary_files;
}

bool aligner::keep_intermediary_files() const {
    return b_keep_intermediary_files;
}

void aligner::set_pileup_file(string pileup_file) {
    this->pileup_file = pileup_file;
}

string aligner::get_pileup_file() const {
    return pileup_file;
}

void aligner::set_gene_coverage_threshold(double gene_coverage_threshold) {
    this->gene_coverage_threshold = gene_coverage_threshold;
}

double aligner::get_gene_coverage_threshold() const {
    return gene_coverage_threshold;
}

void aligner::set_snp_quality_threshold(int snp_quality_threshold) {
    this->snp_quality_threshold = snp_quality_threshold;
}

int aligner::get_snp_quality_threshold() const {
    return snp_quality_threshold;
}

void aligner::set_indel_quality_threshold(int indel_quality_threshold) {
    this->indel_quality_threshold = indel_quality_threshold;
}

int aligner::get_indel_quality_threshold() const {
    return indel_quality_threshold;
}

void aligner::set_max_coverage_threshold(int max_coverage_threshold) {
    this->max_coverage_threshold = max_coverage_threshold;
}

int aligner::get_max_coverage_threshold() const {
    return max_coverage_threshold;
}

void aligner::set_min_coverage_threshold(int min_coverage_threshold) {
    this->min_coverage_threshold = min_coverage_threshold;
}

int aligner::get_min_coverage_threshold() const {
    return min_coverage_threshold;
}

void aligner::set_report_ambiguous_bases(bool b_report_ambiguous_bases) {
    this->b_report_ambiguous_bases = b_report_ambiguous_bases;
}

bool aligner::report_ambiguous_bases() const {
    return b_report_ambiguous_bases;
}

void aligner::set_ignore_unpaired_reads(bool b_ignore_unpaired_reads) {
    this->b_ignore_unpaired_reads = b_ignore_unpaired_reads;
}

bool aligner::ignore_unpaired_reads() const {
    return b_ignore_unpaired_reads;
}

void aligner::set_accepted_deletion_file(string accepted_deletion_file) {
    this->accepted_deletion_file = accepted_deletion_file;
}

string aligner::get_accepted_deletion_file() const {
    return accepted_deletion_file;
}

void aligner::set_accepted_insertion_file(string accepted_insertion_file) {
    this->accepted_insertion_file = accepted_insertion_file;
}

string aligner::get_accepted_insertion_file() const {
    return accepted_insertion_file;
}

void aligner::set_accepted_snp_file(string accepted_snp_file) {
    this->accepted_snp_file = accepted_snp_file;
}

string aligner::get_accepted_snp_file() const {
    return accepted_snp_file;
}

void aligner::set_snp_file(string snp_file) {
    this->snp_file = snp_file;
}

string aligner::get_snp_file() const {
    return snp_file;
}

void aligner::set_log_file(string log_file) {
    this->log_file = log_file;
}

string aligner::get_log_file() const {
    return log_file;
}

void aligner::set_accepted_deletion_count(int accepted_deletion_count) {
    this->accepted_deletion_count = accepted_deletion_count;
}

int aligner::get_accepted_deletion_count() const {
    return accepted_deletion_count;
}

void aligner::set_accepted_insertion_count(int accepted_insertion_count) {
    this->accepted_insertion_count = accepted_insertion_count;
}

int aligner::get_accepted_insertion_count() const {
    return accepted_insertion_count;
}

void aligner::set_accepted_snp_count(int accepted_snp_count) {
    this->accepted_snp_count = accepted_snp_count;
}

int aligner::get_accepted_snp_count() const {
    return accepted_snp_count;
}

void aligner::set_deletion_count(int deletion_count) {
    this->deletion_count = deletion_count;
}

int aligner::get_deletion_count() const {
    return deletion_count;
}

void aligner::set_insertion_count(int insertion_count) {
    this->insertion_count = insertion_count;
}

int aligner::get_insertion_count() const {
    return insertion_count;
}

void aligner::set_indel_count(int indel_count) {
    this->indel_count = indel_count;
}

int aligner::get_indel_count() const {
    return indel_count;
}

void aligner::set_snp_count_n(int snp_count_n) {
    this->snp_count_n = snp_count_n;
}

int aligner::get_snp_count_n() const {
    return snp_count_n;
}

void aligner::set_snp_count(int snp_count) {
    this->snp_count = snp_count;
}

int aligner::get_snp_count() const {
    return snp_count;
}

void aligner::set_deletion_file(string deletion_file) {
    this->deletion_file = deletion_file;
}

string aligner::get_deletion_file() const {
    return deletion_file;
}

void aligner::set_insertion_file(string insertion_file) {
    this->insertion_file = insertion_file;
}

string aligner::get_insertion_file() const {
    return insertion_file;
}

void aligner::set_indel_file(string indel_file) {
    this->indel_file = indel_file;
}

string aligner::get_indel_file() const {
    return indel_file;
}

void aligner::set_filtered_reads(int filtered_reads) {
    this->filtered_reads = filtered_reads;
}

int aligner::get_filtered_reads() const {
    return filtered_reads;
}

void aligner::set_mapped_reads(int mapped_reads) {
    this->mapped_reads = mapped_reads;
}

int aligner::get_mapped_reads() const {
    return mapped_reads;
}

void aligner::set_total_reads(int total_reads) {
    this->total_reads = total_reads;
}

int aligner::get_total_reads() const {
    return total_reads;
}

void aligner::set_max_insert_size(int max_insert_size) {
    this->max_insert_size = max_insert_size;
}

int aligner::get_max_insert_size() const {
    return this->max_insert_size;
}

void aligner::set_accepted_soft_clip_file(string accepted_soft_clip_file) {
    this->accepted_soft_clip_file = accepted_soft_clip_file;
}

string aligner::get_accepted_soft_clip_file() const {
    return this->accepted_soft_clip_file;
}

void aligner::set_accepted_right_soft_clip_count(int accepted_right_soft_clip_count) {
    this->accepted_right_soft_clip_count = accepted_right_soft_clip_count;
}

int aligner::get_accepted_right_soft_clip_count() const {
    return accepted_right_soft_clip_count;
}

void aligner::set_accepted_left_soft_clip_count(int accepted_left_soft_clip_count) {
    this->accepted_left_soft_clip_count = accepted_left_soft_clip_count;
}

int aligner::get_accepted_left_soft_clip_count() const {
    return accepted_left_soft_clip_count;
}


sam* aligner::run(const string& reference_file, const string& left_reads_file, const string& right_reads_file, bool full_run, bool old_indel_format) {
    run(reference_file, left_reads_file, right_reads_file, full_run, old_indel_format, false);
}

sam* aligner::run(const string& reference_file, const string& left_reads_file, const string& right_reads_file, bool full_run, bool old_indel_format, bool process_soft_clips) {
 return run(reference_file, left_reads_file, right_reads_file, full_run, old_indel_format, process_soft_clips, false);
}

sam* aligner::run(const string& reference_file, const string& left_reads_file, const string& right_reads_file, bool full_run, bool old_indel_format, bool process_soft_clips, bool mem) {

    
    cout << "Aligning the reads ..." << endl;
    // align reads the reference
    sam* the_sam = this->the_mapper->align_reads(reference_file, left_reads_file, right_reads_file, aligned_sam_file, unaligned_sam_file, mem);

    // set the mapped read count
    this->mapped_reads = the_sam->get_read_count();
    //    // temporary code
    //    sam* the_sam = new sam(this->aligned_sam_file, true, "XT:A:U");

    // set the appropriate variables in sam object
    the_sam->set_ignore_unpaired_reads(this->b_ignore_unpaired_reads);

    // filter the sam file
    cout << "Filtering alignments ..." << endl;
    //    int filtered_reads_count = the_sam->filter(this->filtered_sam_file, this->max_insert_size, reference_file + ".gff", require_both_unique_reads);
    int filtered_reads_count = the_sam->filter(this->filtered_sam_file, this->max_insert_size, reference_file + ".gff");
    //int filtered_reads_count = the_sam->filter_new(this->filtered_sam_file, this->max_insert_size, reference_file + ".gff");

    // set the filtered read count
    this->filtered_reads = filtered_reads_count;

    if (!this->b_keep_intermediary_files) {
        // compress the unfiltered file        
    }

    // set the filename as filtered_sam_file
    the_sam->set_sam_file(this->filtered_sam_file);

    // convert SAM file to BAM file
    cout << "Sorting alignments and creating BAM file ..." << endl;
    the_sam->to_bam(this->bam_file, full_run);



    // added on 11/9/14 to process soft clips
    if (process_soft_clips) {
        cout << "Processing soft clips ..." << endl;
        string processed_sam_file = the_sam->get_sam_file() + ".sc.processed";
        the_sam->process_soft_clips(processed_sam_file, this->pileup_file, this->bam_file, reference_file, this->min_coverage_threshold, this->max_coverage_threshold, this->gene_coverage_threshold, this->accepted_soft_clip_file);

        // set the new file as the filtered sam file
        utility::copy_file(processed_sam_file, this->filtered_sam_file);
        utility::remove_file(processed_sam_file);

        // convert SAM file to BAM file
        cout << "Sorting alignments and creating BAM file ..." << endl;
        the_sam->to_bam(this->bam_file, full_run);
        
        // Convert the bases into N's after the accepted soft-clip positions. We do this to ensure that we don't get chimeric sequences (particularly when the soft-clip is much before the gene end)
        dna::process_reference_using_soft_clips(reference_file, this->accepted_soft_clip_file, &this->accepted_left_soft_clip_count, &this->accepted_right_soft_clip_count);
    }
    
    if (full_run) {
        // calculate the total read count
        this->total_reads = utility::count_lines(left_reads_file) / 4; // 4 lines per read
        if (!right_reads_file.empty()) {
            this->total_reads *= 2;
        }

        // calculate the base distribution
        cout << "Calculating base distribution ..." << endl;
        the_sam->calculate_base_distribution(this->base_distribution_file);
    }


    // generate variant lists
    cout << "Generating variant lists ..." << endl;
    the_sam->generate_variant_lists(reference_file, this->bam_file, this->b_report_ambiguous_bases,
            this->min_coverage_threshold, this->max_coverage_threshold, this->snp_quality_threshold,
            this->indel_quality_threshold, this->pileup_file, this->snp_file, this->insertion_file, this->deletion_file,
            this->accepted_snp_file, this->accepted_insertion_file, this->accepted_deletion_file,
            &this->snp_count, &this->snp_count_n, &this->insertion_count, &this->deletion_count,
            &this->accepted_snp_count, &this->accepted_insertion_count, &this->accepted_deletion_count, old_indel_format, true);
    //    the_sam->generate_variant_lists_old(reference_file, this->bam_file, this->b_report_ambiguous_bases,
    //            this->min_coverage_threshold, this->max_coverage_threshold, this->snp_quality_threshold,
    //            this->indel_quality_threshold, this->indel_consensus_threshold, this->error_dependency_coefficient,
    //            this->pileup_file, this->snp_file, this->insertion_file, this->deletion_file,
    //            this->accepted_snp_file, this->accepted_insertion_file, this->accepted_deletion_file,
    //            &this->snp_count, &this->snp_count_n, &this->insertion_count, &this->deletion_count,
    //            &this->accepted_snp_count, &this->accepted_insertion_count, &this->accepted_deletion_count);

    if (this->b_report_ambiguous_bases) {
        snp::get_consensus_base(this->accepted_snp_file, this->accepted_snp_file + ".con");
    }
    return the_sam;
}

void aligner::remove_files() {
    utility::remove_file(this->aligned_sam_file);
    utility::remove_file(this->unaligned_sam_file);
    utility::remove_file(this->filtered_sam_file);
    utility::remove_file(this->base_distribution_file);
    utility::remove_file(this->bam_file);
    utility::remove_file(this->bam_file + ".bai");
    utility::remove_file(this->pileup_file);
    utility::remove_file(this->snp_file);
    utility::remove_file(this->indel_file);
    utility::remove_file(this->insertion_file);
    utility::remove_file(this->deletion_file);
    utility::remove_file(this->accepted_snp_file);
    utility::remove_file(this->accepted_insertion_file);
    utility::remove_file(this->accepted_deletion_file);
    utility::remove_file(this->accepted_soft_clip_file);
}

