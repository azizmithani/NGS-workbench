/* 
 * File:   workbench.cpp
 * Author: aziz
 *
 * Created on May 16, 2011, 11:44 AM
 */

#include <cstdlib>
#include <set>
#include <string>
#include <fstream>
#include <iostream>
#include "sam.h"
#include "gff.h"
#include "dna.h"
#include "snp.h"
#include "utility.h"
#include "aligner.h"
#include "mapper.h"
#include "bwa.h"
#include "time.h"
#include <time.h>

using namespace std;

// Function Declarations
int calculate_base_distribution(int argc, char** argv);
int check_diploid_snps(int argc, char** argv);
int compare_snps(int argc, char** argv);
int extract_exonic_sequences(int argc, char** argv);
int extract_subsequence(int argc, char** argv);
int extract_reads_fastq(int argc, char** argv);
int extract_sequences(int argc, char** argv);
int extract_snps(int argc, char** argv);
int fasta_to_multi_fasta(int argc, char** argv);
int filter_genes_with_problematic_snps(int argc, char** argv);
int filter_problematic_snps(int argc, char** argv);
int filter_sam_file(int argc, char** argv);
int filter_sam_file_new(int argc, char** argv);
int generate_variant_lists(int argc, char** argv);
int generate_variant_lists_old(int argc, char** argv);
int get_codons(int argc, char** argv);
int get_codon_frequencies(int argc, char** argv);
int get_consensus_bases(int argc, char** argv);
int identify_missing_polyploid_snps(int argc, char** argv);
int multi_fasta_to_fasta(int argc, char** argv);
int mapping_summary(int argc, char** argv);
int process_soft_clips(int argc, char** argv);
int remove_common_snps(int argc, char** argv);
int remove_transcript_variants_ensembl(int argc, char** argv);
int remove_transcript_variants_trinity(int argc, char** argv);
int remove_uncommon_snps(int argc, char** argv);
int run_bwa(int argc, char** argv);
int snp_table_to_list(int argc, char** argv);
int snps_per_gene(int argc, char** argv);
int soft_clip_summary(int argc, char** argv);
int sort_multi_fasta(int argc, char** argv);
int split_fastq_into_pairs(int argc, char** argv);
int tabulate_sequence_lengths(int argc, char** argv);
int tabulate_snps(int argc, char** argv);
int test(int argc, char** argv);

int main(int argc, char** argv) {

    //    string gff_file(argv[1]);
    //    string sam_file(argv[2]);
    //    string left_reads_file(argv[3]);
    //    string right_reads_file(argv[4]);
    //    string new_left_reads_file(argv[5]);
    //    string new_right_reads_file(argv[6]);
    //
    //    string gff_file = "/Volumes/sann2733/tmp/vtest.gff";
    //    string sam_file = "/Volumes/SAN002/BWA_Analysis/BWA_Data/Wheat/T_aestivum/CS/T_aestivum_1.11_CS_filtered.sam";
    //    string left_reads_file = "/Volumes/SAN001/Sequencing_Data/Wheat/WTCHG/T_aestivum/CS_1.fq";
    //    string right_reads_file = "/Volumes/SAN001/Sequencing_Data/Wheat/WTCHG/T_aestivum/CS_2.fq";
    //    string new_left_reads_file = "/Volumes/sann2733/tmp/vtest_1.fq";
    //    string new_right_reads_file = "/Volumes/sann2733/tmp/vtest_2.fq";


    string command = argv[1];

    if (command.compare("base_distribution") == 0) {
        return calculate_base_distribution(argc, argv);
    } else if (command.compare("bwa") == 0) {
        return run_bwa(argc, argv);
    } else if (command.compare("check_diploid_snps") == 0) {
        return check_diploid_snps(argc, argv);
    } else if (command.compare("codon") == 0) {
        return get_codons(argc, argv);
    } else if (command.compare("codon_frequency") == 0) {
        return get_codon_frequencies(argc, argv);
    } else if (command.compare("compare_snps") == 0) {
        return compare_snps(argc, argv);
    } else if (command.compare("consensus_bases") == 0) {
        return get_consensus_bases(argc, argv);
    } else if (command.compare("extract_reads") == 0) {
        return extract_reads_fastq(argc, argv);
    } else if (command.compare("extract_sequences") == 0) {
        return extract_sequences(argc, argv);
    } else if (command.compare("extract_snps") == 0) {
        return extract_snps(argc, argv);
    } else if (command.compare("fasta_to_multi_fasta") == 0) {
        return fasta_to_multi_fasta(argc, argv);
    } else if (command.compare("filter_genes_with_problematic_snps") == 0) { // unassigned or ambiguous snps(diploids)
        return filter_genes_with_problematic_snps(argc, argv);
    } else if (command.compare("filter_problematic_snps") == 0) { // unassigned or ambiguous snps
        return filter_problematic_snps(argc, argv);
    } else if (command.compare("filter_sam") == 0) { // filter sam file
        return filter_sam_file(argc, argv);
    } else if (command.compare("filter_sam_new") == 0) { // filter sam file
        return filter_sam_file_new(argc, argv);
    } else if (command.compare("identify_missing_polyploid_snps") == 0) { // identify missing polyploid snps
        return identify_missing_polyploid_snps(argc, argv);
    } else if (command.compare("multi_fasta_to_fasta") == 0) {
        return multi_fasta_to_fasta(argc, argv);
    } else if (command.compare("remove_common_snps") == 0) {
        return remove_common_snps(argc, argv);
    } else if (command.compare("remove_uncommon_snps") == 0) {
        return remove_uncommon_snps(argc, argv);
    } else if (command.compare("snp_table_to_list") == 0) {
        return snp_table_to_list(argc, argv);
    } else if (command.compare("snps_per_gene") == 0) {
        return snps_per_gene(argc, argv);
    } else if (command.compare("sort_multi_fasta") == 0) {
        return sort_multi_fasta(argc, argv);
    } else if (command.compare("tabulate_sequence_lengths") == 0) {
        return tabulate_sequence_lengths(argc, argv);
    } else if (command.compare("tabulate_snps") == 0) {
        return tabulate_snps(argc, argv);
    } else if (command.compare("variance") == 0) {
        return generate_variant_lists(argc, argv);
    } else if (command.compare("variance_old") == 0) {
        return generate_variant_lists_old(argc, argv);
    } else if (command.compare("test") == 0) {
        return test(argc, argv);
    } else if (command.compare("exons") == 0) {
        return extract_exonic_sequences(argc, argv);
    } else if (command.compare("subsequence") == 0) {
        return extract_subsequence(argc, argv);
    } else if (command.compare("split_fastq") == 0) {
        return split_fastq_into_pairs(argc, argv);
    } else if (command.compare("mapping_summary") == 0) {
        return mapping_summary(argc, argv);
    } else if (command.compare("soft_clip_summary") == 0) {
        return soft_clip_summary(argc, argv);
    } else if (command.compare("process_soft_clips") == 0) {
        return process_soft_clips(argc, argv);
    } else if (command.compare("remove_transcript_variants_ensembl") == 0) {
        return remove_transcript_variants_ensembl(argc, argv);
    } else if (command.compare("remove_transcript_variants_trinity") == 0) {
        return remove_transcript_variants_trinity(argc, argv);
    } else {
        cout << "Invalid option" << endl;
    }


    return 0;
}

int extract_reads_fastq(int argc, char** argv) {

    string sam_file = "";
    string gff_file = "";
    string reads_file = "";
    string out_file = "";

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-s") == 0) { // sam file
            sam_file = argv[i++];
            continue;
        } else if (parameter.compare("-i") == 0) { // input reads file
            reads_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // output file
            out_file = argv[i++];
            continue;
        } else if (parameter.compare("-g") == 0) { // gff file
            gff_file = argv[i++];
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (sam_file.empty()) {
        cerr << "SAM file must be specified." << endl;
        return 1;
    } else if (gff_file.empty()) {
        cerr << "GFF file must be specified." << endl;
        return 1;
    } else if (reads_file.empty()) {
        cerr << "Reads file must be specified." << endl;
        return 1;
    } else if (out_file.empty()) {
        cerr << "Output file must be specified." << endl;
        return 1;
    }

    // read the GFF file
    gff* gff_ptr = new gff(gff_file);
    //    gff* gff_ptr = new gff();
    //    gff_ptr->read_gff_file(gff_file);

    // extract the read names 
    set<string>* read_names = new set<string>;
    sam::extract_read_names(sam_file, gff_ptr, false, read_names);

    // extract the sequencing reads from the fastq file
    dna::filter_fastq(reads_file, out_file, read_names);

    return 0;
}

int extract_sequences(int argc, char** argv) {

    string multi_fasta_file = "";
    string out_file = "";
    string sequence_names_file = "";

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-i") == 0) { // multifasta file
            multi_fasta_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // output file
            out_file = argv[i++];
            continue;
        } else if (parameter.compare("-n") == 0) { // sequence names file
            sequence_names_file = argv[i++];
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (multi_fasta_file.empty()) {
        cerr << "Multifasta file must be specified." << endl;
        return 1;
    } else if (sequence_names_file.empty()) {
        cerr << "Sequence names file must be specified." << endl;
        return 1;
    }
    if (out_file.empty()) {
        out_file = multi_fasta_file + ".extracted";
    }

    dna::extract_sequences(multi_fasta_file, sequence_names_file, out_file);

    return 0;
}

int filter_sequences(int argc, char** argv) {

    string multi_fasta_file = "";
    string out_file = "";
    string sequence_names_file = "";

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-i") == 0) { // multifasta file
            multi_fasta_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // output file
            out_file = argv[i++];
            continue;
        } else if (parameter.compare("-n") == 0) { // sequence names file
            sequence_names_file = argv[i++];
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (multi_fasta_file.empty()) {
        cerr << "Multifasta file must be specified." << endl;
        return 1;
    } else if (sequence_names_file.empty()) {
        cerr << "Sequence names file must be specified." << endl;
        return 1;
    }
    if (out_file.empty()) {
        out_file = multi_fasta_file + ".extracted";
    }

    dna::filter_sequences(multi_fasta_file, sequence_names_file, out_file);

    return 0;
}

int extract_snps(int argc, char** argv) {

    string snp_file = "";
    string coordinate_file = "";
    string out_file = "";
    bool ignore_Ns_in_the_reference = false;

    //    bool nullitetra_data = false;

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-s") == 0) { // snp file
            snp_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // output file
            out_file = argv[i++];
            continue;
        } else if (parameter.compare("-c") == 0) { // coordinate file
            coordinate_file = argv[i++];
            continue;
        } else if (parameter.compare("-n") == 0) { // ignore Ns in the reference
            string str_ignore_Ns_in_the_reference = argv[i++];
            ignore_Ns_in_the_reference = (str_ignore_Ns_in_the_reference[0] == 'T' || str_ignore_Ns_in_the_reference[0] == 't' ? true : false);
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (snp_file.empty()) {
        cerr << "SNP file must be specified." << endl;
        return 1;
    } else if (coordinate_file.empty()) {
        cerr << "Coordinate file must be specified." << endl;
        return 1;
    }
    if (out_file.empty()) {
        out_file = snp_file + ".filtered";
    }


    snp::extract_snps(snp_file, coordinate_file, out_file, ignore_Ns_in_the_reference);

    return 0;
}

int fasta_to_multi_fasta(int argc, char** argv) {

    string fasta_file = "";
    string out_file = "";

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-i") == 0) { // fasta file
            fasta_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // output file
            out_file = argv[i++];
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (fasta_file.empty()) {
        cerr << "Fasta file must be specified." << endl;
        return 1;
    }
    if (out_file.empty()) {
        out_file = fasta_file + ".mfa";
    }

    dna::fasta_to_multi_fasta(fasta_file, out_file);

    return 0;
}

int multi_fasta_to_fasta(int argc, char** argv) {

    string multi_fasta_file = "";
    string out_file = "";
    string fasta_name = "";
    int gap_size = 200;

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-i") == 0) { // multifasta file
            multi_fasta_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // output file
            out_file = argv[i++];
            continue;
        } else if (parameter.compare("-n") == 0) { // sequence name
            fasta_name = argv[i++];
            continue;
        } else if (parameter.compare("-g") == 0) { // gap size
            gap_size = atoi(argv[i++]);
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (multi_fasta_file.empty()) {
        cerr << "Multifasta file must be specified." << endl;
        return 1;
    }
    if (fasta_name.empty()) {
        fasta_name = "sequence";
    }
    if (out_file.empty()) {
        out_file = multi_fasta_file + ".f-numbera";
    }

    dna::multi_fasta_to_fasta(multi_fasta_file, out_file, fasta_name, gap_size);

    return 0;
}

int generate_variant_lists(int argc, char** argv) {

    string reference_file = "";
    string bam_file = "";
    bool report_ambiguous_bases = false;
    int min_coverage_threshold = 5;
    int max_coverage_threshold = 50000;
    int snp_quality_threshold = 20;
    int indel_quality_threshold = 100;
    double indel_consensus_threshold = 0.70;
    double error_dependency_coefficient = -1.0;
    string out_snp_file = "";
    string out_insertion_file = "";
    string out_deletion_file = "";
    int snp_count = 0;
    int snp_count_n = 0;
    int insertion_count = 0;
    int deletion_count = 0;
    int accepted_snp_count = 0;
    int accepted_insertion_count = 0;
    int accepted_deletion_count = 0;
    bool old_indel_format = false;
    string pileup_file = "";
    bool generate_pileup = true;

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-r") == 0) { // reference file
            reference_file = argv[i++];
            continue;
        } else if (parameter.compare("-b") == 0) { // bam file
            bam_file = argv[i++];
            continue;
        } else if (parameter.compare("-a") == 0) { // report ambiguous base
            string str_report_ambigous_base = argv[i++];
            report_ambiguous_bases = (str_report_ambigous_base[0] == 'T' || str_report_ambigous_base[0] == 't' ? true : false);
            continue;
        } else if (parameter.compare("-d") == 0) { // minimum coverage threshold
            min_coverage_threshold = atoi(argv[i++]);
            continue;
        } else if (parameter.compare("-D") == 0) { // maximum coverage threshold
            max_coverage_threshold = atoi(argv[i++]);
            continue;
        } else if (parameter.compare("-theta") == 0) { // theta for samtools pileup consensus (maq consensus error dependency coefficient)
            error_dependency_coefficient = atof(argv[i++]);
            continue;
        } else if (parameter.compare("-sq") == 0) { // snp quality threshold
            snp_quality_threshold = atoi(argv[i++]);
            continue;
        } else if (parameter.compare("-iq") == 0) { // indel quality threshold
            indel_quality_threshold = atoi(argv[i++]);
            continue;
        } else if (parameter.compare("-ic") == 0) { // indel consensus threshold
            indel_consensus_threshold = atof(argv[i++]);
            continue;
        } else if (parameter.compare("-snp") == 0) { // output file for snps
            out_snp_file = argv[i++];
            continue;
        } else if (parameter.compare("-ins") == 0) { // output file for insertions
            out_insertion_file = argv[i++];
            continue;
        } else if (parameter.compare("-del") == 0) { // output file for deletions
            out_deletion_file = argv[i++];
            continue;
        } else if (parameter.compare("-old") == 0) { // old indel format
            string str_old_indel_format = argv[i++];
            old_indel_format = (str_old_indel_format[0] == 'T' || str_old_indel_format[0] == 't' ? true : false);
            continue;
        } else if (parameter.compare("-p") == 0) { // generate pileup file?
            pileup_file = argv[i++];
            generate_pileup = false;
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    // parameter checking
    if (reference_file.empty()) {
        cerr << "Reference file must be specified." << endl;
        return 1;
    }
    if (bam_file.empty()) {
        cerr << "BAM file must be specified." << endl;
        return 1;
    }
    if (out_snp_file.empty()) {
        out_snp_file = utility::remove_extension(bam_file) + ".snp";
    }
    if (out_insertion_file.empty()) {
        out_insertion_file = utility::remove_extension(bam_file) + ".ins";
    }
    if (out_deletion_file.empty()) {
        out_deletion_file = utility::remove_extension(bam_file) + ".del";
    }
    if (generate_pileup) {
        pileup_file = utility::remove_extension(bam_file) + ".pileup";
    }

    string out_accepted_snp_file = out_snp_file + ".accepted";
    string out_accepted_insertion_file = out_insertion_file + ".accepted";
    string out_accepted_deletion_file = out_deletion_file + ".accepted";


    // create a dummy sam object
    sam* the_sam = new sam(utility::remove_extension(bam_file) + ".sam", utility::extract_dir(bam_file), true, "", "");

    the_sam->generate_variant_lists(reference_file, bam_file, report_ambiguous_bases,
            min_coverage_threshold, max_coverage_threshold, snp_quality_threshold,
            indel_quality_threshold, pileup_file, out_snp_file, out_insertion_file, out_deletion_file,
            out_accepted_snp_file, out_accepted_insertion_file, out_accepted_deletion_file,
            &snp_count, &snp_count_n, &insertion_count, &deletion_count, &accepted_snp_count, &accepted_insertion_count,
            &accepted_deletion_count, old_indel_format, generate_pileup);


    cout << "Substitutions:" << endl;
    cout << "\tReported: " << snp_count << endl;
    cout << "\tAccepted: " << accepted_snp_count << endl;
    cout << "\tN (included): " << snp_count_n << endl;
    cout << "\tInsertions: " << endl;
    cout << "\t\tReported: " << insertion_count << endl;
    cout << "\t\tAccepted: " << accepted_insertion_count << endl;
    cout << "\tDeletions:" << endl;
    cout << "\t\tReported: " << deletion_count << endl;
    cout << "\t\tAccepted: " << accepted_deletion_count << endl;

    // remove unwated files
    //utility::remove_file(pileup_file);

    // free the memory
    delete the_sam;

    return 0;
}

int generate_variant_lists_old(int argc, char** argv) {

    string reference_file = "";
    string bam_file = "";
    bool report_ambiguous_bases = false;
    int min_coverage_threshold = 5;
    int max_coverage_threshold = 50000;
    int snp_quality_threshold = 20;
    int indel_quality_threshold = 100;
    double indel_consensus_threshold = 0.70;
    double error_dependency_coefficient = -1.0;
    string out_snp_file = "";
    string out_insertion_file = "";
    string out_deletion_file = "";
    int snp_count = 0;
    int snp_count_n = 0;
    int insertion_count = 0;
    int deletion_count = 0;
    int accepted_snp_count = 0;
    int accepted_insertion_count = 0;
    int accepted_deletion_count = 0;
    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-r") == 0) { // reference file
            reference_file = argv[i++];
            continue;
        } else if (parameter.compare("-b") == 0) { // bam file
            bam_file = argv[i++];
            continue;
        } else if (parameter.compare("-a") == 0) { // report ambiguous base
            string str_report_ambigous_base = argv[i++];
            report_ambiguous_bases = (str_report_ambigous_base[0] == 'T' || str_report_ambigous_base[0] == 't' ? true : false);
            continue;
        } else if (parameter.compare("-d") == 0) { // minimum coverage threshold
            min_coverage_threshold = atoi(argv[i++]);
            continue;
        } else if (parameter.compare("-D") == 0) { // maximum coverage threshold
            max_coverage_threshold = atoi(argv[i++]);
            continue;
        } else if (parameter.compare("-theta") == 0) { // theta for samtools pileup consensus (maq consensus error dependency coefficient)
            error_dependency_coefficient = atof(argv[i++]);
            continue;
        } else if (parameter.compare("-sq") == 0) { // snp quality threshold
            snp_quality_threshold = atoi(argv[i++]);
            continue;
        } else if (parameter.compare("-iq") == 0) { // indel quality threshold
            indel_quality_threshold = atoi(argv[i++]);
            continue;
        } else if (parameter.compare("-ic") == 0) { // indel consensus threshold
            indel_consensus_threshold = atof(argv[i++]);
            continue;
        } else if (parameter.compare("-snp") == 0) { // output file for snps
            out_snp_file = argv[i++];
            continue;
        } else if (parameter.compare("-ins") == 0) { // output file for insertions
            out_insertion_file = argv[i++];
            continue;
        } else if (parameter.compare("-del") == 0) { // output file for deletions
            out_deletion_file = argv[i++];
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    // parameter checking
    if (reference_file.empty()) {
        cerr << "Reference file must be specified." << endl;
        return 1;
    }
    if (bam_file.empty()) {
        cerr << "BAM file must be specified." << endl;
        return 1;
    }
    if (out_snp_file.empty()) {
        out_snp_file = utility::remove_extension(bam_file) + ".snp";
    }
    if (out_insertion_file.empty()) {
        out_insertion_file = utility::remove_extension(bam_file) + ".ins";
    }
    if (out_deletion_file.empty()) {
        out_deletion_file = utility::remove_extension(bam_file) + ".del";
    }

    string pileup_file = utility::remove_extension(bam_file) + ".pileup";
    string out_accepted_snp_file = out_snp_file + ".accepted";
    string out_accepted_insertion_file = out_insertion_file + ".accepted";
    string out_accepted_deletion_file = out_deletion_file + ".accepted";


    // create a dummy sam object
    sam* the_sam = new sam(utility::remove_extension(bam_file) + ".sam", utility::extract_dir(bam_file), true, "", "");

    the_sam->generate_variant_lists_old(reference_file, bam_file, report_ambiguous_bases,
            min_coverage_threshold, max_coverage_threshold, snp_quality_threshold,
            indel_quality_threshold, indel_consensus_threshold, error_dependency_coefficient,
            pileup_file, out_snp_file, out_insertion_file, out_deletion_file,
            out_accepted_snp_file, out_accepted_insertion_file, out_accepted_deletion_file,
            &snp_count, &snp_count_n, &insertion_count, &deletion_count,
            &accepted_snp_count, &accepted_insertion_count, &accepted_deletion_count);

    cout << "Substitutions:" << endl;
    cout << "\tReported: " << snp_count << endl;
    cout << "\tAccepted: " << accepted_snp_count << endl;
    cout << "\tN (included): " << snp_count_n << endl;
    cout << "Insertions: " << endl;
    cout << "\tReported: " << insertion_count << endl;
    cout << "\tAccepted: " << accepted_insertion_count << endl;
    cout << "Deletions: " << endl;
    cout << "\tReported: " << deletion_count << endl;
    cout << "\tAccepted Deletions: " << accepted_deletion_count << endl;

    // remove unwated files
    utility::remove_file(pileup_file);

    // free the memory
    delete the_sam;

    return 0;
}

int filter_genes_with_problematic_snps(int argc, char** argv) {

    string snp_file = "";
    string gff_file = "";
    string out_file = "";
    bool tabular_data = false;

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-s") == 0) { // snp file
            snp_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // output file
            out_file = argv[i++];
            continue;
        } else if (parameter.compare("-g") == 0) { // gff file
            gff_file = argv[i++];
            continue;
        } else if (parameter.compare("-t") == 0) { // is tabular snp file (snp columns)
            string str_tabular_data = argv[i++];
            tabular_data = (str_tabular_data[0] == 'T' || str_tabular_data[0] == 't' ? true : false);
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (snp_file.empty()) {
        cerr << "SNP file must be specified." << endl;
        return 1;
    } else if (gff_file.empty()) {
        cerr << "GFF file must be specified." << endl;
        return 1;
    }
    if (out_file.empty()) {
        out_file = snp_file + ".filtered";
    }

    // read the GFF file
    gff* gff_ptr = new gff(gff_file);

    if (tabular_data) {
        gff_ptr->filter_genes_with_problematic_snps_in_the_table(snp_file, out_file);
    } else {
        gff_ptr->filter_genes_with_problematic_snps(snp_file, out_file);
    }

    // clean up
    delete gff_ptr;

    return 0;
}

int filter_problematic_snps(int argc, char** argv) {

    string snp_file = "";
    string out_file = "";
    bool allow_extended_bases = false;

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-s") == 0) { // snp file
            snp_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // output file
            out_file = argv[i++];
            continue;
        } else if (parameter.compare("-x") == 0) { // allow extended nucleotide set
            string str_allow_extended_bases = argv[i++];
            allow_extended_bases = (str_allow_extended_bases[0] == 'T' || str_allow_extended_bases[0] == 't' ? true : false);
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (snp_file.empty()) {
        cerr << "SNP file must be specified." << endl;
        return 1;
    }
    if (out_file.empty()) {
        out_file = snp_file + ".filtered";
    }

    dna::filter_problematic_snps(snp_file, out_file, allow_extended_bases);


    return 0;
}

int identify_missing_polyploid_snps(int argc, char** argv) {

    string polyploid_snp_file = "";
    string base_distribution_file = "";
    string control_snp_file = "";
    string reference_file = "";
    string out_file = "";

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-s") == 0) { // snp file
            polyploid_snp_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // output file
            out_file = argv[i++];
            continue;
        } else if (parameter.compare("-c") == 0) { // control snp file
            control_snp_file = argv[i++];
            continue;
        } else if (parameter.compare("-r") == 0) { // reference file (if full snp check is to be performed)
            reference_file = argv[i++];
            continue;
        } else if (parameter.compare("-b") == 0) { // base distribution file
            base_distribution_file = argv[i++];
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (polyploid_snp_file.empty()) {
        cerr << "SNP file must be specified." << endl;
        return 1;
    } else if (base_distribution_file.empty()) {
        cerr << "Base distribution file must be specified." << endl;
        return 1;
    } else if (control_snp_file.empty() && reference_file.empty()) {
        cerr << "Control or reference file must be specified." << endl;
        return 1;
    }
    if (out_file.empty()) {
        out_file = polyploid_snp_file + ".corrected";
    }

    // read the sequence lengths
    dna::SEQUENCE_LENGTH_LIST* sequence_lengths = base_distribution::read_sequence_lengths(base_distribution_file);

    // read SNP list
    dna::BASE_LIST* polyploid_snp_list = dna::read_snp_list(polyploid_snp_file, sequence_lengths);

    // read reference base list 
    dna::BASE_LIST* polyploid_reference_base_list = dna::read_reference_base_list(polyploid_snp_file, sequence_lengths);

    // read control SNP and reference base list if file is provided
    dna::BASE_LIST* control_snp_list = NULL;
    dna::BASE_LIST* control_reference_base_list = NULL;
    if (!control_snp_file.empty()) {
        control_snp_list = dna::read_snp_list(control_snp_file, sequence_lengths);
        control_reference_base_list = dna::read_reference_base_list(control_snp_file, sequence_lengths);
    } else {
        control_reference_base_list = dna::read_reference_as_list(reference_file, sequence_lengths);
    }
    // read the base distribution file
    base_distribution* polyploid_base_distribution_list = base_distribution::read(base_distribution_file, sequence_lengths);

    // identify missing SNPs
    if (!control_snp_file.empty()) {
        snp::identify_missing_polyploid_snps(polyploid_snp_list, polyploid_base_distribution_list, control_snp_list, polyploid_reference_base_list, control_reference_base_list);
    } else {
        snp::identify_missing_polyploid_snps(polyploid_snp_list, polyploid_base_distribution_list, polyploid_reference_base_list, control_reference_base_list);
    }

    // write the updated snp list
    dna::write_snp_list(polyploid_snp_list, polyploid_reference_base_list, sequence_lengths, out_file);

    // clean up
    delete polyploid_snp_list;
    delete control_snp_list;
    delete polyploid_reference_base_list;
    delete control_reference_base_list;
    delete polyploid_base_distribution_list;
    delete sequence_lengths;

    return 0;
}

int tabulate_sequence_lengths(int argc, char** argv) {

    string multi_fasta_file = "";
    string out_file = "";

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-i") == 0) { // input reads file
            multi_fasta_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // output file
            out_file = argv[i++];
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (multi_fasta_file.empty()) {
        cerr << "Multifasta file must be specified." << endl;
        return 1;
    }
    if (out_file.empty()) {
        out_file = multi_fasta_file + ".tab";
    }

    // extract the sequencing reads from the fastq file
    dna::tabulate_sequence_lengths(multi_fasta_file, out_file);

    return 0;
}

int tabulate_snps(int argc, char** argv) {

    string snp_filenames_file = "";
    string snp_file = "";
    string out_file = "";
    string snp_position_file = "";
    bool ignore_Ns_in_the_reference = false;
    bool fill_blanks = false;

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-i") == 0) { // input  file
            snp_filenames_file = argv[i++];
            continue;
        } else if (parameter.compare("-s") == 0) { // SNP file
            snp_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // output file
            out_file = argv[i++];
            continue;
        } else if (parameter.compare("-p") == 0) { // position file
            snp_position_file = argv[i++];
            continue;
        } else if (parameter.compare("-n") == 0) { // ignore Ns in the reference
            string str_ignore_Ns_in_the_reference = argv[i++];
            ignore_Ns_in_the_reference = (str_ignore_Ns_in_the_reference[0] == 'T' || str_ignore_Ns_in_the_reference[0] == 't' ? true : false);
            continue;
        } else if (parameter.compare("-b") == 0) { // fill blanks in the SNP table with the reference base
            string str_fill_blanks = argv[i++];
            fill_blanks = (str_fill_blanks[0] == 'T' || str_fill_blanks[0] == 't' ? true : false);
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (snp_filenames_file.empty() && snp_file.empty()) {
        cerr << "SNP file must be specified." << endl;
        return 1;
    } else if (!snp_file.empty() && snp_position_file.empty()) {
        cerr << "SNP position file must be specified." << endl;
        return 1;
    }
    if (out_file.empty()) {
        out_file = snp_filenames_file + ".tab";
    }

    map<string, string>* snp_file_list = NULL;
    dna::SNP_TABLE* snp_table = NULL;

    if (snp_position_file.empty()) {
        // read snp file list
        snp_file_list = utility::read_file_2(snp_filenames_file);
        // tabulate snps
        snp_table = dna::tabulate_snps(snp_file_list, ignore_Ns_in_the_reference);
    } else {
        // tabulate snps using the snp position file
        snp_table = dna::tabulate_snps(snp_file, snp_position_file, ignore_Ns_in_the_reference);
    }

    // write the snp table to output
    dna::write_snp_table(snp_file_list, snp_table, out_file, fill_blanks);

    // clean up
    delete snp_file_list;
    delete snp_table;

    return 0;
}

int check_diploid_snps(int argc, char** argv) {

    string diploid_snp_file = "";
    string polyploid_snp_file = "";
    string position_file = "";
    string base_distribution_file = "";
    string out_file = "";

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-p") == 0) { // polyploid snp file
            polyploid_snp_file = argv[i++];
            continue;
        } else if (parameter.compare("-pos") == 0) { // Position file containing reference bases
            position_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // output file
            out_file = argv[i++];
            continue;
        } else if (parameter.compare("-s") == 0) { // diploid snp file
            diploid_snp_file = argv[i++];
            continue;
        } else if (parameter.compare("-b") == 0) { // base distribution file
            base_distribution_file = argv[i++];
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (base_distribution_file.empty()) {
        cerr << "Base distribution file must be specified." << endl;
        return 1;
    } else if (diploid_snp_file.empty()) {
        cerr << "Diploid SNP file must be specified." << endl;
        return 1;
    }
    if (out_file.empty()) {
        out_file = polyploid_snp_file + ".corrected";
    }

    // read the sequence lengths
    dna::SEQUENCE_LENGTH_LIST* sequence_lengths = base_distribution::read_sequence_lengths(base_distribution_file);

    // read SNP lists
    dna::BASE_LIST* diploid_snp_list = dna::read_snp_list(diploid_snp_file, sequence_lengths);
    dna::BASE_LIST* polyploid_snp_list = NULL;
    if (!polyploid_snp_file.empty()) {
        polyploid_snp_list = dna::read_snp_list(polyploid_snp_file, sequence_lengths);
    }

    // read reference base lists
    dna::BASE_LIST* diploid_reference_base_list;
    if (position_file.empty()) {
        diploid_reference_base_list = dna::read_reference_base_list(diploid_snp_file, sequence_lengths);
    } else {
        diploid_reference_base_list = dna::read_reference_base_list(position_file, sequence_lengths);
    }
    dna::BASE_LIST* polyploid_reference_base_list = NULL;
    if (!polyploid_snp_file.empty()) {
        polyploid_reference_base_list = dna::read_reference_base_list(polyploid_snp_file, sequence_lengths);
    }
    // read the base distribution file
    base_distribution* diploid_base_distribution_list = base_distribution::read(base_distribution_file, sequence_lengths);

    // identify missing SNPs
    if (!polyploid_snp_file.empty()) {
        snp::check_diploid_snp_list(diploid_snp_list, diploid_base_distribution_list, polyploid_snp_list, diploid_reference_base_list, polyploid_reference_base_list);
    } else {
        snp::check_diploid_snp_list(diploid_snp_list, diploid_base_distribution_list, diploid_reference_base_list);
    }
    // write the updated snp list
    dna::write_snp_list(diploid_snp_list, diploid_reference_base_list, sequence_lengths, out_file);

    // clean up
    if (!polyploid_snp_file.empty()) {
        delete polyploid_snp_list;
    }
    delete diploid_snp_list;
    delete diploid_reference_base_list;
    delete diploid_base_distribution_list;
    delete sequence_lengths;

    return 0;
}

int check_snps(int argc, char** argv) {

    string snp_file = "";
    string position_file = "";
    string base_distribution_file = "";
    string out_file = "";

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-o") == 0) { // output file
            out_file = argv[i++];
            continue;
        } else if (parameter.compare("-s") == 0) { // diploid snp file
            snp_file = argv[i++];
            continue;
        } else if (parameter.compare("-p") == 0) { // Position file containing reference bases
            position_file = argv[i++];
            continue;
        } else if (parameter.compare("-b") == 0) { // base distribution file
            base_distribution_file = argv[i++];
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (base_distribution_file.empty()) {
        cerr << "Base distribution file must be specified." << endl;
        return 1;
    } else if (snp_file.empty()) {
        cerr << "SNP file must be specified." << endl;
        return 1;
    } else if (position_file.empty()) {
        cerr << "Position file must be specified." << endl;
        return 1;
    }
    if (out_file.empty()) {
        out_file = snp_file + ".corrected";
    }

    // read the sequence lengths
    dna::SEQUENCE_LENGTH_LIST* sequence_lengths = base_distribution::read_sequence_lengths(base_distribution_file);

    // read SNP list
    dna::BASE_LIST* snp_list = dna::read_snp_list(snp_file, sequence_lengths);

    // read reference base list
    dna::BASE_LIST* diploid_reference_base_list = dna::read_reference_base_list(position_file, sequence_lengths);

    // read the base distribution file
    base_distribution* diploid_base_distribution_list = base_distribution::read(base_distribution_file, sequence_lengths);

    // identify missing SNPs
    snp::check_diploid_snp_list(snp_list, diploid_base_distribution_list, diploid_reference_base_list);

    // write the updated snp list
    dna::write_snp_list(snp_list, diploid_reference_base_list, sequence_lengths, out_file);

    // clean up
    delete snp_list;
    delete diploid_reference_base_list;
    delete diploid_base_distribution_list;
    delete sequence_lengths;

    return 0;
}

int snp_table_to_list(int argc, char** argv) {

    string str_column = "";
    int column = 0;
    string snp_table_file = "";
    string out_file = "";
    bool ignore_Ns_in_the_reference = false;

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-i") == 0) { // input file
            snp_table_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // output file
            out_file = argv[i++];
            continue;
        } else if (parameter.compare("-c") == 0) { // column to extract
            str_column = argv[i++];
            column = atoi(str_column.c_str());
            continue;
        } else if (parameter.compare("-n") == 0) { // ignore Ns in the reference
            string str_ignore_Ns_in_the_reference = argv[i++];
            ignore_Ns_in_the_reference = (str_ignore_Ns_in_the_reference[0] == 'T' || str_ignore_Ns_in_the_reference[0] == 't' ? true : false);
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (snp_table_file.empty()) {
        cerr << "SNP table file must be specified." << endl;
        return 1;
    } else if (column == 0) {
        cerr << "Invalid column." << endl;
        return 1;
    }
    if (out_file.empty()) {
        out_file = snp_table_file + "." + str_column;
    }

    // extract SNP list
    dna::snp_table_to_list(snp_table_file, column, ignore_Ns_in_the_reference, out_file);


    return 0;
}

int remove_common_snps(int argc, char** argv) {

    string snp_file = "";
    string out_file = "";

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-s") == 0) { // snp file
            snp_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // output file
            out_file = argv[i++];
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (snp_file.empty()) {
        cerr << "SNP file must be specified." << endl;
        return 1;
    }
    if (out_file.empty()) {
        out_file = snp_file + ".filtered";
    }

    dna::remove_common_snps(snp_file, out_file);


    return 0;
}

int remove_uncommon_snps(int argc, char** argv) {

    string snp_file = "";
    string out_file = "";

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-s") == 0) { // snp file
            snp_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // output file
            out_file = argv[i++];
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (snp_file.empty()) {
        cerr << "SNP file must be specified." << endl;
        return 1;
    }
    if (out_file.empty()) {
        out_file = snp_file + ".filtered";
    }

    dna::remove_uncommon_snps(snp_file, out_file);


    return 0;
}

int run_bwa(int argc, char** argv) {

    string id = "";
    string reference_file = "";
    string left_reads_file = "";
    string right_reads_file = "";
    //    string out_file = "";
    string log_file = "";
    string dir_name;
    int max_insert_size = -1;
    bool report_ambiguous_bases = false;
    bool both_reads_to_be_mapped = false;
    int min_coverage_threshold = 5;
    int max_coverage_threshold = 50000;
    int snp_quality_threshold = 20;
    int indel_quality_threshold = 100;
    double gene_coverage_threshold = 0.60;
    int threads = 1;
    bool keep_intermediary_files = true;
    bool process_soft_clips = false;
    //    bool require_both_unique_reads = true;

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-r") == 0) { // reference file
            reference_file = argv[i++];
            continue;
        } else if (parameter.compare("-1") == 0) { // left reads file
            left_reads_file = argv[i++];
            continue;
        } else if (parameter.compare("-2") == 0) { // right reads file
            right_reads_file = argv[i++];
            continue;
        } else if (parameter.compare("-l") == 0) { // log file
            log_file = argv[i++];
            continue;
        } else if (parameter.compare("-x") == 0) { // max insert size
            max_insert_size = atoi(argv[i++]);
            continue;
        } else if (parameter.compare("-a") == 0) { // report ambiguous base
            string str_report_ambigous_base = argv[i++];
            report_ambiguous_bases = (str_report_ambigous_base[0] == 'T' || str_report_ambigous_base[0] == 't' ? true : false);
            continue;
        } else if (parameter.compare("-d") == 0) { // minimum coverage threshold
            min_coverage_threshold = atoi(argv[i++]);
            continue;
        } else if (parameter.compare("-D") == 0) { // maximum coverage threshold
            max_coverage_threshold = atoi(argv[i++]);
            continue;
        } else if (parameter.compare("-sq") == 0) { // snp quality threshold
            snp_quality_threshold = atoi(argv[i++]);
            continue;
        } else if (parameter.compare("-iq") == 0) { // indel quality threshold
            indel_quality_threshold = atoi(argv[i++]);
            continue;
        } else if (parameter.compare("-sc") == 0) { // process soft clips?
            string str_process_soft_clips = argv[i++];
            process_soft_clips = (str_process_soft_clips[0] == 'T' || str_process_soft_clips[0] == 't' ? true : false);
            continue;
        } else if (parameter.compare("-gc") == 0) { // gene coverage threshold (for soft clip processing)
            gene_coverage_threshold = atoi(argv[i++]);
            continue;
        } else if (parameter.compare("-t") == 0) { // threads, where applicable
            threads = atoi(argv[i++]);
            continue;
        } else if (parameter.compare("-id") == 0) { // id
            id = argv[i++];
            continue;
        } else if (parameter.compare("-dir") == 0) { // directory to store alignment data
            dir_name = argv[i++];
            continue;
        } else if (parameter.compare("-b") == 0) { // require both reads in a pair to be mapped
            string str_both_reads_to_be_mapped = argv[i++];
            both_reads_to_be_mapped = (str_both_reads_to_be_mapped[0] == 'T' || str_both_reads_to_be_mapped[0] == 't' ? true : false);
            continue;
            //        } else if (parameter.compare("-UR") == 0) { // allow one read to be repeat if other is unique
            //            string str_one_repeat_read_allowed = argv[i++];
            //            require_both_unique_reads = (str_one_repeat_read_allowed[0] == 'T' || str_one_repeat_read_allowed[0] == 't' ? false : true);
            //            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }

    }

    // parameter checking
    if (reference_file.empty()) {
        cerr << "Reference file must be specified." << endl;
        return 1;
    } else if (left_reads_file.empty()) {
        cout << "Read file must be specified." << endl;
        return 1;
    }

    if (id.empty()) {
        string read_file_name = utility::remove_extension(utility::extract_filename(left_reads_file));
        if (!right_reads_file.empty()) {
            read_file_name = read_file_name.substr(0, read_file_name.length() - 2);
        }
        id = utility::remove_extension(utility::extract_filename(reference_file)) + "_" + read_file_name;
    }
    // dir name
    if (dir_name.empty()) {
        dir_name = ".";
    }
    if (!utility::dir_exists(dir_name)) {
        cout << "Creating data directory ..." << endl;
        utility::make_dir(dir_name);
    }
    if (log_file.empty()) {
        log_file = utility::check_dir(dir_name) + id + ".log";
    }

    // open the log file stream
    ofstream ofs_log(log_file.c_str());
    // write the command to the log file
    ofs_log << "Command:";
    for (int i = 0; i < argc; i++) {
        ofs_log << " " << argv[i];
    }
    ofs_log << endl << endl;


    aligner* the_aligner;
    sam* sam_ptr;

    // create the aligner object
    the_aligner = new aligner(id, dir_name, keep_intermediary_files);
    // set the mapper object
    //the_aligner->set_mapper(new bwa(id, dir_name, max_insert_size, threads));
    the_aligner->set_mapper(new bwa(id, dir_name, threads));
    // set necessary parameters
    the_aligner->set_log_file(log_file);
    the_aligner->set_report_ambiguous_bases(report_ambiguous_bases);
    the_aligner->set_ignore_unpaired_reads(both_reads_to_be_mapped);
    the_aligner->set_min_coverage_threshold(min_coverage_threshold);
    the_aligner->set_max_coverage_threshold(max_coverage_threshold);
    the_aligner->set_snp_quality_threshold(snp_quality_threshold);
    the_aligner->set_indel_quality_threshold(indel_quality_threshold);
    the_aligner->set_gene_coverage_threshold(gene_coverage_threshold);
    the_aligner->set_max_insert_size(max_insert_size);

    // run the aligner (full_run = true .. generate all files)
    //    sam_ptr = the_aligner->run(reference_file, left_reads_file, right_reads_file, true, false, require_both_unique_reads);
    sam_ptr = the_aligner->run(reference_file, left_reads_file, right_reads_file, true, false, process_soft_clips);

    // output the stats
    ofs_log << "Alignment:" << endl;
    ofs_log << "\tTotal reads: " << the_aligner->get_total_reads() << endl;
    ofs_log << "\tMapped reads: " << the_aligner->get_mapped_reads() << " (" << utility::round2((double) the_aligner->get_mapped_reads() / (double) the_aligner->get_total_reads() * 100) << "%)" << endl;
    ofs_log << "\tFiltered reads: " << the_aligner->get_filtered_reads() << " (" << utility::round2((double) the_aligner->get_filtered_reads() / (double) the_aligner->get_total_reads() * 100) << "%)" << endl;
    ofs_log << endl;
    ofs_log << "Mutations:" << endl;
    ofs_log << "\tSubstitutions:" << endl;
    ofs_log << "\t\tReported: " << the_aligner->get_snp_count() << endl;
    ofs_log << "\t\tAccepted: " << the_aligner->get_accepted_snp_count() << endl;
    ofs_log << "\t\tN (included): " << the_aligner->get_snp_count_n() << endl;
    ofs_log << "Insertions/Deletions: " << endl;
    ofs_log << "\tReported: " << the_aligner->get_indel_count() << endl;
    ofs_log << "\tAccepted Insertions: " << the_aligner->get_accepted_insertion_count() << endl;
    ofs_log << "\tAccepted Deletions: " << the_aligner->get_accepted_deletion_count() << endl;
    //    ofs_log << "\tInsertions: " << endl;
    //    ofs_log << "\t\tReported: " << the_aligner->get_insertion_count() << endl;
    //    ofs_log << "\t\tAccepted: " << the_aligner->get_accepted_insertion_count() << endl;
    //    ofs_log << "\tDeletions:" << endl;
    //    ofs_log << "\t\tReported: " << the_aligner->get_deletion_count() << endl;
    //    ofs_log << "\t\tAccepted: " << the_aligner->get_accepted_deletion_count() << endl;
    ofs_log << endl;

    // close the log file
    ofs_log.close();

    // free the memory
    delete sam_ptr;
    delete the_aligner;

    return 0;
}

int filter_sam_file(int argc, char** argv) {

    string sam_file;
    string filtered_sam_file = "";
    string gff_file = "";
    string unique_mapping_tag = "";
    string repeat_mapping_tag = "";
    bool paired_data = true;
    int max_insert_size = -1;
    bool both_reads_to_be_mapped = false;
    //    bool require_both_unique_reads = true;

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-i") == 0) { // SAM file
            sam_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // filtered SAM file
            filtered_sam_file = argv[i++];
            continue;
        } else if (parameter.compare("-g") == 0) { // gff file
            gff_file = argv[i++];
            continue;
        } else if (parameter.compare("-x") == 0) { // max insert size
            max_insert_size = atoi(argv[i++]);
            continue;
        } else if (parameter.compare("-u") == 0) { // unique mapping tag
            unique_mapping_tag = argv[i++];
            continue;
        } else if (parameter.compare("-r") == 0) { // repeat mapping tag
            repeat_mapping_tag = argv[i++];
            continue;
        } else if (parameter.compare("-p") == 0) { // paired data?
            string str_paired_data = argv[i++];
            paired_data = (str_paired_data[0] == 'T' || str_paired_data[0] == 't' ? true : false);
            continue;
        } else if (parameter.compare("-b") == 0) { // require both reads in a pair to be mapped
            string str_both_reads_to_be_mapped = argv[i++];
            both_reads_to_be_mapped = (str_both_reads_to_be_mapped[0] == 'T' || str_both_reads_to_be_mapped[0] == 't' ? true : false);
            continue;
            //        } else if (parameter.compare("-UR") == 0) { // allow one read to be repeat if other is unique
            //            string str_one_repeat_read_allowed = argv[i++];
            //            require_both_unique_reads = (str_one_repeat_read_allowed[0] == 'T' || str_one_repeat_read_allowed[0] == 't' ? false : true);
            //            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }

    }

    // parameter checking
    if (sam_file.empty()) {
        cerr << "SAM file must be specified." << endl;
        return 1;
    }

    if (filtered_sam_file.empty()) {
        filtered_sam_file = sam_file + ".filtered";
    }

    sam* sam_ptr = new sam(sam_file, "", paired_data, unique_mapping_tag, repeat_mapping_tag);
    sam_ptr->set_ignore_unpaired_reads(both_reads_to_be_mapped);
    //    sam_ptr->filter(filtered_sam_file, max_insert_size, gff_file, require_both_unique_reads);
    //    long t1 = time(0);
    sam_ptr->filter(filtered_sam_file, max_insert_size, gff_file);
    //    cout << time(0) - t1 << endl;
    //    long t2 = time(0);
    //    sam_ptr->filter_new(filtered_sam_file, max_insert_size, gff_file);
    //    cout << time(0) - t2 << endl;


    // free the memory
    delete sam_ptr;

    return 0;
}

int filter_sam_file_new(int argc, char** argv) {

    string sam_file;
    string filtered_sam_file = "";
    string gff_file = "";
    string unique_mapping_tag = "";
    string repeat_mapping_tag = "";
    bool paired_data = true;
    int max_insert_size = -1;
    bool both_reads_to_be_mapped = false;
    //    bool require_both_unique_reads = true;

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-i") == 0) { // SAM file
            sam_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // filtered SAM file
            filtered_sam_file = argv[i++];
            continue;
        } else if (parameter.compare("-g") == 0) { // gff file
            gff_file = argv[i++];
            continue;
        } else if (parameter.compare("-x") == 0) { // max insert size
            max_insert_size = atoi(argv[i++]);
            continue;
        } else if (parameter.compare("-u") == 0) { // unique mapping tag
            unique_mapping_tag = argv[i++];
            continue;
        } else if (parameter.compare("-r") == 0) { // repeat mapping tag
            repeat_mapping_tag = argv[i++];
            continue;
        } else if (parameter.compare("-p") == 0) { // paired data?
            string str_paired_data = argv[i++];
            paired_data = (str_paired_data[0] == 'T' || str_paired_data[0] == 't' ? true : false);
            continue;
        } else if (parameter.compare("-b") == 0) { // require both reads in a pair to be mapped
            string str_both_reads_to_be_mapped = argv[i++];
            both_reads_to_be_mapped = (str_both_reads_to_be_mapped[0] == 'T' || str_both_reads_to_be_mapped[0] == 't' ? true : false);
            continue;
            //        } else if (parameter.compare("-UR") == 0) { // allow one read to be repeat if other is unique
            //            string str_one_repeat_read_allowed = argv[i++];
            //            require_both_unique_reads = (str_one_repeat_read_allowed[0] == 'T' || str_one_repeat_read_allowed[0] == 't' ? false : true);
            //            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }

    }

    // parameter checking
    if (sam_file.empty()) {
        cerr << "SAM file must be specified." << endl;
        return 1;
    }

    if (filtered_sam_file.empty()) {
        filtered_sam_file = sam_file + ".filtered";
    }

    sam* sam_ptr = new sam(sam_file, "", paired_data, unique_mapping_tag, repeat_mapping_tag);
    sam_ptr->set_ignore_unpaired_reads(both_reads_to_be_mapped);
    //    sam_ptr->filter(filtered_sam_file, max_insert_size, gff_file, require_both_unique_reads);
    //    long t1 = time(0);
    //    sam_ptr->filter(filtered_sam_file, max_insert_size, gff_file);
    //    cout << time(0) - t1 << endl;
    //    long t2 = time(0);
    sam_ptr->filter_new(filtered_sam_file, max_insert_size, gff_file);
    //    cout << time(0) - t2 << endl;


    // free the memory
    delete sam_ptr;

    return 0;
}

int get_codons(int argc, char** argv) {

    string snp_file = "";
    string reference_file = "";
    string gff_file = "";
    string out_file = "";
    bool is_full_gff = false;

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-i") == 0) { // snp file
            snp_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // output file
            out_file = argv[i++];
            continue;
        } else if (parameter.compare("-r") == 0) { // reference file
            reference_file = argv[i++];
            continue;
        } else if (parameter.compare("-g") == 0) { // gff file, optional
            gff_file = argv[i++];
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (snp_file.empty()) {
        cerr << "SNP file must be specified." << endl;
        return 1;
    } else if (reference_file.empty()) {
        cerr << "Reference file must be specified." << endl;
        return 1;
    }

    if (out_file.empty()) {
        out_file = snp_file + ".codon";
    }
    if (gff_file.empty()) {
        gff_file = reference_file + ".gff";
    }


    dna::get_codon(snp_file, reference_file, gff_file, out_file);

    return 0;
}

int snps_per_gene(int argc, char** argv) {

    string snp_file = "";
    string coordinate_file = "";
    string gff_file = "";
    string out_file = "";
    bool ignore_Ns_in_the_reference = false;
    bool is_snp_table = false;
    bool count_snp_positions = false;

    //    bool nullitetra_data = false;

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-s") == 0) { // snp file
            snp_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // output file
            out_file = argv[i++];
            continue;
        } else if (parameter.compare("-c") == 0) { // coordinate file
            coordinate_file = argv[i++];
            continue;
        } else if (parameter.compare("-g") == 0) { // or gff file
            gff_file = argv[i++];
            continue;
        } else if (parameter.compare("-n") == 0) { // ignore Ns in the reference
            string str_ignore_Ns_in_the_reference = argv[i++];
            ignore_Ns_in_the_reference = (str_ignore_Ns_in_the_reference[0] == 'T' || str_ignore_Ns_in_the_reference[0] == 't' ? true : false);
            continue;
        } else if (parameter.compare("-t") == 0) { // is snp table?
            string str_is_snp_table = argv[i++];
            is_snp_table = (str_is_snp_table[0] == 'T' || str_is_snp_table[0] == 't' ? true : false);
            continue;
        } else if (parameter.compare("-p") == 0) { // count snp positions only (unique snps)? .. meaningful for tabulated data only
            string str_count_snp_positions = argv[i++];
            count_snp_positions = (str_count_snp_positions[0] == 'T' || str_count_snp_positions[0] == 't' ? true : false);
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (snp_file.empty()) {
        cerr << "SNP file must be specified." << endl;
        return 1;
    } else if (coordinate_file.empty() && gff_file.empty()) {
        cerr << "GFF or Coordinate file must be specified." << endl;
        return 1;
    }
    if (out_file.empty()) {
        out_file = snp_file + ".snp.per.gene";
    }

    if (gff_file.empty()) {
        snp::snps_per_gene(snp_file, coordinate_file, out_file, ignore_Ns_in_the_reference, false, is_snp_table, count_snp_positions);
    } else {
        snp::snps_per_gene(snp_file, gff_file, out_file, ignore_Ns_in_the_reference, true, is_snp_table, count_snp_positions);
    }
    return 0;
}

int get_codon_frequencies(int argc, char** argv) {

    string reference_file = "";
    string gff_file = "";
    string out_file = "";

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-o") == 0) { // output file
            out_file = argv[i++];
            continue;
        } else if (parameter.compare("-r") == 0) { // reference file
            reference_file = argv[i++];
            continue;
        } else if (parameter.compare("-g") == 0) { // gff file, optional
            gff_file = argv[i++];
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (reference_file.empty()) {
        cerr << "Reference file must be specified." << endl;
        return 1;
    }

    if (out_file.empty()) {
        out_file = reference_file + ".codon.freq";
    }
    if (gff_file.empty()) {
        gff_file = reference_file + ".gff";
    }


    dna::get_codon_frequencies(reference_file, gff_file, out_file);

    return 0;
}

int calculate_base_distribution(int argc, char** argv) {

    string sam_file;
    string base_dist_file = "";

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-i") == 0) { // SAM file
            sam_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // base dist file
            base_dist_file = argv[i++];
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }

    }

    // parameter checking
    if (sam_file.empty()) {
        cerr << "SAM file must be specified." << endl;
        return 1;
    }

    if (base_dist_file.empty()) {
        base_dist_file = utility::remove_extension(sam_file) + ".base.dist";
    }

    sam* sam_ptr = new sam(sam_file, "", true, "", "");
    sam_ptr->calculate_base_distribution(base_dist_file);


    // free the memory
    delete sam_ptr;

    return 0;
}

int get_consensus_bases(int argc, char** argv) {

    string snp_file = "";
    string out_file = "";
    bool GATK = false;


    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-s") == 0) { // snp file
            snp_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // output file
            out_file = argv[i++];
            continue;
        } else if (parameter.compare("-gatk") == 0) { // GATK
            string str_GATK = argv[i++];
            GATK = (str_GATK[0] == 'T' || str_GATK[0] == 't' ? true : false);
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (snp_file.empty()) {
        cerr << "SNP file must be specified." << endl;
        return 1;
    }
    if (out_file.empty()) {
        out_file = snp_file + ".con";
    }

    if (!GATK) {
        snp::get_consensus_base(snp_file, out_file);
    } else {
        snp::get_consensus_base_GATK(snp_file, out_file);
    }
    return 0;
}

int compare_snps(int argc, char** argv) {

    string snp_table_file = "";
    string polyploid_snp_file = "";
    string out_file = "";


    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-s") == 0) { // snp table file
            snp_table_file = argv[i++];
            continue;
        } else if (parameter.compare("-p") == 0) { // polyploid snp file
            polyploid_snp_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // output file
            out_file = argv[i++];
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (snp_table_file.empty()) {
        cerr << "SNP table file must be specified." << endl;
        return 1;
    } else if (polyploid_snp_file.empty()) {
        cerr << "Polyploid SNP file must be specified." << endl;
        return 1;
    }
    if (out_file.empty()) {
        out_file = snp_table_file + ".summary";
    }


    snp::compare_snps(snp_table_file, polyploid_snp_file);
    return 0;
}

int sort_multi_fasta(int argc, char** argv) {

    string multi_fasta_file = "";
    string out_file = "";
    string ordering_file = "";

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-i") == 0) { // multifasta file
            multi_fasta_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // output file
            out_file = argv[i++];
            continue;
        } else if (parameter.compare("-s") == 0) { // ordering file
            ordering_file = argv[i++];
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (multi_fasta_file.empty()) {
        cerr << "Multifasta file must be specified." << endl;
        return 1;
    } else if (ordering_file.empty()) {
        cerr << "Ordering names file must be specified." << endl;
        return 1;
    }
    if (out_file.empty()) {
        out_file = multi_fasta_file + ".sorted";
    }

    dna::sort_multi_fasta(multi_fasta_file, out_file, ordering_file);

    return 0;
}

int split_fastq_into_pairs(int argc, char** argv) {
    string reads_file = "";
    string out_file_1 = "";
    string out_file_2 = "";
    int read_length = 0;

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-i") == 0) { // input reads file
            reads_file = argv[i++];
            continue;
        } else if (parameter.compare("-o1") == 0) { // output file 1
            out_file_1 = argv[i++];
            continue;
        } else if (parameter.compare("-o2") == 0) { // output file 1
            out_file_2 = argv[i++];
            continue;
        } else if (parameter.compare("-l") == 0) { // gff file
            read_length = atoi(argv[i++]);
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    // split the sequencing reads from the fastq file
    dna::split_fastq_into_pairs(reads_file, out_file_1, out_file_2, read_length);

    return 0;
}

int test(int argc, char** argv) {

    string reference_file = "";
    string new_reference_file = "";
    string snp_file = "";
    string deletion_file = "";
    string insertion_file = "";

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-r") == 0) {
            reference_file = argv[i++];
            continue;
        } else if (parameter.compare("-n") == 0) {
            new_reference_file = argv[i++];
            continue;
        } else if (parameter.compare("-s") == 0) {
            snp_file = argv[i++];
            continue;
        } else if (parameter.compare("-d") == 0) {
            deletion_file = argv[i++];
            continue;
        } else if (parameter.compare("-i") == 0) {
            insertion_file = argv[i++];
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (reference_file.empty()) {
        cerr << "Reference file must be specified." << endl;
        return 1;
    }


    dna::fix_mutations(reference_file, new_reference_file, snp_file, deletion_file, insertion_file, false);


    return 0;
}

int extract_exonic_sequences(int argc, char** argv) {

    string fasta_file = "";
    string gff_file = "";
    string out_file = "";

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-i") == 0) { // fasta file
            fasta_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // output file
            out_file = argv[i++];
            continue;
        } else if (parameter.compare("-g") == 0) { // gff file
            gff_file = argv[i++];
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (fasta_file.empty()) {
        cerr << "Fasta file must be specified." << endl;
        return 1;
    } else if (gff_file.empty()) {
        cerr << "GFF file must be specified." << endl;
        return 1;
    } else if (out_file.empty()) {
        cerr << "Output file must be specified." << endl;
        return 1;
    }


    dna::extract_exome(fasta_file, gff_file, out_file);

    return 0;
}

int extract_subsequence(int argc, char** argv) {

    string fasta_file = "";
    string out_file = "";
    int from = -1;
    int to = -1;

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-i") == 0) { // fasta file
            fasta_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // output file
            out_file = argv[i++];
            continue;
        } else if (parameter.compare("-f") == 0) { // from
            from = atoi(argv[i++]);
            continue;
        } else if (parameter.compare("-t") == 0) { // to
            to = atoi(argv[i++]);
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (fasta_file.empty()) {
        cerr << "Fasta file must be specified." << endl;
        return 1;
    }


    if (out_file.empty()) {
        out_file = fasta_file + ".out";
        return 1;
    }

    if (from < 1 || from > to) {
        cerr << "Invalid range" << endl;
        return 1;
    }

    dna::extract_subsequence(fasta_file, out_file, from, to);

    return 0;
}

int mapping_summary(int argc, char** argv) {

    string header_file = "";
    string pileup_file = "";
    string gff_file = "";
    string out_file = "";
    int read_coverage_threshold = 1;

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-h") == 0) { // header file or sam file containing header
            header_file = argv[i++];
            continue;
        } else if (parameter.compare("-p") == 0) { // pileup file
            pileup_file = argv[i++];
            continue;
        } else if (parameter.compare("-g") == 0) { // gff file
            gff_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // out file
            out_file = argv[i++];
            continue;
        } else if (parameter.compare("-c") == 0) { // read coverage threshold
            read_coverage_threshold = atoi(argv[i++]);
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (header_file.empty()) {
        cerr << "SAM file or a header file must be specified." << endl;
        return 1;
    }
    if (pileup_file.empty()) {
        cerr << "Pileup file must be specified." << endl;
        return 1;
    }
    if (gff_file.empty()) {
        cerr << "GFF file must be specified." << endl;
        return 1;
    }


    if (out_file.empty()) {
        out_file = pileup_file + ".map.sum";
        return 1;
    }

    sam::mapping_summary(header_file, pileup_file, gff_file, out_file, read_coverage_threshold);

    return 0;
}

int soft_clip_summary(int argc, char** argv) {

    string sam_file = "";
    string out_sam_file = "";
    string out_file = "";
    string gff_file = "";
    string mapping_summary_file = "";
    int mapping_qual_threshold = -1;
    int read_coverage_threshold = 3;
    double gene_coverage_threshold = 0.6;

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-s") == 0) { // header file or sam file containing header
            sam_file = argv[i++];
            continue;
        } else if (parameter.compare("-os") == 0) { // out sam file
            out_sam_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // out file
            out_file = argv[i++];
            continue;
        } else if (parameter.compare("-m") == 0) { // mapping summary file
            mapping_summary_file = argv[i++];
            continue;
        } else if (parameter.compare("-q") == 0) { // mapping quality threshold
            mapping_qual_threshold = atoi(argv[i++]);
            continue;
        } else if (parameter.compare("-g") == 0) { // gff file
            gff_file = argv[i++];
            continue;
        } else if (parameter.compare("-rc") == 0) { // read coverage threshold
            read_coverage_threshold = atoi(argv[i++]);
            continue;
        } else if (parameter.compare("-gc") == 0) { // gene coverage threshold
            gene_coverage_threshold = atof(argv[i++]);
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (sam_file.empty()) {
        cerr << "SAM file or a header file must be specified." << endl;
        return 1;
    }
    if (out_file.empty()) {
        out_file = sam_file + ".soft.clip.sum";
    }
    if (mapping_qual_threshold < 0) {
        cerr << "Mapping quality must be specified." << endl;
        return 1;
    }
    if (read_coverage_threshold < 0) {
        cerr << "Read coverage threshold must be non-negative." << endl;
        return 1;
    }
    if (gene_coverage_threshold < 0.0 || gene_coverage_threshold > 1.0) {
        cerr << "Gene coverage threshold must be between 0.0 and 1.0." << endl;
        return 1;
    }
    if (!gff_file.empty() && mapping_summary_file.empty()) {
        cerr << "Mapping Summary File should be specified if GFF file is provided" << endl;
        return 1;
    }


    if (out_file.empty()) {
        out_file = out_sam_file + ".map.sum";
    }

    sam::soft_clip_summary(sam_file, out_file, mapping_qual_threshold, out_sam_file, gff_file);
    if (!gff_file.empty()) {
        string accepted_out_file = out_file + ".accepted";
        sam::identify_soft_clip_at_mapping_boundary(mapping_summary_file, out_file, accepted_out_file, read_coverage_threshold, gene_coverage_threshold);
    }
    return 0;
}

int process_soft_clips(int argc, char** argv) {

    string sam_file = "";
    string out_sam_file = "";
    string soft_clip_file = "";

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-s") == 0) { // header file or sam file containing header
            sam_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // out sam file
            out_sam_file = argv[i++];
            continue;
        } else if (parameter.compare("-sc") == 0) { // soft clip file
            soft_clip_file = argv[i++];
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (sam_file.empty()) {
        cerr << "SAM file must be specified." << endl;
        return 1;
    }
    if (soft_clip_file.empty()) {
        cerr << "Soft clip file must be specified." << endl;
        return 1;
    }

    if (out_sam_file.empty()) {
        out_sam_file = sam_file + ".soft.clip.processed";
    }

    sam::process_soft_clips(sam_file, out_sam_file, soft_clip_file);
    return 0;
}

int remove_transcript_variants_ensembl(int argc, char** argv) {

    string fasta_file = "";
    string out_file = "";

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-i") == 0) { // fasta file
            fasta_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // output file
            out_file = argv[i++];
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (fasta_file.empty()) {
        cerr << "Fasta file must be specified." << endl;
        return 1;
    }


    if (out_file.empty()) {
        out_file = fasta_file + ".out";
    }

    dna::remove_transcript_variants_ensembl(fasta_file, out_file);

    return 0;
}

int remove_transcript_variants_trinity(int argc, char** argv) {

    string fasta_file = "";
    string out_file = "";

    // read input arguments
    int i = 2;
    int number_of_arguments = argc - 1;
    while (i < number_of_arguments) {
        string parameter = argv[i++];
        if (parameter.compare("-i") == 0) { // fasta file
            fasta_file = argv[i++];
            continue;
        } else if (parameter.compare("-o") == 0) { // output file
            out_file = argv[i++];
            continue;
        } else {
            cerr << "Invalid option: " << argv[i - 1] << endl;
            return 1;
        }
    }

    if (fasta_file.empty()) {
        cerr << "Fasta file must be specified." << endl;
        return 1;
    }


    if (out_file.empty()) {
        out_file = fasta_file + ".out";
    }

    dna::remove_transcript_variants_trinity(fasta_file, out_file);

    return 0;
}
