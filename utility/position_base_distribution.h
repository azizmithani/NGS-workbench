/* 
 * File:   position_base_distribution.h
 * Author: Aziz Mithani
 *
 * Created on 08 May 2011, 14:16
 */

#ifndef POSITION_BASE_DISTRIBUTION_H
#define	POSITION_BASE_DISTRIBUTION_H

#include <vector>
#include <map>

using namespace std;

class position_base_distribution {
public:
    position_base_distribution();
    position_base_distribution(int coverage_a, int coverage_t, int coverage_g, int coverage_c, int coverage_n);
    virtual ~position_base_distribution();

    // type definitions
    typedef position_base_distribution* SEQ_BASE_DISTRIBUTION_LIST;
    //typedef map < string, SEQ_BASE_DISTRIBUTION_LIST > BASE_DISTRIBUTION_LIST;

    int get_coverage();
    int get_coverage_A() const;
    int get_coverage_C() const;
    int get_coverage_G() const;
    int get_coverage_T() const;
    void increment_coverage(char base);
    string to_string();
    string to_string(string sep);
    void set_position_base_distribution(vector<string>& entries, int offset);
    //void set_position_base_distribution(vector<string>& entries);
    
private:
    int coverage_a;
    int coverage_c;
    int coverage_g;
    int coverage_t;
    int coverage_n;



};

#endif	/* POSITION_BASE_DISTRIBUTION_H */

