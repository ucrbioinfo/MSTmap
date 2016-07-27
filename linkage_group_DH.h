/*
 *  linkage_group.h
 *  ApproxMap
 *
 *  Created by yonghui on 4/9/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef linkage_group_header
#define linkage_group_header

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <ctime>
#include <cassert>
#include <queue>
#include <utility>
#include <set>
#include <algorithm>
#include "constants.h"
#include "MSTOpt.h"

using namespace std;

bool cmp(pair<double, pair<int,int> > element1, pair<double, pair<int,int> > element2);

class DF{
    public:
        virtual double CM(double rp) const {
            cout << "ERROR, specific DF should be used instead" << endl;
            assert(false);
            return 0;
        };
        virtual double RP(double cm) const {
            cout << "ERROR, specific DF should be used instead" << endl;
            assert(false);
            return 0;
        };
        virtual void print_df_name()  const {
            cout << "generic df" << endl;
        };
};

class DF_Haldane:public DF{
    public:
        virtual double CM(double rp) const { // was rp
	    if (rp >= 0.5) return 50.0; //rp = 0.5 - 0.000001;
            //double rp = r / 2*(1-r);
            return -50.0 * log(1 - 2 * rp);
        };
        virtual double RP(double cm) const {
            return -0.5 * (exp(-cm / 50.0) - 1.0);
        };
        virtual void print_df_name() const {
            cout << "Haldane" << endl;
        };
};

class DF_Kosambi:public DF{
    public:
        virtual double CM(double rp) const { // was rp
            if (rp >= 0.5) return 50.0; //rp = 0.5 - 0.000001;
            //double rp = r / 2*(1-r);
            return 25.0 * log((1 + 2.0 * rp) / (1 - 2.0 * rp));
        };
        virtual double RP(double cm) const {
            return 0.5 * ((exp(cm / 25) - 1) / (exp(cm / 25) + 1));
        };
        virtual void print_df_name() const {
            cout << "Kosambi" << endl;
        };
};

class linkage_group {
    public:
        const vector<vector<double> >& get_pair_wise_distance() const;
        void return_order(vector<int>& out_order, 
                          double & _lowerbound, 
                          double & _upper_bound, 
                          double & _cost_after_initialization, 
                          vector<double> & _distances) const;
        void dump_common() const;
        void bad_genotypes(vector<pair<int,int> >& bad_genotypes) const;
        void dump_distance_matrix();
    protected:
    
        // produce pairwise_distance in cM
        void generate_distance_in_cM(vector<vector<double> >& distance_in_cM);
        void generate_distance_in_ML(vector<vector<double> >& distance_in_ML);
        
        bool detect_bad_data;
        ObjFunc objective_function;
        int number_of_bins;
        int number_of_individuals;
        
        // the distance is normalized to be the expected number of cross-overs per meiosis 
        // given the number of individuals 
        // namely, d_{i,j} = number_of_individuals * r
        // where r is the recombination probability
        vector<vector<double> > pair_wise_distances;
        vector<pair<int,int> > missing_data;
        vector<int> bin_sizes;
        
        // added by yonghui on March 8th, 2008. This data structure is updated 
        // after each iteration
        vector<pair<int,int> > suspicious_data;

        vector<int> current_order; //concatenation of the markers in the order of the path
        vector<int> MST; // the current MST
        double MST_lower_bound;
        double current_upper_bound;
        double cost_after_initialization;
        
        // the wrapper of distance functions (i.e. haldane or kosambi)
        // this is just a reference to an object owned by another object (an genetic_map instance)
        // no need to delete this object
        DF* df;
};

class linkage_group_DH: public linkage_group {
    private:
        vector<vector<float> > raw_data;
        
        // 0: if it is the same as the original input
        // 1: if it is missing
        // 2 and up: if it is deleted 
        vector<vector<int> > data_status;
        int iteration_number;
        
        vector<double> suspicious_data_backup;
        
        /*Calculate the pair_wise distance*/
        void calculate_pair_wise_distance();
        void calculate_pair_wise_distance_initialize();
        int detect_bad_markers();
        void revert_suspicious_data();
        
        void estimate_missing_data();
        // a supportive function to be called by estimate_missing_data
        double estimate_one_missing(const vector<int>& rev_perm, int marker_id, int indi_id);
        
    public:
        linkage_group_DH(int _number_of_bins, 
                         int _number_of_individuals,
                         bool _detect_bad_data,
                         ObjFunc _objective_function,
                         DF* _df,
                         const vector<vector<float> >& _raw_data, 
                         const vector<int>& _current_order, 
                         const vector<pair<int,int> >& _missing_data,
                         const vector<int>& _bin_sizes);
        ~linkage_group_DH();
        
        void dump() const;
    
        void order_markers();
};


#endif
