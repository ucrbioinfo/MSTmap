/*
 *  linkage_group_RIL.h
 *  ApproxMap
 *
 *  Created by yonghui on 12/13/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef LINKAGE_GROUP_RIL_HEADER
#define LINKAGE_GROUP_RIL_HEADER

#include <cassert>
#include <vector>
#include "linkage_group_DH.h"
#include "MSTOpt.h"

using namespace std;

struct allel_state {
    double A;
    double B;
    double AB;
    bool missing;
};

class RIL_dist_cal{
    public:
        RIL_dist_cal(int _generation_index, const vector<allel_state>& _marker1, const vector<allel_state>& _marker2)
            :generation_index_(_generation_index),
            marker1_(_marker1),
            marker2_(_marker2){
            num_of_eff_individuals_ = 0;
            assert(marker1_.size() == marker2_.size());
            num_of_individuals_ = marker1_.size();
            count_class();
            upper_bound();
            lower_bound();
        };
        double Dist() const {
            double opt_delta = find_opt_delta();
            return opt_delta * num_of_individuals_;
        };
    private:
        void count_class();
        void upper_bound();
        void lower_bound();
        
        void expected_CDEFG(double delta,
                            double& C,
                            double& D,
                            double& E,
                            double& F,
                            double& G) const; 
                            
        double squared_error(double delta) const;
        
        double find_opt_delta() const;
        
        // ----------------------------
        // private data members section
        // ----------------------------
        
        // generation_index_ should be set to 1 for F1 generation, 2 to F2 generation and so on...
        int generation_index_;
        
        int num_of_individuals_;
        int num_of_eff_individuals_;
        const vector<allel_state>& marker1_;
        const vector<allel_state>& marker2_;
        
        double CC_;
        double DD_;
        double EE_;
        double FG_;
        double delta_upper_bound_;
        double delta_lower_bound_;
};

class linkage_group_RIL: public linkage_group {
    public:
        linkage_group_RIL(int _number_of_bins, 
                          int _number_of_individuals,
                          int _generation_index,
                          DF* _df,
                          const vector<vector<allel_state> >& _raw_data, 
                          const vector<int>& _current_order, 
                          const vector<pair<int,int> >& _missing_data);
        ~linkage_group_RIL();
        
        void dump() const;
    
        void order_markers();
        
    private:

        /*Calculate the pair_wise distance*/
        void calculate_pair_wise_distance();
        
        void estimate_missing_data();
        
        // ----------------------------
        // private data members section
        // ----------------------------
        vector<vector<allel_state> > raw_data;
        int generation_index_;
};

#endif