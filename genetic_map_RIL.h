/*
 *  genetic_map_RIL.h
 *  ApproxMap
 *
 *  Created by yonghui on 12/13/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef GENETIC_MAP_RIL_HEADER
#define GENETIC_MAP_RIL_HEADER

#include "genetic_map_DH.h"
#include "linkage_group_RIL.h"

class genetic_map_RIL: public genetic_map {
    public:
        genetic_map_RIL():genetic_map(){};
      
        virtual void generate_map();
    private:
        void gen_raw_prob_data();
        void calculate_pair_wise_distance();

        /*generate a linkage group*/
        linkage_group_RIL* construct_linkage_group(int group_id);
        linkage_group_RIL* construct_linkage_group_whole_map();
        
        // ----------------------------
        // private data members section
        // ----------------------------
        int generation_index_;
        vector<vector<allel_state> > raw_prob_data_;
};
#endif