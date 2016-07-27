/*
 *  MSTOpt.h
 *  ApproxMap
 *
 *  Created by yonghui on 12/21/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef MSTOpt_HEADER
#define MSTOpt_HEADER

#include <vector>
#include <limits>

using namespace std;

class MSTOpt{
    public:
        MSTOpt(const vector<vector<double> >& _pair_wise_distances, int num_bins, int nested_level);
        void Opt_Order(vector<int>& opt_order, 
                       vector<int>& _MST,
                       double& _lowerbound, 
                       double& _upper_bound, 
                       double& _cost_after_initialization); 
    private:
        // ----------------------------------
        // Block optimization data structures 
        // ----------------------------------
        struct Block{
            bool orientation; // is the orientation of the block
            vector<int> markers;
            int size;
            int first; // first element in the block
            int last; // last element in the block
            // is the index of the previous block in an array. It is set to -1 if current block is the first 
            // block in the chain
            int p_b; 
            // is the index of the next block in an array. It is set to -1 if the current block is the last 
            // block in the chain
            int n_b; 
        };
        
        struct Block_Chain{
            vector<Block> bs;
            int header;
        };


        // -----------------
        // Private functions 
        // -----------------
        void sanity_check();
        
        void find_opt_order();

        /* Prim's minimum spanning tree algorithm */
        // The function will compute the MST and will store it in the MST vector
        double calculate_MST();
                
        double new_serialization();

        /*Execute the OP2 operation to improve the ordering obtained by the approximation algorithm*/
        void local_improvement();

        vector<int> get_the_longest_path();

        double calculate_crt_upper_bound();
        
        void copy_order(vector<int> & order_from, 
                        vector<int> & order_to, 
                        int from_start_at, 
                        int to_start_at, 
                        int length, 
                        bool change_orientation);
        
        bool dis_locate();

        // a set of functions for block optimization
        Block_Chain break_into_blocks();
        // one iteration of block optimization. The function will return true if the order improves
        bool block_optimize_iteration(Block_Chain& bc); 
        void block_fix_orientation(Block_Chain& bc);
        double block_cost(const Block_Chain& bc);
        // block optimization. The function will return true if it improves the order

        void copy_over_order(const Block_Chain& bc);
        void print_bc(const Block_Chain& bc);
        void contract_blocks(const Block_Chain& bc, vector<vector<double> >& block_distances);
        bool block_optimize();
        // ---------------------
        // private data section:
        // ---------------------
        int number_of_bins;
        const vector<vector<double> >& pair_wise_distances;
        vector<int> current_order; //concatenation of the markers in the order of the path
        
        vector<int> MST;
        double MST_lower_bound;
        double current_upper_bound;
        double cost_after_initialization;
        
        // added by yonghui on Feb 1 
        int nested_level_;
        
        bool verbose_;
};
#endif
