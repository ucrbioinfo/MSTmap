/*
 *  MSTOpt.cpp
 *  ApproxMap
 *
 *  Created by yonghui on 12/21/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "MSTOpt.h"
#include <iostream>
#include <queue>
#include <cassert>
#include "constants.h"

MSTOpt::MSTOpt(const vector<vector<double> > & _pair_wise_distances, int _num_bins, int nested_level) 
    : pair_wise_distances(_pair_wise_distances), number_of_bins(_num_bins) { 
    verbose_ = kMSTVerbose;
    assert(pair_wise_distances.size() == number_of_bins);
    nested_level_ = nested_level;
    find_opt_order();
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double MSTOpt::calculate_MST()
{
    MST.resize(number_of_bins);
    vector<double> estimated_distance(number_of_bins, numeric_limits<double>::max());
    vector<bool> visitted(number_of_bins, false);
    estimated_distance[0] = 0;
    MST[0] = 0;
    MST_lower_bound = 0 ;
    
    for (int ii = 0 ; ii < number_of_bins ; ii++)
    {
        /*get the closest vertex*/
        int closest_vertex = -1;
        double shortest_distance = numeric_limits<double>::max();
        
        /* 
        // print out something for debugging
        for (int kk = 0; kk < number_of_bins; kk++){
            cout << visitted[kk] << ' ';
        }
        cout << endl;
        cout << "number of bins:" << number_of_bins << endl;
        
        for (int kk = 0; kk < number_of_bins; kk++){
            cout << estimated_distance[kk] << endl;
        }
        cout << endl;
        */
        
        for (int jj = 0; jj < number_of_bins; jj ++)
        {
            if ((visitted[jj] == false) and (shortest_distance > estimated_distance[jj]))
            {
                closest_vertex = jj;
                shortest_distance = estimated_distance[jj];
            }
        }
        
        assert(closest_vertex != -1);
        
        visitted[closest_vertex] = true;
        MST_lower_bound = MST_lower_bound + shortest_distance;
        
        /*update the estimated_distance vector*/
        for (int jj = 0 ; jj < number_of_bins; jj++)
        {
            if ((visitted[jj] == false) and (estimated_distance[jj] > pair_wise_distances[closest_vertex][jj]))
            {
                estimated_distance[jj] = pair_wise_distances[closest_vertex][jj];
                MST[jj] = closest_vertex;
            }
        }
    }//end for
    
    return MST_lower_bound;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MSTOpt::find_opt_order()
{
    calculate_MST();
    if (verbose_) {
        cout << "finished calculating MST" << endl;
    }
    new_serialization();
    local_improvement();
    // modified by yonghui on Feb 1st, 2008
    // In order to avoid excessive long running time, 
    // block optimization is only applied to level 1 and level 2
    if (nested_level_ <= 2) {
        bool changed = true;
        while (changed) {
            changed = block_optimize();
        }
    }
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MSTOpt::local_improvement()
{
    /*Run the heuristic only if the current_lowerbound is not equal to the current upperbound*/
    if (current_upper_bound > MST_lower_bound)
    {
        bool changed = true;

        /*a vector used to hold temporary elements*/
        vector<int> tmp_storage(number_of_bins, -1);
        
        while(changed)
        {        
            /*execute the heuristic iteratively until no improvement can be made*/
            while (changed)
            {
                changed = false;
            
                /*Emulate the KL heuristic, consider all the possible 2-ops*/
                /*ii is the id of the first edge to be removed, and jj is the id of the second edge to be removed*/

                for (int ii = 1; ii < number_of_bins; ii ++)
                {
                    for (int jj = ii+1; jj < number_of_bins; jj++)
                    {
                    
                        /*step 1: determine if the new path/paths is shorter than the original path*/
                        double gone_edges_weight = pair_wise_distances[current_order[ii-1]][current_order[ii]] +  
                                                   pair_wise_distances[current_order[jj-1]][current_order[jj]];

                        /*  Block order2 
                            order 1: 1 2 3
                            order 2: 2 1 3
                            order 3: 1 3 2 
                        */                    
                        int block_order_type = 0 ; // the default block order type
                        
                        double lightest_edges_weight = gone_edges_weight;
                        
                        bool flip_1 = false; //check if the orientation of the first block is flipped or not
                        bool flip_2 = false; //check if the orientation of the second block is flipped or not 
                        bool flip_3 = false; //check if the orientation of the third block is flipped or not
                        
                        /* for block order type 1 */
                        for (int ff1 = 0 ; ff1 < 2 ; ff1 ++)
                        {
                            for (int ff2 = 0 ; ff2 < 2; ff2++)
                            {
                                for (int ff3 = 0 ; ff3 < 2 ; ff3++)
                                {
                                
                                    int end1 = current_order[ii-1];
                                    int end2 = current_order[ii];
                                    int end3 = current_order[jj-1];
                                    int end4 = current_order[jj];
                                    if (ff1 == 1) end1 = current_order[0];
                                    if (ff2 == 1) { end2 = current_order[jj-1]; end3 = current_order[ii];}
                                    if (ff3 == 1) end4 = current_order[number_of_bins-1];
                                                                
                                    double current_edges_weight = pair_wise_distances[end1][end2] + pair_wise_distances[end3][end4];
                                    
                                    if (lightest_edges_weight > current_edges_weight)
                                    {
                                        lightest_edges_weight = current_edges_weight;
                                        block_order_type = 1;
                                        if (ff1 == 0) flip_1 = false; else flip_1 = true;
                                        if (ff2 == 0) flip_2 = false; else flip_2 = true;
                                        if (ff3 == 0) flip_3 = false; else flip_3 = true;
                                    }
                                }
                            } 
                        }
                        
                        /* for block order type 2 */
                        for (int ff1 = 0 ; ff1 < 2 ; ff1 ++)
                        {
                            for (int ff2 = 0 ; ff2 < 2; ff2++)
                            {
                                for (int ff3 = 0 ; ff3 < 2 ; ff3++)
                                {
                                
                                    int end1 = current_order[jj-1];
                                    int end2 = current_order[0];
                                    int end3 = current_order[ii-1];
                                    int end4 = current_order[jj];
                                    if (ff2 == 1) end1 = current_order[ii];
                                    if (ff1 == 1) { end2 = current_order[ii-1]; end3 = current_order[0];}
                                    if (ff3 == 1) end4 = current_order[number_of_bins-1];
                                                                
                                    double current_edges_weight = pair_wise_distances[end1][end2] + pair_wise_distances[end3][end4];
                                    
                                    if (lightest_edges_weight > current_edges_weight)
                                    {
                                        lightest_edges_weight = current_edges_weight;
                                        block_order_type = 2;
                                        if (ff1 == 0) flip_1 = false; else flip_1 = true;
                                        if (ff2 == 0) flip_2 = false; else flip_2 = true;
                                        if (ff3 == 0) flip_3 = false; else flip_3 = true;
                                    }
                                }
                            } 
                        }
                                            
                        /* for block order type 3 */
                        for (int ff1 = 0 ; ff1 < 2 ; ff1 ++)
                        {
                            for (int ff2 = 0 ; ff2 < 2; ff2++)
                            {
                                for (int ff3 = 0 ; ff3 < 2 ; ff3++)
                                {
                                
                                    int end1 = current_order[ii-1];
                                    int end2 = current_order[jj];
                                    int end3 = current_order[number_of_bins-1];
                                    int end4 = current_order[ii];
                                    if (ff1 == 1) end1 = current_order[0];
                                    if (ff3 == 1) { end2 = current_order[number_of_bins-1]; end3 = current_order[jj];}
                                    if (ff2 == 1) end4 = current_order[jj-1];
                                                                
                                    double current_edges_weight = pair_wise_distances[end1][end2] + pair_wise_distances[end3][end4];
                                    
                                    if (lightest_edges_weight > current_edges_weight)
                                    {
                                        lightest_edges_weight = current_edges_weight;
                                        block_order_type = 3;
                                        if (ff1 == 0) flip_1 = false; else flip_1 = true;
                                        if (ff2 == 0) flip_2 = false; else flip_2 = true;
                                        if (ff3 == 0) flip_3 = false; else flip_3 = true;
                                    }
                                }
                            } 
                        }
                        
                        if (block_order_type == 1)
                        {
                            if (flip_1)
                            {
                                copy_order(current_order, tmp_storage,0,0,ii, false);
                                copy_order(tmp_storage,current_order,0,0,ii, true);
                            }
                            if (flip_2)
                            {
                                copy_order(current_order,tmp_storage,ii,ii,jj-ii, false);
                                copy_order(tmp_storage, current_order,ii,ii,jj-ii, true);
                            }
                            if (flip_3)
                            {
                                copy_order(current_order, tmp_storage,jj, jj, number_of_bins - jj, false);
                                copy_order(tmp_storage, current_order, jj, jj, number_of_bins - jj, true);
                            }
                        }
                        else if (block_order_type == 2)
                        {
                            copy_order(current_order,tmp_storage, 0 , 0, jj, false);
                            
                            if (flip_1)
                            {
                                copy_order(tmp_storage, current_order,0,jj-ii,ii,true);
                            }
                            else
                            {
                                copy_order(tmp_storage, current_order,0,jj-ii,ii,false);
                            }
                            
                            if (flip_2)
                            {
                                copy_order(tmp_storage, current_order,ii,0,jj-ii,true);
                            }
                            else
                            {
                                copy_order(tmp_storage, current_order,ii,0,jj-ii,false);
                            }
                            
                            if (flip_3)
                            {
                                copy_order(current_order, tmp_storage,jj, jj, number_of_bins - jj, false);
                                copy_order(tmp_storage, current_order, jj, jj, number_of_bins - jj, true);
                            }                    
                        }
                        else if (block_order_type == 3)
                        {
                        
                            if (flip_1)
                            {
                                copy_order(current_order, tmp_storage,0,0,ii,false);
                                copy_order(tmp_storage,current_order,0,0,ii,true);
                            }
                            
                            copy_order(current_order,tmp_storage,ii,ii,number_of_bins-ii,false);
                            
                            if (flip_2)
                            {
                                copy_order(tmp_storage, current_order, ii, number_of_bins - (jj-ii), jj - ii , true);
                            }
                            else
                            {
                                copy_order(tmp_storage, current_order, ii, number_of_bins - (jj-ii), jj - ii , false);
                            }
                            
                            if (flip_3)
                            {
                                copy_order(tmp_storage, current_order, jj, ii, number_of_bins - jj, true);
                            }
                            else
                            {
                                copy_order(tmp_storage, current_order, jj, ii, number_of_bins - jj, false);
                            }
                        }
                        
                        if (block_order_type != 0)
                        {
                            changed = true;
                            current_upper_bound = current_upper_bound + lightest_edges_weight  - gone_edges_weight;
                        }
                        
                    } // end for jj
                } // end for ii
                if(verbose_) {
                    cout << "current upper_bound op2:" << current_upper_bound << endl;
                }
        
            }//end while 

            changed = dis_locate();
            if(verbose_) {
                cout << "current upper_bound dislocation:" << current_upper_bound << endl;        
            }
        }
    } // end if 
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double MSTOpt::calculate_crt_upper_bound(){
    double total_cost = 0.0;
    for (int ii = 1; ii < number_of_bins; ii++) {
        total_cost = total_cost + pair_wise_distances[current_order[ii]][current_order[ii - 1]];
    }
    return total_cost;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MSTOpt::copy_order(vector<int> & order_from, 
                        vector<int> & order_to, 
                        int from_start_at, 
                        int to_start_at, 
                        int length, 
                        bool change_orientation) {
    if ( not change_orientation )
    {
        for (int ii = 0 ; ii < length; ii++)
        {
            order_to[to_start_at + ii] = order_from[from_start_at + ii];
        }
    }
    else
    {
        for (int ii = 0 ; ii < length; ii++)
        {
            order_to[to_start_at + ii] = order_from[from_start_at + length - ii - 1];
        }
    }
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool MSTOpt::dis_locate()
{

    /* first is a pointer to previous element, second is the pointer to the next element*/
    //pair<int,int> crt_order[number_of_bins]; 
    pair<int,int> * crt_order = new pair<int,int>[number_of_bins];
    
    for (int ii = 0 ; ii < number_of_bins -1 ; ii++)
    {
        crt_order[current_order[ii]].second = current_order[ii+1];
    }
    crt_order[current_order[number_of_bins -1 ]].second = -1;
    for (int ii = 1 ; ii < number_of_bins ; ii++)
    {
        crt_order[current_order[ii]].first = current_order[ii-1];
    }
    crt_order[current_order[0]].first = -1;
    int header = current_order[0];
    
    
    bool ever_changed = false;
    bool changed = true;
    while (changed)
    {
        changed = false;
        /*try to relocate the node from everywhere to everywhere else*/
        for (int ii = 0 ; ii < number_of_bins; ii ++)
        {
            double incre_cost_splicing = 0;
            int previous = crt_order[ii].first;
            int next = crt_order[ii].second;
            if (previous != -1) {
                incre_cost_splicing = incre_cost_splicing - pair_wise_distances[ii][previous];
            }
            if (next != -1 ) {
                incre_cost_splicing = incre_cost_splicing - pair_wise_distances[ii][next];
            }
            if ((previous != -1) and (next != -1)) {
                incre_cost_splicing = incre_cost_splicing + pair_wise_distances[previous][next];
            }
            
            int opt_pos = -1;
            double incre_cost_insertion = 0;
            if (header != ii) {
                incre_cost_insertion = incre_cost_insertion + pair_wise_distances[ii][header];
            }
            else {
                incre_cost_insertion = incre_cost_insertion + pair_wise_distances[ii][next];
            }
            for (int jj = 0 ; jj < number_of_bins; jj ++)
            {
                if ((jj != ii) and (jj != previous))
                {
                    double tmp_cost = 0 ;                    
                    int jj_next = crt_order[jj].second;
                    tmp_cost = tmp_cost + pair_wise_distances[jj][ii];
                    if (jj_next != -1)
                    {
                        tmp_cost = tmp_cost + pair_wise_distances[ii][jj_next];
                        tmp_cost = tmp_cost - pair_wise_distances[jj][jj_next];
                    }
                    if (tmp_cost < incre_cost_insertion)
                    {
                        incre_cost_insertion = tmp_cost;
                        opt_pos = jj;
                    }
                }
                
            }
            if ((incre_cost_insertion + incre_cost_splicing) < ZERO_MINUS )
            {
                /*Do the relocation*/
                /*Do the splicing first*/
                ever_changed = true;
                changed = true;
                if (previous != -1)
                {
                    crt_order[previous].second = crt_order[ii].second;
                }
                if (next != -1)
                {
                    crt_order[next].first = crt_order[ii].first;
                }
                if (previous == -1)
                {
                    header = next;
                }
                /*Do the insertion next*/
                if (opt_pos == -1)
                {
                    crt_order[ii].second = header;
                    crt_order[ii].first = -1;
                    crt_order[header].first = ii;
                    header = ii;
                }
                else
                {
                    int jj_next = crt_order[opt_pos].second;
                    crt_order[opt_pos].second = ii;
                    crt_order[ii].first = opt_pos;
                    crt_order[ii].second = jj_next;
                    if (jj_next != -1)
                    {
                        crt_order[jj_next].first = ii;
                    }
                }
                current_upper_bound = current_upper_bound + incre_cost_insertion + incre_cost_splicing;
            }            
            
        }
        
    };
    
    current_order[0] = header;
    int nxt_id = crt_order[header].second;
    for (int ii = 1 ; ii < number_of_bins  ; ii++)
    {
        current_order[ii] = nxt_id;
        nxt_id = crt_order[nxt_id].second;
    }
    
    double crt_cost = calculate_crt_upper_bound();
    assert(crt_cost < current_upper_bound + ZERO_PLUS);
    assert(crt_cost > current_upper_bound - ZERO_PLUS);
    
    delete [] crt_order;
    return ever_changed;
    
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

vector<int> MSTOpt::get_the_longest_path()
{
    vector<vector<int> > graph_MST(number_of_bins);
    for (int ii = 1 ; ii < number_of_bins; ii++)
    {
        int jj = MST[ii];
        graph_MST[ii].push_back(jj);
        graph_MST[jj].push_back(ii);
    }
    vector<int> longest_path;
        
    int longest_distance = -1;

    /*Run BFS starting from every node*/
    for (int ii = 0 ; ii < number_of_bins; ii ++)
    {
        vector<bool> visitted(number_of_bins, false);
        vector<int> depth(number_of_bins, 0);
        visitted[ii] = true;
        vector<int> BFS_tree(number_of_bins);
        BFS_tree[ii] = -1;
        queue<int> fifoqueue;
        int last_out = -1;
        fifoqueue.push(ii);
        while ( not fifoqueue.empty())
        {
            last_out = fifoqueue.front();
            fifoqueue.pop();
            for (vector<int>::iterator iter1 = (graph_MST[last_out]).begin(); 
                 iter1 != (graph_MST[last_out]).end(); 
                 iter1 ++ )
            {
                if (not visitted[*iter1])
                {
                    fifoqueue.push(*iter1);
                    visitted[*iter1] = true;
                    depth[*iter1] = depth[last_out] + 1;
                    BFS_tree[*iter1] = last_out;
                }
            }
        }
        if (depth[last_out] > longest_distance)
        {
            longest_distance = depth[last_out];
            longest_path.clear();
            longest_path.push_back(last_out);
            last_out = BFS_tree[last_out];
            while (last_out != -1)
            {
                longest_path.push_back(last_out);
                last_out = BFS_tree[last_out];
            }
        }
    }
    return longest_path;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double MSTOpt::new_serialization()
{
    vector<int> longest_path = get_the_longest_path();
    
    vector<pair<int,int> > crt_path(number_of_bins,make_pair(-1,-1));
    if(verbose_){
        cout << "the length of the longest path:" << longest_path.size() << endl;
    }
    int header;
    header = longest_path[0];
    crt_path[header].first=-1;
    crt_path[header].second = -1;
    int crt_node = header;
    for ( int ii = 1; ii < longest_path.size(); ii++)
    {
        int jj = longest_path[ii];
        crt_path[jj].first = crt_node;
        crt_path[jj].second = -1;
        crt_path[crt_node].second = jj;
        crt_node = jj;
    }
    
    /*get the set of unvisitted nodes*/

    vector<bool> visitted(number_of_bins, false);
    for (int ii = 0 ; ii < longest_path.size(); ii ++)
    {
        visitted[longest_path[ii]] = true;
    }
    
    /*now, insert the unvisited nodes into the backbone one by one*/
    
    for (int ii = 0 ; ii < number_of_bins; ii++)
    {
        if (visitted[ii] == false)
        {    
            visitted[ii] = true;
            double opt_cost = pair_wise_distances[ii][header] ;
            int opt_posi = -1;
            int crt_posi = header;
            while (crt_posi != -1)
            {
                double tmp_cost;
                tmp_cost = pair_wise_distances[crt_posi][ii];
                int kk = crt_path[crt_posi].second;
                if (kk != -1)
                {
                    tmp_cost = tmp_cost + pair_wise_distances[ii][kk];
                    tmp_cost = tmp_cost - pair_wise_distances[crt_posi][kk];
                }
                if (tmp_cost < opt_cost)
                {
                    opt_cost = tmp_cost;
                    opt_posi = crt_posi;
                }
                crt_posi = (crt_path[crt_posi]).second;
            }
            if (opt_posi == -1)
            {
                crt_path[ii].second = header;
                crt_path[ii].first = -1;
                crt_path[header].first = ii;
                header = ii;
            }
            else
            {
                int kk = (crt_path[opt_posi]).second;
                crt_path[ii].second = kk;
                crt_path[ii].first = opt_posi;
                crt_path[opt_posi].second = ii;
                if (kk != -1)
                {
                    crt_path[kk].first = ii;
                }
            }
            
            
        }
    }

    current_order.resize(number_of_bins);
    current_order[0] = header;
    int crt_node_id = crt_path[header].second;
    for (int ii = 1; ii < number_of_bins; ii ++)
    {
        current_order[ii] = crt_node_id;
        crt_node_id = (crt_path[crt_node_id]).second;
    }

    current_upper_bound = 0 ; 
    for (int ii = 1 ; ii < number_of_bins; ii ++)
    {
        current_upper_bound = current_upper_bound + pair_wise_distances[current_order[ii]][current_order[ii-1]];
    }

    
    cost_after_initialization = current_upper_bound;

    return current_upper_bound;

};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

MSTOpt::Block_Chain MSTOpt::break_into_blocks(){
    int bs = 0;
    Block_Chain bc;
    for (int ii = 1; ii < number_of_bins; ii++) {
        // check the validity of the edge (current_order[ii], current_order[ii-1])
        vector<int> nodes_before;
        vector<int> nodes_after;
        int first_node = current_order[ii - 1];
        int second_node = current_order[ii];
        double cost_of_edge = pair_wise_distances[first_node][second_node];
        for (int jj = 0; jj < ii - 1; jj++) {
            nodes_before.push_back(current_order[jj]);
        }
        for (int jj = ii + 1; jj < number_of_bins; jj++) {
            nodes_after.push_back(current_order[jj]);
        }
        bool valid = true;
        for (vector<int>::iterator iter1 = nodes_before.begin(); iter1 != nodes_before.end(); ++iter1) {
            if (pair_wise_distances[*iter1][second_node] < cost_of_edge) {
                valid = false;
                break;
            }
        }
        for (vector<int>::iterator iter1 = nodes_after.begin(); iter1 != nodes_after.end(); ++iter1) {
            if (pair_wise_distances[first_node][*iter1] < cost_of_edge) {
               valid = false;
               break;
            }
        }
        if (not valid) { // a new block
            vector<int> block_markers;
            for (int jj = bs; jj < ii; jj++) {
                block_markers.push_back(current_order[jj]);
            }
            Block crt_block = {true, block_markers, block_markers.size(), 
                               block_markers[0], block_markers[block_markers.size() - 1],
                               -1, -1};
            bc.bs.push_back(crt_block);
            bs = ii;
        } 
    }
    // assemble the last block
    
    vector<int> last_block_markers;
    for (int jj = bs; jj < number_of_bins; jj++) {
        last_block_markers.push_back(current_order[jj]);
    }
    Block last_block = {true, last_block_markers, 
                        last_block_markers.size(), 
                        last_block_markers[0], last_block_markers[last_block_markers.size() - 1],
                        -1, -1};
    bc.bs.push_back(last_block);
    
    // do some correctness check
    int total_bins = 0;
    for (int ii = 0; ii < bc.bs.size(); ii++) {
        total_bins = total_bins + bc.bs[ii].markers.size();
        bc.bs[ii].p_b = ii - 1;
        bc.bs[ii].n_b = ii + 1;
    }
    assert(total_bins == number_of_bins);
    bc.bs[bc.bs.size() - 1].n_b = -1;
    bc.header = 0;
    return bc;
    
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double MSTOpt::block_cost(const Block_Chain& bc){
    double total_cost = 0.0;
    int pre_block = bc.header;
    int crt_block = bc.bs[pre_block].n_b;
    while (crt_block != -1) {
        int from = -1;
        int to = -1;
        if (bc.bs[pre_block].orientation) {
            from = bc.bs[pre_block].last;
        } else {
            from = bc.bs[pre_block].first;
        }
        
        if (bc.bs[crt_block].orientation) {
            to = bc.bs[crt_block].first;
        } else {
            to = bc.bs[crt_block].last;
        }
        total_cost = total_cost + pair_wise_distances[from][to];
        pre_block = crt_block;
        crt_block = bc.bs[crt_block].n_b;
    }
    return total_cost;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MSTOpt::block_fix_orientation(Block_Chain& bc) {
    bool changed = true;
    double initial_cost = block_cost(bc);
    double total_cost_reduction = 0.0;
    while (changed) {
        double cost_reduction_iter = 0.0;
        int crt_index = bc.header;
        while (crt_index != -1){
            Block& crt_block = bc.bs[crt_index];
            double cost_reduction = 0.0;
            if (crt_block.p_b != -1) {
                const Block& p_block = bc.bs[crt_block.p_b];
                int from = -1;
                if (p_block.orientation) {
                    from = p_block.last;
                } else {
                    from = p_block.first;
                }
                if (crt_block.orientation) {
                    cost_reduction = cost_reduction + pair_wise_distances[from][crt_block.first] -
                                     pair_wise_distances[from][crt_block.last];
                } else {
                    cost_reduction = cost_reduction + pair_wise_distances[from][crt_block.last] -
                                     pair_wise_distances[from][crt_block.first];
                }
            }
            if (crt_block.n_b != -1) {
                const Block& n_block = bc.bs[crt_block.n_b];
                int to = -1;
                if (n_block.orientation) {
                    to = n_block.first;
                } else {
                    to = n_block.last;
                }
                if (crt_block.orientation) {
                    cost_reduction = cost_reduction + pair_wise_distances[to][crt_block.last] -
                                     pair_wise_distances[to][crt_block.first];
                } else {
                    cost_reduction = cost_reduction + pair_wise_distances[to][crt_block.first] -
                                     pair_wise_distances[to][crt_block.last];
                }

            }
            if (cost_reduction > 0.0) {
                cost_reduction_iter = cost_reduction_iter + cost_reduction;
                total_cost_reduction = total_cost_reduction + cost_reduction;
                crt_block.orientation = (not crt_block.orientation);
            }
            crt_index = crt_block.n_b;
        }
        if (cost_reduction_iter > ZERO_PLUS) {
            changed = true;
        } else {
            changed = false;
        }
    }
    double new_cost = block_cost(bc);
    assert(new_cost + total_cost_reduction < initial_cost + ZERO_PLUS);
    assert(new_cost + total_cost_reduction > initial_cost - ZERO_PLUS);
    if(verbose_) {
        cout << "total cost reduction:" << total_cost_reduction << endl;
    }
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool MSTOpt::block_optimize_iteration(Block_Chain& bc){
    // try to relocate every block to everywhere else 
    // return true if this operation improves the order
    int total_blocks = bc.bs.size();
    if(verbose_){
        cout << "number of blocks:" << total_blocks << endl;
    }
    if (total_blocks < 3) {
        return false;
    }
    // print_bc(bc);
    double original_cost = block_cost(bc);
    double total_incremental_cost = 0.0;
    
    vector<Block>& bs = bc.bs;
    for (int ii = 0; ii < total_blocks; ii++) {
        Block& crt_block = bs[ii];
        double cost_reduction = 0.0;
        // splice the block out of the chain
        if (crt_block.n_b != -1) {
            Block& next_block = bs[crt_block.n_b];
            next_block.p_b = crt_block.p_b;
            int from = -1;
            int to = -1;
            if (crt_block.orientation) {
                from = crt_block.last;
            } else {
                from = crt_block.first;
            }
            if (next_block.orientation) {
                to = next_block.first;
            } else {
                to = next_block.last;
            }
            cost_reduction = cost_reduction + pair_wise_distances[from][to];
        }
        if (crt_block.p_b != -1) {
            Block& previous_block = bs[crt_block.p_b];
            previous_block.n_b = crt_block.n_b;
            int from = -1;
            int to = -1;
            if (previous_block.orientation) {
                from = previous_block.last;
            } else {
                from = previous_block.first;
            }
            if (crt_block.orientation) {
                to = crt_block.first;
            } else {
                to = crt_block.last;
            }
            cost_reduction = cost_reduction + pair_wise_distances[from][to];
        }
        if ((crt_block.n_b != -1) and (crt_block.p_b != -1)) {
            Block& previous_block = bs[crt_block.p_b];
            Block& next_block = bs[crt_block.n_b];
            int from = -1;
            int to = -1;
            if (previous_block.orientation) {
                from = previous_block.last;
            } else {
                from = previous_block.first;
            }
            if (next_block.orientation) {
                to = next_block.first;
            } else {
                to = next_block.last;
            }
            cost_reduction = cost_reduction - pair_wise_distances[from][to];
        }
        if (ii == bc.header) {
            bc.header = crt_block.n_b;
        }
        crt_block.p_b = -1;
        crt_block.n_b = -1;

        // cout << "one_block is spliced out" << endl;
        // print_bc(bc);
        
        // now try to find the optimum position for crt_block
        // first try to put the block at the begining
        int opt_pos = -1;
        bool opt_orientation = true;
        double opt_cost = numeric_limits<double>::max();
        // compare two different orientations
        int to = -1;
        if (bs[bc.header].orientation) {
            to = bs[bc.header].first;
        } else {
            to = bs[bc.header].last;
        }
        
        if (pair_wise_distances[crt_block.first][to] < pair_wise_distances[crt_block.last][to]) {
            opt_cost = pair_wise_distances[crt_block.first][to];
            opt_orientation = false;
        } else {
            opt_cost = pair_wise_distances[crt_block.last][to];
            opt_orientation = true;
        }
        
        for (int jj = 0; jj < total_blocks; jj++) {
            if (jj == ii) continue;
            // check the forward direction first
            double additional_cost1 = 0.0;
            double additional_cost2 = 0.0;
            Block& block_jj = bs[jj];
            int from = -1;
            if (block_jj.orientation) {
                from = block_jj.last;
            } else {
                from = block_jj.first;
            }
            additional_cost1 = additional_cost1 + pair_wise_distances[from][crt_block.first];
            additional_cost2 = additional_cost2 + pair_wise_distances[from][crt_block.last];
            if (block_jj.n_b != -1) {
                Block& next_jj = bs[block_jj.n_b];
                int to = -1;
                if (next_jj.orientation) {
                    to = next_jj.first;
                } else {
                    to = next_jj.last;
                }
                additional_cost1 = additional_cost1 - pair_wise_distances[from][to];
                additional_cost2 = additional_cost2 - pair_wise_distances[from][to];
                additional_cost1 = additional_cost1 + pair_wise_distances[crt_block.last][to];
                additional_cost2 = additional_cost2 + pair_wise_distances[crt_block.first][to];
            }
            bool jj_orientation = true;
            double additional_cost;
            if (additional_cost1 < additional_cost2) {
                jj_orientation = true;
                additional_cost = additional_cost1;
            } else {
                additional_cost = additional_cost2;
                jj_orientation = false;
            }
            
            if (additional_cost < opt_cost) {
                opt_cost = additional_cost;
                opt_pos = jj;
                opt_orientation = jj_orientation;
            }
        }
        // now place crt_block to it's position
        crt_block.orientation = opt_orientation;
        if (opt_pos == -1) { // the block should be placed at the header of the chain
             crt_block.n_b = bc.header;
             bs[bc.header].p_b = ii;
             bc.header = ii;
        } else {
            int next_pos = bs[opt_pos].n_b;
            crt_block.n_b = next_pos;
            bs[opt_pos].n_b = ii;
            crt_block.p_b = opt_pos;
            if (next_pos != -1) {
                bs[next_pos].p_b = ii;
            }
        }
        total_incremental_cost = total_incremental_cost + opt_cost - cost_reduction;
        // cout << "the block has been put back" << endl;
        // print_bc(bc);
    }
    
    double new_cost = block_cost(bc);
    assert(new_cost < original_cost + total_incremental_cost + ZERO_PLUS);
    assert(new_cost > original_cost + total_incremental_cost - ZERO_PLUS);
    if(verbose_){
        cout << "total incremental cost" << total_incremental_cost << endl;
    }
    if (total_incremental_cost < ZERO_MINUS) {
        return true;
    } else {
        return false;
    }
}; 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MSTOpt::copy_over_order(const Block_Chain& bc){
    int crt_block_id = bc.header;
    int crt_marker_id = 0;
    while (crt_block_id != -1){
        const Block& crt_block = bc.bs[crt_block_id];
        if (crt_block.orientation) {
            for (int ii = 0; ii < crt_block.size; ii++){
                current_order[crt_marker_id] = crt_block.markers[ii];
                crt_marker_id++;
            }
        } else {
            for (int ii = crt_block.size - 1; ii >= 0; ii--) {
                current_order[crt_marker_id] = crt_block.markers[ii];
                crt_marker_id++;
            }
        }
        crt_block_id = crt_block.n_b;
    }
    assert(crt_marker_id == number_of_bins);
    current_upper_bound = calculate_crt_upper_bound();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MSTOpt::contract_blocks(const Block_Chain& bc, vector<vector<double> >& block_distances) {
    int number_of_blocks = bc.bs.size();
    block_distances.resize(number_of_blocks);
    for (int ii = 0; ii < number_of_blocks; ii++) {
        block_distances[ii].resize(number_of_blocks);
    }
    // compute the pair-wise distance
    for (int ii = 0; ii < number_of_blocks; ii++) {
        for (int jj = ii + 1; jj < number_of_blocks; jj++) {
            const Block& b_ii = bc.bs[ii];
            int ii_size = b_ii.size;
            const Block& b_jj = bc.bs[jj];
            int jj_size = b_jj.size;
            double shortest_dist = numeric_limits<double>::max();
            for (int kk = 0; kk < ii_size; kk++) {
                for (int ll = 0; ll < jj_size; ll++) {
                    if (pair_wise_distances[b_ii.markers[kk]][b_jj.markers[ll]] < shortest_dist) {
                        shortest_dist = pair_wise_distances[b_ii.markers[kk]][b_jj.markers[ll]];
                    }
                }
            }
            block_distances[ii][jj] = shortest_dist;
            block_distances[jj][ii] = shortest_dist;
        }
    }
    for (int ii = 0; ii < number_of_blocks; ii++) {
        block_distances[ii][ii] = 0;
    }
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MSTOpt::print_bc(const Block_Chain& bc){
    int crt_ind = bc.header;
    while (crt_ind != -1) {
        const Block& crt_block = bc.bs[crt_ind];
        cout << crt_block.size << '\t';
        crt_ind = bc.bs[crt_ind].n_b;
    }
    cout << endl;
    // print the distances in the order of the chain
    int crt_ind1 = bc.header;
    while (crt_ind1 != -1) {
        const Block& b1 = bc.bs[crt_ind1];
        int from = -1;
        if (b1.orientation) {
            from = b1.last;
        } else {
            from = b1.first;
        }
        int crt_ind2 = bc.header;
        while (crt_ind2 != -1) {
            const Block& b2 = bc.bs[crt_ind2];
            int to = -1;
            if (b2.orientation) {
                to = b2.first;
            } else {
                to = b2.last;
            }
            if (crt_ind1 != crt_ind2) {
                cout << pair_wise_distances[from][to] << '\t';
            } else {
                cout << "0.0\t";
            }
            
            crt_ind2 = bc.bs[crt_ind2].n_b;
        }
        cout << endl;
        crt_ind1 = bc.bs[crt_ind1].n_b;
    }
    
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool MSTOpt::block_optimize() {
    Block_Chain bc = break_into_blocks();
    int number_of_blocks = bc.bs.size();
    
    // if there is no reduction of the problem size. return false
    if (number_of_blocks == number_of_bins) {
        return false;
    }
    
    vector<vector<double> > block_distance;
    contract_blocks(bc, block_distance);
    
    /*
    // print out some debugging information
    cout << "print out debugging information" << endl;
    print_bc(bc);
    for (int ii = 0; ii < number_of_blocks; ii++) {
        for (int jj = 0; jj < number_of_blocks; jj++) {
            cout << block_distance[ii][jj] << "\t";
        }
        cout << endl;
    }
    double ttl_b_dist = 0.0;
    for (int ii = 1; ii < number_of_blocks; ii++) {
        ttl_b_dist = ttl_b_dist + block_distance[ii - 1][ii];
    }
    cout << "initial distance:" << ttl_b_dist << endl;
    */

    MSTOpt b_opt(block_distance, number_of_blocks, nested_level_ + 1);
    vector<int> b_opt_order;
    vector<int> b_mst;
    double b_lowerbound;
    double b_upperbound;
    double b_initial_cost;
    b_opt.Opt_Order(b_opt_order, b_mst, b_lowerbound, b_upperbound, b_initial_cost);
    
    /*
    cout << "print out the order after optimization " << endl;
    for (int ii = 0; ii < number_of_blocks; ii++) {
        for (int jj = 0; jj < number_of_blocks; jj++) {
            cout << block_distance[b_opt_order[ii]][b_opt_order[jj]] << "\t";
        }
        cout << endl;
    }
    double ttl_b_dist2 = 0.0;
    for (int ii = 1; ii < number_of_blocks; ii++) {
        ttl_b_dist2 = ttl_b_dist2 + block_distance[b_opt_order[ii - 1]][b_opt_order[ii]];
    }
    cout << "total distance after optimization:" << ttl_b_dist2 << endl;
    */

    
    // re-organize the blocks
    for (int ii = 1; ii < number_of_blocks; ii++) {
        bc.bs[b_opt_order[ii]].p_b = b_opt_order[ii - 1];
    }
    for (int ii = 0; ii < number_of_blocks - 1; ii++) {
        bc.bs[b_opt_order[ii]].n_b = b_opt_order[ii + 1];
    }
    bc.bs[b_opt_order[0]].p_b = -1;
    bc.bs[b_opt_order[number_of_blocks - 1]].n_b = -1;
    bc.header = b_opt_order[0];
    
    block_fix_orientation(bc);
    
    bool changed = true;
    while (changed) {
        changed = block_optimize_iteration(bc);
    }
    
    vector<int> backup_order;
    backup_order = current_order;
    double backup_upper_bound = current_upper_bound;
    copy_over_order(bc);
    
    local_improvement();
    
    double new_upper_bound = calculate_crt_upper_bound();
    
    if(verbose_){
        cout << "backup_upper_bound: " << backup_upper_bound << "new_upper_bound: " << new_upper_bound << endl; 
    }
    
    bool improved = false;
    if (new_upper_bound < backup_upper_bound - ZERO_PLUS) {
        improved = true;
    } else {
        improved = false;
        current_upper_bound = backup_upper_bound;
        current_order = backup_order;
    }
    return improved;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MSTOpt::Opt_Order(vector<int>& out_order,
                       vector<int>& _MST,
                       double & _lowerbound, 
                       double & _upper_bound, 
                       double & _cost_after_initialization) {
    out_order = current_order;
    _MST = MST;
    _lowerbound = MST_lower_bound;
    _upper_bound = current_upper_bound;
    _cost_after_initialization = cost_after_initialization;
    sanity_check();
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// perform a sanity check of the result obtained.
void MSTOpt::sanity_check(){
    /*Consistency checking*/
    
    /*Step 1: check that the current order is valid*/
    vector<int> tmp(number_of_bins,0);
    for (int ii = 0 ; ii < number_of_bins; ii++)
    {
        tmp[current_order[ii]] = tmp[current_order[ii]] + 1;
    }
    for (int ii = 0; ii < number_of_bins; ii++)
    {
        if (tmp[ii] != 1 )
        {
            cout << "ERROR, not a valid order. " << ii << " appeared " << tmp[ii] << "times " << endl;
            assert(false); // crash the program on error
        }
    }
    
    double epsilon = 0.000001;
    /*Step 2: check that the current upper_bound is also correct*/
    double tmp_upper_bound = 0 ;
    for (int ii = 1 ; ii < number_of_bins; ii ++)
    {
        tmp_upper_bound = tmp_upper_bound + pair_wise_distances[current_order[ii]][current_order[ii-1]];
    }
    
    if (not ((tmp_upper_bound - current_upper_bound < epsilon) and (tmp_upper_bound - current_upper_bound > -epsilon))){
        cout << "ERROR, current_upper_bound is not correct" << endl;
        cout << "tmp_upper_bound:" << tmp_upper_bound << " upper_bound:" << current_upper_bound << endl;
        assert(false); // crash the program on error
    }
    
    /*Step 3: make sure that the lowerbound is calculated correctly*/
    double tmp_lower_bound = 0 ; 
    for (int ii = 1 ; ii < number_of_bins; ii ++)
    {
        tmp_lower_bound = tmp_lower_bound + pair_wise_distances[MST[ii]][ii];
    }
    if (not ((tmp_lower_bound - MST_lower_bound < epsilon) and (tmp_lower_bound - MST_lower_bound > -epsilon))) {
        cout << "ERROR, the lowerbound is not correct" << endl;
        cout << "tmp_lower_bound:" << tmp_lower_bound << " lower_bound:" << MST_lower_bound << endl;
        assert(false); // crash the program on error
    }
   
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
