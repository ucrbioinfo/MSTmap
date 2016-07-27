/*
 *  genetic_map_RIL.cpp
 *  ApproxMap
 *
 *  Created by yonghui on 12/13/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "genetic_map_RIL.h"
#include <cstdlib>

void genetic_map_RIL::gen_raw_prob_data(){
    raw_prob_data_.resize(number_of_loci);
    for (int ii = 0; ii < number_of_loci; ii++) {
        raw_prob_data_[ii].resize(number_of_individual);
    }
    for (int ii = 0; ii < number_of_loci; ii++) {
        for (int jj = 0; jj < number_of_individual; jj++) {
            if (raw_mapping_data[ii][jj] == 'A') {
                raw_prob_data_[ii][jj].A = 1.0;
                raw_prob_data_[ii][jj].B = 0.0;
                raw_prob_data_[ii][jj].AB = 0.0;
                raw_prob_data_[ii][jj].missing = false;
            } else if (raw_mapping_data[ii][jj] == 'B') {
                raw_prob_data_[ii][jj].A = 0.0;
                raw_prob_data_[ii][jj].B = 1.0;
                raw_prob_data_[ii][jj].AB = 0.0;
                raw_prob_data_[ii][jj].missing = false;
            } else if (raw_mapping_data[ii][jj] == 'X') {
                raw_prob_data_[ii][jj].A = 0.0;
                raw_prob_data_[ii][jj].B = 0.0;
                raw_prob_data_[ii][jj].AB = 1.0;
                raw_prob_data_[ii][jj].missing = false;                
            } else if (raw_mapping_data[ii][jj] == '-') {
                raw_prob_data_[ii][jj].A = 1.0 / 3.0;
                raw_prob_data_[ii][jj].B = 1.0 / 3.0;
                raw_prob_data_[ii][jj].AB = 1.0 / 3.0;
                raw_prob_data_[ii][jj].missing = true;                
            } else {
                cout << "ERROR! invalid genotype" << endl;
                assert(false);
            }
        }
    }
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void genetic_map_RIL::calculate_pair_wise_distance(){
    pair_wise_distances.resize(number_of_loci);
    for (int ii = 0 ; ii < number_of_loci; ii++)
    {
        pair_wise_distances[ii].resize(number_of_loci);
    }

    cout << "start calculating pair-wise distance" << time(NULL);
    for (int ii = 0 ; ii < number_of_loci; ii++)
    {
        cout << "finished one marker" << endl;
        for (int jj = ii + 1 ; jj < number_of_loci; jj++)
        {
            RIL_dist_cal ril_dist(generation_index_, raw_prob_data_[ii], raw_prob_data_[jj]);
            double distance_ii_jj = ril_dist.Dist();
            pair_wise_distances[ii][jj] = distance_ii_jj;
            pair_wise_distances[jj][ii] = distance_ii_jj;
        }
    }    
    cout << "finished calculating the pair-wise distance" << endl;
    for (int ii = 0; ii < number_of_loci; ii++) {
        pair_wise_distances[ii][ii] = 0.0;
    }
    cout << time(NULL) << endl;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

linkage_group_RIL* genetic_map_RIL::construct_linkage_group(int group_id){
    int _number_of_bins = (linkage_group_bins[group_id]).size();
    int _number_of_individuals = number_of_individual;

    /*Store the probability for each allele to be A*/
    vector<vector<allel_state> > _raw_data ;
    vector<pair<int, int> >  _missing_data;
    vector<int> _current_order;
    
    _raw_data.resize(_number_of_bins);
    for (int ii = 0 ; ii < _number_of_bins; ii++) {
        _raw_data[ii] = raw_prob_data_[linkage_group_bins[group_id][ii][0]];
        for (int jj = 0; jj < _number_of_individuals; jj++) {
            if (raw_mapping_data[linkage_group_bins[group_id][ii][0]][jj] == '-') {
                /*ii is the id for the marker, and jj is the id for the individual*/
                _missing_data.push_back(make_pair(ii,jj)); 
            }
        }
    }
    for (int ii = 0 ; ii < _number_of_bins; ii ++) {
        _current_order.push_back(ii);
    }
    linkage_group_RIL * to_be_returned = new linkage_group_RIL(_number_of_bins, 
                                                               _number_of_individuals, 
                                                               generation_index_,
                                                               df_,
                                                               _raw_data, 
                                                               _current_order, 
                                                               _missing_data);
    return to_be_returned;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

linkage_group_RIL* genetic_map_RIL::construct_linkage_group_whole_map(){
    int _number_of_bins = number_of_loci;
    int _number_of_individuals = number_of_individual;
    
    /*Store the probability for each allele to be A*/
    vector<vector<allel_state> > _raw_data;
    vector<pair<int, int> >  _missing_data;
    vector<int> _current_order;
    
    _raw_data.resize(_number_of_bins);
    for (int ii = 0 ; ii < _number_of_bins; ii++) {
        _raw_data[ii] = raw_prob_data_[ii];
        for (int jj = 0; jj < _number_of_individuals; jj++) {
            if (raw_mapping_data[ii][jj] == '-') {
                _missing_data.push_back(make_pair(ii,jj)); /*ii is the id for the marker, and jj is the id for the individual*/
            }
        }
    }
    for (int ii = 0 ; ii < _number_of_bins; ii ++) {
        _current_order.push_back(ii);
    }
    linkage_group_RIL* to_be_returned = new linkage_group_RIL(_number_of_bins, 
                                                             _number_of_individuals, 
                                                             generation_index_,
                                                             df_,
                                                             _raw_data, 
                                                             _current_order, 
                                                             _missing_data);
    return to_be_returned;
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void genetic_map_RIL::generate_map()
{
   
    gen_raw_prob_data();
    
    const char* ppl_type = population_type.c_str();
    assert(population_type.length() >= 4);
    assert(ppl_type[0] == 'R');
    assert(ppl_type[1] == 'I');
    assert(ppl_type[2] == 'L');    
    
    generation_index_ = atoi(&(ppl_type[3]));

    pair_wise_distances.resize(number_of_loci);
    for (int ii = 0 ; ii < number_of_loci; ii++)
    {
        pair_wise_distances[ii].resize(number_of_loci, 0.0);
    }    
    
    /*
      if the total number of missing observations exceeds certain threshold, 
      we need to estimate the missing data before clustering
    */
    if ((total_number_of_missing_obs >= ESTIMATION_BEFORE_CLUSTERING * number_of_loci * number_of_individual) and
        (estimation_before_clustering)) {
        linkage_group_RIL * linkage_group_whole_map = construct_linkage_group_whole_map();
        linkage_group_whole_map->order_markers();
        const vector<vector<double> > & new_dist = linkage_group_whole_map-> get_pair_wise_distance();
        for (int ii = 0 ; ii < number_of_loci; ii++)
        {
            for (int jj = 0 ; jj < number_of_loci; jj++)
            {
                pair_wise_distances[ii][jj] = new_dist[ii][jj];
            }
        }
        delete linkage_group_whole_map;
    } else {
        cout << "calculating the pair-wise hamming distance" << endl;
        calculate_pair_wise_distance();
        cout << "finished calculating the pair-wise hamming distance" << endl;
    }

    cluster();
    cout << "found " << number_of_connected_components << " connected components" << endl;
    
    condense_markers_into_bins();
    
    orders.resize(number_of_connected_components);
    upperbounds.resize(number_of_connected_components);
    lowerbounds.resize(number_of_connected_components);
    approx_bounds.resize(number_of_connected_components);
    distance_between_adjacent_pairs.resize(number_of_connected_components);
    
    for (int ii = 0 ; ii < number_of_connected_components; ii++)
    {
        linkage_group_RIL * current_linkage_group = construct_linkage_group(ii);

        current_linkage_group->order_markers();
        current_linkage_group->return_order(orders[ii], 
                                            lowerbounds[ii], 
                                            upperbounds[ii], 
                                            approx_bounds[ii], 
                                            distance_between_adjacent_pairs[ii]);
        current_linkage_group->dump();
        delete current_linkage_group;
        cout <<"finished the " << ii+1 << " linkage group" << endl; 
    }
    
    // Added by Yonghui on Oct 20, 2007 
    // The last step is to condense adjacent bins if they are too close to each other
    condense_bin();
    
    dump_connected_components_edges();
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
