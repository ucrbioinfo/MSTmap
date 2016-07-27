/*
 *  single_mapping_population_raw_data.h
 *  ApproxMap
 *
 *  Created by yonghui on 4/7/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef single_mapping_population_raw_data_header
#define single_mapping_population_raw_data_header

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <ctime>
#include "constants.h"
#include <queue>
#include "linkage_group_DH.h"

using namespace std;


class genetic_map{
    protected:
        /*The following members are supposed to be constant after initialization*/
        int number_of_loci;
        int number_of_individual;
        int total_number_of_missing_obs;
        
        string population_name;
        string distance_function;
        string population_type;
        vector<vector<char> > raw_mapping_data;
        vector<string> marker_names;
        vector<string> individual_names;

        // The function used to convert cm to rp 
        DF* df_;

        double clustering_prob_cut_off;
        bool estimation_before_clustering;
        bool detect_bad_data;
        ObjFunc objective_function;

        // if a small set of markers (by default 3) is more than no_map_threshold away from the rest of the 
        // markers, it is not mapped.
        double no_map_dist; 
        int no_map_size;
        double missing_threshold;
        
        vector<vector<double> > pair_wise_distances;
        
        //the number of linkage groups in the map
        int number_of_connected_components;
        
        //the markers for each linkage group
        vector<vector<int> > connected_components;
        
        //level of reference: linkage group, bin, markers
        vector<vector<vector<int> > > linkage_group_bins;
        vector<vector<int> > orders;
        vector<vector<double> > distance_between_adjacent_pairs; // this is the distance is cM
        
        vector<double>  upperbounds;
        vector<double>  lowerbounds;
        vector<double>  approx_bounds;

        /*calculate the chernoff_bound*/
        double calculate_hoeffding_bound(double prob_cut_off);
    
        /*condense the markers into bins*/
        void condense_markers_into_bins();

        /*cluster the markers into linkage groups*/
        int cluster();
        
        // Added by Yonghui on Oct 20, 2007
        // When the input genotype data contains a lot of missing observations, the number of bins could be 
        // unnecessarily too many. As a result, the distance between adjacent bins is almost 0. 
        // When this happens, we really need to combine the adjacent bins into one single bin. 
        // The following function does what I have just described above. 
        // The function should be called by the end of the generate_map function
        // level of reference: linkage group, bin, markers
        vector<vector<vector<int> > > lg_bins_condensed;
        vector<vector<double> > dist_condensed; //this is the distance in cM
        void condense_bin();
        
     public:
        genetic_map();
        ~genetic_map();
                
        /*return 0 if the function finishes successfully*/
        int read_raw_mapping_data(string file_path);
        
        /*used for debugging purposes*/
        void dump();
        
        void dump_distance_matrix();
        
        // this function is supposed to be called after the order of each lg has been figured out
        void dump_connected_components_edges();    
        
        virtual void generate_map() = 0;
        
        void write_output(ostream& _output);
};

class genetic_map_DH: public genetic_map {
    private:
    	vector<pair<string, string> > suspicious_data;
        /*calculate the pair-wise distance*/
        void calculate_pair_wise_distance();
 
        /*generate a linkage group*/
        linkage_group_DH* construct_linkage_group(int group_id);
        linkage_group_DH* construct_linkage_group_whole_map();
        void print_suspicious_data();
    public:
        genetic_map_DH():genetic_map(){};
        ~genetic_map_DH();
        virtual void generate_map();
        void print_double_cross_overs();
};

#endif
