/*
 *  single_mapping_population_raw_data.cpp
 *  ApproxMap
 *
 *  Created by yonghui on 4/7/07.linkage_group_DH
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "genetic_map_DH.h"

genetic_map::genetic_map() {
    clustering_prob_cut_off = PROB_HOEFFDING_CUT_OFF;
    number_of_loci = 0;
    number_of_individual = 0; 
    population_name = "";
    population_type = "";
    number_of_connected_components = 0;
    total_number_of_missing_obs = 0;
    objective_function = OBJF_COUNT;
    detect_bad_data = false;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

genetic_map::~genetic_map() {
    delete df_;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int genetic_map::read_raw_mapping_data(string file_path) {
    total_number_of_missing_obs = 0 ; 
    ifstream raw_mapping_data_file(file_path.c_str());
    string tmp_str;

    raw_mapping_data_file >> tmp_str;
    if (tmp_str != "population_type")
    {
        cout << "ERROR in the file format" << endl;
        assert(false); // crash the program on error
        return -1;
    }
    raw_mapping_data_file >> population_type;

    raw_mapping_data_file >> tmp_str;
    if (tmp_str != "population_name")
    {
        cout << "ERROR in the file format" << endl;
        assert(false); // crash the program on error
        return -1;
    }
    raw_mapping_data_file >> population_name;

    raw_mapping_data_file >> tmp_str;
    if (tmp_str != "distance_function")
    {
        cout << "ERROR in the file format" << endl;
        assert(false); // crash the program on error
        return -1;
    }
    raw_mapping_data_file >> distance_function;

    if (distance_function == HALDANE) {
        df_ = new DF_Haldane();
    } else if (distance_function == KOSAMBI) {
        df_ = new DF_Kosambi();
    } else {
        cout << "un-recognized distance function name" << endl;
        assert(false); // crash the program on error
    }

    raw_mapping_data_file >> tmp_str;
    if (tmp_str != "cut_off_p_value")
    {
        cout << "ERROR in the file format" << endl;
        assert(false); // crash the program on error
        return -1;
    }
    raw_mapping_data_file >> clustering_prob_cut_off;

    raw_mapping_data_file >> tmp_str;
    if (tmp_str != "no_map_dist")
    {
        cout << "ERROR in the file format" << endl;
        assert(false); // crash the program on error
        return -1;
    }
    raw_mapping_data_file >> no_map_dist;

    raw_mapping_data_file >> tmp_str;
    if (tmp_str != "no_map_size")
    {
        cout << "ERROR in the file format" << endl;
        assert(false); // crash the program on error
        return -1;
    }
    raw_mapping_data_file >> no_map_size;

    raw_mapping_data_file >> tmp_str;
    if (tmp_str != "missing_threshold")
    {
        cout << "ERROR in the file format" << endl;
        assert(false); // crash the program on error
        return -1;
    }
    raw_mapping_data_file >> missing_threshold;    

    raw_mapping_data_file >> tmp_str;
    if (tmp_str != "estimation_before_clustering")
    {
        cout << "ERROR in the file format" << endl;
        assert(false); // crash the program on error
        return -1;
    }
    raw_mapping_data_file >> tmp_str;
    if (tmp_str == "yes") {
        estimation_before_clustering = true;
    } else if (tmp_str == "no") {
        estimation_before_clustering = false;
    } else {
        cout << "unrecognized file format" << endl;
        assert(false); // crash the program on error
    }

    raw_mapping_data_file >> tmp_str;
    if (tmp_str != "detect_bad_data")
    {
        cout << "ERROR in the file format" << endl;
        assert(false); // crash the program on error
        return -1;
    }
    raw_mapping_data_file >> tmp_str;
    if (tmp_str == "yes") {
        detect_bad_data = true;
    } else if (tmp_str == "no") {
        detect_bad_data = false;
    } else {
        cout << "unrecognized file format" << endl;
        assert(false); // crash the program on error
    }
    
    // objective function
    raw_mapping_data_file >> tmp_str;
    if (tmp_str != "objective_function")
    {
        cout << "ERROR in the file format" << endl;
        assert(false); // crash the program on error
        return -1;
    }
    raw_mapping_data_file >> tmp_str;
    if (tmp_str == "ML") {
        objective_function = OBJF_ML;
    } else if (tmp_str == "COUNT") {
        objective_function = OBJF_COUNT;
    } else if (tmp_str == "CM") {
        objective_function = OBJF_CM;
    } else {
        cout << "unrecognized file format" << endl;
        assert(false); // crash the program on error
    }

    raw_mapping_data_file >> tmp_str;
    if (tmp_str != "number_of_loci")
    {
        cout << "ERROR in the file format" << endl;
        assert(false); // crash the program on error
        return -1;
    }
    raw_mapping_data_file >> number_of_loci;

    raw_mapping_data_file >> tmp_str;
    if (tmp_str != "number_of_individual")
    {
        cout << "ERROR in the file format" << endl;
        return -1;
    }
    raw_mapping_data_file >> number_of_individual;
    
    // added by yonghui on Mar 7th
    // read in the individual names
    individual_names.resize(number_of_individual);
    string name_ii;
    raw_mapping_data_file >> name_ii; // read off a dummy column
    for (int ii = 0; ii < number_of_individual; ii++) {
        raw_mapping_data_file >> name_ii;
        individual_names[ii] = name_ii;
    }

    int killed_markers = 0;
    for (int ii = 0 ; ii < number_of_loci; ii ++)
    {   
        int num_missing = 0; 
        vector<char> marker_data;
        marker_data.resize(number_of_individual);
        string marker_name_ii;
        raw_mapping_data_file >> marker_name_ii;
        for (int jj = 0 ; jj < number_of_individual ; jj++)
        {
            string SNP_jj;
            raw_mapping_data_file >> SNP_jj;
            if ((SNP_jj == "a") or (SNP_jj == "A")) { // homozygous A
                marker_data[jj] = 'A';
            } else if ((SNP_jj == "b") or (SNP_jj == "B")) { // homozygous B
                marker_data[jj] = 'B';
            } else if ((SNP_jj == "X") or (SNP_jj == "x") or (SNP_jj == "AB") or (SNP_jj == "ab")) { // heterozygous 
                marker_data[jj] = 'X';
            } else if ((SNP_jj == "U") or (SNP_jj == "-") or (SNP_jj == "_") or (SNP_jj == "u")) { // missing allele
                num_missing = num_missing + 1;
                marker_data[jj] = '-';
            } else {
                cout << "unrecognzed marker at line " << ii+1 << " marker:" << marker_name_ii << " column " << jj + 1 << endl;
                assert(false); // crash the program on error
                return -1;
            }
        }
        if (num_missing < missing_threshold * number_of_individual) {
            raw_mapping_data.push_back(marker_data);
            marker_names.push_back(marker_name_ii);
            total_number_of_missing_obs = total_number_of_missing_obs + num_missing;                            
        } else {
            killed_markers = killed_markers + 1;
            cout << "caution! marker:" << marker_name_ii << "was killed due to too many missing genotype calls" << endl;
        }
    }
    assert(number_of_loci == raw_mapping_data.size() + killed_markers);
    assert(raw_mapping_data.size() == marker_names.size());
    number_of_loci = raw_mapping_data.size();
    raw_mapping_data_file.close();
    return 0;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void genetic_map::condense_markers_into_bins()
{
    vector<bool> visitted(number_of_loci, false);
    for (int ii = 0 ; ii < number_of_connected_components; ii++) {
        vector<vector<int> > bins_ii;
        for (int jj = 0 ; jj < (connected_components[ii]).size(); jj ++) {
            if (not visitted[connected_components[ii][jj]]) {
                vector<int> bin_ii_jj_markers;
                int kk = connected_components[ii][jj];
                bin_ii_jj_markers.push_back(kk);
                for (vector<int>::iterator iter1 = (connected_components[ii]).begin(); 
                     iter1 != (connected_components[ii]).end(); 
                     iter1++) {
                    if ((pair_wise_distances[kk][*iter1] <= ZERO_PLUS ) and (*iter1 != kk) and (not visitted[*iter1])) {
                        bin_ii_jj_markers.push_back(*iter1);
                    }
                }
                for (vector<int>::iterator iter1 = bin_ii_jj_markers.begin(); 
                     iter1 != bin_ii_jj_markers.end(); 
                     iter1++) {
                    visitted[*iter1] = true;
                }
                bins_ii.push_back(bin_ii_jj_markers);
            }
        }
        linkage_group_bins.push_back(bins_ii);
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void genetic_map::dump_distance_matrix() {
    char buffer[10];
    cout << "matrix dimension:" << pair_wise_distances.size() << endl;
    for (int ii = 0; ii < pair_wise_distances.size(); ii++) {
        for (int jj = 0; jj < pair_wise_distances[ii].size(); jj++) {
            sprintf(buffer, "%.2f ", pair_wise_distances[ii][jj]);
            cout << buffer;
        }
        cout << endl;
    } 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void genetic_map::dump()
{
    cout << population_name << endl;
    cout << population_type << endl;
    cout << distance_function << endl;
    cout << number_of_loci << endl;
    cout << number_of_individual << endl;
    for (int ii = 0 ; ii < number_of_loci; ii ++) {
        for (int jj = 0 ; jj < number_of_individual; jj ++) {
            if (raw_mapping_data[ii][jj] == 'A') { 
                cout << '.';
            } else if (raw_mapping_data[ii][jj] == 'B') {
                cout << '#';
            } else if (raw_mapping_data[ii][jj] == 'X') { 
                cout << 'X'; // heterozygous
            } else {
                cout << '-'; // heterozygous
            }
            
        }
        cout << endl;
    }
    cout << "the number of connected components " << number_of_connected_components << endl;
    for (int ii = 0 ; ii < number_of_connected_components; ii++)
    {
        cout << (connected_components[ii]).size() << ',';
    }
    cout << endl;
    
    /*perform some consistency check*/
    
    /*1. Make sure the total number of markers in the linkage groups add up to number_of_loci*/
    int tmp_total = 0;
    for (int ii = 0 ; ii < number_of_connected_components; ii++)
    {
        tmp_total = tmp_total + (connected_components[ii]).size();
        int tmp_total_ii = 0 ;
        for (int jj = 0 ; jj < (linkage_group_bins[ii]).size(); jj++) {
            tmp_total_ii = tmp_total_ii + (linkage_group_bins[ii][jj]).size();
        }
        if (tmp_total_ii != (connected_components[ii]).size()) {
            cout << "ERROR, tmp_total_ii NOT equal to connected_components[ii]" << endl;
        }

    }
    if (tmp_total != number_of_loci)
    {
        cout << "ERROR, tmp_total NOT equal to number_of_loci" << endl;
    }
    
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void genetic_map::dump_connected_components_edges() {
    cout << "dump edges" << endl;
    double threshold = calculate_hoeffding_bound(clustering_prob_cut_off);
    cout << "calculate_hoeffding_bound:" << threshold << endl;
    for (int ii = 0; ii < number_of_connected_components; ii++) {
        cout << "==============================================" << endl;
        cout << "\t";
        vector<int> markers;
        for (int jj = 0; jj < linkage_group_bins[ii].size(); jj++) {
            markers.insert(markers.end(), 
                           linkage_group_bins[ii][orders[ii][jj]].begin(), 
                           linkage_group_bins[ii][orders[ii][jj]].end());
        }
        assert(markers.size() == connected_components[ii].size()); 
        for (int jj = 0; jj < markers.size(); jj++) {
            cout << marker_names[markers[jj]] << "\t";
        }
        cout << endl;
        for (int jj = 0; jj < markers.size(); jj++) {
            cout << marker_names[markers[jj]] << "\t";
            for (int kk = 0; kk < markers.size(); kk++) {
                if (pair_wise_distances[markers[jj]][markers[kk]] < threshold) {
                    // cout << pair_wise_distances[markers[jj]][markers[kk]];
                    cout << "#";
                } else {
                    cout << ".";
                }
                cout << "\t";
            }
            cout << endl;
        }
    }

    for (int ii = 0; ii < number_of_connected_components; ii++) {
        cout << "==============================================" << endl;
        cout << "\t";
        vector<int> markers;
        for (int jj = 0; jj < linkage_group_bins[ii].size(); jj++) {
            markers.insert(markers.end(), 
                           linkage_group_bins[ii][orders[ii][jj]].begin(), 
                           linkage_group_bins[ii][orders[ii][jj]].end());
        }
        assert(markers.size() == connected_components[ii].size()); 
        for (int jj = 0; jj < markers.size(); jj++) {
            cout << marker_names[markers[jj]] << "\t";
        }
        cout << endl;
        for (int jj = 0; jj < markers.size(); jj++) {
            cout << marker_names[markers[jj]] << "\t";
            for (int kk = 0; kk < markers.size(); kk++) {
                if (pair_wise_distances[markers[jj]][markers[kk]] < threshold) {
                    cout << pair_wise_distances[markers[jj]][markers[kk]];
                } else {
                    cout << ".";
                }
                cout << "\t";
            }
            cout << endl;
        }
    }

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double genetic_map::calculate_hoeffding_bound(double prob_cut_off) {
    double t;
    if (prob_cut_off >= 1)
    {
        t = numeric_limits<double>::max();
        return t;
    }
    else
    {
        t = sqrt((log(prob_cut_off)) / (-2*number_of_individual)); //according to the hoeffding bound
    }
    return number_of_individual*(0.5-t);
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int genetic_map::cluster() {
    double number_of_recombinations_cut_off  = calculate_hoeffding_bound(clustering_prob_cut_off);
    cout << "number_of_recombinations_cut_off: " << number_of_recombinations_cut_off << endl;
    
    double no_map_threshold = number_of_individual * df_->RP(no_map_dist);
    set<int> un_mapped_markers;
    
    if (no_map_threshold < number_of_recombinations_cut_off) {
        // splice our those suspicious markers
        vector<bool> visitted; 
        visitted.resize(number_of_loci);
        for (int ii = 0; ii < number_of_loci ; ii++) {
            visitted[ii] = false;
        }
        
        for (int ii = 0 ; ii < number_of_loci ; ii ++) {
            if (visitted[ii] == false) {//start a new connected component
                vector<int> crt_cc;
                queue<int> fifo_queue;
                fifo_queue.push(ii);
                while (not fifo_queue.empty()) {
                    int head = fifo_queue.front();
                    fifo_queue.pop();
                    if (visitted[head] == false) {
                        visitted[head] = true;
                        crt_cc.push_back(head);
                        for (int jj = 0 ; jj < number_of_loci; jj ++) {
                            if (pair_wise_distances[head][jj] < no_map_threshold) {
                                fifo_queue.push(jj);
                            }
                        }
                    }
                }//end while
                if (crt_cc.size() <= no_map_size) {
                    connected_components.push_back(crt_cc);
                    un_mapped_markers.insert(crt_cc.begin(), crt_cc.end());
                }
            }
        }
    }
    
    vector<bool> visitted; 
    visitted.resize(number_of_loci);
    for (int ii = 0; ii < number_of_loci ; ii++)
    {
        visitted[ii] = false;
    }
    for (set<int>::iterator iter1 = un_mapped_markers.begin(); iter1 != un_mapped_markers.end(); ++iter1) {
        visitted[*iter1] = true;
    }
    
    for (int ii = 0 ; ii < number_of_loci ; ii ++)
    {
        if (visitted[ii] == false) //start a new connected component
        {
            queue<int> fifo_queue;
            connected_components.push_back(vector<int>());
            int last_cc_id = connected_components.size() - 1;
            fifo_queue.push(ii);
            while (not fifo_queue.empty()) {
                int head = fifo_queue.front();
                fifo_queue.pop();
                if (visitted[head] == false) {
                    visitted[head] = true;
                    connected_components[last_cc_id].push_back(head);
                    for (int jj = 0 ; jj < number_of_loci; jj ++) {
                        if ((pair_wise_distances[head][jj] < number_of_recombinations_cut_off) and
                            (visitted[jj] == false)) { 
                            fifo_queue.push(jj);
                        }
                    }
                }
            }//end while
        }
    }
    number_of_connected_components = connected_components.size();
    return connected_components.size();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void genetic_map::condense_bin(){
    lg_bins_condensed.resize(linkage_group_bins.size());
    dist_condensed.resize(linkage_group_bins.size());
    // for each linkage group condense two bins if they are too close to each other
    for (int ii = 0; ii < linkage_group_bins.size(); ii++) {
        lg_bins_condensed[ii].push_back(linkage_group_bins[ii][orders[ii][0]]);
        int crt_bin_id = 0;
        for (int jj = 1; jj < orders[ii].size(); jj++) {
            // determine whether or not to create a new bin
            if (distance_between_adjacent_pairs[ii][jj-1] < COMBINE_BINS_TH) { 
                // condense the next bin into the current bin
                lg_bins_condensed[ii][crt_bin_id].insert(lg_bins_condensed[ii][crt_bin_id].end(),
                                                         linkage_group_bins[ii][orders[ii][jj]].begin(),
                                                         linkage_group_bins[ii][orders[ii][jj]].end());
            } else {
                crt_bin_id = crt_bin_id + 1;
                lg_bins_condensed[ii].push_back(linkage_group_bins[ii][orders[ii][jj]]);
                dist_condensed[ii].push_back(distance_between_adjacent_pairs[ii][jj-1]);
            }
        }
    }
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void genetic_map::write_output(ostream& _output)
{
    _output << ";number of linkage groups: " << number_of_connected_components << endl;
    _output << ";The size of the linkage groups are: " << endl << ";\t\t";
    for (int ii = 0 ; ii < number_of_connected_components ; ii++)
    {
        _output << (connected_components[ii]).size() << '\t';
    }
    _output << endl;
    
    _output << ";The number of bins in each linkage group: " << endl << ";\t\t";
    for (int ii = 0 ; ii < number_of_connected_components; ii++)
    {
        _output << (lg_bins_condensed[ii]).size() << '\t';
    }
    _output << endl << endl << endl;
    for (int ii = 0 ; ii < number_of_connected_components; ii++)
    {
        char buffer1[100];
        char buffer2[100];
        char buffer3[100];
        sprintf(buffer1,"%.3f",lowerbounds[ii]);
        sprintf(buffer2,"%.3f",upperbounds[ii]);
        sprintf(buffer3,"%.3f",approx_bounds[ii]);
        _output << ";lowerbound:" << buffer1 << " upperbound: " << buffer2;
        _output << " cost after initialization:" << buffer3 << endl;
        _output << "group lg" << ii << endl;
        _output << ";BEGINOFGROUP" << endl;
        assert(lg_bins_condensed[ii].size() > 0);
        for (vector<int>::iterator iter2 = (lg_bins_condensed[ii][0]).begin(); 
            iter2 != (lg_bins_condensed[ii][0]).end(); 
            iter2++) {
            _output << marker_names[*iter2] << '\t' << "0.000" << endl;
        }
        double cum_dist = 0.0;
        assert(lg_bins_condensed[ii].size() == dist_condensed[ii].size() + 1);
        for (int jj = 1; jj < lg_bins_condensed[ii].size(); jj++) {
            cum_dist = cum_dist + dist_condensed[ii][jj-1];
            for (vector<int>::iterator iter2 = (lg_bins_condensed[ii][jj]).begin(); 
                 iter2 != (lg_bins_condensed[ii][jj]).end(); 
                 iter2++) {
                char buffer[100];
                sprintf(buffer, "%.3f", cum_dist); 
                _output << marker_names[*iter2] << '\t' << buffer << endl;
            }
        }
        _output << ";ENDOFGROUP" << endl << endl << endl;
    }
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

genetic_map_DH::~genetic_map_DH(){
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void genetic_map_DH::calculate_pair_wise_distance()
{
    pair_wise_distances.resize(number_of_loci);
    for (int ii = 0 ; ii < number_of_loci; ii++)
    {
        pair_wise_distances[ii].resize(number_of_loci, 0.0);
    }

    cout << "start calculating pair-wise distance" << time(NULL) << endl;
    for (int ii = 0; ii < number_of_loci; ii++)
    {
        for (int jj = ii ; jj < number_of_loci; jj++)
        {
            double distance_ii_jj = 0;
            double none_missing = 0;
            for (int kk = 0 ; kk < number_of_individual; kk++)
            {
                
                if ((raw_mapping_data[ii][kk] != '-') and 
                    (raw_mapping_data[jj][kk] != '-')) {
                    none_missing = none_missing + 1.0;
                    if (raw_mapping_data[ii][kk] != raw_mapping_data[jj][kk]) {
                        distance_ii_jj = distance_ii_jj + 1.0;
                    }
                }
            }
            if (none_missing < 0.5 * number_of_individual) { 
                cout << "caution: too many missing for pair:(" 
                     << marker_names[ii] << "," << marker_names[jj] << ")" << endl;
            }
            
            if (none_missing < 0.25 * number_of_individual) { // almost everything is missing, adjust the estimate 
                distance_ii_jj = 0.5 * number_of_individual;
                none_missing = number_of_individual;
            }
            pair_wise_distances[ii][jj] = (distance_ii_jj / none_missing) * number_of_individual;
            pair_wise_distances[jj][ii] = pair_wise_distances[ii][jj];
        }
    }    
    cout << "finished calculating pair-wise distance:" << time(NULL) << endl;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

linkage_group_DH* genetic_map_DH::construct_linkage_group(int group_id)
{
    int _number_of_bins = (linkage_group_bins[group_id]).size();
    int _number_of_individuals = number_of_individual;
    
    /*Store the probability for each allele to be A*/
    vector<vector< float > > _raw_data ;
    
    vector<pair<int, int> >  _missing_data;
    
    vector<int> _current_order;
    
    _raw_data.resize(_number_of_bins);
    for (int ii = 0 ; ii < _number_of_bins; ii++)
    {
        _raw_data[ii].resize(_number_of_individuals);
        for (int jj = 0; jj < _number_of_individuals; jj++)
        {
            if (raw_mapping_data[linkage_group_bins[group_id][ii][0]][jj] == 'A') 
            {
                /*If an allele is A, then its probability being A is A*/
                _raw_data[ii][jj] = 1.0;
            }
            else if (raw_mapping_data[linkage_group_bins[group_id][ii][0]][jj] == 'B') 
            {
                /*If an allele is B, then its probability of being A is 0*/
                _raw_data[ii][jj] = 0.0;
            }
            else 
            {
                /*If an allele is missing, then, assign probability 0.5 for it to be A*/
                _raw_data[ii][jj] = 0.5;
                _missing_data.push_back(make_pair(ii,jj)); /*ii is the id for the marker, and jj is the id for the individual*/
            }
        }
    }
    for (int ii = 0 ; ii < _number_of_bins; ii ++)
    {
        _current_order.push_back(ii);
    }
    vector<int> bin_sizes;
    for (int ii = 0; ii < _number_of_bins; ii++) {
        bin_sizes.push_back(linkage_group_bins[group_id][ii].size());
    }
    linkage_group_DH * to_be_returned = new linkage_group_DH(_number_of_bins, 
                                                             _number_of_individuals,
                                                             detect_bad_data, 
                                                             objective_function,
                                                             df_,
                                                             _raw_data, 
                                                             _current_order, 
                                                             _missing_data,
                                                             bin_sizes);
    return to_be_returned;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

linkage_group_DH* genetic_map_DH::construct_linkage_group_whole_map()
{
    int _number_of_bins = number_of_loci;
    int _number_of_individuals = number_of_individual;
    
    /*Store the probability for each allele to be A*/
    vector<vector< float > > _raw_data ;
    
    vector<pair<int, int> >  _missing_data;
    
    vector<int> _current_order;
    
    _raw_data.resize(_number_of_bins);
    for (int ii = 0 ; ii < _number_of_bins; ii++)
    {
        _raw_data[ii].resize(_number_of_individuals);
        for (int jj = 0; jj < _number_of_individuals; jj++)
        {
            if (raw_mapping_data[ii][jj] == 'A') 
            {
                /*If an allele is A, then its probability being A is A*/
                _raw_data[ii][jj] = 1.0;
            }
            else if (raw_mapping_data[ii][jj] == 'B') 
            {
                /*If an allele is B, then its probability of being A is 0*/
                _raw_data[ii][jj] = 0.0;
            }
            else 
            {
                /*If an allele is missing, then, assign probability 0.5 for it to be A*/
                _raw_data[ii][jj] = 0.5;
                _missing_data.push_back(make_pair(ii,jj)); /*ii is the id for the marker, and jj is the id for the individual*/
            }
        }
    }
    for (int ii = 0 ; ii < _number_of_bins; ii ++)
    {
        _current_order.push_back(ii);
    }
    vector<int> bin_sizes;
    for (int ii = 0; ii < _number_of_bins; ii++) {
        bin_sizes.push_back(1);
    }
    linkage_group_DH* to_be_returned = new linkage_group_DH(_number_of_bins, 
                                                             _number_of_individuals,
                                                             false, // this is fixed to be false for whole map 
                                                             OBJF_COUNT,
                                                             df_,
                                                             _raw_data, 
                                                             _current_order, 
                                                             _missing_data,
                                                             bin_sizes);
    return to_be_returned;
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void genetic_map_DH::print_suspicious_data(){
    cout << endl;
    for (int ii = 0; ii < suspicious_data.size(); ii++) {
        cout << suspicious_data[ii].first;
        cout << '\t';
        cout << suspicious_data[ii].second;
        cout << endl;
    }
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void genetic_map_DH::generate_map()
{
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
        linkage_group_DH * linkage_group_whole_map = construct_linkage_group_whole_map();
        linkage_group_whole_map->order_markers();
        const vector<vector<double> > & new_dist = linkage_group_whole_map-> get_pair_wise_distance();
        for (int ii = 0 ; ii < number_of_loci; ii++)
        {
            for (int jj = 0 ; jj < number_of_loci; jj++)
            {
                pair_wise_distances[ii][jj] = new_dist[ii][jj];
            }
        }
        // linkage_group_whole_map->dump_distance_matrix();
        delete linkage_group_whole_map;
    } else {
        cout << "calculating the pair-wise hamming distance" << endl;
        calculate_pair_wise_distance();
        cout << "finished calculating the pair-wise hamming distance" << endl;
    }
    // dump_distance_matrix(); // STELO
    cluster();
    // dump_distance_matrix(); // STELO
    cout << "found " << number_of_connected_components << " connected components" << endl;
    
    condense_markers_into_bins(); 
    // dump_distance_matrix(); // STELO
    
    orders.resize(number_of_connected_components);
    upperbounds.resize(number_of_connected_components);
    lowerbounds.resize(number_of_connected_components);
    approx_bounds.resize(number_of_connected_components);
    distance_between_adjacent_pairs.resize(number_of_connected_components);
    
    for (int ii = 0 ; ii < number_of_connected_components; ii++)
    {
        linkage_group_DH * current_linkage_group = construct_linkage_group(ii);

        current_linkage_group->order_markers();
        current_linkage_group->return_order(orders[ii], 
                                            lowerbounds[ii], 
                                            upperbounds[ii], 
                                            approx_bounds[ii], 
                                            distance_between_adjacent_pairs[ii]);
        vector<pair<int, int> > bad_data_ii;                                    
        current_linkage_group->bad_genotypes(bad_data_ii);
        for (int jj = 0; jj < bad_data_ii.size(); jj++) {
            int bin_id = bad_data_ii[jj].first;
            int indi_id = bad_data_ii[jj].second;
            for (int kk = 0; kk < linkage_group_bins[ii][bin_id].size(); kk++) {
                string marker_name = marker_names[linkage_group_bins[ii][bin_id][kk]];
                string indi_name = individual_names[indi_id];
                suspicious_data.push_back(make_pair(marker_name, indi_name));
            }
        }
        current_linkage_group->dump();
        delete current_linkage_group;
        cout <<"finished the " << ii+1 << " linkage group" << endl; 
    }
    
    // Added by Yonghui on Oct 20, 2007 
    // The last step is to condense adjacent bins if they are too close to each other
    condense_bin(); 
    cout << "suspicious data detected by our algorithm" << endl;
    print_suspicious_data();
    cout << "double cross overs based on the current order" << endl;
    print_double_cross_overs();
    // dump the distance matrix
    // dump_connected_components_edges(); // STELO
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void genetic_map_DH::print_double_cross_overs(){
    for (int ii = 0; ii < lg_bins_condensed.size(); ii++) {
        if (lg_bins_condensed[ii].size() < 3) {
            continue;
        }
        for (int jj = 0; jj < lg_bins_condensed[ii].size(); jj++){
            if (lg_bins_condensed[ii][jj].size() > 1) {
                continue;
            }
            int marker_id = lg_bins_condensed[ii][jj][0];           
            int pre_marker_id = -1;
            if (jj == 0){ // it is the first bin
                pre_marker_id = lg_bins_condensed[ii][1][0];
            } else {
                pre_marker_id = lg_bins_condensed[ii][jj - 1][0];
            }
            
            int next_marker_id = -1;
            if (jj == lg_bins_condensed[ii].size() - 1) { // it is the last bin
                next_marker_id = lg_bins_condensed[ii][lg_bins_condensed[ii].size() - 2][0];
            } else {
                next_marker_id = lg_bins_condensed[ii][jj + 1][0];
            }
            
            for (int kk = 0; kk < number_of_individual; kk++) {
                if (raw_mapping_data[marker_id][kk] == '-') { // ignore missing 
                    continue;
                }
                if ((raw_mapping_data[marker_id][kk] != raw_mapping_data[pre_marker_id][kk]) and 
                    (raw_mapping_data[marker_id][kk] != raw_mapping_data[next_marker_id][kk])) {
                    // this is a double cross-over
                    cout << marker_names[marker_id] << "," << individual_names[kk] << endl;
                }
            }
        }
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

