/*
 *  linkage_group_RIL.cpp
 *  ApproxMap
 *
 *  Created by yonghui on 12/13/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "linkage_group_RIL.h"
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void RIL_dist_cal::count_class(){
    CC_ = 0.0; DD_ = 0.0; EE_ = 0.0; FG_ = 0.0;
    double epsilon = 0.000001;
    for (int ii = 0; ii < num_of_individuals_; ii++) {
        double marker1_total = marker1_[ii].A + marker1_[ii].B + marker1_[ii].AB;
        double marker2_total = marker2_[ii].A + marker2_[ii].B + marker2_[ii].AB;
        assert(marker1_total < 1 + epsilon);
        assert(marker2_total < 1 + epsilon);
        assert(marker1_total > 1 - epsilon);
        assert(marker2_total > 1 - epsilon);
    }
    for (int ii = 0; ii < num_of_individuals_; ii++) {
        if ((not marker1_[ii].missing) and (not marker2_[ii].missing)) {
            num_of_eff_individuals_++;
            CC_ = CC_ + marker1_[ii].A * marker2_[ii].A + marker1_[ii].B * marker2_[ii].B;
            DD_ = DD_ + marker1_[ii].A * marker2_[ii].B + marker1_[ii].B * marker2_[ii].A;
            EE_ = EE_ + marker1_[ii].AB * marker2_[ii].A + marker1_[ii].AB * marker2_[ii].B 
                    + marker2_[ii].AB * marker1_[ii].A + marker2_[ii].AB * marker1_[ii].B;
            FG_ = FG_ + marker1_[ii].AB * marker2_[ii].AB;
        }
    }
    double total = CC_ + DD_ + EE_ + FG_;

    assert(total < num_of_eff_individuals_ + epsilon);
    assert(total > num_of_eff_individuals_ - epsilon);
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void RIL_dist_cal::upper_bound(){
    double R = (DD_ + EE_ + FG_) / (CC_ + DD_ + EE_ + FG_);
    delta_upper_bound_ = R / (2 - 2 * R);
    if (delta_upper_bound_ > 0.5) delta_upper_bound_ = 0.5;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void RIL_dist_cal::lower_bound(){
    double R = DD_ / (CC_ + DD_ + EE_ + FG_);
    delta_lower_bound_ = R / (2 - 2 * R);
    if (delta_lower_bound_ > 0.5) delta_lower_bound_ = 0.5;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void RIL_dist_cal::expected_CDEFG(double delta,
                                  double& C,
                                  double& D,
                                  double& E,
                                  double& F,
                                  double& G) const {
    assert(generation_index_ > 1);
    double CC = 0.0;
    double DD = 0.0;
    double EE = 0.0;
    double FF = 1.0;
    double GG = 0.0;
    for (int ii = 0; ii < generation_index_ - 1; ii++) {
        double CC_prime = CC + 0.5 * EE + 0.25 * (1 - delta) * (1 - delta) * FF + 0.25 * delta * delta * GG;
        double DD_prime = DD + 0.5 * EE + 0.25 * delta * delta * FF + 0.25 * (1 - delta) * (1 - delta) * GG; 
        double EE_prime = 0.5 * EE + 0.5 * delta * (1 - delta) * (FF + GG);
        double FF_prime = 0.5 * (1 - delta) * (1 - delta) * FF + 0.5 * delta * delta * GG;
        double GG_prime = 0.5 * delta * delta * FF + 0.5 * (1 - delta) * (1 - delta) * GG;
        CC = CC_prime; DD = DD_prime; EE = EE_prime; FF = FF_prime; GG = GG_prime;
    }
    double epsilon = 0.000001;
    double total = 2 * CC + 2 * DD + 4 * EE + FF + GG;
    assert(total < 1 + epsilon);
    assert(total > 1 - epsilon);
    C = CC; D = DD; E = EE; F = FF; G = GG;

};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double RIL_dist_cal::squared_error(double delta) const {
    double CC1, DD1, EE1, FF1, GG1, FG1;
    expected_CDEFG(delta, CC1, DD1, EE1, FF1, GG1);
    // to normalize CC1, DD1, EE1, FG1 such that they add up to 1
    CC1 = 2 * CC1;
    DD1 = 2 * DD1;
    EE1 = 4 * EE1;
    FG1 = FF1 + GG1;
    double epsilon = 0.000001;
    
    assert(CC1 + DD1 + EE1 + FG1 < 1 + epsilon);
    assert(CC1 + DD1 + EE1 + FG1 > 1 - epsilon);
    
    double CC2 = CC_ / num_of_eff_individuals_;
    double DD2 = DD_ / num_of_eff_individuals_;
    double EE2 = EE_ / num_of_eff_individuals_;
    double FG2 = FG_ / num_of_eff_individuals_;
    double loss = 0.0;
    loss = loss + (CC1 - CC2) * (CC1 - CC2)
                + (DD1 - DD2) * (DD1 - DD2)
                + (EE1 - EE2) * (EE1 - EE2)
                + (FG1 - FG2) * (FG1 - FG2);
    // cout << "CC1:" << CC1 << " DD1:" << DD1 << " EE1:" << EE1 << " FG1:" << FG1 << endl;
    // cout << "CC2:" << CC2 << " DD2:" << DD2 << " EE2:" << EE2 << " FG2:" << FG2 << endl;
    return loss;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double RIL_dist_cal::find_opt_delta() const {
    double bin_width = 0.1 / num_of_individuals_;
    double opt_delta = delta_lower_bound_;
    double crt_delta = delta_lower_bound_;
    double opt_loss = squared_error(delta_lower_bound_);
    while (crt_delta <= delta_upper_bound_){
        double crt_loss = squared_error(crt_delta);
        if (crt_loss < opt_loss) {
            opt_loss = crt_loss;
            opt_delta = crt_delta;
        }
        crt_delta = crt_delta + bin_width;
    }
    // cout << opt_delta << delta_lower_bound_ << delta_upper_bound_ << endl;
    return opt_delta;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

linkage_group_RIL::linkage_group_RIL(int _number_of_bins, 
                                     int _number_of_individuals, 
                                     int _generation_index,
                                     DF* _df,
                                     const vector<vector<allel_state> >& _raw_data, 
                                     const vector<int>& _current_order, 
                                     const vector<pair<int,int> >& _missing_data){
    number_of_bins = _number_of_bins;
    number_of_individuals = _number_of_individuals;
    raw_data = _raw_data;
    generation_index_ = _generation_index;
    current_order = _current_order;
    missing_data = _missing_data;
    df = _df;

    /*perform some consistency check*/
    if (raw_data.size() != number_of_bins) {
        cout << "BAD DATA" << endl;
        assert(false); // crash the problem on error
    }
    
    pair_wise_distances.resize(number_of_bins);
    for (int ii = 0; ii < number_of_bins; ii++) {
        (pair_wise_distances[ii]).resize(number_of_bins);
    }
    
    /*Calculate the pair-wise distance*/
    calculate_pair_wise_distance();
    
    current_upper_bound = 0 ; 
    for (int ii = 1 ; ii < number_of_bins; ii++) {
        current_upper_bound = current_upper_bound + pair_wise_distances[current_order[ii-1]][current_order[ii]];
    }
    cost_after_initialization = 0 ; 
    MST_lower_bound = 0;
    
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

linkage_group_RIL::~linkage_group_RIL(){
    return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void linkage_group_RIL::calculate_pair_wise_distance(){
    pair_wise_distances.resize(number_of_bins);
    for (int ii = 0 ; ii < number_of_bins; ii++)
    {
        pair_wise_distances[ii].resize(number_of_bins);
    }
    for (int ii = 0 ; ii < number_of_bins; ii++)
    {
        for (int jj = ii + 1 ; jj < number_of_bins; jj++)
        {
            RIL_dist_cal ril_dist(generation_index_, raw_data[ii], raw_data[jj]);
            double distance_ii_jj = ril_dist.Dist();
            pair_wise_distances[ii][jj] = distance_ii_jj;
            pair_wise_distances[jj][ii] = distance_ii_jj;
        }
    }
    for (int ii = 0; ii < number_of_bins; ii++) {
        pair_wise_distances[ii][ii] = 0.0;
    }    

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void linkage_group_RIL::estimate_missing_data() {
    // need to do something here
    return;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void linkage_group_RIL::dump() const {
    dump_common();
    cout << "generation_index:" << generation_index_ << endl;
    cout << "The raw data ordered" << endl;
    for (int ii = 0 ; ii < number_of_bins; ii++) {
        int jj = current_order[ii];
        for (int kk = 0 ; kk < number_of_individuals ; kk++) {
            if ((raw_data[jj][kk].A > raw_data[jj][kk].B) and
                (raw_data[jj][kk].A > raw_data[jj][kk].AB))  {
                cout << '.';
            } else if ((raw_data[jj][kk].B > raw_data[jj][kk].A) and
                       (raw_data[jj][kk].B > raw_data[jj][kk].AB))  {
                cout << '#';
            } else {
                cout << '+';
            }
        }
        cout << endl;
    };
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void linkage_group_RIL::order_markers() {
    vector<int> last_order = current_order;
    double last_upper_bound = current_upper_bound;
    bool to_continue = true;
    while (to_continue) {
        to_continue = false;
        vector<vector<double> >  distance_in_cM;
        generate_distance_in_cM(distance_in_cM);
//      MSTOpt opt_iter(pair_wise_distances, number_of_bins);
        MSTOpt opt_iter(distance_in_cM, number_of_bins, 1);
        opt_iter.Opt_Order(current_order, MST, MST_lower_bound, current_upper_bound, cost_after_initialization);
        
        bool order_changed = false;
        for (int ii = 0 ; ii < number_of_bins; ii ++) {
            if (last_order[ii] != current_order[ii]) {
                order_changed = true;
            }
        }
        if ((order_changed) and (missing_data.size() > 0) and (current_upper_bound <  last_upper_bound + ZERO_MINUS )) {
            to_continue = true;
            last_order = current_order;
            last_upper_bound = current_upper_bound;
            estimate_missing_data();
            cout << "start a new iteration" << endl;
        }
    }
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
