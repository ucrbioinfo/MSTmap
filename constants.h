/*
 *  constants.h
 *  ApproxMap
 *
 *  Created by yonghui on 4/17/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef CONSTANTS_HEADER
#define CONSTANTS_HEADER

#include <string>

using namespace std;
const double PROB_HOEFFDING_CUT_OFF = 0.000001;
const double ZERO_MINUS = -0.0001;
const double ESTIMATION_BEFORE_CLUSTERING = 0.01;
const double ZERO_PLUS = 0.0001;
const double Missing_Threshold = 0.30; // a marker will be removed if more than 40% of its genotype calls are missing
const double COMBINE_BINS_TH = 0.1; // 0.000000000000001;
const string HALDANE = "haldane";
const string KOSAMBI = "kosambi";
const double kMaskThreshold = 0.75;
const double kMinMaskThreshold = 0.75;
const double kMaskDecrement = 0.02;
enum ObjFunc{OBJF_ML, OBJF_COUNT, OBJF_CM};
const int kMaxErrorDectionRounds = 20;
const int kMaxMissingEstRounds = 10;
const bool kMSTVerbose = false;
const int kBadDetMaxNum = 8;
#endif
