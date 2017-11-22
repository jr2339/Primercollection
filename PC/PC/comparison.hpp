//
//  comparison.hpp
//  PC
//
//  Created by jr2339 on 10/30/17.
//  Copyright © 2017 jr2339. All rights reserved.
//

#ifndef comparison_hpp
#define comparison_hpp

#include <stdio.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>
#include <string>
#include <string>

using namespace std;

typedef struct TOOL{
    int longest_run;
    int most_matches;
    int jackpot;
    int hairpin_score;
    int* matches;
}TOOL;


typedef struct SCORE{
    int longest_run;
    int most_matches;
    int r_jackpot;
    int l_jackpot;
    int max_hairpin_score;
}SCORE;

string reverse_compliment(string sequence);
int get_primer_size(string p1, string p2);
int* generate_match_list(string p1, string p2);

int get_hairpin(string p1, string p2);

TOOL* get_tools(string p1, string p2);
SCORE* get_score(string primer_1, string primer_2);
SCORE* self_hybridization(string primer);
SCORE* cross_hybridization(string p1, string p2);



#endif /* comparison_hpp */
