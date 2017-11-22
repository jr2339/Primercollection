//
//  comparison.cpp
//  PC
//
//  Created by jr2339 on 10/30/17.
//  Copyright © 2017 jr2339. All rights reserved.
//

#include "comparison.hpp"


string reverse_compliment(string sequence){
    
    reverse(sequence.begin(), sequence.end());
    
    for (int i=0; i<sequence.length(); i++) {
        if (sequence[i] == 'A') {
            sequence[i] = 'T';
        }
        else if (sequence[i] == 'C'){
            sequence[i] = 'G';
        }
        else if (sequence[i] == 'G'){
            sequence[i] = 'C';
        }
        else if (sequence[i] == 'T'){
            sequence[i] = 'A';
        }
    }
    
    return sequence;
}

int get_primer_size(string p1, string p2){
    int size = 0;
    
    if (p1.length() > p2.length()) {
        size = (int)p1.length();
    }
    else{
        size = (int)p2.length();
    }
    return size;
}

int* generate_match_list(string p1, string p2){
    
    int size = get_primer_size(p1,p2);
    
    int* matches = (int*)malloc(size * sizeof(int));
    
    for (int i=0; i<(int)size; i++) {
        if (p1[i] == p2[i]) {
            matches[i] = 1;
        }
        else{
            matches[i] = 0;
        }
    }
    
    
    return matches;
    
}


int get_hairpin(string p1, string p2){
    int hairpin = 0;
    int* matches = generate_match_list(p1,p2);
    int size = get_primer_size(p1, p2);
    
    
    if (matches[size-1] == 0) {
        hairpin = 1;
    }
    else{
        hairpin = 0;
    }
    
    
    return hairpin;
}


TOOL* get_tools(string p1, string p2){
    TOOL* tool = new TOOL;
    tool->longest_run = 0;
    tool->hairpin_score = 0;
    tool->jackpot = 0;
    tool->most_matches = 0;
    
    int group = 1; //at least we have 1 group
    int sum_of_a_group = 0;
    int index = 0;
    int size = get_primer_size(p1, p2);
    int hairpin = get_hairpin(p1, p2);
    char store[2*size];
    tool->matches = generate_match_list(p1,p2);
    
    for (int i=0; i<(int)size; i++){
        tool->most_matches += tool->matches[i];
    }
    
    if ( hairpin == 1) {
        for (int i=0; i<(int)size; i++){
            tool->hairpin_score+= tool->matches[i];
        }
    }
    else{
        tool->hairpin_score = 0;
    }
    
    if (size == 1) {
        store[index] = (tool->matches[0]+'0');
        index++;
        store[index] = ' ';
    }
    else{
        for (int i= 0; i<size; i++) {
            if (i != size-1) {
                
                if (tool->matches[i] != tool->matches[i+1]) {
                    store[index] = (tool->matches[i]+'0');
                    
                    index++;
                    store[index]  = ' ';
                    
                    index++;
                }
                else{
                    store[index] = (tool->matches[i]+'0');
                    
                    index++;
                }
            }
            else{
                
                if (tool->matches[i-1] != tool->matches[i]) {
                    store[index] = ' ';
                    
                    index++;
                    store[index] = (tool->matches[i]+'0');
                    
                    index++;
                    store[index] = ' ';
                }
                else{
                    store[index] = (tool->matches[i]+'0');
                    
                    index++;
                    store[index] = ' ';
                    
                }
            }
        }
        
    }
    
    
    
    //check how many group do we have
    
    for(int i=0;i<index;i++){
        if (store[i] == ' ') {
            group++;
        }
        
    }
    
    if (group == 1) {
        
        for(int i=0;i<size;i++){
            
            tool->jackpot += tool->matches[i];
        }
    }
    
    vector<int> groups;
    groups.reserve(group);
    
    
    if (size == 1) {
        sum_of_a_group = tool->matches[0];
        groups.push_back(sum_of_a_group);
        
    }
    else{
        for(int i=0;i<index+1;i++){
            
            if (store[i] != ' ') {
                sum_of_a_group += (store[i]-'0');
            }
            else{
                groups.push_back(sum_of_a_group);
                sum_of_a_group = 0;
                
            }
            
        }
    }
    
    tool->longest_run = *max_element(groups.begin(),groups.end());
    
    
    return tool;
}



SCORE* get_score(string primer_1, string primer_2){
    
    SCORE* score = new SCORE;
    
    score->longest_run = 0;
    score->most_matches = 0;
    score->r_jackpot = 0;
    score->l_jackpot = 0;
    score->max_hairpin_score = 0;
    
    int l1 = (int)primer_1.length();
    int l2 = (int)primer_2.length();
    
    
    
    for (int i = 0; i < l1; i++) {
        string p1 = primer_1.substr(l1-i-1,l1);
        string p2 = primer_2.substr(0,i+1);
        
        
        
        TOOL* tool = get_tools(p1,p2);
        
        
        
        
        if (tool->longest_run > score->longest_run) {
            score->longest_run = tool->longest_run;
        }
        if (tool->most_matches > score->most_matches ) {
            score->most_matches = tool->most_matches;
        }
        if (tool->jackpot > score->r_jackpot ) {
            score->r_jackpot  = tool->jackpot;
        }
        
        if (tool->hairpin_score > score->max_hairpin_score) {
            score->max_hairpin_score = tool->hairpin_score;
        }
        
    }
    
    for (int i=0; i < l2; i++) {
        
        string p1 = primer_1.substr(0,i+1);
        string p2 = primer_2.substr(l2-i-1,l2);
        
        
        
        
        TOOL* tool = get_tools(p1,p2);
        
        
        
        
        
        if (tool->longest_run > score->longest_run) {
            score->longest_run = tool->longest_run;
            
        }
        if (tool->most_matches > score->most_matches ) {
            score->most_matches = tool->most_matches;
        }
        if (tool->jackpot > score->l_jackpot ) {
            score->l_jackpot  = tool->jackpot;
        }
        if (tool->hairpin_score > score->max_hairpin_score) {
            score->max_hairpin_score = tool->hairpin_score;
        }
        
    }
    
    
    return score;
    
}


SCORE* self_hybridization(string primer){
    SCORE* score = new SCORE;
    string primer_rc = reverse_compliment(primer);
    
    score = get_score(primer, primer_rc);
    
    return score;
}


SCORE* cross_hybridization(string p1, string p2){
    SCORE* score = new SCORE;
    string short_primer, long_primer;
    if (p1.length() < p2.length()) {
        short_primer = reverse_compliment(p1);
        long_primer = p2;
    }
    else if (p2.length() < p1.length()){
        short_primer = reverse_compliment(p2);
        long_primer = p1;
    }
    p1 = long_primer;
    p2 = short_primer;
    
    score = get_score(p1,p2);
    
    
    
    return score;
    
}

