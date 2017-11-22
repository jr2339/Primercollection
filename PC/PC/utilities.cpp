//
//  utilities.cpp
//  PC
//
//  Created by jr2339 on 10/30/17.
//  Copyright © 2017 jr2339. All rights reserved.
//

#include "utilities.hpp"

int count_headers(string fasta_file_path){
    int header_lines = 0;
    string line;
    ifstream input;
    input.open(fasta_file_path);
   
    if (!input) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }
    
    while (getline(input,line)) {
        if (line[0] == '>') {
            header_lines++;
        }
    }
    
    
    return header_lines;
}

int count_sequence(string fasta_file_path){
    int sequence_lines = 0;
    string line;
    ifstream input;
    input.open(fasta_file_path);
    
    while (getline(input, line)) {
        if (line[0] != '>') {
            sequence_lines++;
        }
    }
    
    return sequence_lines;
}


FASTA* read_fasta_file(string fasta_file_path){
    FASTA* fasta = new FASTA;
    
    string line;
   
    ifstream input;
    input.open(fasta_file_path);
    
    while (getline(input, line)){
        if (line[0] == '>') {
            size_t position_1 = line.find('>');
            size_t position_2 = line.find_last_of(' ');
            fasta->fasta_headers += line.substr(position_1+1,position_2);
        }
        else{
            size_t position = line.find('\n');
            fasta->fasta_sequences += line.substr(0,position);
        }
    }
    
    
    
    while (1) {
        if (fasta->fasta_sequences[0] == ' ') {
            fasta->fasta_sequences =  fasta->fasta_sequences.substr(1,fasta->fasta_sequences.length());
        }
        else{
            break;
        }
    }
    
    //cout<<fasta->fasta_headers<<endl;
    //cout<<fasta->fasta_sequences<<endl;
    return fasta;
}



COEFFICIENT* get_melting_temp_coeffients(string sequence){
    string name[] = {"AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"};
    int counts[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    
    float H_value[] = {-7.9,-8.4,-7.8,-7.2,-8.5,-8.0,-10.6,-7.8,-8.2,-9.8,-8.0,-8.4,-7.2,-8.2,-8.5,-7.9};
    float S_value[] = {-22.2,-22.4,-21.0,-20.4,-22.7,-19.9,-27.2,-21.0,-22.2,-24.4,-19.9,-22.4,-21.3,-22.2,-22.7,-22.2};
    
    COEFFICIENT* coefficient = new COEFFICIENT ;
    
    
    
    vector<string>pair(sequence.length()-1);
    
    
    
    for (int i=0; i<sequence.length()-1; i++) {
        char*  temp = (char*)malloc(2);
        
        temp[0] = sequence[i];
        temp[1] = sequence[i+1];
        pair[i].assign(temp, 2);
        
    }
    
    
    
    for (int i=0; i<sequence.length()-1; i++) {
        
        for (int j=0; j<16; j++) {
            if (pair[i] == name[j]) {
                counts[j] = counts[j] + 1;
                
            }
        }
        
    }
    
    
    
    for (int i=0; i<16; i++) {
        
        
        coefficient->delta_h_co += counts[i]*H_value[i];
        
        coefficient->delta_s_co += counts[i]*S_value[i];
    }

    return coefficient;
}

DELTA* get_melting_temp_deltas(string sequence){
    DELTA* delta = new DELTA;
    string sequence_reversed = sequence;
    reverse(sequence_reversed.begin(), sequence_reversed.end());
    if (sequence == sequence_reversed) {
        delta->delta_s = -1.4;
    }
    else{
        delta->delta_s = 0.0;
    }
    
    
    char start = sequence[0];
    char end = sequence[sequence.length()-1];
    if (start == 'C' or start == 'G') {
        delta->delta_h += 0.1;
        delta->delta_s -= 2.8;
    }
    else if (start == 'A' or start == 'T'){
        delta->delta_h += 2.3;
        delta->delta_s += 4.1;
    }
    if (end == 'C' or end == 'G') {
        delta->delta_h += 0.1;
        delta->delta_s -= 2.8;
    }
    else if (end == 'A' or end =='T'){
        delta->delta_h += 2.3;
        delta->delta_s += 4.1;
    }
    
    return delta;
}

string compliment(string sequence){
    char* compliment = (char*)malloc((sequence.length()-1));
    
    for (int i=0; i<sequence.length(); i++) {
        if (sequence[i]=='A') {
            compliment[i] = 'T';
        }
        else if (sequence[i]=='T'){
            compliment[i] = 'A';
        }
        else if (sequence[i]=='G'){
            compliment[i] = 'C';
        }
        else if (sequence[i]=='C'){
            compliment[i] = 'G';
        }
        else{
            compliment[i] = sequence[i];
        }
    }
    
    
    
    return compliment;
}

float get_gc_content(string sequence){
    int gc_count = 0;
    for (int i=0; i < sequence.length(); i++) {
        if (sequence[i] == 'G' or sequence[i] == 'C') {
            gc_count++;
        }
    }
    float gc_content = (float)gc_count/(sequence.length());
    
    return gc_content;
}

float get_melting_temp(string sequence){
    
    float melting_temp =0;
    float cation_ratio,denominator,a,d,g;
    float DNA_c = 5000.0 * 1e-9;
    float Na_c = 50.0 * 1e-3;
    float Mg_c =  0.0 * 1e-3;
    float dNTPs_c = 0.0 * 1e-3;
    float R = 1.987;
    
    DELTA* delta = get_melting_temp_deltas(sequence);
    
    
    COEFFICIENT* coefficient = get_melting_temp_coeffients(sequence);
    
    
    delta->delta_h += coefficient->delta_h_co;
    delta->delta_s += coefficient->delta_s_co;
    

    
    //melting temperature is calculated
    melting_temp = (1000 * delta->delta_h)/(delta->delta_s + (R*log(DNA_c)));
    //cout<<melting_temp<<endl;

    
    
    float gc_content = get_gc_content(sequence) ;
    
    float Ka = 3e4;
    
    
    float D = pow((Ka * dNTPs_c - Ka * Mg_c + 1),2) + (4 * Ka * Mg_c);
    
    float Fmg = (-(Ka * dNTPs_c - Ka * Mg_c + 1) + sqrt(D))/(2 * Ka);
    
    if (Na_c > 0) {
        cation_ratio = sqrt(Fmg)/Na_c;
    }
    else{
        cation_ratio = 7.0;
    }
    
    if (cation_ratio < 0.22) {
        denominator = ((1/melting_temp) + ((4.29 * gc_content - 3.95) * log(Na_c) +  0.94 * pow(log(Na_c),2)) * 1e-5);
        melting_temp = 1/denominator;
    }
    else{
        a = 3.92;
        d = 1.42;
        g = 8.31;
        Fmg = Mg_c;
        if (cation_ratio < 6.0) {
            a = a * (0.843 - 0.352 * sqrt(Na_c) * log(Na_c));
            d =d * (1.279 - 4.03 * log(Na_c) *
                    1e-3 - 8.03 * pow(log(Na_c),2) * 1e-3);
            g =g * (0.486 - 0.258 * log(Na_c) + 5.25 * pow(log(Na_c),3) * 1e-3);
        }
        melting_temp = 1.0 / (
                              (1.0 / melting_temp) +
                              (a - 0.911 * log(Fmg) + gc_content * (6.26 + d * log(Fmg)) +
                               1.0 / (2.0 * (sequence.length() - 1.0)) *
                               (-48.2 + 52.5 * log(Fmg) + g * pow(log(Fmg),2))) * 1e-5);
    }
    
    
    melting_temp -= 273.15;
    return melting_temp;
}













