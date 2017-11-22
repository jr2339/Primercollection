//
//  utilities.hpp
//  PC
//
//  Created by jr2339 on 10/30/17.
//  Copyright © 2017 jr2339. All rights reserved.
//

#ifndef utilities_hpp
#define utilities_hpp

#include <stdio.h>
#include <stdlib.h>
#include "record.hpp"
#include <stdio.h>
#include <string>
#include <math.h>
#include <iostream>
#include <fstream>
#include<vector>

typedef struct FASTA{
    string fasta_headers;
    string fasta_sequences;
}FASTA;


typedef struct COEFFICIENT{
    float delta_h_co;
    float delta_s_co;
}COEFFICIENT;


typedef struct DELTA{
    float delta_h;
    float delta_s;
}DELTA;

int count_headers(string fasta_file_path);
int count_sequence(string fasta_file_path);
FASTA* read_fasta_file(string fasta_file_path);
COEFFICIENT* get_melting_temp_coeffients(string sequence);
DELTA* get_melting_temp_deltas(string sequence);
string compliment(string sequence);
float get_gc_content(string sequence);
float get_melting_temp(string sequence);

#endif /* utilities_hpp */
