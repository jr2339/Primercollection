//
//  collection.hpp
//  PC
//
//  Created by jr2339 on 10/30/17.
//  Copyright © 2017 jr2339. All rights reserved.
//

#ifndef collection_hpp
#define collection_hpp

#include <stdio.h>
#include <stdlib.h>
#include "utilities.hpp"
#include "comparison.hpp"
#include "record.hpp"
#include <cstring>
#include <string>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/types.h>
#include <sys/wait.h>

typedef struct ITEM{
    string  identifier;
    int start;
    int end;
    string reference_fp;
    string fp;
    int poi_record;
    string sequence;
}ITEM;

typedef struct MATCH{
    int* perfect_match;
    int** near_match;
}MATCH;

typedef struct STAFF{
    GenomeRecord** genome_stats;
    RegionRecord** region_stats;
    PrimerStatsRecord** primer_stats;
    PrimerCharRecord** primer_char_stats;
}STAFF;



int count_lines(string infile_fp);
int count_item(string infile_fp);
vector<string> get_lines(string infile_fp);
ITEM** parse_input(string infile_fp);
int total_primer(string infile_fp, int primer_length_min, int primer_length_max);
vector<string> generate_chunk_ids(string infile_fp, int primer_length_min, int primer_length_max);
GenomeRecord** get_genome(string infile_fp, string working_dirp);
RegionRecord** get_region(string infile_fp, int primer_length_min, int primer_length_max,string working_dirp);
PrimerStatsRecord** get_primer(string infile_fp, int primer_length_min, int primer_length_max,string working_dirp);
Homopolymers* get_homopolymers(string sequence);
Dimerzation* get_dimerzation(string sequence);
string run_blast(string db_location,string query_fp,string working_dirp,int word_size);
string make_blast_db(string subject_fp, string working_dirp);
string generate_query_fp(string infile_fp, int primer_length_min, int primer_length_max,string working_dirp);
string generate_subject_fp(string infile_fp, int primer_length_min, int primer_length_max,string working_dirp);
MATCH* get_match(string infile_fp, int primer_length_min, int primer_length_max,string working_dirp);
PrimerCharRecord** get_primer_char(string infile_fp, int primer_length_min, int primer_length_max,string working_dirp);












#endif /* collection_hpp */
