//
//  record.hpp
//  PC
//
//  Created by jr2339 on 10/30/17.
//  Copyright © 2017 jr2339. All rights reserved.
//

#ifndef record_h
#define record_h

#include <iostream>
#include <fstream>
#include <stdint.h>

#include <cstdlib>

#include<fstream>
#include<cstring>

using namespace std;


typedef struct GenomeRecord {
    string identifier;
    string fasta_fp;
}GenomeRecord;


typedef struct RegionRecord {
    string identifier;
    int start;
    string sequence;
    int poi_record;
    string headers;
    string chunk_ids;
}RegionRecord;


typedef struct PrimerStatsRecord{
    string identifier;
    int start;
    string sequence;
    string amp_direction;
    string orientation;
    string region_id;
}PrimerStatsRecord;


typedef struct Homopolymers{
    int base;
    string longest;
    float percent;
    
}Homopolymers;

typedef struct Dimerzation{
    int base;
    int longest;
    float percent;
    int hairpin;
    string jackpots;
}Dimerzation;

typedef struct Specificity{
    int occurrences;
    int* similar;
}Specificity;


typedef struct PrimerCharRecord{
    string identifier;
    float melting_temp;
    int gc_content;
    int length;
    Homopolymers* Hom;
    Dimerzation* Dim;
    Specificity* Spec;
}PrimerCharRecord;

typedef struct WarningRecord{
    float melting_temp;
    int gc_content;
    int length;
    Homopolymers* Hom;
    Dimerzation* Dim;
    Specificity* Spec;
}WarningRecord;






#endif /* record_h */
