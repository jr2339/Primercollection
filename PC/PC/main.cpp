//
//  main.cpp
//  PC
//
//  Created by jr2339 on 10/30/17.
//  Copyright © 2017 jr2339. All rights reserved.
//

#include <iostream>
#include "utilities.hpp"
#include "comparison.hpp"
#include "collection.hpp"
int main(int argc, const char * argv[]) {
    if (argc !=5) {
        perror("we need pass four arguments");
        exit(1);//if the nunber is not 0, not access to error
    }
    const char* inputNmae = argv[1];
    int primer_min_length = atoi(argv[2]);
    int primer_max_length = atoi(argv[3]);
    const char* working_dirp = argv[4];
    //get_melting_temp("TATGAATCCAGGTGGCATTCG");
    //parse_input(inputNmae);
    //total_primer(inputNmae, primer_min_length, primer_max_length);
    //generate_chunk_ids(inputNmae, primer_min_length, primer_max_length);
    //get_genome(inputNmae,working_dirp);
    //get_region(inputNmae,primer_min_length,primer_max_length,working_dirp);
    //get_primer(inputNmae,primer_min_length,primer_max_length,working_dirp);
    //get_match(inputNmae, primer_min_length,primer_max_length,working_dirp);
    get_primer_char(inputNmae,primer_min_length,primer_max_length,working_dirp);
    return 0;
}
