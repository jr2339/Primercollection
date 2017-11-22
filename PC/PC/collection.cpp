//
//  collection.cpp
//  PC
//
//  Created by jr2339 on 10/30/17.
//  Copyright © 2017 jr2339. All rights reserved.
//

#include "collection.hpp"
int count_lines(string infile_fp){
    int lines = 0;
    string line;
    ifstream input;
    input.open(infile_fp);
    
    if (!input) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }
    
    while (getline(input,line)) {
            lines++;
    }
    return lines;
}


int count_item(string infile_fp){
    int count =1;
    string line;
    ifstream input;
    input.open(infile_fp);
    getline(input, line);
    
    
    for (int i=0; i<line.length(); i++) {
        if (line[i]== '\t') {
            count++;
        }
    }
    return count;
}


vector<string> get_lines(string infile_fp){
    int count =  count_lines(infile_fp);
    vector<string> lines;
    lines.reserve(count);
    string line;
    ifstream input;
    input.open(infile_fp);
    while (getline(input, line)) {
        lines.push_back(line);
    }
    return lines;
}

ITEM** parse_input(string infile_fp){
    int count =  count_lines(infile_fp);
    vector<string> lines = get_lines(infile_fp);
    ITEM** item = (ITEM**)malloc(sizeof(ITEM*)*count);
    for (int i=0; i<count; i++) {
        item[i] = new ITEM;
        
        item[i]->identifier = lines[i].substr(0,lines[i].find_first_of("\t"));
        lines[i] = lines[i].substr(lines[i].find_first_of('\t')+1,lines[i].size());
        
        item[i]->start = stoi(lines[i].substr(0,lines[i].find_first_of("\t")));
        lines[i] = lines[i].substr(lines[i].find_first_of('\t')+1,lines[i].size());
        
        item[i]->reference_fp = lines[i].substr(0,lines[i].find_first_of('\t'));
        item[i]->fp = item[i]->reference_fp.substr(lines[i].find_first_of('/')+1);
        for (int j=0; j<4; j++) {
            item[i]->fp = item[i]->fp.substr(item[i]->fp.find_first_of('/')+1);
            item[i]->fp = item[i]->fp.substr(item[i]->fp.find_first_of('\t')+1,item[i]->fp.size());
        }
        item[i]->fp = "./" + item[i]->fp;
        
        lines[i] = lines[i].substr(lines[i].find_first_of('\t')+1,lines[i].size());
        
        item[i]->end  = stoi(lines[i].substr(0,lines[i].find_first_of("\t")));

    
        item[i]->poi_record =  item[i]->end - item[i]->start;
 
        lines[i] = lines[i].substr(lines[i].find_first_of('\t')+1,lines[i].size());
        
        item[i]->sequence = lines[i].substr(0,lines[i].find_first_of("\t"));
    }
    

    
    return item;
    
}


int total_primer(string infile_fp, int primer_length_min, int primer_length_max){
    int total = 0;
    int count =  count_lines(infile_fp);
    ITEM** item = parse_input(infile_fp);
    int middle_point[count];
    for (int i=0; i<count; i++) {
        middle_point[i] = (2*item[i]->poi_record - primer_length_min)/2;

    }
    
    for (int i=0; i<count; i++) {
        for (int index = 0; index<2*item[i]->poi_record; index++) {
            for (int span =0; span < primer_length_max - primer_length_min; span++) {
                if (index <= middle_point[i]) {
                    if (index + primer_length_min + span <= item[i]->poi_record-primer_length_min/2){

                        total++;
                    }
                    }
                if (index >= middle_point[i]+primer_length_min+1) {
                    if (index + primer_length_min + span <= 2*item[i]->poi_record) {
 
                        total++;
                    }
                }
                
                }
            }

    }
    
    return total;
}

vector<string> generate_chunk_ids(string infile_fp, int primer_length_min, int primer_length_max){
    int count =  count_lines(infile_fp);
    vector<string> chunk_ids(count);
    ITEM** item = parse_input(infile_fp);
    int middle_point[count];
    for (int i=0; i<count; i++) {
        middle_point[i] = (2*item[i]->poi_record - primer_length_min)/2;
        
    }
    
    for (int i=0; i<count; i++) {
        for (int index = 0; index<2*item[i]->poi_record; index++) {
            for (int span =0; span < primer_length_max - primer_length_min; span++) {
                if (index <= middle_point[i]) {
                    if (index + primer_length_min + span <= item[i]->poi_record-primer_length_min/2){
                        chunk_ids[i].append(item[i]->identifier);
                        chunk_ids[i].append("_");
                        chunk_ids[i].append(to_string(index));
                        chunk_ids[i].append("_");
                        chunk_ids[i].append(to_string(primer_length_min+span));
                        chunk_ids[i].append(",");
                        chunk_ids[i].append(" ");
                       
                    }
                }
                if (index >= middle_point[i]+primer_length_min+1) {
                    if (index + primer_length_min + span <= 2*item[i]->poi_record) {
                        chunk_ids[i].append(item[i]->identifier);
                        chunk_ids[i].append("_");
                        chunk_ids[i].append(to_string(index));
                        chunk_ids[i].append("_");
                        chunk_ids[i].append(to_string(primer_length_min+span));
                        chunk_ids[i].append(",");
                        chunk_ids[i].append(" ");
                       
                    }
                }
                
            }
        }
        
    }

    return chunk_ids;
}


GenomeRecord** get_genome(string infile_fp,string working_dirp){
    string genome_stats_fp = working_dirp + "genome_stats.txt";
    ofstream outputFile;
    outputFile.open(genome_stats_fp);
    int count =  count_lines(infile_fp);
    GenomeRecord** genome_stats = (GenomeRecord**)malloc(sizeof(GenomeRecord*)*count);
    ITEM** item = parse_input(infile_fp);
    vector<FASTA*> fasta(sizeof(FASTA*)*count);
    for (int i = 0; i < count; i++) {
        fasta[i] = read_fasta_file(item[i]->reference_fp);
    }
    for (int i=0; i<count; i++) {
        genome_stats[i] = new GenomeRecord;
        genome_stats[i]->identifier = fasta[i]->fasta_headers;
        genome_stats[i]->fasta_fp = item[i]->reference_fp;
        outputFile<<genome_stats[i]->identifier<<":"<<genome_stats[i]->identifier<<" "<<genome_stats[i]->fasta_fp<<endl;
    }
    return genome_stats;
}


RegionRecord** get_region(string infile_fp, int primer_length_min, int primer_length_max,string working_dirp){
    string region_stats_fp = working_dirp + "region_stats.txt";
    ofstream outputFile;
    outputFile.open(region_stats_fp);
    int count =  count_lines(infile_fp);
    RegionRecord** region_stats = (RegionRecord**)malloc(sizeof(RegionRecord*)*count);
    ITEM** item = parse_input(infile_fp);
    vector<FASTA*> fasta(sizeof(FASTA*)*count);
    vector<string> chunk_ids = generate_chunk_ids(infile_fp, primer_length_min, primer_length_max);
    
    for (int i = 0; i < count; i++) {
        fasta[i] = read_fasta_file(item[i]->reference_fp);
    }
    
    for (int i=0; i<count; i++) {
        region_stats[i] = new RegionRecord;
        region_stats[i]->identifier= item[i]->identifier;
        region_stats[i]->start= item[i]->start;
        region_stats[i]->sequence = item[i]->sequence;
        region_stats[i]->poi_record = item[i]->poi_record;
        region_stats[i]->headers = fasta[i]->fasta_headers;
        region_stats[i]->chunk_ids = chunk_ids[i];

        outputFile<<region_stats[i]->identifier<<":"<<region_stats[i]->identifier<<" "<<region_stats[i]->start<<" "<<region_stats[i]->sequence<<" "<<"["<< region_stats[i]->poi_record<<"]"<<" "<<region_stats[i]->headers<<" "<<region_stats[i]->chunk_ids<<endl;

    }
    
    return region_stats;
}

PrimerStatsRecord** get_primer(string infile_fp, int primer_length_min, int primer_length_max,string working_dirp){
    string primer_stats_fp = working_dirp + "primer_stats.txt";
    ofstream outputFile;
    outputFile.open(primer_stats_fp);
    int count =  count_lines(infile_fp);
    int total = total_primer(infile_fp, primer_length_min,primer_length_max);
     ITEM** item = parse_input(infile_fp);
    int middle_point[count];
    for (int i=0; i<count; i++) {
        middle_point[i] = (2*item[i]->poi_record - primer_length_min)/2;
        
    }
    
    string s;
    s.reserve(primer_length_max);
    vector<string>identifiers(total,s);
    vector<int>starts(total);
    vector<string>sequences(total,s);
    vector<string>amp_directions(total,s);
    vector<string>region_ids(total,s);
 
    string orientation ="5";
    PrimerStatsRecord** primer_stats = (PrimerStatsRecord**)malloc(sizeof(PrimerStatsRecord*)*total);
    
    int track = 0;
    for (int i=0; i<count; i++) {
        for (int index = 0; index<2*item[i]->poi_record; index++) {
            for (int span =0; span < primer_length_max - primer_length_min; span++) {
                if (index <= middle_point[i]) {
                    if (index + primer_length_min + span <= item[i]->poi_record-primer_length_min/2){
                        identifiers[track]+=item[i]->identifier;
                        identifiers[track]+="_";
                        identifiers[track]+=to_string(index);
                        identifiers[track]+="_";
                        identifiers[track]+=to_string(primer_length_min+span);
                        starts[track] = index;
                        sequences[track] += item[i]->sequence.substr(index,(primer_length_min+span));
                        amp_directions[track]+="r";
                        amp_directions[track]+=to_string(i+1);
                        region_ids[track] +=item[i]->identifier;
                        track++;
                    }
                }
                if (index >= middle_point[i]+primer_length_min+1) {
                    if (index + primer_length_min + span <= 2*item[i]->poi_record) {
                         identifiers[track]+=item[i]->identifier;
                        identifiers[track]+="_";
                        identifiers[track]+=to_string(index);
                        identifiers[track]+="_";
                        identifiers[track]+=to_string(primer_length_min+span);
                        starts[track] = index;
                        sequences[track] += item[i]->sequence.substr(index,(primer_length_min+span));
                        amp_directions[track]+="r";
                        amp_directions[track]+=to_string(i+1);
                        region_ids[track] +=item[i]->identifier;
                        track++;
                    }
                }
                
            }
        }
        
    }
 
    for (int i=0; i<total; i++) {
       primer_stats[i] = new PrimerStatsRecord;
       primer_stats[i]->identifier = identifiers[i];
       primer_stats[i]->start = starts[i];
       primer_stats[i]->sequence = sequences[i];
       primer_stats[i]->amp_direction = amp_directions[i];
       primer_stats[i]->orientation = orientation;
       primer_stats[i]->region_id = region_ids[i];
       outputFile<<primer_stats[i]->identifier<<":"<<primer_stats[i]->identifier<<"  "<<primer_stats[i]->start<<"  "<<primer_stats[i]->sequence<<"  "<<primer_stats[i]->amp_direction<<"  "<<primer_stats[i]->orientation<<"  "<<primer_stats[i]->region_id<<endl;
    }
 
    return primer_stats;
}



Homopolymers* get_homopolymers(string sequence){
    Homopolymers* homopolymers = new Homopolymers;
    int index = 0;
    int group = 0;//how many groups do we have
    int length = 0;
    int size = (int)sequence.size();
    char store[2*size];
    //group by
    for (int i =0; i < size; i++) {
        if (i != size-1) {
            if (sequence[i] != sequence[i+1]) {
                store[index] = sequence[i];
                index++;
                store[index] = ' ';
                index++;
            }
            else{
                store[index] = sequence[i];
                index++;
            }
        }
        else{
            if (sequence[i-1] != sequence[i]) {
                store[index] = ' ';
                index++;
                store[index] = sequence[i];
                index++;
                store[index] = ' ';
            }
            else{
                store[index] = sequence[i];
                index++;
                store[index] = ' ';
                index++;
            }
        }
    }
    
    for (int i = 0; i<index; i++) {
        if (store[i] == ' ') {
            group++;
        }
    }
    
    vector<int>lengths;
    lengths.reserve(group);
    for (int i=0; i<index+1; i++) {
        if (store[i] != ' ') {
            length++;
        }
        else{
            lengths.push_back(length);
            length = 0;
        }
    }
    homopolymers->base = 0;
    
    for (int i = 0; i < group; i++) {
        if (homopolymers->base <= lengths[i] && lengths[i] >=3) {
            homopolymers->base += lengths[i];
        }
    }
    int flg = 0;
    int position = 0;
    int max_length = *max_element(lengths.begin(),lengths.end());
    while (lengths[flg] != max_length) {
        position += lengths[flg];
        flg++;
    }
    if (max_length >=3) {
        homopolymers->longest.append("[");
        for (int m =0; m<max_length-1; m++) {
            homopolymers->longest.append(sequence.substr(position+m,1));
            homopolymers->longest.append(",");
            
        }
        homopolymers->longest.append(sequence.substr(position+max_length-1,1));
        homopolymers->longest.append("]");
    }
    else{
        homopolymers->longest.append("[");
        homopolymers->longest.append(" ");
        homopolymers->longest.append("]");
    }
    homopolymers->percent = (float)homopolymers->base/sequence.length();
    
    return homopolymers;
}


Dimerzation* get_dimerzation(string sequence){
    Dimerzation* dimerzation = new Dimerzation;
    SCORE* score = self_hybridization(sequence);
    dimerzation->base = score->most_matches;
    dimerzation->longest = score->longest_run;
    dimerzation->hairpin = score->max_hairpin_score;
    dimerzation->jackpots += ("[");
    dimerzation->jackpots += (to_string(score->l_jackpot));
    dimerzation->jackpots += (",");
    dimerzation->jackpots += (to_string(score->r_jackpot));
    dimerzation->jackpots += ("]");
    dimerzation->percent = (float)dimerzation->base/sequence.length();
    return dimerzation;
}



string run_blast(string db_location,string query_fp,string working_dirp,int word_size){
    string blast_out_fp = working_dirp + "blast.out";
    string command = "blastn -query "+ query_fp +" -out "+ blast_out_fp + " -db "+ db_location + " -outfmt 6 -num_threads 10 -word_size "+ to_string(word_size) +" -max_target_seqs 1  | sort -u > " + blast_out_fp;
 
    int pid = fork();
    if (pid == 0) {
        execlp("sh", "sh", "-c",command.c_str(),(char*)NULL);
    }
    wait(&pid);
    return blast_out_fp;
}

string make_blast_db(string subject_fp, string working_dirp){
    string db_location = working_dirp + "blast_db/" + "blast_db";
    string command = "makeblastdb -in " + subject_fp + " -out " + db_location + " -dbtype nucl ";
    int pid = fork();
    if (pid == 0) {
        execlp("sh", "sh", "-c",command.c_str(),(char*)NULL);
    }
    
    wait(&pid);
   
    return db_location;
}

string generate_query_fp(string infile_fp, int primer_length_min, int primer_length_max,string working_dirp){
    int total = total_primer(infile_fp, primer_length_min,primer_length_max);
    PrimerStatsRecord** psr = get_primer(infile_fp,primer_length_min,primer_length_max,working_dirp);
    ofstream outputFile;
    string query_fp = working_dirp + "query_fp.fasta";
    outputFile.open(query_fp);
    for (int i=0; i<total; i++) {
        outputFile<<">"<<psr[i]->identifier<<endl;
        outputFile<<psr[i]->sequence<<endl;
    }
    return query_fp;
}

string generate_subject_fp(string infile_fp, int primer_length_min, int primer_length_max,string working_dirp){
    string subject_fp = working_dirp + "subject_fp.fasta";
    int count =  count_lines(infile_fp);
    GenomeRecord** gr = get_genome(infile_fp,working_dirp);
     vector<int>unique_index(count);
    int index =1;
    for (int i=0; i<count-1; i++) {
        if (gr[i]->identifier == gr[i+1]->identifier) {
            unique_index[index-1] = i+1;
        }
        if (gr[i]->identifier != gr[i+1]->identifier) {
            index++;
            unique_index[index-1] = i+1;
        }
    }
    ofstream outputFile;
    outputFile.open(subject_fp);
    vector<FASTA*> fasta(sizeof(FASTA*)*index);
    for (int i = 0; i <index; i++) {
        fasta[i] = read_fasta_file(gr[unique_index[i]]->fasta_fp);
        outputFile<<">"<<fasta[i]->fasta_headers<<endl;
        outputFile<<fasta[i]->fasta_sequences<<endl;
    }
    return subject_fp;
}



MATCH* get_match(string infile_fp, int primer_length_min, int primer_length_max,string working_dirp){
    int total = total_primer(infile_fp, primer_length_min,primer_length_max);
    int count = count_item(infile_fp);
    PrimerStatsRecord**  primer_stats = get_primer(infile_fp, primer_length_min,primer_length_max,working_dirp);
    MATCH* match = new MATCH;
    match->perfect_match = (int*)malloc(sizeof(int)*total);
    match->near_match = (int**)malloc(sizeof(int*)*total);
    for (int i=0; i<total; i++) {
        match->perfect_match[i] = 0;
        match->near_match[i] =(int*)malloc(sizeof(int)*4);
        for (int j =0; j<4; j++) {
            match->near_match[i][j] = 0;
        }
       
    }
 
    int word_size = 7;
    string subject_fp = generate_subject_fp(infile_fp,primer_length_min,primer_length_max,working_dirp);
    string db_location = make_blast_db(subject_fp, working_dirp);
    string query_fp = generate_query_fp(infile_fp,primer_length_min,primer_length_max,working_dirp);
    string blast_out = run_blast(db_location, query_fp, working_dirp, word_size);
    //string blast_out = working_dirp + "blast.out";
    int lines = count_lines(blast_out);


    
    
    ifstream file;
    file.open(blast_out,ios::in);
    string line;
    line.reserve(primer_length_max*count);
    vector<string> blast_lines(lines,line);
 
    string s;
    s.reserve(primer_length_max);
    vector<string>query(lines,s);
    vector<string>alen(lines,s);
    vector<string>mm(lines,s);
    vector<string>gaps(lines,s);

    int index = 0;
    while (getline(file, line)) {
    
        blast_lines[index] = line;
        index++;
    }
   
    for (int i=0; i<lines; i++) {
        
        query[i] = blast_lines[i].substr(0,blast_lines[i].find_first_of('\t'));
        blast_lines[i] = blast_lines[i].substr(blast_lines[i].find_first_of('\t')+1,blast_lines[i].size());
        blast_lines[i] = blast_lines[i].substr(blast_lines[i].find_first_of('\t')+1,blast_lines[i].size());
        blast_lines[i] = blast_lines[i].substr(blast_lines[i].find_first_of('\t')+1,blast_lines[i].size());
        alen[i] =blast_lines[i].substr(0,blast_lines[i].find_first_of('\t'));
        blast_lines[i] = blast_lines[i].substr(blast_lines[i].find_first_of('\t')+1,blast_lines[i].size());
        mm[i] =blast_lines[i].substr(0,blast_lines[i].find_first_of('\t'));
       
        blast_lines[i] = blast_lines[i].substr(blast_lines[i].find_first_of('\t')+1,blast_lines[i].size());
        gaps[i] = blast_lines[i].substr(0,blast_lines[i].find_first_of('\t'));
        
    }
    int j =0;
    int track = 0;
    for (int i =0; i<lines; i++) {
        
        for (j =0; j<total; j++) {
            if (primer_stats[j]->identifier == query[i]) {
                track++;
                break;
            }
        }
        int edit_distance =  (int)primer_stats[j]->sequence.length() - stoi(alen[i]) + stoi(mm[i]) + stoi(gaps[i]);
        
        if (edit_distance == 0) {
            match->perfect_match[j] += 1;
        }
        else if (edit_distance <= 3 ){
            match->near_match[j][edit_distance-1] += 1;
        }
        else {
            match->near_match[j][3] += 1;
        }
        
    }
       
    
    return match;
}



PrimerCharRecord** get_primer_char(string infile_fp, int primer_length_min, int primer_length_max,string working_dirp){
   
    int total = total_primer(infile_fp, primer_length_min,primer_length_max);
    PrimerStatsRecord** psr = get_primer(infile_fp,primer_length_min,primer_length_max,working_dirp);
    PrimerCharRecord** primer_char = (PrimerCharRecord**)malloc(sizeof(PrimerCharRecord*)*total);
    MATCH* match = get_match(infile_fp,primer_length_min,primer_length_max,working_dirp);

    string primer_chars = working_dirp + "primer_chars.txt";
    ofstream outputFile;
    outputFile.open(primer_chars);
    for (int i=0; i<total; i++) {
        primer_char[i] = new PrimerCharRecord;
        primer_char[i]->identifier =psr[i]->identifier;

        primer_char[i]->melting_temp = get_melting_temp(psr[i]->sequence);

        primer_char[i]->gc_content = get_gc_content(psr[i]->sequence);

        primer_char[i]->length =(int) psr[i]->sequence.length();

        primer_char[i]->Hom = get_homopolymers(psr[i]->sequence);
 
        primer_char[i]->Dim = get_dimerzation(psr[i]->sequence);

        primer_char[i]->Spec = new Specificity;
        primer_char[i]->Spec->occurrences = match->perfect_match[i]-1;
       
       
        primer_char[i]->Spec->similar = (int*)malloc(sizeof(int)*4);
        primer_char[i]->Spec->similar = match->near_match[i];
        
        outputFile<<primer_char[i]->identifier<<"   "<< primer_char[i]->melting_temp<<"   "<<primer_char[i]->gc_content<<"  "<<primer_char[i]->length<<"   "<<primer_char[i]->Hom->base<<"   "<<primer_char[i]->Hom->longest<<"   "<<primer_char[i]->Hom->percent<<"  "<<primer_char[i]->Dim->base<<"   "<<primer_char[i]->Dim->longest<<"   "<<primer_char[i]->Dim->percent<<"   "<<primer_char[i]->Dim->hairpin<<"   "<<primer_char[i]->Dim->jackpots<<"   "<<primer_char[i]->Spec->occurrences<<"   "<<"["<<primer_char[i]->Spec->similar[0]<<" "<<primer_char[i]->Spec->similar[1]<<" "<<primer_char[i]->Spec->similar[2]<<" "<<primer_char[i]->Spec->similar[3]<<"]"<<endl;
        
    }
    
    return primer_char;
        
}



