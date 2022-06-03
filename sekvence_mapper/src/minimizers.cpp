#include <iostream>
#include <vector>
#include <tuple>
#include <cstring>
#include <string.h>
#include <list>
#include <iomanip>
#include <sstream>
using namespace std;

tuple<unsigned int, unsigned int> find_min(vector<tuple<unsigned int, unsigned int>> kmer_vector){
	tuple<unsigned int, unsigned int> minimizer;
	for (auto i : kmer_vector ) {
		
    		if (i == kmer_vector.front()){
    			minimizer=i;
    		}
    		else{
    			if (i < minimizer){
    				minimizer=i;
    			}
    		}
    		
    		//cout<< "Kmer " <<get<0>(i) << "poz " << get<1>(i) <<'\n';
    		
	}
	
	return minimizer;
}

string reverse_complement(const char* sequence, unsigned int sequence_len){
	string rev;
	for(int i = sequence_len; i >= 0; i--){
		switch (sequence[i]) {
		    case 'A':
		        rev+='T';
		        break;
		    case 'C':
		        rev+='G';
		        break;
		    case 'G':
		        rev+='C';
		        break;
		    case 'T':
		        rev+='A';
		        break;
		}
	}
	
	return rev;
}

unsigned int kmer_mask(unsigned int window_len){
	unsigned int mask=0;
	for (int i = 0; i < window_len; i++){
		if (i==0){
			mask=(mask << 2) + 0;
		}
		else{
			mask = (mask << 2) + 3;
		}
	}
	return mask;
}

unsigned int get_kmer_v(const char* kmer, unsigned int kmer_len, unsigned int value) {
    unsigned int i = 0;
	
    while (i < kmer_len) {
        switch (kmer[i]) {
            case 'A':
                value = (value << 2) + 0;
                break;
            case 'C':
                value = (value << 2) + 1;
                break;
            case 'G':
                value = (value << 2) + 2;
                break;
            case 'T':
                value = (value << 2) + 3;
                break;
        }
        i++;
    }
    return value;
}


vector<tuple<unsigned int, unsigned int, bool>>
Minimize(const char *sequence, unsigned int sequence_len, unsigned int kmer_len, unsigned int window_len, bool origin)
{
    string tr_seq;
    if(origin == false){
    
    	tr_seq = reverse_complement((sequence), sequence_len).c_str();
    	
    }
    else{
    	tr_seq = (sequence);
    }
    
    vector<tuple<unsigned int, unsigned int, bool>> minimizers{};
    unsigned int fract_len;
    fract_len = window_len + kmer_len -1;
    vector<tuple<unsigned int, unsigned int>> kmer_vector;
    unsigned int kmer_value;
    tuple<unsigned int, unsigned int> minimizer;
    int mask = kmer_mask(kmer_len);
    
    for(unsigned int i=0; i<= sequence_len-1; i++){
    	int z=0;
	int pomak = 0;
    	if (i==0){
    		
        	kmer_value = get_kmer_v(&(tr_seq.c_str()[i]), kmer_len, 0);
        	kmer_vector.push_back(make_tuple(kmer_value, i+1));
    		minimizer=make_tuple(kmer_value, i+1);		
    		z=1;
    	}
    	else if (i < window_len){
    		
    		int prev_kmer = get<0>(kmer_vector.back());
    		prev_kmer = prev_kmer & mask;
    		kmer_value = get_kmer_v(&(tr_seq.c_str()[i+kmer_len-1]) , 1, prev_kmer);
    		kmer_vector.push_back(make_tuple(kmer_value, i+1));
    		if (kmer_value < get<0>(minimizer)){
    			minimizer= make_tuple(kmer_value, i+1);
    			z=1;
    		}
    		//cout<<&(sequence[i+window_len-1])<<' ';
    		//cout<<minimizer<<'\n';
    		
    	}
    	else if (i <= sequence_len-kmer_len){
    		tuple<unsigned int, unsigned int> pot_min = kmer_vector.front();
    		int prev_kmer = get<0>(kmer_vector.back());
    		
    		prev_kmer = prev_kmer & mask;
    		kmer_value = get_kmer_v(&(tr_seq.c_str()[i+kmer_len-1]) , 1, prev_kmer);
    		kmer_vector.push_back(make_tuple(kmer_value, i+1));
    		kmer_vector.erase(kmer_vector.begin());
    		
    		if (pot_min == minimizer){
    			
    			minimizer = find_min(kmer_vector);
    			z=1;
    			
    			
    		}
    		else{
	    		if (kmer_value < get<0>(minimizer)){
	    			minimizer= make_tuple(kmer_value, i+1);
	    			z=1;
	    		}
		}	
    		
    	}
    	else{
    		if (!kmer_vector.empty()){
	    		tuple<unsigned int, unsigned int> pot_min = kmer_vector.front();
	    		
	    		kmer_vector.erase(kmer_vector.begin());
	    		
	    		if (pot_min == minimizer){
	    			minimizer = find_min(kmer_vector);
	    			
	    			z=1;
	    			
	    		}
    		}
    		//cout<<&(sequence[i+window_len-1])<<' ';
    		//cout<<minimizer<<'\n';
    		
    	}
    	
    	if(z==1){
    		minimizers.push_back(make_tuple(get<0>(minimizer), get<1>(minimizer), origin));
    	}
    	
    	
    } 
    
    return minimizers;
}

