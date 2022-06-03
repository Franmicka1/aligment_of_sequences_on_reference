#include <iostream>
#include <getopt.h>

#include <stdlib.h>
#include <string.h>
#include <memory>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <stdio.h>
#include <fstream>

#include "src/alignment.h"
#include "src/minimizers.h"
#include <bits/stdc++.h>

#include "bioparser/fasta_parser.hpp"
#include "thread_pool/thread_pool.hpp"

using namespace std;

//INICIJALIZACIJA


bool cigar_flag;
int kmer_len = 15;
int window_len = 5;
double minimizer_freq = 0.001;
int match_cost = 1;
int mismatch_cost = - 1;
int gap_cost = - 1;
int thread_num = 2;
AlignmentType type = GLOBAL;
int algorithm;
bool optimize = true;

static std::string HELP = "-h or --help for help\n"
                          "-v or --version for version\n"
                          "-m value for matching\n"
                          "-n value for mismatching\n"
                          "-g value for gap\n"
                          "-a alignment type (0 - GLOBAL, 1 - LOCAL, 2 - SEMI-GLOBAL)\n"
                          "-c cigar string enabled\n"
                          "-k k-mer length\n"
                          "-w window length\n"
                          "-f top f frequent minimizers that will not be taken in account\n"
                          "-o turn the optimization algorithm on (only for GLOBAL align)\n"
                          "query name of the file with redings in FASTQ format\n"
                          "target name of the reference file in FASTA format\n";

class Sequence
{
public:
    Sequence(const char *name, std::uint32_t name_len,
             const char *data, std::uint32_t data_len) : names(name, name_len),
                                                         all_data(data, data_len)
    {
    }
    std::string names, all_data;
};


//POCETAK FUNCKIJA

vector<tuple<bool, int,unsigned int>> longest_increasing_subsequence(vector<tuple<bool, unsigned int, unsigned int>>& matches) {

    int N = matches.size();
    int X[N];   
    auto itr = 0;
    for (auto e : matches){
    	X[itr] = get<2>(e);
    	++itr;
    	
    }
    int P[N];
    int M[N + 1];
    int L = 0;
    int lo, hi, mid, newL;
    for (int i=0; i < N; i++){
    	lo = 1;
    	hi = L + 1;
    	while (lo < hi){
    	    mid = lo + floor((hi-lo)/2);
    	    if (X[M[mid]] < X[i]){
    	    	lo = mid + 1;
    	    }
    	    else{
    	    	hi = mid;
    	    }
    	}
    	newL = lo;
    	
    	P[i] = M[newL - 1];
    	M[newL] = i;
    	
    	if (newL > L){
    		L = newL;
    	}
    
    }
    vector<tuple<bool, int,unsigned int>> result;
    int k = M[L];
    for (int i=L-1; i > 0; i--){
    	if (result.empty()){
    		result.push_back(matches[k]);
    	}
    	
    	if ((int(get<1>(result.back()))- int(get<1>(matches[k]))) < 5000){
    		result.push_back(matches[k]);
    	}
    	else{
    		result.pop_back();
    	}
    	
    	k = P[k];
    }
    
    int cnt1=0, cnt2 = 0;
    if(!result.empty()){
	     for(auto i:result){
	     	if (get<0>(i)==1){
	     		cnt1+=1;
	     	}
	     	else cnt2+=1;
	     }
	     if (cnt1 > cnt2){
	     	for (int i=0; i <= cnt2; i++){
	     		result.pop_back();
	     	}
	     }
	     else {
	     	for (int i=0; i <= cnt1; i++){
	     		result.erase(result.begin());
	     	}
	     }
	     //cout << cnt1 << "   " <<cnt2 << "\n";
    }
    return result;

}


/*string reverse_complement(const char* sequence, unsigned int sequence_len){
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
*/
void make_minimizer_index(const std::unique_ptr<Sequence> &sequence,
                          std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>> &index,
                          bool origin)
{
    
    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers = 
    Minimize(sequence->all_data.c_str(),sequence->all_data.size(), kmer_len, window_len, origin);
    
    for (auto minimizer : minimizers)
    { 
        index[std::get<0>(minimizer)].emplace_back(std::make_pair(std::get<1>(minimizer), std::get<2>(minimizer)));
        
    }        
}

void make_reference_index(const std::unique_ptr<Sequence> &sequence,
                          std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>> &index)
{
    cerr << "Making reference index...\n"; 
    bool origin = true;
    make_minimizer_index(sequence, index, origin);
    origin = false;
    make_minimizer_index(sequence, index, origin);
    std::vector<std::pair<unsigned int, unsigned int>> min_occ;
    min_occ.reserve(index.size());
    int num_singl = 0;
    for (auto &i : index)
    {
        if (i.second.size() == 1)
        {
            num_singl += 1;
        }
        min_occ.emplace_back(std::make_pair(i.second.size(), i.first));
    }
    cerr << "Sorting minimizers by occurence...\n";
    sort(min_occ.begin(), min_occ.end(), std::greater<std::pair<unsigned int, unsigned int>>());
    size_t skip = ceil(index.size() * minimizer_freq);
    
    int index_size = index.size();
    if (skip >= index_size)
    {
        skip = index_size - 1;
    }

    cerr << "Number of distinct minimizers: " << index.size() << "\n";
    cerr << "Fraction of singletons: " << ((double) num_singl / index.size()) << "\n"; 
    cerr << "NUmber of occurences of the most frequent minimizers with top " << minimizer_freq * 100 << "% most frequent ignored: " << skip << "\n";
    
    for (int i = 0; i < skip; i++) {
        index.erase(min_occ[i].second);
    }
    
    
}

void best_match_cluster(unordered_map<unsigned int, vector<pair<unsigned int, bool>>>& fragment_index,
                        unordered_map<unsigned int, vector<pair<unsigned int, bool>>>& reference_index,
                        vector<tuple<bool, int, unsigned int>>& result) {

        
    vector<tuple<bool, unsigned int, unsigned int>> matches;
    for (auto e : fragment_index) {
        if (reference_index.count(e.first) != 0) {
            for (auto location : e.second) {
                for (auto r_location : reference_index[e.first]) {
                   
                    matches.emplace_back(make_tuple(r_location.second , r_location.first, location.first));
                }
            }
        }
    }
    if (!matches.empty()) {
        sort(matches.begin(), matches.end());
        result = longest_increasing_subsequence(matches);
    	/*for (auto i: result){
    	    
    	    cout << get<1>(i)<< " : " << get<2>(i) << " - "<< get<0>(i) <<'\t';
    	}
    	*/	
    }
    
}

string paf_format(const unique_ptr<Sequence>& reference, 
const unique_ptr<Sequence>& fragment, vector<tuple<bool, int, unsigned int>> result, string* cigar){
    
    string paf = "";
    
    int tr_broj = 0;
    int ukuZb = 0; 
    int matchZb = 0;
    
    
    /*
    for (char& i : *cigar){
        
        if (isdigit(i)){
            tr_broj = (tr_broj* 10) + (int(i)- 48); 
        }
        else{
            ukuZb+= tr_broj;
            if (i == 'M'){
                matchZb+=tr_broj;
            }
            tr_broj=0;
        }
    }*/
    
    if (result.size()<=5){
    	paf += "Skipping fragment: ";
    	paf += (fragment->names.c_str());
    	paf += '\t';
    	paf += to_string(fragment->all_data.size()); 
    	paf += " because not enough matches were found.\n";
    	return paf;
    }
    
    paf += fragment->names.c_str();
    paf += '\t';
    paf += to_string(fragment->all_data.size());
    paf += '\t';
    paf += to_string(get<2>(result[result.size()-1]));
    paf += '\t';
    paf += to_string(get<2>(result[0]));
    paf += '\t';
    if(get<0>(result[1]) == 0){
    	paf += "-";
    	paf += '\t';
    	paf += reference->names.c_str();
	paf += '\t';
	paf += to_string(reference->all_data.size());
	paf += '\t';
	paf += to_string(reference->all_data.size() - get<1>(result[0]));
	paf += '\t';
	paf += to_string(reference->all_data.size() - get<1>(result[result.size()-1]));
	
	
    }
    else{
    	paf += "+";
    	paf += '\t';
    	paf += reference->names.c_str();
    	paf += '\t';
    	paf += to_string(reference->all_data.size());
    	paf += '\t';
    	paf += to_string(get<1>(result[result.size()-1]));
    	paf += '\t';
    	paf += to_string(get<1>(result[0]));
    }
    
    
    paf += '\n';
    
    
    
    /*
    result += to_string(matchZb);
    result += '\t';
    result += to_string(ukuZb);
    result += '\t';
    result += '\n';
    if (cigar_flag) {
        result += "cg:Z";
        result += *cigar;
    }
    */

    return paf;
}


void ispis_opcenito(std::vector<unique_ptr<Sequence>> &fragments, std::vector<unique_ptr<Sequence>> &reference){
    int frag_size = (int)fragments.size();
    int ref_size = (int)reference.size();
    cerr << "Reference genome sequences: \n";
    for (int i = 0; i < int(reference.size()); i++) {
        std::cerr << "\t" << reference[i]->names.c_str() << ", length = " << reference[i]->all_data.size() << "\n";
        cerr << "\n";
    }
    cerr << "Number of fragments: " << fragments.size() << "\n";
    uint64_t len_sum = 0;
    vector<size_t> lengths(fragments.size());
    for (int i = 0; i < int(fragments.size()); i++) {
    
        lengths[i] = fragments[i]->all_data.size();
        len_sum += lengths[i];
    }
    sort(lengths.begin(), lengths.end(), greater<size_t>());

    uint64_t N50 = -1, tmp_sum = 0;
    for (int i = 0; i < int(fragments.size()); i++) {
    
        tmp_sum += lengths[i];
        if (tmp_sum * 2 >= len_sum) {
            N50 = lengths[i];
            break;
        }
    }

    cerr << "Average length: " << len_sum * 1.0 / fragments.size() << "\n";
    cerr << "N50 length: " << N50 << "\n";
    cerr << "Minimal length: " << lengths.back() << "\n";
    cerr << "Maximal length: " << lengths.front() << "\n";
    
    int sum = 0;
    for (int i = 0; i < frag_size; i++) {
        sum += fragments[i]->all_data.size();
    } 
    float avg_size = sum / frag_size;
    std::vector<size_t> fragment_vector;
    for (int i = 0; i < frag_size; i++) {
        fragment_vector.push_back(fragments[i]->all_data.size());
    }
    std::sort(fragment_vector.begin(), fragment_vector.end());
    fragment_vector.clear();

}

void mapping(unordered_map<unsigned int, vector<pair<unsigned int, bool>>> &reference_index, vector<unique_ptr<Sequence>> &fragments, vector<unique_ptr<Sequence>> &reference, int frag_begin, int frag_per_thread){
    
    string cigar;
    string target = "";
    ofstream MyFile("mapping.paf");
    int target_len = 0;
    int query_len = 0;
    int align_score = 0;
    unsigned int target_begin;
    
    for (int i = frag_begin; i < frag_per_thread; i++) {
    //for (int i= 10; i < 20; i++){
        std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>> fragment_index;
        //cout << fragments[i]->all_data.size() << "\n";
        make_minimizer_index(fragments[i], fragment_index, true);
        vector<tuple<bool, int, unsigned int>> result;
        
        best_match_cluster(fragment_index, reference_index, result);
        

        if (result.size() < 10){
        	continue;
        }

        target_len =int(get<1>(result[0])) - int(get<1>(result[result.size()-1]));
        query_len =int(get<2>(result[0])) - int(get<2>(result[result.size()-1]));  
        //cout << query_len <<  "    "<< target_len << "\n";
        
        
        auto query = fragments[i]->all_data.substr(int(get<2>(result[result.size()-1])), query_len);
        
        if (!(get<0>(result[0]))){
            auto target_r = reference.front()->all_data.substr((reference.front()->all_data.size()) - int(get<1>(result[0])), target_len);
            target = reverse_complement(target_r.c_str(), target_len);
        }
        else {
            auto target_r = reference.front()->all_data.substr(int(get<1>(result[result.size()-1])), target_len);
            target = target_r;
        }
        //cout << query <<  "\n"<< target << "\n";
        
        if (optimize && type == GLOBAL) {

        	align_score = LinearnaSloz(query.c_str(), query_len, target.c_str(), target_len, match_cost, mismatch_cost, gap_cost, &cigar);
    	} 
    	else {
        	align_score = Align(query.c_str(), query_len, target.c_str(), target_len, type, match_cost, mismatch_cost, gap_cost, &cigar, &target_begin);
    	}
    
	cout << align_score << "\n";
        
        string res = paf_format(reference.front(), fragments[i], result, &cigar);

        MyFile << res;    
    }
    MyFile.close();
    
}





void display_version()
{
    std::cout << "v0.1.0"
              << "\n";
}

void display_help()
{
    std::cout << HELP
              << "\n";
}

//MAIN PROGRAM 

int main(int argc, char *argv[])
{
    int option;
    const char *optstring = "m:g:n:a:k:w:f:t:hvco";

    while ((option = getopt(argc, argv, optstring)) != -1)
    {
        switch (option)
        {
        case 'v':
            display_version();
            return 0;
            break;
        case 'h':
            display_help();
            return 0;
            break;
        case 'm':
            match_cost = atoi(optarg);
            break;
        case 'n':
            mismatch_cost = atoi(optarg);
            break;
        case 'g':
            gap_cost = atoi(optarg);
            break;
        case 'a':
            algorithm = atoi(optarg);
            break;
        case 'k':
            kmer_len = atoi(optarg);
            break;
        case 'w':
            window_len = atoi(optarg);
            break;
        case 'f':
            minimizer_freq = atof(optarg);
            break;
        case 't':
            thread_num = atoi(optarg);
            break;
        case 'c':
            cigar_flag = true;
            break;
        case 'o':
            optimize = true;
            break;

        default:
            std::cout<<HELP;
                     
            exit(1);
        }
    }
    switch (algorithm)
    {
        case 0:
            type = GLOBAL;
            break;
        case 1:
            type = LOCAL;
            break;
        case 2:
            type = SEMI_GLOBAL;
            break;
        default:
            break;
    }

    if (optind < argc)
    {
        std::string path = argv[optind];   
        
    }
    
    
    //ucitavanje i parsiranje reference
    auto ref = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(argv[argc - 2]);
    auto reference = ref->Parse(-1);
    
    //stvaranje indeksa reference
    unordered_map<unsigned int, vector<pair<unsigned int, bool>>> reference_index;
    make_reference_index(reference.front(), reference_index);
    
    //ucitavanje i parsiranje fragmenata
    auto frag = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(argv[argc - 1]);
    auto fragments = frag->Parse(-1);
    
    //ispis opcenitih podataka o referenci i fragmentima 
    ispis_opcenito( fragments, reference);
    
    
    int frag_per_thread = int(fragments.size() / (double)thread_num);
    int frag_begin = 0;
    thread_pool::ThreadPool tp{};
    std::vector<std::future<void>> f;
    for (std::size_t i = 0; i < thread_num-1; ++i) {
    	f.emplace_back(tp.Submit(mapping, std::ref(reference_index), std::ref(fragments), std::ref(reference), frag_begin, frag_per_thread));
	frag_begin += frag_per_thread;
    }
    f.emplace_back(tp.Submit(mapping, std::ref(reference_index), std::ref(fragments), std::ref(reference), frag_begin, int(fragments.size()))); 
    for (auto& it : f) {
        it.get();
    }
    

    

    return 0;
}
