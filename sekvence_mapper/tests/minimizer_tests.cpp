#include <gtest/gtest.h>
#include "../src/minimizers.h"
#include <iostream>

using namespace std;

TEST(kmer_tests, test1){
	const char* target = "GTCATGCACGTTCAC";
	unsigned int target_len = 15;
	unsigned int kmer_len = 3;
	unsigned int window_len = 4;
	vector<tuple<unsigned int, unsigned int, bool>> result = Minimize(target, target_len, kmer_len, window_len,1);
	
	vector<tuple<unsigned int, unsigned int, bool>> test;
	test.push_back(make_tuple(45, 1, true));
	test.push_back(make_tuple(19, 3, true));
	test.push_back(make_tuple(14, 4, true));
	test.push_back(make_tuple(6, 8, true));
	test.push_back(make_tuple(27, 9, true));
	test.push_back(make_tuple(17, 13, true));
	EXPECT_EQ(kmer_len, 3);
	EXPECT_EQ(test, result);
}

TEST(kmer_tests, test2){
	const char* target = "GTCATGCACGTTCAC";
	unsigned int target_len = 15;
	unsigned int kmer_len = 3;
	unsigned int window_len = 3;
	vector<tuple<unsigned int, unsigned int, bool>> result{};
	result = Minimize(target, target_len, kmer_len, window_len,1);
	
	vector<tuple<unsigned int, unsigned int, bool>> test;
	test.push_back(make_tuple(45, 1, true));
	test.push_back(make_tuple(19, 3, true));
	test.push_back(make_tuple(14, 4, true));
	test.push_back(make_tuple(17, 7, true));
	test.push_back(make_tuple(6, 8, true));
	test.push_back(make_tuple(27, 9, true));
	test.push_back(make_tuple(47, 10, true));
	test.push_back(make_tuple(17, 13, true));
	
	ASSERT_EQ(result.size(), test.size()) << "Vectors size and test are of unequal length";
	EXPECT_EQ(result, test);
	
}

TEST(kmer_tests, test3){
	const char* target = "TACGCATAAGCGCCAAAAGCACGCCGGGCGACCATAATGAGAAGATGCTCACCGCCAGCGTCAACAGTAAGTAGACGATGAGAGCGCCAGGGTAAACAGTTTGCACGGGCCGCGACGAAGGCAGTCGATAATGCCGCCATTGCTGGGATGCTCGCCCCCAGACGCGCATAGGCATAACCGGAAAACATCGCCACAATACCGCCAAAGCAAAGGCGACCCCAGGTCGAGGCTTCCATTAGCAATGCAGCCTGCCCCAGCAGCGCGAAGATCGCCGCCCCCACCATTGCCCCGGATTATCGATGGAAACGACGTTCCAGACCGCGGGTTTGTTACCGTTATTACCTTCCGTGTTTCATCATGTAATAGCAGCCCTTAGTAAACACGTTATAGCCGAAAAATTGCTTAGCGACCGATGCCACCATTTTTGACCGCGAACACTGTTGCCATCTTCGTCCAGCGGTGGTGAGAGAACGCGGCAATTCCCATCACTCCAGGGGACGACCACGCCAGAATACCGCCACCTACACCCGGTTTACGCCCGGTAAACCAACACGATACGCCCAGTCACCGGGAACATACCCAAGACTACAGCCCTTCCATCATCATTTCGTAAGAATGTACGGCACGTTGTCGGCCTGAAGAACGCGTTCCTGCGTCAACGATTCGCACACCACCTGCGCCGCGCTTGTCGCGCCAAGCGTTGCCAGTTCGCTACGTGGTCGAGGAGCCCTTGGAGCACTGACGGGTATACACGTCACAGGCTTCCCATTGCATCACAATAGAGATATCCGGCGGAGTAAGCAGCCAGGCTATGGCCCGGTTATGGAGACTTTATTTGTTTGTTCCGACTGGTTAACTTCGTCAGAGAGCGCTACCTGCTCGCCAGCAGTTTGCTTTGGATATGTAAATTCGCTGCCAGCGTTGTTCAACAGTTTTCAGCGTTAATCAGGCTATCTTAGGCAATAGCGCCAGGATAGATTCACCAGTGGCGGGAAATCGTTTGCCGCGGATGCAACTCTAAGGCGATAACTGAGTTACGGGGCAAGATCCGGTCGGGTCCCCACTAGGGATCTTTGTGTCGTACCGCCTGCGGGCCGCACATCTTCCAACGCAAGGGCTAACGTACAGAGACTTTCGAGATGGATTCCAGTGCAAAGCGTAATCACTGTCACCGCACTACGACGTTGCCATCGCGCAGTACACGATAGCCACTGCCGCCAGTTGACCTGGTACATTCGCCGAAGGGAATGTAATCGGGCATTTTGTCCGCCGTTAAGTGAGTGAGAAATTGGGTGTAAGCCTGATCTGTTACGGAATTGCCTGCTGTAATTTGTTTGCATCTAACAATATTTTCTTGTTAACTCCTTTTATAAGTCTCGGGAGGTAATTCCTCACCGGCTGGTGCCGATTTCAGGCATCCTGATTTAACTTAGCACCGGAGACTTAACTACAGGAAAACACAAAGAGATAAATGTCTAATCCTGATGCAAAAACTGCATCAGCAAATTTTTAATCTTTACGGACTTTTACCGCCTGGTTTATTAATTTCTTGACCTTCCCCTTGCTGAAGGTTTAACCTTTATCACAGCCAGTCAAAACCGTGTGTAAAGGGATGTTTATGTCAAACTATCGACCTGACCCTCTGGACGGCCTGTCTGCGGCGGTCACTGCGTTAAACGCGTGAAAGAAAGTCTTGAACAGCGTCGTCCGGGATGTTGAGCAGGCGGATGTGTCTATCACTGAAGCACGTTACCGGGACTGCCGAACAGGCTCGCTAAATTGAAACCATCAAACAAGCGGGTTATGACGCATCTGTAAGCCACCCAAAGGCTAAACCGCCACAGCCGTGACTGTTATGTT";
	unsigned int target_len = strlen(target);
	unsigned int kmer_len = 15;
	unsigned int window_len = 5;
	vector<tuple<unsigned int, unsigned int, bool>> result{};
	result = Minimize(target, target_len, kmer_len, window_len, 1);
	
	vector<tuple<unsigned int, unsigned int, bool>> test;
	EXPECT_EQ(result.size(), 0);
	EXPECT_EQ(result, test);
}

TEST(kmer_tests, test4){
	const char* target = "TACGCATAAGCGCCAAAAGCACGCCGGGCGACCATAATGAGAAGATGCTCACCGCCAGCGTCAACAGTAAGTAGACGATGAGAGCGCCAGGGTAAACAGTTTGCACGGGCCGCGACGAAGGCAGTCGATAATGCCGCCATTGCTGGGATGCTCGCCCCCAGACGCGCATAGGCATAACCGGAAAACATCGCCACAATACCGCCAAAGCAAAGGCGACCCCAGGTCGAGGCTTCCATTAGCAATGCAGCCTGCCCCAGCAGCGCGAAGATCGCCGCCCCCACCATTGCCCCGGATTATCGATGGAAACGACGTTCCAGACCGCGGGTTTGTTACCGTTATTACCTTCCGTGTTTCATCATGTAATAGCAGCCCTTAGTAAACACGTTATAGCCGAAAAATTGCTTAGCGACCGATGCCACCATTTTTGACCGCGAACACTGTTGCCATCTTCGTCCAGCGGTGGTGAGAGAACGCGGCAATTCCCATCACTCCAGGGGACGACCACGCCAGAATACCGCCACCTACACCCGGTTTACGCCCGGTAAACCAACACGATACGCCCAGTCACCGGGAACATACCCAAGACTACAGCCCTTCCATCATCATTTCGTAAGAATGTACGGCACGTTGTCGGCCTGAAGAACGCGTTCCTGCGTCAACGATTCGCACACCACCTGCGCCGCGCTTGTCGCGCCAAGCGTTGCCAGTTCGCTACGTGGTCGAGGAGCCCTTGGAGCACTGACGGGTATACACGTCACAGGCTTCCCATTGCATCACAATAGAGATATCCGGCGGAGTAAGCAGCCAGGCTATGGCCCGGTTATGGAGACTTTATTTGTTTGTTCCGACTGGTTAACTTCGTCAGAGAGCGCTACCTGCTCGCCAGCAGTTTGCTTTGGATATGTAAATTCGCTGCCAGCGTTGTTCAACAGTTTTCAGCGTTAATCAGGCTATCTTAGGCAATAGCGCCAGGATAGATTCACCAGTGGCGGGAAATCGTTTGCCGCGGATGCAACTCTAAGGCGATAACTGAGTTACGGGGCAAGATCCGGTCGGGTCCCCACTAGGGATCTTTGTGTCGTACCGCCTGCGGGCCGCACATCTTCCAACGCAAGGGCTAACGTACAGAGACTTTCGAGATGGATTCCAGTGCAAAGCGTAATCACTGTCACCGCACTACGACGTTGCCATCGCGCAGTACACGATAGCCACTGCCGCCAGTTGACCTGGTACATTCGCCGAAGGGAATGTAATCGGGCATTTTGTCCGCCGTTAAGTGAGTGAGAAATTGGGTGTAAGCCTGATCTGTTACGGAATTGCCTGCTGTAATTTGTTTGCATCTAACAATATTTTCTTGTTAACTCCTTTTATAAGTCTCGGGAGGTAATTCCTCACCGGCTGGTGCCGATTTCAGGCATCCTGATTTAACTTAGCACCGGAGACTTAACTACAGGAAAACACAAAGAGATAAATGTCTAATCCTGATGCAAAAACTGCATCAGCAAATTTTTAATCTTTACGGACTTTTACCGCCTGGTTTATTAATTTCTTGACCTTCCCCTTGCTGAAGGTTTAACCTTTATCACAGCCAGTCAAAACCGTGTGTAAAGGGATGTTTATGTCAAACTATCGACCTGACCCTCTGGACGGCCTGTCTGCGGCGGTCACTGCGTTAAACGCGTGAAAGAAAGTCTTGAACAGCGTCGTCCGGGATGTTGAGCAGGCGGATGTGTCTATCACTGAAGCACGTTACCGGGACTGCCGAACAGGCTCGCTAAATTGAAACCATCAAACAAGCGGGTTATGACGCATCTGTAAGCCACCCAAAGGCTAAACCGCCACAGCCGTGACTGTTATGTT";
	unsigned int target_len = strlen(target);
	unsigned int kmer_len = 15;
	unsigned int window_len = 5;
	vector<tuple<unsigned int, unsigned int, bool>> result{};
	result = Minimize(target, target_len, kmer_len, window_len, 0);
	
	vector<tuple<unsigned int, unsigned int, bool>> test;
	EXPECT_EQ(result.size(), 0);
	EXPECT_EQ(result, test);
}

TEST(kmer_tests, test5){
	const char* target = "CGGGGACGATGCGTCCAGAGATACTGCACCGTGGTAAACCCAGTTTATATTTTCGATTTCAACCCCGCGATCCGTGAAAGCCGTCGTTTGCA";
	
	unsigned int target_len = strlen(target);
	unsigned int kmer_len = 15;
	unsigned int window_len = 5;
	vector<tuple<unsigned int, unsigned int, bool>> result{};
	result = Minimize(target, target_len, kmer_len, window_len, 1);
	
	vector<tuple<unsigned int, unsigned int, bool>> test;
	EXPECT_EQ(result.size(), 0);
	EXPECT_EQ(result, test);
	
}

