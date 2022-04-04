#include <gtest/gtest.h>
#include "../src/minimizers.h"


TEST(kmer_tests, test1){
	std::vector<std::tuple<unsigned long int, unsigned int, bool>> result = new std::vector<std::tuple<unsigned long int, unsigned int, bool>>;
	const char* target = "GTCATGCACGTTCAC"
	unsigned int target_len = 15
	unsigned int kmer_len = 3
	unsigned int window_len = 3
	result = Minimize(target, target_len, kmer_len, window_len)
	EXPECT_EQ(result, 0);
}
