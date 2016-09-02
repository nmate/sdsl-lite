#include "sdsl/bit_vectors.hpp" // for r3d3_vector
#include "gtest/gtest.h"
#include <string>

using namespace sdsl;
using namespace std;

string test_file;

namespace
{

template<class T>
class bit_vector_test : public ::testing::Test { };

using testing::Types;

typedef Types<
rrri_vector<64>,
rrri_vector<256>,
rrri_vector<129>,
rrri_vector<192>,
rrri_vector<255>,
rrri_vector<15>,
rrri_vector<31>,
rrri_vector<63>,
rrri_vector<83>,
rrri_vector<127>,
rrri_vector<128>,
r3d3_vector<64>,
r3d3_vector<256>,
r3d3_vector<129>,
r3d3_vector<192>,
r3d3_vector<255>,
r3d3_vector<15>,
r3d3_vector<31>,
r3d3_vector<63>,
r3d3_vector<83>,
r3d3_vector<127>,
r3d3_vector<128>,
r3d3i_vector<64>,
r3d3i_vector<256>,
r3d3i_vector<129>,
r3d3i_vector<192>,
r3d3i_vector<255>,
r3d3i_vector<15>,
r3d3i_vector<31>,
r3d3i_vector<63>,
r3d3i_vector<83>,
r3d3i_vector<127>,
r3d3i_vector<128>
> Implementations;



TYPED_TEST_CASE(bit_vector_test, Implementations);


//! Test operator[]
TYPED_TEST(bit_vector_test, access)
{
    static_assert(sdsl::util::is_regular<TypeParam>::value, "Type is not regular");
    bit_vector bv;
    ASSERT_TRUE(load_from_file(bv, test_file));
    TypeParam c_bv(bv);
    ASSERT_EQ(bv.size(), c_bv.size());
    for (uint64_t j=0; j < bv.size(); ++j) {
        ASSERT_EQ((bool)(bv[j]), (bool)(c_bv[j]));
    }
    TypeParam mo_bv = TypeParam(bv);
    ASSERT_EQ(bv.size(), mo_bv.size());
    for (uint64_t j=0; j < bv.size(); ++j) {
        ASSERT_EQ((bool)(bv[j]), (bool)(c_bv[j]));
    }
}

//Todo-1: implement get_int in r3d3_vector
//TYPED_TEST(bit_vector_test, get_int)
//{
//    bit_vector bv;
//    ASSERT_TRUE(load_from_file(bv, test_file));
//    TypeParam c_bv(bv);
//    ASSERT_EQ(bv.size(), c_bv.size());
//    uint8_t len = 63;
//    for (uint64_t j=0; j+len < bv.size(); j+=len) {
//        ASSERT_EQ(bv.get_int(j, len), c_bv.get_int(j, len));
//    }
//}

//TYPED_TEST(bit_vector_test, get_int_all_block_sizes)
//{
//    bit_vector bv(10000, 0);
//    std::mt19937_64 rng;
//    std::uniform_int_distribution<uint64_t> distribution(0, 9);
//    auto dice = bind(distribution, rng);
//    for (size_t i=1001; i < bv.size(); ++i) {
//        if (0 == dice())
//            bv[i] = 1;
//    }
// 
//    TypeParam c_bv(bv);
//    for (uint8_t len=1; len<=64; ++len) {
//        for (size_t i=0; i+len <= bv.size(); ++i) {
//            ASSERT_EQ(bv.get_int(i,len), c_bv.get_int(i,len))
//                    << "i="<<i<<" len="<<(int)len<<endl;
//        }
//    }
//}

TYPED_TEST(bit_vector_test, swap)
{
    bit_vector bv;
    ASSERT_TRUE(load_from_file(bv, test_file));
    TypeParam c_bv(bv);
    ASSERT_EQ(bv.size(), c_bv.size());
    TypeParam bv_empty;
    ASSERT_EQ((uint64_t)0, bv_empty.size());
    bv_empty.swap(c_bv);
    ASSERT_EQ((uint64_t)0, c_bv.size());
    ASSERT_EQ(bv.size(), bv_empty.size());
    for (uint64_t j=0; j < bv.size(); ++j) {
        ASSERT_EQ((bool)(bv[j]), (bool)(bv_empty[j]));
    }
}

}// end namespace

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc < 2) {
        // LCOV_EXCL_START
        cout << "Usage: " << argv[0] << " FILE " << endl;
        cout << "  Reads a bitvector from FILE and executes tests." << endl;
        return 1;
        // LCOV_EXCL_STOP
    }
    test_file = argv[1];
    return RUN_ALL_TESTS();
}

