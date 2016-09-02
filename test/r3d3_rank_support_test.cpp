#include "sdsl/bit_vectors.hpp"
#include "sdsl/rank_support.hpp"
#include "gtest/gtest.h"
#include <string>

using namespace sdsl;
using namespace std;

string test_file;

namespace
{

template<class T>
class rank_support_test : public ::testing::Test { };

using testing::Types;

typedef Types<
rank_support_rrri<1, 64>,
rank_support_rrri<1, 192>,
rank_support_rrri<1, 256>,
rank_support_rrri<1, 255>,
rank_support_rrri<1, 15>,
rank_support_rrri<1, 31>,
rank_support_rrri<1, 63>,
rank_support_rrri<1, 83>,
rank_support_rrri<1, 127>,
rank_support_rrri<1, 128>,
rank_support_rrri<1, 129>,
rank_support_r3d3<1, 64>,
rank_support_r3d3<1, 192>,
rank_support_r3d3<1, 256>,
rank_support_r3d3<1, 255>,
rank_support_r3d3<1, 15>,
rank_support_r3d3<1, 31>,
rank_support_r3d3<1, 63>,
rank_support_r3d3<1, 83>,
rank_support_r3d3<1, 127>,
rank_support_r3d3<1, 128>,
rank_support_r3d3<1, 129>,
rank_support_r3d3i<1, 64>,
rank_support_r3d3i<1, 192>,
rank_support_r3d3i<1, 256>,
rank_support_r3d3i<1, 255>,
rank_support_r3d3i<1, 15>,
rank_support_r3d3i<1, 31>,
rank_support_r3d3i<1, 63>,
rank_support_r3d3i<1, 83>,
rank_support_r3d3i<1, 127>,
rank_support_r3d3i<1, 128>,
rank_support_r3d3i<1, 129>
> Implementations;

TYPED_TEST_CASE(rank_support_test, Implementations);

//! Test the rank method
TYPED_TEST(rank_support_test, rank_method)
{
    static_assert(sdsl::util::is_regular<TypeParam>::value, "Type is not regular");
    bit_vector bvec;
    ASSERT_TRUE(load_from_file(bvec, test_file));
    typename TypeParam::bit_vector_type bv(bvec);
    TypeParam rs(&bv);
    uint64_t rank=0;
    for (uint64_t j=0; j < bvec.size(); ++j) {
        ASSERT_EQ(rank, rs.rank(j));
        bool found = (j >= TypeParam::bit_pat_len-1);
        for (uint8_t k=0; found and k < TypeParam::bit_pat_len; ++k) {
            found &= bvec[j-k] == ((TypeParam::bit_pat>>k)&1);
        }
        rank += found;
    }
    EXPECT_EQ(rank, rs.rank(bvec.size()));
}

}// end namespace

int main(int argc, char** argv)
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

