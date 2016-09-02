#include "sdsl/bit_vectors.hpp"
#include "sdsl/select_support.hpp"
#include "gtest/gtest.h"
#include <string>

using namespace sdsl;
using namespace std;

string test_file;

namespace
{

template<class T>
class select_support_test : public ::testing::Test { };

using testing::Types;

typedef Types<
select_support_rrri<1, 256>,
select_support_rrri<1, 129>,
select_support_rrri<1, 192>,
select_support_rrri<1, 255>,
select_support_rrri<1, 15>,
select_support_rrri<1, 31>,
select_support_rrri<1, 63>,
select_support_rrri<1, 127>,
select_support_rrri<1, 128>,
select_support_r3d3<1, 256>,
select_support_r3d3<1, 129>,
select_support_r3d3<1, 192>,
select_support_r3d3<1, 255>,
select_support_r3d3<1, 15>,
select_support_r3d3<1, 31>,
select_support_r3d3<1, 63>,
select_support_r3d3<1, 127>,
select_support_r3d3<1, 128>,
select_support_r3d3i<1, 256>,
select_support_r3d3i<1, 129>,
select_support_r3d3i<1, 192>,
select_support_r3d3i<1, 255>,
select_support_r3d3i<1, 15>,
select_support_r3d3i<1, 31>,
select_support_r3d3i<1, 63>,
select_support_r3d3i<1, 127>,
select_support_r3d3i<1, 128>
> Implementations;

TYPED_TEST_CASE(select_support_test, Implementations);

//! Test the select method
TYPED_TEST(select_support_test, select_method)
{
    static_assert(sdsl::util::is_regular<TypeParam>::value, "Type is not regular");
    bit_vector bvec;
    ASSERT_TRUE(load_from_file(bvec, test_file));
    typename TypeParam::bit_vector_type bv(bvec);
    TypeParam ss(&bv);
    for (uint64_t j=0, select=0; j < bvec.size(); ++j) {
        bool found = (j >= TypeParam::bit_pat_len-1);
        for (uint8_t k=0; found and k < TypeParam::bit_pat_len; ++k) {
            found &= bvec[j-k] == ((TypeParam::bit_pat>>k)&1);
        }
        if (found) {
            ++select;
            ASSERT_EQ(j, ss.select(select));
        }
//        if (bvec[j] == TypeParam::bit_pat) {
//            ++select;
//            ASSERT_EQ(j, ss.select(select));
//        }
    }
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
