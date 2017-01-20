#include "sdsl/bit_vectors.hpp" // for ef_pure
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
ef_pure<>
> Implementations;

TYPED_TEST_CASE(bit_vector_test, Implementations);


//! Test operator[]
TYPED_TEST(bit_vector_test, access)
{
  //static_assert(sdsl::util::is_regular<TypeParam>::value, "Type is not regular");
    bit_vector bv;
    ASSERT_TRUE(load_from_file(bv, test_file));
    TypeParam c_bv(bv);
    ASSERT_EQ(bv.size(), c_bv.size());
    uint64_t st = 0;
    for (uint64_t j=0; j < bv.size(); ++j) {
        ASSERT_EQ((bool)(bv[j]), (bool)(c_bv[j]))<<"[mismatch at]: "<<j;
        if (j % 10000 == 0){
          cout<<"Stage: "<<st++<<endl;
        }
    }
    TypeParam mo_bv = TypeParam(bv);
    ASSERT_EQ(bv.size(), mo_bv.size());
    for (uint64_t j=0; j < bv.size(); ++j) {
        ASSERT_EQ((bool)(bv[j]), (bool)(c_bv[j]));
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

