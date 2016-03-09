#include "sdsl/suffix_arrays.hpp"
#include "sdsl/coder.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <string>

namespace
{

using namespace sdsl;
using namespace std;

typedef int_vector<>::size_type size_type;

string test_file;
string temp_file;
string temp_dir;
bool in_memory;

class csa_byte_test_base : public ::testing::Test
{
protected:
    static cache_config s_config;

public:
    static void SetupTestCase()
    {
    }
    
    static void TearDownTestCase()
    {
        // csa_rao handles multiple zero symbols but other classes may not.
        util::delete_all_files(s_config.file_map);
    }
};

template<class T>
class csa_byte_test : public csa_byte_test_base
{
protected:
    using csa_byte_test_base::s_config;
};

cache_config csa_byte_test_base::s_config = {};

using testing::Types;

typedef Types<
       csa_wt<>,
       csa_sada<>,
       csa_rao<>,
       csa_sada<enc_vector<coder::fibonacci>>,
       csa_sada<enc_vector<coder::elias_gamma>>,
       csa_rao<csa_rao_spec<2, 2>>,
       csa_rao<csa_rao_spec<2, 4>>,
       csa_rao<csa_rao_spec<4, 2>>,
       csa_rao<csa_rao_spec<5, 2>>,
       csa_rao<csa_rao_spec<4, 4>>,
       csa_wt<wt_huff<>, 8, 16, text_order_sa_sampling<>>,
       csa_wt<wt_huff<>,32,32,fuzzy_sa_sampling<>>,
       csa_wt<wt_huff<>,32,32,fuzzy_sa_sampling<bit_vector, bit_vector>, fuzzy_isa_sampling_support<>>,
       csa_wt<wt_huff<>,32,32,fuzzy_sa_sampling<>, fuzzy_isa_sampling_support<>>,
       csa_wt<wt_huff<>,32,32,text_order_sa_sampling<>,isa_sampling<>>,
       csa_wt<wt_huff<>,32,32,text_order_sa_sampling<>,text_order_isa_sampling_support<>>,
       csa_sada<enc_vector<>, 32,32,text_order_sa_sampling<>,isa_sampling<>>,
       csa_sada<enc_vector<>, 32,32,text_order_sa_sampling<>,text_order_isa_sampling_support<>>,
       csa_wt<wt_huff<>, 8, 16, sa_order_sa_sampling<>>,
       csa_wt<wt_huff<>, 8, 16, sa_order_sa_sampling<>, isa_sampling<>, succinct_byte_alphabet<bit_vector, rank_support_v<>, select_support_mcl<>>>,
       csa_wt<wt_huff<>, 8, 16, sa_order_sa_sampling<>, isa_sampling<>, succinct_byte_alphabet<>>,
       csa_bitcompressed<>
       > Implementations;

TYPED_TEST_CASE(csa_byte_test, Implementations);

TYPED_TEST(csa_byte_test, create_and_store_test)
{
    static_assert(sdsl::util::is_regular<TypeParam>::value, "Type is not regular");
    TypeParam csa;
    csa_byte_test_base::s_config = cache_config(false, temp_dir, util::basename(test_file));
    construct(csa, test_file, csa_byte_test_base::s_config, 1);
    ASSERT_TRUE(store_to_file(csa, temp_file));
}

//! Test backward_search
TYPED_TEST(csa_byte_test, backward_search)
{
    TypeParam csa;
    ASSERT_TRUE(load_from_file(csa, temp_file));
    int_vector<8> text;
    ASSERT_TRUE(load_vector_from_file(text, test_file, 1));
    auto expected_interval_member = csa.psi[0];
    size_type count, l_res, r_res;
    // search for full text
    count = backward_search(csa, 0, csa.size() - 1, text.begin(), text.end(), l_res, r_res);
    ASSERT_EQ((size_type)1, count);
    ASSERT_EQ(l_res, expected_interval_member);
    ASSERT_EQ(r_res, expected_interval_member);
    // search for short phrase
    text.resize(min((int_vector<8>::size_type)4, text.size()));
    count = backward_search(csa, 0, csa.size() - 1, text.begin(), text.end(), l_res, r_res);
    ASSERT_LE((size_type)1, count);
    ASSERT_LE(l_res, expected_interval_member);
    ASSERT_GE(r_res, expected_interval_member);
    // search for empty phrase
    text.resize(0);
    count = backward_search(csa, 0, csa.size() - 1, text.begin(), text.end(), l_res, r_res);
    ASSERT_EQ(csa.size(), count);
    ASSERT_EQ(l_res, (size_type)0);
    ASSERT_EQ(r_res, (size_type)(csa.size() - 1));
}

//! Test forward_search
TYPED_TEST(csa_byte_test, forward_search)
{
    TypeParam csa;
    ASSERT_TRUE(load_from_file(csa, temp_file));
    int_vector<8> text;
    ASSERT_TRUE(load_vector_from_file(text, test_file, 1));
    auto expected_interval_member = csa.psi[0];
    size_type count, l_res, r_res;
    // search for full text
    count = forward_search(csa, 0, csa.size() - 1, text.begin(), text.end(), l_res, r_res);
    ASSERT_EQ((size_type)1, count);
    ASSERT_EQ(l_res, expected_interval_member);
    ASSERT_EQ(r_res, expected_interval_member);
    // search for short phrase
    text.resize(min((int_vector<8>::size_type)4, text.size()));
    count = forward_search(csa, 0, csa.size() - 1, text.begin(), text.end(), l_res, r_res);
    ASSERT_LE((size_type)1, count);
    ASSERT_LE(l_res, expected_interval_member);
    ASSERT_GE(r_res, expected_interval_member);
    // search for empty phrase
    text.resize(0);
    count = forward_search(csa, 0, csa.size() - 1, text.begin(), text.end(), l_res, r_res);
    ASSERT_EQ(csa.size(), count);
    ASSERT_EQ(l_res, (size_type)0);
    ASSERT_EQ(r_res, (size_type)(csa.size() - 1));
}

//! Test sigma member
TYPED_TEST(csa_byte_test, sigma)
{
    TypeParam csa;
    ASSERT_TRUE(load_from_file(csa, temp_file));
    int_vector<8> text;
    ASSERT_TRUE(load_vector_from_file(text, test_file, 1));
    text.resize(text.size()+1);
    text[text.size()-1] = 0; // add 0-character to the end
    ASSERT_EQ(text.size(), csa.size());
    bit_vector occur(256, 0);
    uint16_t sigma = 0;
    for (size_type j=0; j<text.size(); ++j) {
        if (!occur[text[j]]) {
            occur[text[j]] = 1;
            ++sigma;
        }
    }
    ASSERT_EQ(sigma, csa.sigma);
}

//! Test suffix array access methods
TYPED_TEST(csa_byte_test, sa_access)
{
    TypeParam csa;
    ASSERT_TRUE(load_from_file(csa, temp_file));
    int_vector<> sa;
    load_from_file(sa, csa_byte_test_base::s_config.file_map[conf::KEY_SA]);
    size_type n = 1 + util::file_size(test_file);
    ASSERT_EQ(n, csa.size());
    // csa_rao may have padding, remove its effect.
    auto const padding(csa.padding());
    for (size_type j=0; j<n; ++j) {
        ASSERT_EQ(sa[j + padding], csa[j])<<" j="<<j;
    }
}

//! Test inverse suffix access methods
TYPED_TEST(csa_byte_test, isa_access)
{
    TypeParam csa;
    ASSERT_TRUE(load_from_file(csa, temp_file));
    int_vector<> isa;
    size_type n = 0;
    {
        int_vector<> sa;
        load_from_file(sa, csa_byte_test_base::s_config.file_map[conf::KEY_SA]);
        n = 1 + util::file_size(test_file);
        ASSERT_EQ(n, csa.size());
        isa = sa;
        // csa_rao may have padding, remove its effect.
        auto const padding(csa.padding());
        for (size_type j=0; j<n; ++j) {
            isa[sa[j + padding]] = j;
        }
    }
    
    for (size_type j=0; j<n; ++j) {
        ASSERT_EQ(isa[j], csa.isa[j])<<" j="<<j;
    }
}

//! Test text access methods
TYPED_TEST(csa_byte_test, text_access)
{
    if (csa_byte_test_base::s_config.file_map.find(conf::KEY_TEXT) != csa_byte_test_base::s_config.file_map.end()) {
        TypeParam csa;
        ASSERT_TRUE(load_from_file(csa, temp_file));
        int_vector<8> text;
        ASSERT_TRUE(load_vector_from_file(text, test_file, 1));
        text.resize(text.size()+1);
        text[text.size()-1] = 0; // add 0-character to the end
        size_type n = text.size();
        ASSERT_EQ(n, csa.size());
        for (size_type j=0; j<n; ++j) {
            ASSERT_EQ(text[j], csa.text[j])<<" j="<<j;
        }
        auto len = std::min(csa.size(),
                            std::max(csa.size()/10, (decltype(csa.size()))20));
        auto ex_text = extract(csa, 0, len-1);
        for (size_type j=0; j<len; ++j) {
            ASSERT_EQ(text[j], (decltype(text[j]))ex_text[j])<<" j="<<j;
        }
    }
}

//! Test Burrows-Wheeler access methods
TYPED_TEST(csa_byte_test, bwt_access)
{
    if (csa_byte_test_base::s_config.file_map.find(conf::KEY_BWT) != csa_byte_test_base::s_config.file_map.end()) {
        TypeParam csa;
        ASSERT_TRUE(load_from_file(csa, temp_file));
        int_vector<8> bwt;
        load_from_file(bwt, csa_byte_test_base::s_config.file_map[conf::KEY_BWT]);
        size_type n = 1 + util::file_size(test_file);
        ASSERT_EQ(n, csa.size());
        // csa_rao may have padding, remove its effect.
        auto const padding(csa.padding());
        for (size_type j=0; j<n; ++j) {
            ASSERT_EQ(bwt[j + padding], csa.bwt[j])<<" j="<<j;
        }
    }
}

TYPED_TEST(csa_byte_test, f_access)
{
    if (csa_byte_test_base::s_config.file_map.find(conf::KEY_TEXT) != csa_byte_test_base::s_config.file_map.end()) {
        TypeParam csa;
        ASSERT_TRUE(load_from_file(csa, temp_file));
        int_vector<8> text;
        ASSERT_TRUE(load_vector_from_file(text, test_file, 1));
        text.resize(text.size()+1);
        text[text.size()-1] = 0; // add 0-character to the end
        std::sort(begin(text),end(text));
        size_type n = text.size();
        ASSERT_EQ(n, csa.size());
        for (size_type j=0; j<n; j+=200) {
            ASSERT_EQ(text[j], csa.F[j])<<" j="<<j;
        }
    }
}


//! Test Psi access methods
TYPED_TEST(csa_byte_test, psi_access)
{
    if (csa_byte_test_base::s_config.file_map.find(conf::KEY_PSI) != csa_byte_test_base::s_config.file_map.end()) {
        TypeParam csa;
        ASSERT_TRUE(load_from_file(csa, temp_file));
        int_vector<> psi;
        load_from_file(psi, csa_byte_test_base::s_config.file_map[conf::KEY_PSI]);
        size_type n = psi.size();
        ASSERT_EQ(n, csa.size());
        for (size_type j=0; j<n; ++j) {
            ASSERT_EQ(psi[j], csa.psi[j])<<" j="<<j;
        }
    }
}

//! Test if Psi[LF[i]]=i
TYPED_TEST(csa_byte_test, psi_lf_access)
{
    TypeParam csa;
    ASSERT_TRUE(load_from_file(csa, temp_file));
    for (size_type j=0; j<csa.size(); ++j) {
        size_type lf = csa.lf[j];
        ASSERT_TRUE(lf < csa.size());
        ASSERT_EQ(j, csa.psi[lf])<<" j="<<j;
    }
}


//! Test access after swap
TYPED_TEST(csa_byte_test, swap)
{
    TypeParam csa1;
    ASSERT_TRUE(load_from_file(csa1, temp_file));
    TypeParam csa2;
    csa1.swap(csa2);
    int_vector<> sa;
    load_from_file(sa, csa_byte_test_base::s_config.file_map[conf::KEY_SA]);
    size_type n = 1 + util::file_size(test_file);
    ASSERT_EQ(n, csa2.size());
    // csa_rao may have padding, remove its effect.
    auto const padding(csa2.padding());
    for (size_type j=0; j<n; ++j) {
        ASSERT_EQ((typename TypeParam::value_type)sa[j + padding], csa2[j]);
    }
}


TYPED_TEST(csa_byte_test, delete_)
{
    sdsl::remove(temp_file);
    util::delete_all_files(csa_byte_test_base::s_config.file_map);
}

}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc < 4) {
        // LCOV_EXCL_START
        cout << "Usage: " << argv[0] << " test_file temp_file tmp_dir [in-memory]" << endl;
        cout << " (1) Generates a CSA out of test_file; stores it in temp_file." << endl;
        cout << "     Temporary files (SA/BWT/TEXT) are stored in tmp_dir." << endl;
        cout << "     If `in-memory` is specified, the in-memory construction is tested." << endl;
        cout << " (2) Performs tests." << endl;
        cout << " (3) Deletes temp_file." << endl;
        return 1;
        // LCOV_EXCL_STOP
    }
    test_file = argv[1];
    temp_file = argv[2];
    temp_dir  = argv[3];
    in_memory    = argc > 4;
    if (in_memory) {
        temp_dir = "@";
        int_vector<8> data;
        load_vector_from_file(data, test_file, 1);
        test_file = ram_file_name(test_file);
        store_to_plain_array<uint8_t>(data, test_file);
        temp_file = ram_file_name(temp_file);
    }
    return RUN_ALL_TESTS();
}
