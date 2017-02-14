/* sdsl - succinct data structures library
   Copyright (C) 2011-2013 Simon Gog

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file r3d3i_vector.hpp
  \brief r3d3i_vector.hpp contains the sdsl::r3d3i_vector class, and
  classes which support rank and select for r3d3i_vector.
  \author Mate Nagy
*/
#ifndef INCLUDED_SDSL_R3D3I_VECTOR
#define INCLUDED_SDSL_R3D3I_VECTOR

#include "int_vector.hpp"
#include "util.hpp"
#include "r3d3_helper.hpp" // for r3d3 helper class
#include "iterators.hpp"
#include <vector>
#include <algorithm> // for next_permutation
#include <iostream>

//! Namespace for the succinct data structure library
namespace sdsl
{

  // forward declaration needed for friend declaration
  template<uint8_t t_b=1, uint16_t t_bs=15, class t_rac=int_vector<>, uint16_t t_k=32>
  class rank_support_r3d3i;                // in r3d3i_vector

  // forward declaration needed for friend declaration
  template<uint8_t t_b=1, uint16_t t_bs=15, class t_rac=int_vector<>, uint16_t t_k=32>
  class select_support_r3d3i;             // in r3d3i_vector

  //! An indexed  \f$Elias-Fanof$-compressed bitvector representation.
  /*!
   *   \tparam t_bs   Size of a basic block.
   *   \tparam t_rac  Random access integer vector. Use to store the block types.
   *                  It is possible to use WTs for t_rac.
   *   \tparam t_k    A rank sample value is stored before every t_k-th basic block.
   *
   *   References:
   *    - Rasmus Pagh
   *      Low redundancy in dictionaries with O(1) worst case lookup time
   *      Technical Report 1998.
   *      ftp://ftp.cs.au.dk/BRICS/Reports/RS/98/28/BRICS-RS-98-28.pdf, Section 2.
   *    - Rajeev Raman, V. Raman and S. Srinivasa Rao
   *      Succinct Indexable Dictionaries with Applications to representations
   *      of k-ary trees and multi-sets.
   *      SODA 2002.
   *    - On the fly-decoding and encoding was discovered in;
   *      Gonzalo Navarro, Eliana Providel:
   *      Fast, Small, Simple Rank/Select on Bitmaps.
   *      SEA 2012
   *    - INFOCOM paper
   *
   *    In this version the block size can be adjust by the template parameter t_bs!
   */
  template<uint16_t t_bs=63, class t_rac=int_vector<>, uint16_t t_k=32>
  class r3d3i_vector
  {
    static_assert(t_bs >= 3 and t_bs <= 256 , "r3d3i_vector_ef: block size t_bs must be 3 <= t_bs <= 256.");
    static_assert(t_k > 1, "r3d3i_vector_ef: t_k must be > 0.");
  public:
    typedef bit_vector::size_type                       size_type;
    typedef bit_vector::value_type                      value_type;
    typedef bit_vector::difference_type                 difference_type;
    typedef t_rac                                       rac_type;
    typedef random_access_const_iterator<r3d3i_vector> iterator;
    typedef bv_tag                                      index_category;

    typedef rank_support_r3d3i<1, t_bs, t_rac, t_k>   rank_1_type;
    typedef rank_support_r3d3i<0, t_bs, t_rac, t_k>   rank_0_type;
    typedef select_support_r3d3i<1, t_bs, t_rac, t_k> select_1_type;
    typedef select_support_r3d3i<0, t_bs, t_rac, t_k> select_0_type;

    friend class rank_support_r3d3i<0, t_bs, t_rac, t_k>;
    friend class rank_support_r3d3i<1, t_bs, t_rac, t_k>;
    friend class select_support_r3d3i<0, t_bs, t_rac, t_k>;
    friend class select_support_r3d3i<1, t_bs, t_rac, t_k>;

    typedef r3d3_helper<t_bs> r3d3_helper_type;
    typedef typename r3d3_helper_type::number_type number_type;

    enum { block_size = t_bs };
  private:
    size_type    m_size = 0;  // Size of the original bit_vector.
    rac_type     m_bt;     // Vector for the block types (bt). bt equals the
    // number of set bits in the block.
    bit_vector   m_btnr;   // Compressed block type numbers.
    int_vector<> m_sbtnrp;  // Sample pointers for superblocks
    int_vector<> m_bbtnrp;  // Sample (relative) pointers for each block
    int_vector<> m_srank;   // Sample rank values (for superblocks)
    int_vector<> m_brank;   // Sample rank values (for blocks-relative)
    bit_vector   m_invert; // Specifies if a superblock (i.e. t_k blocks)
    // have to be considered as inverted i.e. 1 and
    // 0 are swapped
    uint32_t number_type_size;

    void copy(const r3d3i_vector& r3d3i)
    {
      m_size = r3d3i.m_size;
      m_bt = r3d3i.m_bt;
      m_btnr = r3d3i.m_btnr;
      m_sbtnrp = r3d3i.m_sbtnrp;
      m_bbtnrp = r3d3i.m_bbtnrp;
      m_srank = r3d3i.m_srank;
      m_brank = r3d3i.m_brank;
      m_invert = r3d3i.m_invert;
    }

  public:
    const rac_type& bt     = m_bt;
    const bit_vector& btnr = m_btnr;

    //! Default constructor
    r3d3i_vector() {};

    //! Copy constructor
    r3d3i_vector(const r3d3i_vector& r3d3i)
    {
      copy(r3d3i);
    }

    //! Move constructor
    r3d3i_vector(r3d3i_vector&& r3d3i) : m_size(std::move(r3d3i.m_size)),
                                         m_bt(std::move(r3d3i.m_bt)),
                                         m_btnr(std::move(r3d3i.m_btnr)),
                                         m_sbtnrp(std::move(r3d3i.m_sbtnrp)),
                                         m_bbtnrp(std::move(r3d3i.m_bbtnrp)),
                                         m_srank(std::move(r3d3i.m_srank)),
                                         m_brank(std::move(r3d3i.m_brank)),
                                         m_invert(std::move(r3d3i.m_invert)) {}

    //! Constructor
    /*!
     *  \param bv  Uncompressed bitvector.
     *  \param k   Store rank samples and pointers each k-th blocks.
     */
    r3d3i_vector(const bit_vector& bv)
    {
      m_size = bv.size();
      number_type_size = 8*sizeof(number_type);
      int_vector<> bt_array;
      bt_array.width(bits::hi(t_bs)+1);
      bt_array.resize((m_size+t_bs)/((size_type)t_bs)); // blocks for the bt_array + a dummy block at the end,
      // if m_size%t_bs == 0

      // (1) calculate the block types and store them in m_bt
      size_type pos = 0, i = 0, x = 0;
      size_type btnr_pos = 0;
      size_type sum_rank = 0;

      // Going through all the superblocks and decide whether will be inverted or not
      // This precalculates m_invert. We need this information for being able to
      // store inverted offset in EF format.
      size_type gt_half_t_bs = 0; //greater_than_half_of_the_bt
      size_type superblock_width = t_k*t_bs;
      size_type pos_s = superblock_width, pos_b = 0, s_id = 0;
      size_type num_of_superblocks = m_size/superblock_width+1;
      m_invert = bit_vector((bt_array.size()+t_k-1)/t_k, 0);
      while (pos_s + superblock_width <= m_size){
        gt_half_t_bs = 0;
        while (pos_b + t_bs <= pos_s){
          x = r3d3_helper_type::get_bt(bv, pos_b, t_bs);
          if (x > t_bs/2) {gt_half_t_bs++;}
          pos_b += t_bs;
        }
        //set invert bit (and below invert bt value)
        if (gt_half_t_bs > t_k/2){
          m_invert[s_id] = 1;
        }
        s_id++;
        pos_s += superblock_width;
      }

      //btnr_pos stores the length of a block
      uint16_t msb_bit_pos = 0;
      while (pos + t_bs <= m_size) { // handle all blocks full blocks
        number_type bin;
        bt_array[ i++ ] = x = r3d3_helper_type::get_bt(bv, pos, t_bs);
        bin = r3d3_helper_type::get_int(bv, pos, t_bs); //bin is t_bs long slice of input
        sum_rank += x;

        size_type s_id = pos / superblock_width;
        if (m_invert[s_id]){
          bin = ~bin;
          if (number_type_size > t_bs){
            number_type mask = ((number_type)1<<t_bs)-(number_type)1;//t_bs long mask
            bin = bin & mask;
          }
          msb_bit_pos = r3d3_helper_type::hi(bin);
#ifndef R3D3_C1_NO_OPT
          // When C=1 (we have exactly one piece of 0 in
          // the inverted block) we do not use Elias-Fano
          if (t_bs-x == 1){
            btnr_pos += r3d3_helper_type::space_for_class_1(msb_bit_pos);
          }
          else {
            btnr_pos += r3d3_helper_type::space_for_bt_i(msb_bit_pos, t_bs-x);
          }
#else
          btnr_pos += r3d3_helper_type::space_for_bt_i(msb_bit_pos, t_bs-x);
#endif
        }
        else{
          msb_bit_pos = r3d3_helper_type::hi(bin);
#ifndef R3D3_C1_NO_OPT
          // When C=1 we do not use Elias-Fano
          if (x == 1){
            btnr_pos += r3d3_helper_type::space_for_class_1(msb_bit_pos);
          }
          else {
            btnr_pos += r3d3_helper_type::space_for_bt_i(msb_bit_pos, x);
          }
#else
          btnr_pos += r3d3_helper_type::space_for_bt_i(msb_bit_pos, x);
#endif
        }
        pos += t_bs;
      }
      if (pos < m_size) { // handle last not full block
        number_type bin;
        bt_array[ i++ ] = x = r3d3_helper_type::get_bt(bv, pos, m_size - pos);
        bin = r3d3_helper_type::get_int(bv, pos, m_size-pos); //bin is t_bs long slice of input
        msb_bit_pos = r3d3_helper_type::hi(bin);
        sum_rank += x;
#ifndef R3D3_C1_NO_OPT
        // When C=1 we do not use Elias-Fano
        if (x == 1){
          btnr_pos += r3d3_helper_type::space_for_class_1(msb_bit_pos);
        }
        else {
          btnr_pos += r3d3_helper_type::space_for_bt_i(msb_bit_pos, x); //last block is not inverted
        }
#else
        btnr_pos += r3d3_helper_type::space_for_bt_i(msb_bit_pos, x); //last block is not inverted
#endif
      }

      //use btnr_pos to allocate memory for below structures
      m_btnr  = bit_vector(std::max(btnr_pos, (size_type)64), 0);      // max necessary for case: t_bs == 1
      m_sbtnrp = int_vector<>((bt_array.size()+t_k-1)/t_k, 0,  bits::hi(btnr_pos)+1);

      //m_bbtnr length: #of blocks, width: log2(max len of a compressed offset)=2*t_bs
      m_bbtnrp = int_vector<>(bt_array.size(), 0,  bits::hi(t_k*2*t_bs)+1);
      m_srank  = int_vector<>((bt_array.size()+t_k-1)/t_k + ((m_size % (t_k*t_bs))>0), 0, bits::hi(sum_rank)+1);

      //m_brank stores the relative cumulative rank values in a superblock
      //if superblock is inverted, still m_brank stores the original rank values
      m_brank  = int_vector<>(bt_array.size(), 0, bits::hi(t_k*t_bs)+1);
      //                                                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      //   only add a finishing block, if the last block of the superblock is not a dummy block

      // Saving the size of the offset (m_btnr) for error handling (see below)
      size_type btnr_pos_final = btnr_pos;

      // (2) calculate block type numbers and pointers into btnr and rank samples
      pos = 0; i = 0;
      btnr_pos= 0, sum_rank = 0; pos_b=0;
      bool invert = false;
      uint16_t rank_b = 0;
      while (pos + t_bs <= m_size) {  // handle all full blocks
        if ((i % t_k) == (size_type)0) {
          m_sbtnrp[ i/t_k ] = btnr_pos;
          m_srank[ i/t_k ] = sum_rank;

          // calculate invert bt for that superblock
          m_invert[i/t_k] == 1 ? invert = true : invert = false;
          if (i+t_k <= bt_array.size()) {
            if (invert){
              for (size_type j=i; j < i+t_k; ++j) {
                bt_array[j] = t_bs - bt_array[j];
              }
            }
          }
        }
        //obtaining offset's length
        number_type bin = r3d3_helper_type::get_int(bv, pos, t_bs);
        //do block inversion
        if (invert){
          bin = ~bin;
          if (number_type_size > t_bs){
            number_type mask = ((number_type)1<<t_bs)-(number_type)1;//t_bs long mask
            bin = bin & mask; //cutting off the invalid part of the inverted block
          }//e.g.129 is stored as uint256_t (256-129 part are rubbish after inversion)
        }
        msb_bit_pos = r3d3_helper_type::hi(bin);
        uint16_t space_for_bt;
        x = bt_array[i];
#ifndef R3D3_C1_NO_OPT
        // space_for_bt depends on the coding (EF/straight)
        if (x == 1){
          space_for_bt = r3d3_helper_type::space_for_class_1(msb_bit_pos);
        } else {
          space_for_bt = r3d3_helper_type::space_for_bt_i(msb_bit_pos, x);
        }
#else
        space_for_bt = r3d3_helper_type::space_for_bt_i(msb_bit_pos, x);
#endif

        sum_rank += (invert ? (t_bs - x) : x);

        //calculating 'l' value for this block's EF compression
        uint16_t l = (x > 0 ? log2(msb_bit_pos/x) : 0);

        if (space_for_bt) { //filling m_btnr (offset container) - compress block
          // Error handling
          if (btnr_pos >= btnr_pos_final) {
            std::cout<<"Error: trying to overwrite m_btnr!"<<std::endl;
            std::cout<<" btnr_pos value: "<<btnr_pos<<" allocated space for m_btnr: "
                     <<btnr_pos_final<<std::endl;
            std::cout<<"Exit now!"<<std::endl;
            exit(-1);
          }
#ifndef R3D3_C1_NO_OPT
          // When C=1, use straight-coding
          if (x == 1){
            r3d3_helper_type::trait::set_int(m_btnr, btnr_pos, msb_bit_pos, space_for_bt);
          } else {
            r3d3_helper_type::compress_ef(m_btnr, btnr_pos, bin, l);
          }
#else
          r3d3_helper_type::compress_ef(m_btnr, btnr_pos, bin, l);
#endif
        }
        btnr_pos += space_for_bt;

        //if at superblock border
        if ((i % t_k) == (size_type)0) {
          pos_b = space_for_bt;
          rank_b = (invert ? (t_bs - x) : x);
        }
        else {
          pos_b += space_for_bt;
          rank_b += (invert ? (t_bs - x) : x);
        }

        //if no superblock border, save relative position
        //and relative rank
        m_bbtnrp[i] = pos_b;
        m_brank[i]  = rank_b;

        pos += t_bs;
        i++;
      }
      if (pos < m_size) { // handle last not full block
        if ((i % t_k) == (size_type)0) {
          m_sbtnrp[ i/t_k ] = btnr_pos;
          m_srank[ i/t_k ] = sum_rank;
          m_invert[ i/t_k ] = 0; // default: set last block to not inverted

          invert = false;
        }

        number_type bin = r3d3_helper_type::get_int(bv, pos, m_size-pos);
        //if last block is in an inverted superblock
        if (invert){
          bin = ~bin; //invert block
          if (number_type_size > t_bs){ //filter the rubbish
            number_type mask = ((number_type)1<<t_bs)-(number_type)1;//t_bs long mask
            bin = bin & mask;
          }
        }
        msb_bit_pos = r3d3_helper_type::hi(bin);
        uint16_t space_for_bt;
#ifndef R3D3_C1_NO_OPT
        if (x == 1){
          space_for_bt = r3d3_helper_type::space_for_class_1(msb_bit_pos);
        }
        else {
          space_for_bt = r3d3_helper_type::space_for_bt_i(msb_bit_pos, x = bt_array[i]);
        }
#else
        space_for_bt = r3d3_helper_type::space_for_bt_i(msb_bit_pos, x=bt_array[i]);
#endif

        // no extra dummy block added to bt_array, therefore this condition should hold
        assert(i+1 == bt_array.size());
        sum_rank += invert ? (t_bs - x) : x;

        //calculating 'l' value for this block's EF compression
        uint16_t l = (x > 0 ? log2(msb_bit_pos/x) : 0);
        if (space_for_bt) {
#ifndef R3D3_C1_NO_OPT
          // When C=1, use straight-coding
          if (x == 1){
            r3d3_helper_type::trait::set_int(m_btnr, btnr_pos, msb_bit_pos, space_for_bt);
          } else {
            r3d3_helper_type::compress_ef(m_btnr, btnr_pos, bin, l);
          }
#else
          r3d3_helper_type::compress_ef(m_btnr, btnr_pos, bin, l);
#endif
        }
        btnr_pos += space_for_bt;

        if ((i % t_k) == (size_type)0) {
          pos_b = space_for_bt;
          rank_b = (invert ? (t_bs - x) : x);
        }
        else {
          pos_b += space_for_bt;
          rank_b += (invert ? (t_bs - x) : x);
        }
        m_bbtnrp[i] = pos_b;
        m_brank[i]  = rank_b;
        i++;

        //assert(m_srank.size()-1 == ((i+1+t_k-1)/t_k));//mate
      } else { // handle last empty full block
        assert(m_srank.size()-1 == ((i+t_k-1)/t_k));
      }
      // for technical reasons we add a last element to m_rank
      m_srank[ m_srank.size()-1 ] = sum_rank; // sum_rank contains the total number of set bits in bv
      m_bt = bt_array;
    }

    //! Swap method
    void swap(r3d3i_vector& r3d3i)
    {
      if (this != &r3d3i) {
        std::swap(m_size, r3d3i.m_size);
        m_bt.swap(r3d3i.m_bt);
        m_btnr.swap(r3d3i.m_btnr);
        m_bbtnrp.swap(r3d3i.m_bbtnrp);
        m_sbtnrp.swap(r3d3i.m_sbtnrp);
        m_srank.swap(r3d3i.m_srank);
        m_brank.swap(r3d3i.m_brank);
        m_invert.swap(r3d3i.m_invert);
      }
    }

    void printUncompressedBlock(const bit_vector& bv, const size_type &pos){
      number_type bin = r3d3_helper_type::get_int(bv, pos, t_bs);
      std::cout<<" Block: "<<bin<<std::endl;
    }

    void printCompressedBlock(const size_type &id){
      size_type bt_idx = id; //block id
      uint16_t bt = m_bt[bt_idx];
      size_type s_id = bt_idx/t_k; //superblock's id
      size_type btnrp = m_sbtnrp[ s_id ]; //pos until superblock
      btnrp += m_bbtnrp[bt_idx]; //relative index

      //length of the block
      uint16_t next_btnrp = m_bbtnrp[bt_idx+1];
      uint16_t btnrlen = abs(next_btnrp - btnrp);

      std::cout<<" Compr Block "<<id<<" : "<<r3d3_helper_type::get_int(m_btnr, btnrp, btnrlen)<<std::endl;
    }

    void printSizes(){
      std::cout<<"  R3D3I structure element [Mbyte]:"<<std::endl;
      std::cout <<"   m_bt     : "<<size_in_mega_bytes(m_bt)<<std::endl;
      std::cout <<"   m_btnr   : "<<size_in_mega_bytes(m_btnr)<<std::endl;
      std::cout <<"   m_sbtnrp : "<<size_in_mega_bytes(m_sbtnrp)<<std::endl;
      std::cout <<"   m_bbtnrp : "<<size_in_mega_bytes(m_bbtnrp)<<std::endl;
      std::cout <<"   m_srank   : "<<size_in_mega_bytes(m_srank)<<std::endl;
      std::cout <<"   m_brank   : "<<size_in_mega_bytes(m_brank)<<std::endl;
      std::cout <<"   m_invert : "<<size_in_mega_bytes(m_invert)<<std::endl;
    }

    //! Accessing the i-th element of the original bit_vector
    value_type operator[](size_type i)const
    {
      size_type bt_idx = i/t_bs; //start indexing with 0
      uint16_t bt = m_bt[bt_idx];
      size_type s_id = bt_idx/t_k; //superblock's id

#ifndef RRR_NO_OPT
      if (bt == 0 or bt == t_bs) { // very effective optimization
        if (m_invert[s_id]) bt = t_bs - bt;
        return bt>0;
      }
#endif
      uint16_t off = i % t_bs; //i - bt_idx*t_bs;

      size_type btnrp = m_sbtnrp[ s_id ];

      uint16_t btnrlen;
      //if at superblock border, first block's m_bbtnrp holds its length
      if (bt_idx % t_k == 0) {
        btnrlen = m_bbtnrp[bt_idx];
      }
      else { //if not at superblock border set btnrlen and step btnrp
        btnrlen = abs(m_bbtnrp[bt_idx]-m_bbtnrp[bt_idx-1]);
        btnrp += m_bbtnrp[bt_idx-1];
      }

      //obtaining l by btnrlen and class
      uint16_t l = (bt > 0 ? r3d3_helper_type::size_of_l(btnrlen, bt) : 0);

      bool is_bit_set = r3d3_helper_type::decode_bit_ef_i(m_btnr, btnrp, btnrlen, l, bt, off);
      if (m_invert[s_id]) is_bit_set = !is_bit_set;
      return is_bit_set;
    }

    //! Get the integer value of the binary string of length len starting at position idx.
    /*! \param idx Starting index of the binary representation of the integer.
     *  \param len Length of the binary representation of the integer. Default value is 64.
     *   \returns The integer value of the binary string of length len starting at position idx.
     *
     *  \pre idx+len-1 in [0..size()-1]
     *  \pre len in [1..64]
     */
    //uint64_t get_int(size_type idx, uint8_t len=64)const {Todo-1}


    //! Assignment operator
    r3d3i_vector& operator=(const r3d3i_vector& r3d3i)
    {
      if (this != &r3d3i) {
        copy(r3d3i);
      }
      return *this;
    }

    //! Move assignment operator
    r3d3i_vector& operator=(r3d3i_vector&& r3d3i)
    {
      swap(r3d3i);
      return *this;
    }

    //! Returns the size of the original bit vector.
    size_type size()const
    {
      return m_size;
    }

    //! Answers select queries
    //! Serializes the data structure into the given ostream
    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
    {
      structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
      size_type written_bytes = 0;
      written_bytes += write_member(m_size, out, child, "size");
      written_bytes += m_bt.serialize(out, child, "bt");
      written_bytes += m_btnr.serialize(out, child, "btnr");
      written_bytes += m_sbtnrp.serialize(out, child, "sbtnrp");
      written_bytes += m_bbtnrp.serialize(out, child, "bbtnrp");
      written_bytes += m_srank.serialize(out, child, "superblock rank_samples");
      written_bytes += m_brank.serialize(out, child, "block rank_samples");
      written_bytes += m_invert.serialize(out, child, "invert");
      structure_tree::add_size(child, written_bytes);
      return written_bytes;
    }

    //! Loads the data structure from the given istream.
    void load(std::istream& in)
    {
      read_member(m_size, in);
      m_bt.load(in);
      m_btnr.load(in);
      m_sbtnrp.load(in);
      m_srank.load(in);
      m_brank.load(in);
      m_invert.load(in);
    }

    iterator begin() const
    {
      return iterator(this, 0);
    }

    iterator end() const
    {
      return iterator(this, size());
    }
  };

  template<uint8_t t_bit_pattern>
  struct rank_support_r3d3i_trait {
    typedef bit_vector::size_type size_type;
    static size_type adjust_rank(size_type r, SDSL_UNUSED size_type n)
    {
      return r;
    }
  };

  template<>
  struct rank_support_r3d3i_trait<0> {
    typedef bit_vector::size_type size_type;
    static size_type adjust_rank(size_type r, size_type n)
    {
      return n - r;
    }
  };

  //! rank_support for the r3d3i_vector class
  /*!
   * \tparam t_b   The bit pattern of size one. (so `0` or `1`)
   * \tparam t_bs  The block size of the corresponding r3d3i_vector
   * \tparam t_rac Type used to store the block type in the corresponding r3d3i_vector.
   *  TODO: Test if the binary search can be speed up by
   *        saving the (n/2)-th rank value in T[0], the (n/4)-th in T[1],
   *        the (3n/4)-th in T[2],... for small number of rank values
   *    is this called hinted binary search???
   *    or is this called
   */
  template<uint8_t t_b, uint16_t t_bs, class t_rac, uint16_t t_k>
  class rank_support_r3d3i
  {
    static_assert(t_b == 1u or t_b == 0u , "rank_support_r3d3i: bit pattern must be `0` or `1`");
  public:
    typedef r3d3i_vector<t_bs, t_rac, t_k> bit_vector_type;
    typedef typename bit_vector_type::size_type size_type;
    typedef typename bit_vector_type::r3d3_helper_type r3d3_helper_type;
    typedef typename r3d3_helper_type::number_type number_type;
    enum { bit_pat = t_b };
    enum { bit_pat_len = (uint8_t)1 };

  private:
    const bit_vector_type* m_v; //!< Pointer to the rank supported r3d3i_vector

  public:
    //! Standard constructor
    /*! \param v Pointer to the r3d3i_vector, which should be supported
     */
    explicit rank_support_r3d3i(const bit_vector_type* v=nullptr)
    {
      set_vector(v);
    }

    //! Answers rank queries
    /*! \param i Argument for the length of the prefix v[0..i-1], with \f$0\leq i \leq size()\f$.
      \returns Number of 1-bits in the prefix [0..i-1] of the original bit_vector.
      \par Time complexity
      \f$ \Order{ sample\_rate of the r3d3i\_vector} \f$
    */
    const size_type rank(size_type i)const
    {
      assert(m_v != nullptr);
      assert(i <= m_v->size()); //start indexing with 0
      size_type bt_idx = i/t_bs;//start indexing with 0
      size_type s_id = bt_idx/t_k;
      size_type btnrp = m_v->m_sbtnrp[ s_id ];
      size_type srank  = m_v->m_srank[ s_id ];
      if (s_id+1 < m_v->m_srank.size()) {
        size_type diff_rank  = m_v->m_srank[ s_id+1 ] - srank;
#ifndef RRR_NO_OPT
        if (diff_rank == (size_type)0) {
          return  rank_support_r3d3i_trait<t_b>::adjust_rank(srank, i);
        } else if (diff_rank == (size_type)t_bs*t_k) {
          size_type adj_rank = srank + i - s_id*t_k*t_bs;
          return  rank_support_r3d3i_trait<t_b>::adjust_rank(adj_rank, i);
        }
#endif
      }
      uint16_t off = i % t_bs;

      if (!off && i == m_v->size() && bt_idx > 0) {   // needed for special case: if i=size() is a multiple of t_bs
        // the access to m_bt would cause a invalid memory access
        if (bt_idx % t_k == 0) {
          return rank_support_r3d3i_trait<t_b>::adjust_rank(srank, i);
        }
        else {//not at superblock border
          return rank_support_r3d3i_trait<t_b>::adjust_rank(srank + m_v->m_brank[bt_idx-1], i);
        }
      }
      uint16_t bt = m_v->m_bt[ bt_idx ];

      uint16_t btnrlen;
      //if at superblock border, first block's m_bbtnrp holds its length
      if (bt_idx % t_k == 0) {
        btnrlen = m_v->m_bbtnrp[bt_idx];
      }
      else { //if not at superblock border set btnrlen and step btnrp, srank
        btnrlen = abs(m_v->m_bbtnrp[bt_idx]-m_v->m_bbtnrp[bt_idx-1]);
        btnrp += m_v->m_bbtnrp[bt_idx-1];
        srank += m_v->m_brank[bt_idx-1];
      }
      uint16_t l = (bt > 0 ? r3d3_helper_type::size_of_l(btnrlen, bt) : 0);

      //returns rank in a block
      uint16_t rank_b = r3d3_helper_type::rank_on_ef_block_i(m_v->m_btnr, btnrp, btnrlen, l, bt, off);
      if (m_v->m_invert[s_id]) rank_b = off - rank_b;

      return rank_support_r3d3i_trait<t_b>::adjust_rank(srank + rank_b, i);
    }

    //! Short hand for rank(i)
    const size_type operator()(size_type i)const
    {
      return rank(i);
    }

    //! Returns the size of the original vector
    const size_type size()const
    {
      return m_v->size();
    }

    //! Set the supported vector.
    void set_vector(const bit_vector_type* v=nullptr)
    {
      m_v = v;
    }

    rank_support_r3d3i& operator=(const rank_support_r3d3i& rs)
    {
      if (this != &rs) {
        set_vector(rs.m_v);
      }
      return *this;
    }

    void swap(rank_support_r3d3i&) { }

    //! Load the data structure from a stream and set the supported vector.
    void load(std::istream&, const bit_vector_type* v=nullptr)
    {
      set_vector(v);
    }

    //! Serializes the data structure into a stream.
    size_type serialize(std::ostream&, structure_tree_node* v=nullptr, std::string name="")const
    {
      structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
      structure_tree::add_size(child, 0);
      return 0;
    }
  };


  //! Select support for the r3d3i_vector class.
  /*
   * \tparam t_b   The bit pattern of size one. (so `0` or `1`)
   * \tparam t_bs  The block size of the corresponding r3d3i_vector
   * \tparam t_rac Type used to store the block type in the corresponding r3d3i_vector.
   *
   * Possible TODO: Add heap which contains the 10 first items of
   * each binary search could increase performance.
   * Experiments on select_support_interleaved showed about
   * 25%.
   */
  template<uint8_t t_b, uint16_t t_bs, class t_rac, uint16_t t_k>
  class select_support_r3d3i
  {
    static_assert(t_b == 1u or t_b == 0u , "select_support_r3d3i: bit pattern must be `0` or `1`");
  public:
    typedef r3d3i_vector<t_bs, t_rac, t_k> bit_vector_type;
    typedef typename bit_vector_type::size_type size_type;
    typedef typename bit_vector_type::r3d3_helper_type r3d3_helper_type;
    typedef typename r3d3_helper_type::number_type number_type;
    enum { bit_pat = t_b };
    enum { bit_pat_len = (uint8_t)1 };

  private:
    const bit_vector_type* m_v; //!< Pointer to the rank supported r3d3i_vector

    size_type select1(size_type i)const
    {
      if (m_v->m_srank[m_v->m_srank.size()-1] < i)
        return size();

      //if ( i <= 0 ) return -1;

      //  (1) binary search for the answer in the rank_samples (superblocks)
      size_type begin_sb=0, end_sb=m_v->m_srank.size()-1; // min included, max excluded
      size_type idx, rank;
      // invariant:  m_rank[end]   >= i
      //             m_rank[begin]  < i
      while (end_sb-begin_sb > 1) {
        idx  = (begin_sb+end_sb) >> 1; // idx in [0..m_rank.size()-1]
        rank = m_v->m_srank[idx];
        if (rank >= i)
          end_sb = idx;
        else { // rank < i
          begin_sb = idx;
        }
      }
      //   (2) linear search within the samples (superblock)
      // Modified to binary search within the samples
      rank = m_v->m_srank[begin_sb]; // now i>rank
      idx = begin_sb * t_k; // initialize idx for select result
      size_type diff_rank  = m_v->m_srank[end_sb] - rank;
#ifndef RRR_NO_OPT
      if (diff_rank == (size_type)t_bs*t_k) {// optimisation for select<1>
        return idx*t_bs + i-rank -1;
      }
#endif
      const bool inv = m_v->m_invert[ begin_sb ];
      size_type btnrp = m_v->m_sbtnrp[ begin_sb ];
      uint16_t bt = 0, btnrlen = 0;
      //begin/end_block of the superblock
      size_type begin_b = begin_sb*t_k;
      size_type end_b;
      //if not last superblock
      if ( end_sb != m_v->m_srank.size()-1 )
        end_b = end_sb*t_k;
      else
        end_b = m_v->m_brank.size()-1;

      size_type rank_b = rank; //holds cummulative rank
      while (end_b-begin_b > 1){
        idx  = (begin_b+end_b) >> 1;
        rank_b = m_v->m_brank[idx] + rank;

        if (rank_b >= i){
          end_b = idx;
        }
        else {
          begin_b = idx;
        }
      }

      // Correction because m_brank[0] != 0
      // but holds the rank of the first block
      // if not at superblock border or at
      // sb border but end_b holds the required block
      if (begin_b % t_k != 0 || i > rank + m_v->m_brank[begin_b]){
        idx = end_b;
        btnrlen = abs(m_v->m_bbtnrp[idx]-m_v->m_bbtnrp[idx-1]);
        btnrp += m_v->m_bbtnrp[idx-1];
      }
      else {
        idx = begin_b;
        btnrlen = m_v->m_bbtnrp[idx];
      }

      bt = m_v->m_bt[idx];
      rank_b = rank + m_v->m_brank[idx]-(inv ? t_bs - bt: bt);

      uint16_t l = (bt > 0 ? r3d3_helper_type::size_of_l(btnrlen, bt) : 0);
      int pos_in_b = r3d3_helper_type::select_on_ef_block(m_v->m_btnr, btnrp, btnrlen, l, bt, i-rank_b, inv);

      return idx * t_bs + pos_in_b;
    }

   /* size_type select0(size_type i)const {Todo} */

  public:
    explicit select_support_r3d3i(const bit_vector_type* v=nullptr)
    {
      set_vector(v);
    }

    //! Answers select queries
    size_type select(size_type i)const
    {
      //return  t_b ? select1(i) : select0(i);
      return  select1(i);
    }

    const size_type operator()(size_type i)const
    {
      return select(i);
    }

    const size_type size()const
    {
      return m_v->size();
    }

    void set_vector(const bit_vector_type* v=nullptr)
    {
      m_v = v;
    }

    select_support_r3d3i& operator=(const select_support_r3d3i& rs)
    {
      if (this != &rs) {
        set_vector(rs.m_v);
      }
      return *this;
    }

    void swap(select_support_r3d3i&) { }

    void load(std::istream&, const bit_vector_type* v=nullptr)
    {
      set_vector(v);
    }

    size_type serialize(std::ostream&, structure_tree_node* v=nullptr, std::string name="")const
    {
      structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
      structure_tree::add_size(child, 0);
      return 0;
    }
  };

}// end namespace sdsl

#endif
