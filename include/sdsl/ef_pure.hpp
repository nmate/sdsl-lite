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
/*! \file ef_pure.hpp
  \brief ef_pure.hpp contains the sdsl::ef_pure class, and
  classes which support rank and select for ef_pure.
  \author Mate Nagy
*/
#ifndef INCLUDED_SDSL_EF_PURE
#define INCLUDED_SDSL_EF_PURE

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
  class rank_support_ef_pure;               

  // forward declaration needed for friend declaration
  template<uint8_t t_b=1, uint16_t t_bs=15, class t_rac=int_vector<>, uint16_t t_k=32>
  class select_support_ef_pure;

  //! A \f$Elias-Fanof$ compressed bitvector representation.
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
   *    \sa sdsl::ef_pure for a specialized version for block_size=15
   */
  template<class t_rac=int_vector<> >
  class ef_pure
  {
  public:
    typedef bit_vector::size_type                       size_type;
    typedef bit_vector::value_type                      value_type;
    typedef t_rac                                       rac_type;
    typedef random_access_const_iterator<ef_pure> iterator;

    //typedef rank_support_ef_pure<1, t_rac>   rank_1_type;
    //typedef rank_support_ef_pure<0, t_rac>   rank_0_type;
    //typedef select_support_ef_pure<1, t_rac> select_1_type;
    //typedef select_support_ef_pure<0, t_rac> select_0_type;
    // 
    //friend class rank_support_ef_pure<0, t_rac>;
    //friend class rank_support_ef_pure<1, t_rac>;
    //friend class select_support_ef_pure<0, t_rac>;
    //friend class select_support_ef_pure<1, t_rac>;

    //could use ef_pure_helper...
    //typedef r3d3_helper<t_bs> r3d3_helper_type;
    //typedef typename r3d3_helper_type::number_type number_type;

  private:
    size_type    m_size = 0;  // Size of the original bit_vector.
    //rac_type     m_bt;     // Vector for the block types (bt). bt equals the
    // number of set bits in the block.
    size_type    m_popcnt; // popcount of the original bv
    bit_vector   m_lba; //lower bits
    bit_vector   m_uba; //upper bits

    //int_vector<> m_btnrp;  // Sample pointers into m_btnr.
    //int_vector<> m_rank;   // Sample rank values.
    //bit_vector   m_invert; // Specifies if a superblock (i.e. t_k blocks)
    // have to be considered as inverted i.e. 1 and
    // 0 are swapped
    //uint32_t number_type_size;

    void copy(const ef_pure& ef_pure)
    {
      m_size = ef_pure.m_size;
      m_popcnt = ef_pure.m_popcnt;
      m_lba = ef_pure.m_lba;
      m_uba = ef_pure.m_uba;
    }

  public:
    const size_type& popcnt = m_popcnt;
    const bit_vector& lba   = m_lba;
    const bit_vector& uba   = m_uba;

    //! Default constructor
    ef_pure() {};

    //! Copy constructor
    ef_pure(const ef_pure& ef_pure)
    {
      copy(ef_pure);
    }

    //! Move constructor
    ef_pure(ef_pure&& ef_pure) : m_size(std::move(ef_pure.m_size)),
				 m_popcnt(std::move(ef_pure.m_popcnt)),
				 m_uba(std::move(ef_pure.m_uba)), 
                                 m_lba(std::move(ef_pure.m_lba)) {}

    //! Constructor
    /*!
     *  \param bv  Uncompressed bitvector.
     *  \param k   Store rank samples and pointers each k-th blocks.
     */
    ef_pure(const bit_vector& bv)
    {
      m_size = bv.size();
      bit_vector::const_iterator bvIt;

      // Calculate popcnt and MSB position
      int msb_pos = m_size;
      int popcnt = 0;
      for (bvIt = bv.begin(); bvIt != bv.end(); ++bvIt){
	if (*bvIt) {
         popcnt++;
	}
      }

      bvIt = bv.end();
      while (*bvIt != 1){
       	  --bvIt;
          --msb_pos;
      }
      std::cout<<"msb_pos: "<<msb_pos<<std::endl;
      std::cout<<"popcnt : "<<popcnt<<std::endl;

      // Calculate popcount
      //m_popcnt = Todo

      //msb bit of bv

      //calculate l

    }

    //! Swap method
    //void swap(r3d3_vector& r3d3)
    //{
    //  if (this != &r3d3) {
    //    std::swap(m_size, r3d3.m_size);
    //    m_bt.swap(r3d3.m_bt);
    //    m_btnr.swap(r3d3.m_btnr);
    //    m_btnrp.swap(r3d3.m_btnrp);
    //    m_rank.swap(r3d3.m_rank);
    //    m_invert.swap(r3d3.m_invert);
    //  }
    //}
    // 
    //void printCompressedBlock(const size_type &id){
    //  size_type bt_idx = id; //block id
    //  uint16_t bt = m_bt[bt_idx];
    //  size_type s_id = bt_idx/t_k; //superblock's id
    //  size_type btnrp = m_btnrp[ s_id ]; //pos until superblock
    // 
    //  for (size_type j = s_id*t_k; j < bt_idx; ++j) {
    //    btnrp += r3d3_helper_type::space_for_bt(m_bt[j]);
    //  }
    //  uint16_t btnrlen = r3d3_helper_type::space_for_bt(bt);
    // 
    //  std::cout<<" bt_id: "<<bt_idx<<" ; bt: "<<bt<<" ; s_id: "<<s_id<<std::endl;
    //  std::cout<<" btnrp: "<<btnrp<<" ; btnrlen: "<<btnrlen<<std::endl;
    //  std::cout<<" Compr Block "<<id<<" : "<<r3d3_helper_type::get_int(m_btnr, btnrp, btnrlen)<<std::endl;
    //}

    //! Accessing the i-th element of the original bit_vector
    //value_type operator[](size_type i)const
    //{
    //  size_type bt_idx = i/t_bs; //start indexing with 0
    //  uint16_t bt = m_bt[bt_idx];
    //  size_type s_id = bt_idx/t_k; //superblock's id
    // 
    //#ifndef RRR_NO_OPT
    //  if (bt == 0 or bt == t_bs) { // very effective optimization
    //    if (m_invert[s_id]) bt = t_bs - bt;
    //    return bt>0;
    //  }
    //#endif
    //  uint16_t off = i % t_bs; //i - bt_idx*t_bs;
    // 
    //  uint16_t l = (bt > 0 ? log2(t_bs/bt) : 0);
    // 
    //  size_type btnrp = m_btnrp[ s_id ];
    //  for (size_type j = s_id*t_k; j < bt_idx; ++j) {
    //    btnrp += r3d3_helper_type::space_for_bt(m_bt[j]);
    //  }
    //  uint16_t btnrlen = r3d3_helper_type::space_for_bt(bt);
    // 
    //  // and give it to block decoding function than
    //  bool is_bit_set = r3d3_helper_type::decode_bit_ef(m_btnr, btnrp, btnrlen, l, bt, off);
    // 
    //  if (m_invert[s_id]) is_bit_set = !is_bit_set;
    //  return is_bit_set;
    //}

    /*! \param idx Starting index of the binary representation of the integer.
     *  \param len Length of the binary representation of the integer. Default value is 64.
     *   \returns The integer value of the binary string of length len starting at position idx.
     *
     *  \pre idx+len-1 in [0..size()-1]
     *  \pre len in [1..64]
     *
     *  Todo-1: this should be implemented for R3D3 for testing purposes
     *          Anyhow bit exactness is already ensured, so this seems to be an overkill
     * 
     */
    //uint64_t get_int(size_type idx, uint8_t len=64)const { Todo-1 }

    //! Assignment operator
    //r3d3_vector& operator=(const r3d3_vector& r3d3)
    //{
    //  if (this != &r3d3) {
    //    copy(r3d3);
    //  }
    //  return *this;
    //}
    // 
    ////! Move assignment operator
    //r3d3_vector& operator=(r3d3_vector&& r3d3)
    //{
    //  swap(r3d3);
    //  return *this;
    //}

    //! Returns the size of the original bit vector.
    size_type size()const
    {
      return m_size;
    }

    //! Answers select queries
    //! Serializes the data structure into the given ostream
    //size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
    //{
    //  structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    //  size_type written_bytes = 0;
    //  written_bytes += write_member(m_size, out, child, "size");
    //  written_bytes += m_bt.serialize(out, child, "bt");
    //  written_bytes += m_btnr.serialize(out, child, "btnr");
    //  written_bytes += m_btnrp.serialize(out, child, "btnrp");
    //  written_bytes += m_rank.serialize(out, child, "rank_samples");
    //  written_bytes += m_invert.serialize(out, child, "invert");
    //  structure_tree::add_size(child, written_bytes);
    //  return written_bytes;
    //}

    //! Loads the data structure from the given istream.
    //void load(std::istream& in)
    //{
    //  read_member(m_size, in);
    //  m_bt.load(in);
    //  m_btnr.load(in);
    //  m_btnrp.load(in);
    //  m_rank.load(in);
    //  m_invert.load(in);
    //}

    //ezmi?
    //iterator begin() const
    //{
    //  return iterator(this, 0);
    //}
    // 
    //iterator end() const
    //{
    //  return iterator(this, size());
    //}
  };
  //todo put back rank and select support stuff

}// end namespace sdsl

#endif
