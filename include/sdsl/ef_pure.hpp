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
    size_type    m_popcnt; // popcount of the original bv
    bit_vector   m_lba; //lower bits
    bit_vector   m_uba; //upper bits
    size_type    m_lba_length;
    size_type    m_uba_length;
    size_type    m_lower_width;
    size_type    m_msb_pos;

    void copy(const ef_pure& ef_pure)
    {
      m_size = ef_pure.m_size;
      m_popcnt = ef_pure.m_popcnt;
      m_lba = ef_pure.m_lba;
      m_uba = ef_pure.m_uba;
    }

  public:
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

      assert(m_size > 0);
      // Calculate popcnt and MSB position
      uint32_t msb_pos = m_size;
      m_popcnt = 0;
      for (bvIt = bv.begin(); bvIt != bv.end(); ++bvIt){
	if (*bvIt) {
          m_popcnt++;
	}
      }
      assert(m_popcnt > 0);

      bvIt = bv.end();
      while (*bvIt != 1){
       	  --bvIt;
          --msb_pos;
      }
      m_msb_pos = msb_pos;

      // Calculate l
      m_lower_width = log2(m_msb_pos/m_popcnt);

      // Allocate memory for m_lba, m_uba
      //m_lba = bit_vector(std::max(m_popcnt*m_lower_width, (size_type)64), 0);
      m_lba = bit_vector(m_popcnt*m_lower_width, 0);
      m_uba = bit_vector(m_popcnt+(m_msb_pos/2), 0);

      // Compression:
      uint32_t bit_pos = 0;
      uint32_t cntOnes = 0;
      uint32_t preUpVal = 0; //initially
      uint64_t upPos = 0;
      for (bvIt = bv.begin(); bvIt != bv.end(); ++bvIt){
        if (*bvIt == 1){
          cntOnes++;

          // LBA
          if (m_lower_width){
	    uint32_t maskLBA = (1 << m_lower_width) - 1;
            uint32_t lowVal = bit_pos & maskLBA;
            m_lba.set_int((cntOnes-1)*m_lower_width, lowVal, m_lower_width);
          }

          // UBA
          uint32_t upVal = bit_pos >> m_lower_width;
          uint16_t gap = upVal - preUpVal;
          uint16_t numOfZeros = 0;
          while (gap--){
            numOfZeros++;
          }
          preUpVal = upVal; //save
          uint64_t upValUnary = 1 << numOfZeros;
          m_uba.set_int(upPos, upValUnary, numOfZeros+1); 
          upPos += numOfZeros + 1;
        }
        bit_pos++;
      }
      m_lba_length = m_popcnt*m_lower_width;
      m_uba_length = upPos;
    }

    void printCompressedData(){
      std::cout<<"=== EF pure compression ==="<<std::endl;
      std::cout<<" m_msb_pos     : "<<m_msb_pos<<std::endl;
      std::cout<<" m_popcnt      : "<<m_popcnt<<std::endl;
      std::cout<<" m_lower_width : "<<m_lower_width<<std::endl;
      std::cout<<" m_lba_length  : "<<m_lba_length<<std::endl;
      std::cout<<" m_lba.size()  : "<<m_lba.size()<<std::endl;
      std::cout<<" m_uba_length  : "<<m_uba_length<<std::endl;
      std::cout<<" m_uba.size()   : "<<m_uba.size()<<std::endl;

      bit_vector::const_iterator bvIt;
      std::cout<<std::endl;
      std::cout<<"LBA: "<<std::endl;
      for (bvIt = m_lba.begin(); bvIt != m_lba.end(); ++bvIt){
	if (*bvIt == 1) {
          std::cout<<"1 ";
	}
        else {
	  std::cout<<"0 ";
        }
      }
      std::cout<<std::endl;
      std::cout<<"UBA: "<<std::endl;
      for (bvIt = m_uba.begin(); bvIt != m_uba.end(); ++bvIt){
	if (*bvIt == 1) {
          std::cout<<"1 ";
	}
        else {
	  std::cout<<"0 ";
        }
      }
      std::cout<<std::endl;
      std::cout<<std::endl;
    }

    //! Accessing the i-th element of the original bit_vector
    value_type operator[](size_type i)const
    {
      
      // LBA part
      uint64_t lowVal = m_lba.get_int(i*m_lower_width, m_lower_width);

      // UBA part
      uint64_t upVal;
      uint16_t gap = 1 << m_lower_width;
      uint16_t gap_id = i / gap;
      uint16_t nr_gaps = m_uba_length - m_popcnt;

      // nothing is stored, value is 0
      if ( gap_id > nr_gaps ) return false;

      uint16_t gap_pos = 0;
      if ( gap_id != 0 ){
        select_support_mcl<> sb(m_uba);//special select!!, thats a cheat!
        gap_pos = sb(gap_id) + 1;
      }
      //innen todo!!!

      return 0;
    }

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
