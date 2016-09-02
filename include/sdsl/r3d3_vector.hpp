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
/*! \file r3d3_vector.hpp
  \brief r3d3_vector.hpp contains the sdsl::r3d3_vector class, and
  classes which support rank and select for r3d3_vector.
  \author Mate Nagy
*/
#ifndef INCLUDED_SDSL_R3D3_VECTOR
#define INCLUDED_SDSL_R3D3_VECTOR

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
  class rank_support_r3d3;                // in r3d3_vector

  // forward declaration needed for friend declaration
  template<uint8_t t_b=1, uint16_t t_bs=15, class t_rac=int_vector<>, uint16_t t_k=32>
  class select_support_r3d3;                // in r3d3_vector

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
   *    \sa sdsl::r3d3_vector for a specialized version for block_size=15
   */
  template<uint16_t t_bs=63, class t_rac=int_vector<>, uint16_t t_k=32>
  class r3d3_vector
  {
    static_assert(t_bs >= 3 and t_bs <= 256 , "r3d3_vector: block size t_bs must be 3 <= t_bs <= 256.");
    static_assert(t_k > 1, "r3d3_vector: t_k must be > 0.");
  public:
    typedef bit_vector::size_type                       size_type;
    typedef bit_vector::value_type                      value_type;
    typedef bit_vector::difference_type                 difference_type;
    typedef t_rac                                       rac_type;
    typedef random_access_const_iterator<r3d3_vector> iterator;
    typedef bv_tag                                      index_category;

    typedef rank_support_r3d3<1, t_bs, t_rac, t_k>   rank_1_type;
    typedef rank_support_r3d3<0, t_bs, t_rac, t_k>   rank_0_type;
    typedef select_support_r3d3<1, t_bs, t_rac, t_k> select_1_type;
    typedef select_support_r3d3<0, t_bs, t_rac, t_k> select_0_type;

    friend class rank_support_r3d3<0, t_bs, t_rac, t_k>;
    friend class rank_support_r3d3<1, t_bs, t_rac, t_k>;
    friend class select_support_r3d3<0, t_bs, t_rac, t_k>;
    friend class select_support_r3d3<1, t_bs, t_rac, t_k>;

    typedef r3d3_helper<t_bs> r3d3_helper_type;
    typedef typename r3d3_helper_type::number_type number_type;

    enum { block_size = t_bs };
  private:
    size_type    m_size = 0;  // Size of the original bit_vector.
    rac_type     m_bt;     // Vector for the block types (bt). bt equals the
    // number of set bits in the block.
    bit_vector   m_btnr;   // Compressed block type numbers.
    int_vector<> m_btnrp;  // Sample pointers into m_btnr.
    int_vector<> m_rank;   // Sample rank values.
    bit_vector   m_invert; // Specifies if a superblock (i.e. t_k blocks)
    // have to be considered as inverted i.e. 1 and
    // 0 are swapped
    uint32_t number_type_size;

    void copy(const r3d3_vector& r3d3)
    {
      m_size = r3d3.m_size;
      m_bt = r3d3.m_bt;
      m_btnr = r3d3.m_btnr;
      m_btnrp = r3d3.m_btnrp;
      m_rank = r3d3.m_rank;
      m_invert = r3d3.m_invert;
    }

  public:
    const rac_type& bt     = m_bt;
    const bit_vector& btnr = m_btnr;

    //! Default constructor
    r3d3_vector() {};

    //! Copy constructor
    r3d3_vector(const r3d3_vector& r3d3)
    {
      copy(r3d3);
    }

    //! Move constructor
    r3d3_vector(r3d3_vector&& r3d3) : m_size(std::move(r3d3.m_size)),
				      m_bt(std::move(r3d3.m_bt)),
				      m_btnr(std::move(r3d3.m_btnr)), 
                                      m_btnrp(std::move(r3d3.m_btnrp)),
				      m_rank(std::move(r3d3.m_rank)), 
                                      m_invert(std::move(r3d3.m_invert)) {}

    //! Constructor
    /*!
     *  \param bv  Uncompressed bitvector.
     *  \param k   Store rank samples and pointers each k-th blocks.
     */
    r3d3_vector(const bit_vector& bv)
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
        //set invert bit and invert bt value
        if (gt_half_t_bs > t_k/2){
          m_invert[s_id] = 1;
        }
        s_id++;
        pos_s += superblock_width;
      }

      while (pos + t_bs <= m_size) { // handle all full blocks
        bt_array[ i++ ] = x = r3d3_helper_type::get_bt(bv, pos, t_bs);
        sum_rank += x;

        size_type s_id = pos / superblock_width;
        // If inverted superblock
        if (m_invert[s_id]){
          btnr_pos += r3d3_helper_type::space_for_bt(t_bs-x);
        }
        // No inversion
        else{
          btnr_pos += r3d3_helper_type::space_for_bt(x);
        }
        pos += t_bs;
      }//end full blocks
      if (pos < m_size) { // handle last not full block
        bt_array[ i++ ] = x = r3d3_helper_type::get_bt(bv, pos, m_size - pos);
        sum_rank += x;

        btnr_pos += r3d3_helper_type::space_for_bt(x); //last block is not inverted

      } //end last not full block

      m_btnr  = bit_vector(std::max(btnr_pos, (size_type)64), 0);      // max necessary for case: t_bs == 1
      m_btnrp = int_vector<>((bt_array.size()+t_k-1)/t_k, 0,  bits::hi(btnr_pos)+1);
      m_rank  = int_vector<>((bt_array.size()+t_k-1)/t_k + ((m_size % (t_k*t_bs))>0), 0, bits::hi(sum_rank)+1);
      //                                                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      // (2) calculate block type numbers (class) and pointers into btnr and rank samples
      pos = 0; i = 0;
      btnr_pos= 0, sum_rank = 0;
      bool invert = false;
      while (pos + t_bs <= m_size) {  // handle all full blocks
        if ((i % t_k) == (size_type)0) {
          m_btnrp[ i/t_k ] = btnr_pos;
          m_rank[ i/t_k ] = sum_rank;
          // calculate invert bt for that superblock
          m_invert[i/t_k] == 1 ? invert = true : invert = false;

          if (i+t_k <= bt_array.size() && invert) {
            for (size_type j=i; j < i+t_k; ++j) {
              bt_array[j] = t_bs - bt_array[j];
            }
          }
        }

        uint16_t space_for_bt = r3d3_helper_type::space_for_bt(x = bt_array[i++]);
        sum_rank += (invert ? (t_bs - x) : x);

        //calculating 'l' value for this block's EF compression
        uint16_t l = (x > 0 ? log2(t_bs/x) : 0);
        if (space_for_bt) { //filling m_btnr (offset container)
          number_type bin = r3d3_helper_type::get_int(bv, pos, t_bs);

          if (invert) { //store the negate of the block
            bin = ~bin;

            if (number_type_size > t_bs){
              number_type mask = ((number_type)1<<t_bs)-(number_type)1; //t_bs long mask
              bin = bin & mask; //cutting off the invalid part of the inverted block
              //e.g.129 is stored as uint256_t (256-129 part are rubbish after inversion)
            }
          }
#ifndef R3D3_C1_NO_OPT
          if (x == 1){
            uint16_t msb_bit_pos = r3d3_helper_type::hi(bin);
            r3d3_helper_type::trait::set_int(m_btnr, btnr_pos, msb_bit_pos, space_for_bt);
          }
          else {
            r3d3_helper_type::compress_ef(m_btnr, btnr_pos, bin, l);
          }
#else
          r3d3_helper_type::compress_ef(m_btnr, btnr_pos, bin, l);
#endif
        }
        btnr_pos += space_for_bt;
        pos += t_bs;
      }//end full blocks
      if (pos < m_size) { // handle last not full block
        if ((i % t_k) == (size_type)0) {
          m_btnrp[ i/t_k ] = btnr_pos;
          m_rank[ i/t_k ] = sum_rank;
          m_invert[ i/t_k ] = 0; // default: set last superblock to not inverted
          invert = false;
        }

        uint16_t space_for_bt = r3d3_helper_type::space_for_bt(x=bt_array[i++]);
        assert(i == bt_array.size());

        sum_rank += (invert ? (t_bs - x) : x);

        //calculating 'l' value for this block's EF compression
        uint16_t l = (x > 0 ? log2(t_bs/x) : 0);
        if (space_for_bt) {
          number_type bin = r3d3_helper_type::get_int(bv, pos, m_size-pos);
          if (invert) bin = ~bin; //store the negate of offset
#ifndef R3D3_C1_NO_OPT
          if (x == 1){
            uint16_t msb_bit_pos = r3d3_helper_type::hi(bin);
            r3d3_helper_type::trait::set_int(m_btnr, btnr_pos, msb_bit_pos, space_for_bt);
          } else {
            r3d3_helper_type::compress_ef(m_btnr, btnr_pos, bin, l);
          }
#else
          r3d3_helper_type::compress_ef(m_btnr, btnr_pos, bin, l);
#endif
        }

        assert(m_rank.size()-1 == ((i+t_k-1)/t_k));
      } else { // handle last empty full block
        assert(m_rank.size()-1 == ((i+t_k-1)/t_k));
      }
      // for technical reasons we add a last element to m_rank
      m_rank[ m_rank.size()-1 ] = sum_rank; // sum_rank contains the total number of set bits in bv
      m_bt = bt_array;
    }

    //! Swap method
    void swap(r3d3_vector& r3d3)
    {
      if (this != &r3d3) {
        std::swap(m_size, r3d3.m_size);
        m_bt.swap(r3d3.m_bt);
        m_btnr.swap(r3d3.m_btnr);
        m_btnrp.swap(r3d3.m_btnrp);
        m_rank.swap(r3d3.m_rank);
        m_invert.swap(r3d3.m_invert);
      }
    }

    void printCompressedBlock(const size_type &id){
      size_type bt_idx = id; //block id
      uint16_t bt = m_bt[bt_idx];
      size_type s_id = bt_idx/t_k; //superblock's id
      size_type btnrp = m_btnrp[ s_id ]; //pos until superblock

      for (size_type j = s_id*t_k; j < bt_idx; ++j) {
        btnrp += r3d3_helper_type::space_for_bt(m_bt[j]);
      }
      uint16_t btnrlen = r3d3_helper_type::space_for_bt(bt);

      std::cout<<" bt_id: "<<bt_idx<<" ; bt: "<<bt<<" ; s_id: "<<s_id<<std::endl;
      std::cout<<" btnrp: "<<btnrp<<" ; btnrlen: "<<btnrlen<<std::endl;
      std::cout<<" Compr Block "<<id<<" : "<<r3d3_helper_type::get_int(m_btnr, btnrp, btnrlen)<<std::endl;
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

      uint16_t l = (bt > 0 ? log2(t_bs/bt) : 0);

      size_type btnrp = m_btnrp[ s_id ];
      for (size_type j = s_id*t_k; j < bt_idx; ++j) {
        btnrp += r3d3_helper_type::space_for_bt(m_bt[j]);
      }
      uint16_t btnrlen = r3d3_helper_type::space_for_bt(bt);

      // and give it to block decoding function than
      bool is_bit_set = r3d3_helper_type::decode_bit_ef(m_btnr, btnrp, btnrlen, l, bt, off);

      if (m_invert[s_id]) is_bit_set = !is_bit_set;
      return is_bit_set;
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
    r3d3_vector& operator=(const r3d3_vector& r3d3)
    {
      if (this != &r3d3) {
        copy(r3d3);
      }
      return *this;
    }

    //! Move assignment operator
    r3d3_vector& operator=(r3d3_vector&& r3d3)
    {
      swap(r3d3);
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
      written_bytes += m_btnrp.serialize(out, child, "btnrp");
      written_bytes += m_rank.serialize(out, child, "rank_samples");
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
      m_btnrp.load(in);
      m_rank.load(in);
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
  struct rank_support_r3d3_trait {
    typedef bit_vector::size_type size_type;
    static size_type adjust_rank(size_type r, SDSL_UNUSED size_type n)
    {
      return r;
    }
  };

  template<>
  struct rank_support_r3d3_trait<0> {
    typedef bit_vector::size_type size_type;
    static size_type adjust_rank(size_type r, size_type n)
    {
      return n - r;
    }
  };

  //! rank_support for the r3d3_vector class
  /*!
   * \tparam t_b   The bit pattern of size one. (so `0` or `1`)
   * \tparam t_bs  The block size of the corresponding rrr_vector
   * \tparam t_rac Type used to store the block type in the corresponding rrr_vector.
   *  TODO: Test if the binary search can be speed up by
   *        saving the (n/2)-th rank value in T[0], the (n/4)-th in T[1],
   *        the (3n/4)-th in T[2],... for small number of rank values
   *    is this called hinted binary search???
   *    or is this called
   */
  template<uint8_t t_b, uint16_t t_bs, class t_rac, uint16_t t_k>
  class rank_support_r3d3
  {
    static_assert(t_b == 1u or t_b == 0u , "rank_support_r3d3: bit pattern must be `0` or `1`");
  public:
    typedef r3d3_vector<t_bs, t_rac, t_k> bit_vector_type;
    typedef typename bit_vector_type::size_type size_type;
    typedef typename bit_vector_type::r3d3_helper_type r3d3_helper_type;
    typedef typename r3d3_helper_type::number_type number_type;
    enum { bit_pat = t_b };
    enum { bit_pat_len = (uint8_t)1 };

  private:
    const bit_vector_type* m_v; //!< Pointer to the rank supported rrr_vector

  public:
    //! Standard constructor
    /*! \param v Pointer to the rrr_vector, which should be supported
     */
    explicit rank_support_r3d3(const bit_vector_type* v=nullptr)
    {
      set_vector(v);
    }

    //! Answers rank queries
    const size_type rank(size_type i)const
    {
      assert(m_v != nullptr);
      assert(i <= m_v->size()); //start indexing with 0
      size_type bt_idx = i/t_bs;//start indexing with 0
      size_type s_id = bt_idx/t_k;
      size_type btnrp = m_v->m_btnrp[ s_id ];
      size_type rank  = m_v->m_rank[ s_id ];
      const bool inv = m_v->m_invert[ s_id ];
      if (s_id+1 < m_v->m_rank.size()) {
        size_type diff_rank  = m_v->m_rank[ s_id+1 ] - rank;
#ifndef RRR_NO_OPT
        if (diff_rank == (size_type)0) {
          return  rank_support_r3d3_trait<t_b>::adjust_rank(rank, i);
        } else if (diff_rank == (size_type)t_bs*t_k) {
          size_type adj_rank = rank + i - s_id*t_k*t_bs;
          return  rank_support_r3d3_trait<t_b>::adjust_rank(
							    adj_rank, i);
        }
#endif
      }
      for (size_type j = s_id*t_k; j < bt_idx; ++j) {
        uint16_t r = m_v->m_bt[j];
        rank  += (inv ? t_bs - r: r);
        btnrp += r3d3_helper_type::space_for_bt(r);
      }
      uint16_t off = i % t_bs;

      uint16_t bt = m_v->m_bt[ bt_idx ];

      uint16_t l = (bt > 0 ? log2(t_bs/bt) : 0);
      uint16_t btnrlen = r3d3_helper_type::space_for_bt(bt);

      //returns rank in a block
      uint16_t rank_b = r3d3_helper_type::rank_on_ef_block(m_v->m_btnr, btnrp, btnrlen, l, bt, off);
      if (m_v->m_invert[s_id]) rank_b = off - rank_b;

      return rank_support_r3d3_trait<t_b>::adjust_rank(rank + rank_b, i);
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

    rank_support_r3d3& operator=(const rank_support_r3d3& rs)
    {
      if (this != &rs) {
        set_vector(rs.m_v);
      }
      return *this;
    }

    void swap(rank_support_r3d3&) { }

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


  //! Select support for the r3d3_vector class.
  /*
   * \tparam t_b   The bit pattern of size one. (so `0` or `1`)
   * \tparam t_bs  The block size of the corresponding rrr_vector
   * \tparam t_rac Type used to store the block type in the corresponding rrr_vector.
   *
   */
  template<uint8_t t_b, uint16_t t_bs, class t_rac, uint16_t t_k>
  class select_support_r3d3
  {
    static_assert(t_b == 1u or t_b == 0u , "select_support_rrr_ef: bit pattern must be `0` or `1`");
  public:
    typedef r3d3_vector<t_bs, t_rac, t_k> bit_vector_type;
    typedef typename bit_vector_type::size_type size_type;
    typedef typename bit_vector_type::r3d3_helper_type r3d3_helper_type;
    typedef typename r3d3_helper_type::number_type number_type;
    enum { bit_pat = t_b };
    enum { bit_pat_len = (uint8_t)1 };

  private:
    const bit_vector_type* m_v; //!< Pointer to the rank supported rrr_vector

    size_type select1(size_type i)const
    {
      if (m_v->m_rank[m_v->m_rank.size()-1] < i)
        return size();

      //  (1) binary search for the answer in the rank_samples (superblocks)
      size_type begin=0, end=m_v->m_rank.size()-1; // min included, max excluded
      size_type idx, rank;
      // invariant:  m_rank[end]   >= i
      //             m_rank[begin]  < i
      while (end-begin > 1) {
        idx  = (begin+end) >> 1; // idx in [0..m_rank.size()-1]
        rank = m_v->m_rank[idx];
        if (rank >= i)
          end = idx;
        else { // rank < i
          begin = idx;
        }
      }
      //   (2) linear search within the samples (superblock)
      rank = m_v->m_rank[begin]; // now i>rank
      idx = begin * t_k; // initialize idx for select result
      size_type diff_rank  = m_v->m_rank[end] - rank;
#ifndef RRR_NO_OPT
      if (diff_rank == (size_type)t_bs*t_k) {// optimisation for select<1>
        return idx*t_bs + i-rank -1;
      }
#endif
      const bool inv = m_v->m_invert[ begin ];
      size_type btnrp = m_v->m_btnrp[ begin ];
      uint16_t bt = 0, btnrlen = 0; // temp variables for block_type and space for block type
      while (i > rank) {
        bt = m_v->m_bt[idx++]; //bt = inv ? t_bs-bt : bt;
        rank += (inv ? t_bs - bt: bt);
        btnrp += (btnrlen=r3d3_helper_type::space_for_bt(bt));
      }
      rank -= (inv ? t_bs - bt: bt);

      uint16_t l = (bt > 0 ? log2(t_bs/bt) : 0);

      int pos_in_b = r3d3_helper_type::select_on_ef_block(m_v->m_btnr, btnrp-btnrlen, btnrlen, l, bt, i-rank, inv);
      return (idx-1) * t_bs + pos_in_b;
    }

    /*size_type select0(size_type i)const {Todo} */ 


  public:
    explicit select_support_r3d3(const bit_vector_type* v=nullptr)
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

    select_support_r3d3& operator=(const select_support_r3d3& rs)
    {
      if (this != &rs) {
        set_vector(rs.m_v);
      }
      return *this;
    }

    void swap(select_support_r3d3&) { }

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
