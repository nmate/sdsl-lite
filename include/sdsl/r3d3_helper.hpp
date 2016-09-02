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
/*! \file r3d3_helper.hpp
  \author Mate Nagy
*/
#ifndef SDSL_R3D3_HELPER
#define SDSL_R3D3_HELPER

#ifdef RRR_NO_OPT
#ifndef RRR_NO_BS
#define RRR_NO_BS
#endif
#endif

#include <algorithm> // for next permutation
#include <iostream>
#include "bits.hpp"
#include "uint128_t.hpp"
#include "uint256_t.hpp"

namespace sdsl
{

  //! Trait struct for the binomial coefficient struct to handle different type of integers.
  /*! This generic implementation works for 64-bit integers.
   */
  template<uint16_t log_n>
  struct r3d3_trait {
    typedef uint64_t number_type;
    static inline uint16_t hi(number_type x) {
      return bits::hi(x);
    }

    //! Read a \f$len\f$-bit integer of type number_type from a bitvector.
    /*!
     * \param bv   A bit_vector of int_vector from which we extract the integer.
     * \param pos  Position of the least significant bit of the integer which should be read.
     * \param len  bit-width of the integer which should be read.
     * \return The len-bit integer.
     */
    template<class bit_vector_type>
    static inline number_type get_int(const bit_vector_type& bv,
                                      typename bit_vector_type::size_type pos,
                                      uint16_t len) {
      return bv.get_int(pos, len);
    }

    //! Write a \f$len\f$-bit integer x of type number_type to a bitvector.
    /*!
     * \param bv     A bit_vecor or int_vector in which we write the integer.
     * \param pos    Position of the least significant bit of the integer which should be written.
     * \param x    The integer x which should be written.
     * \param len  Bit-width of x.
     */
    template<class bit_vector_type>
    static void set_int(bit_vector_type& bv, typename bit_vector_type::size_type pos,
                        number_type x, uint16_t len) {
      bv.set_int(pos, x, len);
    }

    //! Count the number of set bits in x.
    /*!
     *  \param x The integer x.
     */
    static inline uint16_t popcount(number_type x) {
      return bits::cnt(x);
    }

    static inline uint16_t select(number_type &word, uint32_t index){
      return bits::sel(word, index);
    }

    static inline bool isBitSet(number_type &word, uint32_t index){
      number_type mask = (uint64_t)1 << index;
      return word & mask;
    }
  };

  //! Specialization of binomial_coefficients_trait for 128-bit integers.
  template<>
  struct r3d3_trait<7> {
    typedef uint128_t number_type;
    static inline uint16_t hi(number_type x) {
      if ((x >> 64)) {
        return bits::hi(x >> 64) + 64;
      } else {
        return bits::hi(x);
      }
    }

    template<class bit_vector_type>
    static inline number_type get_int(const bit_vector_type& bv,
                                      typename bit_vector_type::size_type pos,
                                      uint16_t len) {
      if (len <= 64) {
        return bv.get_int(pos, len);
      } else {
        return ((((number_type) bv.get_int(pos+64, len-64))<<64) + bv.get_int(pos, 64));
      }
    }

    template<class bit_vector_type>
    static void set_int(bit_vector_type& bv,
                        typename bit_vector_type::size_type pos,
                        number_type x, uint16_t len) {
      if (len <= 64) {
        bv.set_int(pos, x, len);
      } else {
        bv.set_int(pos, (uint64_t)x, 64); bv.set_int(pos+64, x>>64, len-64);
      }
    }

    static inline uint16_t popcount(number_type x) {
      return bits::cnt(x >> 64) + bits::cnt(x);
    }

    static inline uint16_t select(number_type &word, uint32_t index){
      return word.select(index);
    }

    static inline bool isBitSet(number_type &word, uint32_t index){
      return word.isBitSet(index);
    }
  };

  //! Specialization of r3d3_trait for 256-bit integers.
  template<>
  struct r3d3_trait<8> {
    typedef uint256_t number_type;
    static inline uint16_t hi(number_type x) {
      return x.hi();
    }

    template<class bit_vector_type>
    static inline number_type get_int(const bit_vector_type& bv,
                                      typename bit_vector_type::size_type pos,
                                      uint16_t len) {
      if (len <= 64) {
        return number_type(bv.get_int(pos, len));
      } else if (len <= 128) {
        return number_type(bv.get_int(pos, 64), bv.get_int(pos+64, len-64));
      } else if (len <= 192) {
        return number_type(bv.get_int(pos, 64), bv.get_int(pos + 64, 64),
                           (uint128_t)bv.get_int(pos + 128, len-128));
      } else { // > 192
        return number_type(bv.get_int(pos, 64), bv.get_int(pos+64, 64),
                           (((uint128_t)bv.get_int(pos+192, len-192))<<64) | bv.get_int(pos+128, 64));
      }
    }

    template<class bit_vector_type>
    static void set_int(bit_vector_type& bv,
                        typename bit_vector_type::size_type pos,
                        number_type x,
                        uint16_t len) {
      if (len <= 64) {
        bv.set_int(pos, x, len);
      } else if (len <= 128) {
        bv.set_int(pos, x, 64); bv.set_int(pos+64, x>>64, len-64);
      } else if (len <= 192) {
        bv.set_int(pos, x, 64); bv.set_int(pos+64, x>>64, 64);
        bv.set_int(pos+128, x>>128, len-128);
      } else { // > 192
        bv.set_int(pos, x, 64); bv.set_int(pos+64, x>>64, 64);
        bv.set_int(pos+128, x>>128, 64); bv.set_int(pos+192, x>>192, len-192);
      }
    }

    static inline uint16_t popcount(number_type x) {
      return x.popcount();
    }

    static inline uint16_t select(number_type &word, uint32_t index){
      return word.select(index);
    }

    static inline bool isBitSet(number_type &word, uint32_t index){
      return word.isBitSet(index);
    }

  };

  // For fixed length EF compression (R3D3). 
  // Length only depends on 'c' but not on 
  // 'B', however it consumes more space
  template<uint16_t n> //n is blocksize
  struct r3d3_len_ef_table {
    static struct impl {
      uint16_t table[n+1];              

      impl() {
        //initialize with 0s
        for (uint16_t k=0; k < n+1; ++k) {
          table[k] = 0;
        }

        //fill in with class dependent 'length'
        for (uint16_t c=1; c <= n; ++c) {
          int l = 0;

          l = log2(n/c);
          table[c] = n/(1<<l) + c + c*l;
        }
      }
    } data;
  };

  // For variable length EliasFano blocks (R3D3i)
  // 'b' and 'c' defines the length of the compressed offset
  template<uint16_t n> //n is blocksize
  struct r3d3i_len_ef_table {
    static struct impl {
      uint16_t table[n+1][n+1];

      impl() {
        //initialize with 0s
        for (uint16_t j=0; j<n+1; ++j)
          for (uint16_t k=0; k < n+1; ++k) {
            table[j][k] = 0;
          }

        //fill in with B and c (class) dependent 'length'
        // 0 <= B <= n-1 ; msb bit position starts with 0
        // 1 <= c <= B+1 ; class valid interval
        for (uint16_t B = 0; B < n; ++B)
          for (uint16_t c=1; c <= B+1; ++c) {
            int l = 0;

            l = (B > 0 ? log2(B/c) : 0);
            table[B][c] = B/(1<<l) + c + c*l;
          }
      }
    } data;
  };

  // For obtaining 'l' values. In case of using
  // variable length offsets we do not know the
  // MSB bit position 'B' which is needed for
  // the calculation of 'l'. Instead we use this
  // offline table: [size, class]-> l
  // 0 <= size <= 2*t_bs, 0 <= class <= t_bs
  template<uint16_t n> //n is blocksize
  struct r3d3i_l_ef_table {
    static struct impl {
      uint16_t table[2*(n+1)][n+1];

      impl() {
        //initialize with 0s
        for (uint16_t j=0; j<2*(n+1); ++j)
          for (uint16_t k=0; k < n+1; ++k) {
            table[j][k] = 0;
          }

        //B: msb_bit_pos goes 0..t_bs-1, e.g. 0..15
        //but if B=0, l=0
        for (uint16_t B = 1; B <= n-1; ++B)
          for (uint16_t c = 1; c <= n; ++c) {
            int l = 0;

            //fill up extended table
            l = log2(B/c);
            uint16_t offset_size =  B/(1<<l) + c + c*l;
            table[offset_size][c] = l;
          }
      }
    } data;
  };

  // Datastructure that stores number of bits that is
  // necessary for storing a single bit position value
  // (c=1 optimization).
  template<uint16_t n> //n is blocksize
  struct c1_len_table {
    static struct impl {
      uint16_t table[n+1];

      impl(){
        for (uint16_t bp = 0; bp < n+1; ++bp){
          uint16_t space;
          space = (bp > 0) ? log2(bp) + 1 : 0;
          table[bp] = space;
        }
      }
    } data;
  };

  template<uint16_t n>
  typename r3d3_len_ef_table<n>::impl r3d3_len_ef_table<n>::data;

  template<uint16_t n>
  typename r3d3i_len_ef_table<n>::impl r3d3i_len_ef_table<n>::data;

  template<uint16_t n>
  typename r3d3i_l_ef_table<n>::impl r3d3i_l_ef_table<n>::data;

  template<uint16_t n>
  typename c1_len_table<n>::impl c1_len_table<n>::data;

  //! R3D3 helper for EliasFano compressed offset
  template<uint16_t n>
  struct r3d3_helper{
    typedef r3d3_len_ef_table<n>    r3d3_len;
    typedef r3d3i_len_ef_table<n>   r3d3i_len;
    typedef r3d3i_l_ef_table<n>     r3d3i_l;
    typedef c1_len_table<n>         c1_len;

    //let's use predefined functions
    enum {MAX_LOG = (n>128 ? 8 : (n > 64 ? 7 : 6))};
    typedef r3d3_trait<MAX_LOG> trait;
    typedef typename trait::number_type number_type;

    static inline uint16_t space_for_bt(const uint16_t &c) {
      return r3d3_len::data.table[c];
    }

    static inline uint16_t space_for_bt_i(const uint16_t &b, const uint16_t &c) {
      return r3d3i_len::data.table[b][c];
    }

    static inline uint16_t size_of_l(const uint16_t &offset_size, const uint16_t &c){
      return r3d3i_l::data.table[offset_size][c];
    }

    static inline uint16_t space_for_class_1(const uint16_t &bp) {
      return c1_len::data.table[bp];
    }

    //! createLBA(bv, pos, bin, l):
    //
    //  Builds up the Lower Bits Array
    //
    // @Params:
    //  bv : offset, this is the output
    //  pos: put compressed result at pos
    //  bin: number to be compressed, input
    //  l  : length of an LBA item
    //
    template<class bit_vector_type>
    static inline void createLBA(bit_vector_type& btnr, 
                                 typename bit_vector_type::size_type pos,
                                 const number_type &bin, const uint16_t &l){
      uint32_t lowVal = 0;
      number_type num = bin;
      int maskLBA = (1 << l) - 1;

      if (l == 0) return; //nothing to do

      //leftmost mask
      uint16_t bitPos   = 0; //starting position with 0!
      number_type probe = 1;
      uint32_t posInLBA = pos;
      while ( num != (number_type)0 ){
        if ( num & probe ){
          //insert LBA value
          lowVal = bitPos & maskLBA;

          //insert number into the offset
          trait::set_int(btnr, posInLBA, lowVal, l);
          posInLBA += l;
        }
        bitPos++;
        num = num >> 1;
      }
    }

    //! createUBA(bv, pos, bin, l):
    //
    // Builds up the Upper-Bits Array
    //
    // @Params:
    //  bv : the bitvector in which we store the output
    //  pos: starting position in the bitvector
    //  bin: number to be compressed, input
    //  l  : length of an LBA item
    //
    //  We use negated unary code 
    template<class bit_vector_type>
    static inline void createUBA(bit_vector_type& btnr, 
                                 typename bit_vector_type::size_type pos,
                                 const number_type &bin, const uint16_t &l){
      number_type num     = bin;
      number_type upVal   = 0;
      uint16_t bitPos     = 0; //starting position with 0!
      number_type probe   = 1;
      uint16_t popCnt     = trait::popcount(num);
      number_type prevUpVal  = 0;
      number_type posInUBA   = pos + l*popCnt; //LBA comes first
      while ( num != (number_type)0 ){
        if ( num & probe ){
          upVal = bitPos;
          upVal = upVal >> l; //cut off the lower part
          uint16_t gap = upVal - prevUpVal; //calc gap
          uint16_t numOfZeros = 0;
          while ( gap-- ){
            numOfZeros++;
          }
          prevUpVal = upVal;
          //reuse upVal
          upVal = 1;
          upVal = upVal << numOfZeros;

          //negate upVal
          upVal = ~upVal;

          //insert upVal into m_btnr
          trait::set_int(btnr, posInUBA, upVal, numOfZeros+1);
          posInUBA += numOfZeros + 1;
        }
        //move forward
        bitPos++;
        num = num >> 1;
      }
    }

    //! Elias-Fano compression
    //
    // @Params
    //  bv : bitvector that stores the output
    //  bin: input that is to be compressed
    template<class bit_vector_type>
    static inline void compress_ef(bit_vector_type& bv, typename bit_vector_type::size_type pos,
                                   const number_type &bin, const uint16_t &l){
      createLBA(bv, pos, bin, l);
      createUBA(bv, pos, bin, l);
    }

    //! Decompress the index'th bit position value from EF offset.
    template<class bit_vector_type>
    static inline uint16_t decode_value_ef(bit_vector_type& bv, typename bit_vector_type::size_type btnrp,
                                           const uint16_t &btnrlen, const uint16_t &l, const uint16_t& bt, const uint16_t &index){

      //to avoid return of invalid data
      if (index == 0 || index > bt){
        return -1; //hope it will be visible
      }

#ifndef RRR_NO_OPT
      if (bt == n) {  // if bt==n, then the encoded block consists only of ones
        return index-1;
      } else if (bt == 0) { // if bt==0 then the encoded block consists only of zeros
        return 0;
      }
#endif

#ifndef R3D3_C1_NO_OPT
      if (bt == 1){
        uint16_t bit_pos = trait::get_int(bv, btnrp, btnrlen);
        return bit_pos;
      }
#endif

      uint16_t lba_width = bt*l;
      uint16_t uba_width = btnrlen-lba_width;
      uint16_t chunk_size = 8*sizeof(number_type);
      uint16_t chunk_id = 0, _index=index;
      number_type uba[2];
      uint16_t popcnt_0;
      uint16_t nr_gaps_0 = 0;

      //read the first UBA chunk, and 2nd it it exists
      if (uba_width <= chunk_size){
        uba[0] = trait::get_int(bv, btnrp + lba_width, uba_width);
      }
      else{
        uba[0] = trait::get_int(bv, btnrp + lba_width, chunk_size);
        uba[1] = trait::get_int(bv, btnrp + lba_width + chunk_size, uba_width-chunk_size);
        
        nr_gaps_0 = trait::popcount(uba[0]);
        popcnt_0 = chunk_size - nr_gaps_0;

        //For Select operation choose the right chunk
        if (popcnt_0 < index){
          _index -= popcnt_0;
          chunk_id++;
        }
      }

      // UBA is stored in negated unary code (0: value, 1:gap)
      // so now we need to switch it because we do a search
      // on values
      uba[chunk_id] = ~uba[chunk_id]; // value: 1
      uint16_t uba_pos = trait::select(uba[chunk_id], _index) + 1;

      if (chunk_id != 0)  uba_pos += chunk_size;
      uint16_t uba_val = uba_pos - index;

      //extracting lbaValue
      uint16_t lba_val = trait::get_int(bv, btnrp+(index-1)*l, l);

      //calculating final value
      return (1<<l)*uba_val + lba_val;
    }

    // Do rank on a block for r3d3_vector datastructure
    template<class bit_vector_type>
    static inline uint16_t rank_on_ef_block(bit_vector_type& bv, typename bit_vector_type::size_type btnrp,
                                            uint16_t btnrlen, const uint16_t &l, const uint16_t& bt,
                                            const uint16_t &r_i){
      //optimization
#ifndef RRR_NO_OPT
      if (bt == 0) {  // if bt==0, then the encoded block consists only zeros
        return 0;
      } else if (r_i == 0){
        return 0;
      } else if (bt == n) { // if bt==n then the encoded block consists only of ones     
        return r_i;
      }
#endif

#ifndef R3D3_C1_NO_OPT
      if (bt == 1){
        uint16_t bit_pos = trait::get_int(bv, btnrp, btnrlen);
        return (r_i > bit_pos) ? 1 : 0;
      }
#endif
     
      //choose proper uba chunk
      uint16_t lba_width = bt*l;
      uint16_t uba_width = btnrlen-lba_width;
      uint16_t uba_width_eff; //effective uba_width
      uint16_t chunk_size = 8*sizeof(number_type);
      uint16_t chunk_id = 0;
      number_type uba[2];
      uint16_t gap = 1<<l;
      uint16_t gap_id = r_i / gap;
      uint16_t nr_gaps_0 = 0;//only used if uba[1] is valid
      uint16_t nr_gaps;
      uint16_t rank = 0;

      // Read the first chunk, and 2nd if exists 
      if (uba_width <= chunk_size){
        uba[0] = trait::get_int(bv, btnrp + lba_width, uba_width);
        nr_gaps = trait::popcount(uba[0]);
        uba_width_eff = nr_gaps + bt;
        //requested gap is greater than what we have stored
        if (gap_id > nr_gaps) return bt; 
      }
      // If UBA fits into 2 chunks, we need to
      // find the chunk in which the gap resides
      // (belonging to r_i)
      else {
        uba[0] = trait::get_int(bv, btnrp + lba_width, chunk_size);
        uba[1] = trait::get_int(bv, btnrp + lba_width + chunk_size, uba_width-chunk_size);
        nr_gaps_0 = trait::popcount(uba[0]);
        uba_width_eff = nr_gaps_0 + trait::popcount(uba[1]) + bt;
        nr_gaps = uba_width_eff - bt;        
        
        // Check if gap is in the 2nd chunk
        // if it is, normalize gap_id and 
        // add popcnt(uba[0]) to rank
        if (nr_gaps_0 < gap_id){          
          if (gap_id > nr_gaps) return bt;
          chunk_id = 1;
        }
      }

      // Find the position of the value next to the
      // gap_id. In the UBA we store negated unary
      // code, where 0 indicates value and 1 means gap.
      // Val_pos points to the pos of gap_id. The next
      // bit has meaning for us. 
      uint16_t gap_pos = 0;
      uint16_t bit_pos_value = 0, lba_i;
      uint16_t curr_gap_id = chunk_id ? gap_id-nr_gaps_0 : gap_id;

      if (gap_id != 0){ //select(0) is invalid!
        gap_pos = trait::select(uba[chunk_id], curr_gap_id) + 1;
        rank = gap_pos - curr_gap_id;
        if ( chunk_id ) rank += chunk_size - nr_gaps_0;

        // If gap is the last element of the first chunk
        // then start iteration on the second chunk
        if (gap_pos == chunk_size){
          chunk_id = 1;
          gap_pos = 0;
          curr_gap_id = gap_id-nr_gaps_0;
        }
      }

      // If there are values (0) stored where gap_pos
      // points at, extract them and add to rank if
      // those bit position values are smaller than r_i
      while( gap_pos < uba_width_eff &&
             !trait::isBitSet(uba[chunk_id], gap_pos) ){
        lba_i = gap_pos - curr_gap_id;
        bit_pos_value = gap_id*gap + trait::get_int(bv, btnrp+(lba_i)*l, l);
        if (bit_pos_value < r_i) {
          rank = lba_i + 1;
        }
        else { break; }
        gap_pos++;
        // We switch to uba[1]
        if (gap_pos == chunk_size) {
          chunk_id = 1;
          gap_pos = 0;
          curr_gap_id = gap_id-nr_gaps_0;
        }
      }
      return rank;
    }


    // Do rank on a block for r3d3i_vector datastructure
    template<class bit_vector_type>
    static inline uint16_t rank_on_ef_block_i(bit_vector_type& bv, typename bit_vector_type::size_type btnrp,
                                              uint16_t btnrlen, const uint16_t &l, const uint16_t& bt,
                                              const uint16_t &r_i){
      //optimization
#ifndef RRR_NO_OPT
      if (bt == 0) {  // if bt==0, then the encoded block consists only zeros
        return 0;
      } else if (r_i == 0){
        return 0;
      } else if (bt == n) { // if bt==n then the encoded block consists only of ones     
        return r_i;
      }
#endif

#ifndef R3D3_C1_NO_OPT
      if (bt == 1){
        uint16_t bit_pos = trait::get_int(bv, btnrp, btnrlen);
        return (r_i > bit_pos) ? 1 : 0;
      }
#endif
     
      //choose proper uba chunk
      uint16_t lba_width = bt*l;
      uint16_t uba_width = btnrlen-lba_width;
      uint16_t chunk_size = 8*sizeof(number_type);
      uint16_t chunk_id = 0;
      number_type uba[2];
      uint16_t gap = 1<<l;
      uint16_t gap_id = r_i / gap;
      uint16_t nr_gaps_0 = 0;
      uint16_t nr_gaps = uba_width - bt;
      uint16_t rank = 0;
 
      // If requested gap_id is greater what 
      // we store, rank is bt 
      if (gap_id > nr_gaps) return bt;

      // Read the first chunk, and 2nd if exists
      if (uba_width <= chunk_size){
        uba[0] = trait::get_int(bv, btnrp + lba_width, uba_width);
      }
      // If UBA fits into 2 chunks, we need to
      // find the chunk in which the gap resides
      // (belonging to r_i)
      else {
        uba[0] = trait::get_int(bv, btnrp + lba_width, chunk_size);
        uba[1] = trait::get_int(bv, btnrp + lba_width + chunk_size, uba_width-chunk_size);
        nr_gaps_0 = trait::popcount(uba[0]);

        // Check if gap is in the 2nd chunk
        // if it is, normalize gap_id and 
        // add popcnt(uba[0]) to rank
        if (nr_gaps_0 < gap_id){
          chunk_id = 1;
        }
      }

      // Find the position of the value next to the
      // gap_id. In the UBA we store negated unary
      // code, where 0 indicates value and 1 means gap.
      // Val_pos points to the pos of gap_id. The next
      // bit has meaning for us. 
      uint16_t gap_pos = 0;
      uint16_t bit_pos_value = 0, lba_i;
      uint16_t curr_gap_id = chunk_id ? gap_id-nr_gaps_0 : gap_id;

      if (gap_id != 0){ //select(0) is invalid!
        gap_pos = trait::select(uba[chunk_id], curr_gap_id) + 1;
        rank = gap_pos - curr_gap_id;
        if ( chunk_id ) rank += chunk_size - nr_gaps_0;

        // If gap is the last element of the first chunk
        // then start iteration on the second chunk
        if (gap_pos == chunk_size){
          chunk_id = 1;
          gap_pos = 0;
          curr_gap_id = gap_id-nr_gaps_0;
        }
      }

      // If there are values (0) stored where gap_pos
      // points at, extract them and add to rank if
      // those bit position values are smaller than r_i
      while( gap_pos < uba_width &&
             !trait::isBitSet(uba[chunk_id], gap_pos) ){
        lba_i = gap_pos - curr_gap_id;
        bit_pos_value = gap_id*gap + trait::get_int(bv, btnrp+(lba_i)*l, l);
        if (bit_pos_value < r_i) {
          rank = lba_i + 1;
        }
        else { break; }
        gap_pos++;
        // We switch to uba[1]
        if (gap_pos == chunk_size) {
          chunk_id = 1;
          gap_pos = 0;
          curr_gap_id = gap_id-nr_gaps_0;
        }
      }

      return rank;
    }


    //! Select on an EF compressed block
    template<class bit_vector_type>
    static inline uint16_t select_on_ef_block(bit_vector_type& bv, typename bit_vector_type::size_type btnrp,
                                              const uint16_t &btnrlen, const uint16_t &l, const uint16_t& bt,
                                              const uint16_t &s_i, const uint16_t &inv){

      //if not inverted, return stored index
      if (!inv) return decode_value_ef(bv, btnrp, btnrlen, l, bt, s_i);

      // ------- for inverted offset we need below part -------
      //         this part should not be called much 
      //         since p < 0.5 is our use case
      //else
      // 1.) decompress 0 bit positions and define from those
      //     values the num_of_1s until that bit position
      uint16_t num_of_1s_high = 0, num_of_1s_low = 0;
      uint16_t ef_index = 1, bit_pos_0_high = 0, bit_pos_0_low = 0, last_ef_index = 1;

      while ( num_of_1s_high < s_i && ef_index <= bt ){
        //save lower values
        num_of_1s_low = num_of_1s_high;
        bit_pos_0_low = bit_pos_0_high;
        last_ef_index = ef_index;
        //try new value
        bit_pos_0_high = decode_value_ef(bv, btnrp, btnrlen, l, bt, ef_index);
        num_of_1s_high = bit_pos_0_high + 1 - ef_index;
        ef_index++;
      }

      // 2.) now we know that between bit_pos_0_low and bit_pos_0_high
      //     there are only 1s, count them until s_i and step the index too.

      // if query is less than first index in EF s_i = s_i-1
      if ((last_ef_index == 1 && num_of_1s_high >= s_i) || bt == 0)
        return s_i-1;

      uint16_t bit_pos;
      uint16_t num_of_1s;
      // if I reached the end of block and still do not have
      // the proper num of 1's. Query is more than the class.
      if (last_ef_index == bt && num_of_1s_high < s_i){
        bit_pos = bit_pos_0_high;
        num_of_1s = num_of_1s_high;
      }
      else{
        bit_pos = bit_pos_0_low;
        num_of_1s = num_of_1s_low;
      }

      while (num_of_1s != s_i){
        num_of_1s++;
        bit_pos++;
      }

      return bit_pos;
    }

    template<class bit_vector_type>
    static inline number_type get_int(const bit_vector_type& bv,
                                      typename bit_vector_type::size_type btnrp,
                                      uint16_t btnrlen) {
      return trait::get_int(bv, btnrp, btnrlen);
    }

    //returns the popcnt of a block_size long part of the
    //original bitmap
    template<class bit_vector_type>
    static inline uint16_t get_bt(const bit_vector_type& bv, typename bit_vector_type::size_type pos,
                                  uint16_t block_size) {
      return trait::popcount(trait::get_int(bv, pos, block_size));
    }

    //returns most significant bit
    static inline uint16_t hi(number_type x) {
      return trait::hi(x);
    }

    //! For the fixed length EF compressed offset. In this case there is a padding
    //  in the offset and we need to separate it from the useful part (see uba_width_eff).
    //  Gives back if pos'th bit was 1 or 0 (access operation)
    template<class bit_vector_type>
    static inline bool decode_bit_ef(bit_vector_type& bv, typename bit_vector_type::size_type btnrp,
                                     const uint16_t &btnrlen, const uint16_t &l, const uint16_t& bt,
                                     const uint16_t &pos_query){
      //optimization
#ifndef RRR_NO_OPT
      if (bt == 0) {  // if bt==0, then the encoded block consists only of zeros
        return false;
      } else if (bt == n) { // if bt==n then the encoded block consists only of ones
        return true;
      }
#endif
      // c=1 optimization
#ifndef R3D3_C1_NO_OPT
      if (bt == 1){
        uint16_t bit_pos = trait::get_int(bv, btnrp, btnrlen);
        return (bit_pos == pos_query) ? true : false;
      }
#endif

      //choose proper uba chunk
      uint16_t lba_width = bt*l;
      uint16_t uba_width = btnrlen-lba_width; //holds padding
      uint16_t uba_width_eff; //effective uba_width
      uint16_t chunk_size = 8*sizeof(number_type);
      uint16_t chunk_id = 0;
      number_type uba[2];
      uint16_t gap = 1<<l;
      uint16_t gap_id = pos_query / gap;
      uint16_t nr_gaps_0 = 0; //only used if uba[1] valid
      uint16_t nr_gaps;

      // Read the first chunk, and 2nd if exists
      if (uba_width <= chunk_size){
        uba[0] = trait::get_int(bv, btnrp + lba_width, uba_width);
        nr_gaps = trait::popcount(uba[0]);
        uba_width_eff = nr_gaps + bt;
        if (gap_id > nr_gaps) return false;
      }
      // If UBA fits into 2 chunks, we need to
      // find the chunk in which the gap resides
      // (belonging to pos_query)
      else {
        uba[0] = trait::get_int(bv, btnrp + lba_width, chunk_size);
        uba[1] = trait::get_int(bv, btnrp + lba_width + chunk_size, uba_width-chunk_size);
        nr_gaps_0 = trait::popcount(uba[0]);
        uba_width_eff = nr_gaps_0 + trait::popcount(uba[1]) + bt;
        nr_gaps = uba_width_eff - bt;

        // Check if gap is in the 2nd chunk
        // if it is, normalize gap_id
        if (nr_gaps_0 < gap_id){
          if (gap_id > nr_gaps) return false;
          chunk_id = 1;
        }
      }

      // Find the position of the value, next to the
      // gap_id. In the UBA we store negated unary
      // code, where 0 indicates value and 1 means gap
      uint16_t gap_pos = 0;
      uint16_t bit_pos_value, lba_i;
      uint16_t curr_gap_id = chunk_id ? gap_id-nr_gaps_0 : gap_id;

      if (gap_id != 0){ //select(0) is invalid!
        gap_pos = trait::select(uba[chunk_id], curr_gap_id) + 1;

        // If gap is the last element of the first chunk
        // then start iteration on the second chunk
        if (gap_pos == chunk_size){
          chunk_id = 1;
          gap_pos = 0;
          curr_gap_id = gap_id-nr_gaps_0;
        }
      }

      // If there is no gap between values, extract them!
      while( gap_pos < uba_width_eff &&
             !trait::isBitSet(uba[chunk_id], gap_pos) ){
        // Extract the bit pos value
        lba_i = gap_pos - curr_gap_id;
        bit_pos_value = gap_id*gap + trait::get_int(bv, btnrp+(lba_i)*l, l);
        if (bit_pos_value == pos_query) return true;
        gap_pos++;
        // We switch to uba[1]
        if (gap_pos == chunk_size) {
          chunk_id = 1;
          gap_pos = 0;
          curr_gap_id = gap_id-nr_gaps_0;
        }
      }
      return false;
    }

    //! For variable length EF offset. In this case there is no
    //  padding added to the offset, we only store useful bits.
    //  Gives back if pos'th bit was 1 or 0 (access operation)
    template<class bit_vector_type>
    static inline bool decode_bit_ef_i(bit_vector_type& bv, typename bit_vector_type::size_type btnrp,
                                       const uint16_t &btnrlen, const uint16_t &l, const uint16_t& bt,
                                       const uint16_t &pos_query){
      //optimization
#ifndef RRR_NO_OPT
      if (bt == 0) {  // if bt==0, then the encoded block consists only of zeros
        return false;
      } else if (bt == n) { // if bt==n then the encoded block consists only of ones
        return true;
      }
#endif
      // c=1 optimization
#ifndef R3D3_C1_NO_OPT
      if (bt == 1){
        uint16_t bit_pos = trait::get_int(bv, btnrp, btnrlen);
        return (bit_pos == pos_query) ? true : false;
      }
#endif

      //choose proper uba chunk
      uint16_t lba_width = bt*l;
      uint16_t uba_width = btnrlen-lba_width;
      uint16_t chunk_size = 8*sizeof(number_type);
      uint16_t chunk_id = 0;
      number_type uba[2];
      uint16_t gap = 1<<l;
      uint16_t gap_id = pos_query / gap;
      uint16_t nr_gaps_0 = 0;
      uint16_t nr_gaps = uba_width - bt;

      // If gap_id is greater than nr of gaps we have in UBA
      // then bit position is not stored, thus it's value is 0
      if (gap_id > nr_gaps) return false;

      // Read the first chunk, and 2nd if exists
      if (uba_width <= chunk_size){
        uba[0] = trait::get_int(bv, btnrp + lba_width, uba_width);
      }
      // If UBA fits into 2 chunks, we need to
      // find the chunk in which the gap resides
      // (belonging to pos_query)
      else {
        uba[0] = trait::get_int(bv, btnrp + lba_width, chunk_size);
        uba[1] = trait::get_int(bv, btnrp + lba_width + chunk_size, uba_width-chunk_size);
        nr_gaps_0 = trait::popcount(uba[0]);

        // Check if gap is in the 2nd chunk
        // if it is, normalize gap_id
        if (nr_gaps_0 < gap_id){
          chunk_id = 1;
        }
      }

      // Find the position of the value, next to the 
      // gap_id. In the UBA we store negated unary
      // code, where 0 indicates value and 1 means gap
      uint16_t gap_pos = 0;
      uint16_t bit_pos_value, lba_i;

      // curr_gap_id refers to the gap_id within a chunk
      uint16_t curr_gap_id = chunk_id ? gap_id-nr_gaps_0 : gap_id;

      if (gap_id != 0){ //select(0) is invalid!
        gap_pos = trait::select(uba[chunk_id], curr_gap_id) + 1;

        // If gap is the last element of the first chunk
        // then start iteration on the second chunk
        if (gap_pos == chunk_size){
          chunk_id = 1;
          gap_pos = 0;
          curr_gap_id = gap_id-nr_gaps_0;
        }
      }

      // If there is no gap between values, extract them!
      while( gap_pos < uba_width &&
             !trait::isBitSet(uba[chunk_id], gap_pos) ){
        // Extract the bit pos value
        lba_i = gap_pos - curr_gap_id;
        bit_pos_value = gap_id*gap + trait::get_int(bv, btnrp+(lba_i)*l, l);
        if (bit_pos_value == pos_query) return true;
        gap_pos++;
        // We switch to uba[1]
        if (gap_pos == chunk_size) {
          chunk_id = 1;
          gap_pos = 0;
          curr_gap_id = gap_id-nr_gaps_0;
        }
      }
      return false;
    }

    static inline void print_R3d3_LenTable(){
      std::cout<<"Table: "<<std::endl;
      for (uint16_t c=0; c < n + 1; ++c) {
        std::cout<<" ["<<n<<", "<<c<<"]: "<<r3d3_len::data.table[c];
      }
      std::cout<<std::endl;
    }

    static inline void print_R3d3i_LowerValueTable(){
      std::cout<<"Table_l: "<<std::endl;
      for (uint16_t off_size=0; off_size < 2*(n + 1); ++off_size){
        for (uint16_t c=0; c <= n; ++c) {
          int l = ((off_size > 0 && c>0) ? log2(off_size/c) : 0);
          std::cout<<" ["<<off_size<<", "<<c<<"]: "<<r3d3i_l::data.table(off_size, c);
        }
      }
      std::cout<<std::endl;
    }

    static inline void print_R3d3i_C1_LenTable(){
      std::cout<<"Table_class_1: "<<std::endl;
      for (uint16_t c=0; c < n + 1; ++c) {
        std::cout<<" ["<<c<<"]: "<<c1_len::data.table_ef_bp_space[c];
      }
      std::cout<<std::endl;
    }

  };

} // end namespace
#endif
