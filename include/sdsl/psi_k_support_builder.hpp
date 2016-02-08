/*
 Copyright (c) 2015 Tuukka Norri
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

#ifndef INCLUDED_SDSL_PSI_K_SUPPORT_BUILDER
#define INCLUDED_SDSL_PSI_K_SUPPORT_BUILDER

#include <cstdint>
#include <sdsl/elias_inventory.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/psi_k_index.hpp>


namespace sdsl
{
	//! Provides boilerplate for instantiating psi_k_support.
	/*! \tparam	t_csa		The CSA class.
	 *  \tparam t_text_buf	Text buffer class.
	 *  \tparam t_sa_buf	SA buffer class.
	 *  \tparam t_alphabet	Policy for alphabet representation.
	 *  \tparam t_psi_k_fn	A class that provides the Ψ_k values with the function call operator.
	 *  \sa sdsa::psi_k_support
	 *  
	 *  \par Reference
	 *  S. Srinivasa Rao:
	 *  Time-space trade-offs for compressed suffix arrays.
	 *  Information Processing Letters 82(6): 307–311 (2002)
	 *  
	 */
	// The comments refer to Rao's article.
	/*
		template <class t_builder>
		class delegate
		{
		public:
			int_vector<0>::size_type stored_count(t_builder &builder, uint32_t partition) { return 0; }
	 		uint8_t stored_width(t_builder &builder, uint32_t partition) { return 0; }
			bool psi_k(t_builder &builder, uint32_t partition, uint64_t i, typename psi_k_index<t_sa_buf>::value_type &psi_k, uint64_t &j) { return false; }
		};
	*/
	// TODO: verify time and space complexity.
	template<
		class t_csa,
		class t_text_buf,
		class t_sa_buf,
		class t_alphabet,
		class t_psi_k_fn
	>
	class psi_k_support_builder
	{
	public:
		typedef psi_k_support_builder<t_csa, t_text_buf, t_sa_buf, t_alphabet, t_psi_k_fn> builder_type;
		
	public:
		t_csa const &m_csa;
		t_text_buf &m_text_buf;
		t_sa_buf &m_sa_buf;
		t_alphabet const &m_alphabet;
		t_psi_k_fn &m_psi_k_fn;
		
	public:
		t_csa const			&csa() const { return m_csa; }
		t_text_buf			&text_buf() { return m_text_buf; }
		t_sa_buf			&sa_buf() { return m_sa_buf; }
		t_alphabet const	&alphabet() { return m_alphabet; }
		t_psi_k_fn			&psi_k_fn() { return m_psi_k_fn; }
		
		psi_k_support_builder() = delete;
		psi_k_support_builder(psi_k_support_builder const &) = delete;
		psi_k_support_builder(psi_k_support_builder &&) = default;
		psi_k_support_builder &operator=(psi_k_support_builder const &) = delete;
		psi_k_support_builder &operator=(psi_k_support_builder &&) = default;
		
		
		psi_k_support_builder(
			t_csa const &csa,
			t_text_buf &text_buf,
			t_sa_buf &sa_buf,
			t_alphabet const &alphabet,
			t_psi_k_fn &psi_k_fn
		):
			m_csa(csa),
			m_text_buf(text_buf),
			m_sa_buf(sa_buf),
			m_alphabet(alphabet),
			m_psi_k_fn(psi_k_fn)
		{
		}
		
		
		template <class t_psi_k_support, class t_delegate>
		void build(
			t_psi_k_support &psi_k_support,
			uint32_t const partition,
			t_delegate &delegate
		);
	};
		
	
	template<
		class t_csa,
		class t_text_buf,
		class t_sa_buf,
		class t_alphabet,
		class t_psi_k_fn
	>
	template <
		class t_psi_k_support,
		class t_delegate
	>
	void psi_k_support_builder<t_csa, t_text_buf, t_sa_buf, t_alphabet, t_psi_k_fn>::build(
		t_psi_k_support &psi_k_support,
		uint32_t const partition,
		t_delegate &delegate
	)
	{
		// (3 (4.2)) Create the L^k_j lists and L_k.
		// Use 0-based indexing.
		// The length of each L^k_j isn't known beforehand (although the sum of the lengths is).
		// FIXME: could the process be streamlined?
		
		int_vector<0>::size_type const stored_count(delegate.stored_count(*this, partition));
		
		// (3 (4.1)) Create V_k. Use 0-based indexing for i.
		bit_vector v_values(m_sa_buf.size(), 0);
		
		// The list number to which the jth element in the subsequence belongs to (3 (4.2)), 0-based.
		int_vector<0> l_k_values(stored_count, 0);

		// Cumulative sums of the counts of elements in the L^k lists, 1-based (3 (4.3)).
		int_vector<64> c_k_values;
		
		// Contents of each L^k_j (Ψ_k) with the order of j (but not j itself) preserved.
		// Temporary, we combine the values into psi_k_values.
		std::vector<std::vector<uint64_t>> l_vec;
		
		{
			// Create the lists for (3 (4.2)).
			// Assume that the numbering of the lists (j) is sparse.
			// This could affect the time complexity aimed for, though,
			// which is O(n(log n)^(1/2) log(sigma)), since the last factor becomes log(n/l).
			// Solve the problem by re-mapping the indices to consecutive unsigned integers.
			
			// The L lists by list index j.
			std::map<uint64_t, std::vector<uint64_t>> l_map_tmp;
			
			// The i values by list index j, i.e. which L list stores the Ψ_k(i).
			std::unordered_multimap<uint64_t, uint64_t> l_k_values_tmp;
			
			for (uint64_t i(0), count(m_sa_buf.size()); i < count; ++i)
			{
				typename psi_k_index<t_sa_buf>::value_type psi_k_i(0);
				uint64_t j(0);
				if (delegate.psi_k(*this, partition, i, psi_k_i, j))
				{
					v_values[i] = 1;
					
					// Insert Ψ_k(i) into L^k_j. The lists end up sorted (Lemma 3).
					// Also store i into an inverted index.
					l_map_tmp[j].push_back(psi_k_i);
					l_k_values_tmp.emplace(j, i);
				}
			}
			
			// Re-map and invert the l_k_values_tmp index.
			// Since the lists are numbered non-consequtively, we also re-map their indices to be able to use C_k
			// (as only the sizes of the non-empty lists are stored).
			// XXX: this remapping is not in Rao's paper but should be obvious.
			uint64_t i(0);
			uint64_t ii(0);
			uint64_t prev_j(0);
			uint64_t rep_j(0);
			for (auto l_it(l_map_tmp.begin()), l_end(l_map_tmp.end()); l_it != l_end; ++l_it)
			{
				// l_it->first gives the list index j.
				auto const j_val(l_it->first);
				
				// Replace j_val. Handle the case in which j_val = 0.
				if (1 + j_val != prev_j)
				{
					prev_j = 1 + j_val;
					++rep_j;
				}
				
				l_vec.emplace_back(std::move(l_it->second)); // l_it->second gives a vector of values of Ψ_k.
				auto range(l_k_values_tmp.equal_range(j_val)); // Returns a range containing all elements with the given key in the container.
				for (auto kv_it(range.first); kv_it != range.second; ++kv_it)
				{
					// kv_it->second gives the value of i. Consequtive indices are needed, though.
					l_k_values[ii] = rep_j - 1; // rep_j starts from 1.
					++ii;
				}
				
				// Prevent re-inserting the values into l_k_values.
				// (Erase removes all values for a given key.)
				l_k_values_tmp.erase(j_val);
				++i;
			}
		}
		
		// (3 (4.3)) Create C_k. Use 1-based indexing for i.
		c_k_values = int_vector<64>(1 + l_vec.size(), 0);
		
		{
			uint64_t sum(0);
			for (uint64_t i(0), count(l_vec.size()); i < count; ++i)
			{
				sum += l_vec[i].size();
				c_k_values[1 + i] = sum;
			}
		}
		
		// (3 (4.4) and corollary 2) Compact representation of Ψ_k.
		elias_inventory<typename t_psi_k_support::s_bit_vector> psi_k_values(l_vec, c_k_values);
		
		t_psi_k_support tmp_support(v_values, l_k_values, c_k_values, psi_k_values);
		psi_k_support = std::move(tmp_support);
	}
	
	
	template<
		class t_csa,
		class t_text_buf,
		class t_sa_buf,
		class t_alphabet,
		class t_psi_k_fn
	>
	auto construct_psi_k_support_builder(
		t_csa const &csa,
		t_text_buf &text_buf,
		t_sa_buf &sa_buf,
		t_alphabet const &alphabet,
		t_psi_k_fn &psi_k_fn
	) -> psi_k_support_builder<t_csa, t_text_buf, t_sa_buf, t_alphabet, t_psi_k_fn>
	{
		psi_k_support_builder<
			t_csa,
			t_text_buf,
			t_sa_buf,
			t_alphabet,
			t_psi_k_fn
		> retval(csa, text_buf, sa_buf, alphabet, psi_k_fn);
		return retval;
	}
}

#endif
