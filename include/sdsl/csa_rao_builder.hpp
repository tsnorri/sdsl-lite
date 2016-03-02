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

#ifndef INCLUDED_SDSL_CSA_RAO_BUILDER
#define INCLUDED_SDSL_CSA_RAO_BUILDER

#include <sdsl/csa_alphabet_strategy.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include <sdsl/psi_k_index.hpp>


namespace sdsl
{
	//! A builder class for sdsl::csa_rao.
	/*! \tparam	t_csa_rao	The CSA class.
	 *  \sa sdsa::csa_rao
	 *  
	 *  \par Reference
	 *  S. Srinivasa Rao:
	 *  Time-space trade-offs for compressed suffix arrays.
	 *  Information Processing Letters 82(6): 307–311 (2002)
	 *  
	 */
	// The comments refer to Rao's article.
	// TODO: verify time and space complexity.
	template<class t_csa_rao>
	class csa_rao_builder
	{
	public:
		typedef csa_rao_builder<t_csa_rao> builder_type;
		typedef typename t_csa_rao::spec_type::delegate_type delegate_type;
		typedef psi_k_support<typename t_csa_rao::spec_type::r_bit_vector, typename t_csa_rao::spec_type::s_bit_vector> psi_k_support_type;
		
	protected:
		friend typename t_csa_rao::spec_type::delegate_type;
		
		template <class t_builder, class t_sa_buf>
		class psi_k_support_builder_delegate
		{
		protected:
			int_vector<0> const &m_d_values;
			uint8_t const m_ll{0};
			
			// Maximum value for the k * m_partition_count^current_level symbols in the pattern in base-σ.
			uint64_t m_kc_max{0};
			
		public:
			psi_k_support_builder_delegate(int_vector<0> const &d_values, uint8_t ll):
				m_d_values(d_values),
				m_ll(ll)
			{
			}
			
			void set_kc_max(uint64_t kc_max) { m_kc_max = kc_max; }
			int_vector<0>::size_type stored_count(t_builder &builder, uint32_t partition);
	 		uint8_t stored_width(t_builder &builder, uint32_t partition);
			bool psi_k(t_builder &builder, uint32_t partition, uint64_t i, typename psi_k_index<t_sa_buf>::value_type &psi_k, uint64_t &j);
		};
	
	protected:
		delegate_type m_delegate;
		
		// Builder lifetime is shorter than any of the two below.
		t_csa_rao &m_csa;
		cache_config &m_config;
		
		int_vector_buffer<t_csa_rao::alphabet_type::int_width> m_text_buf;
		int_vector_buffer<> m_sa_buf;
		uint64_t m_d_size{0};
		
	public:
		csa_rao_builder(t_csa_rao &csa, cache_config &config):
			m_csa(csa),
			m_config(config)
		{
			m_delegate.setup(m_csa, *this);
		}
		
		csa_rao_builder(builder_type const &other) = delete;
		csa_rao_builder(builder_type &&other) = default;
		csa_rao_builder &operator=(csa_rao_builder const &) & = delete;
		csa_rao_builder &operator=(csa_rao_builder &&) & = default;
	
		void build();
		
		// Called by build.
		void create_alphabet();
		void read_sa();
		void check_parameters();
		void compress_sa();
		
		template<class t_sa_buf_type>
		void compress_level(uint8_t const ll, t_sa_buf_type &sa_buf);
	};
	

	template<class t_csa_rao>
	void csa_rao_builder<t_csa_rao>::build()
	{
		{
			auto event(memory_monitor::event("construct csa-alphabet"));
			create_alphabet();
		}
		
		{
			auto event(memory_monitor::event("compress SA"));
			read_sa();
			check_parameters();
			compress_sa();
		}
	}

	
	template<class t_csa_rao>
	void csa_rao_builder<t_csa_rao>::create_alphabet()
	{
		auto const KEY_BWT(key_trait<t_csa_rao::alphabet_type::int_width>::KEY_BWT);
		assert(cache_file_exists(KEY_BWT, m_config));
	
		int_vector_buffer<t_csa_rao::alphabet_type::int_width> bwt_buf(cache_file_name(KEY_BWT, m_config));
		typename t_csa_rao::size_type n(bwt_buf.size());
		typename t_csa_rao::alphabet_type tmp_alphabet(bwt_buf, n);
		m_csa.m_alphabet.swap(tmp_alphabet);
	}
	
	
	template<class t_csa_rao>
	void csa_rao_builder<t_csa_rao>::read_sa()
	{
		auto const KEY_SA(conf::KEY_SA);
		auto const KEY_TEXT(key_text_trait<t_csa_rao::alphabet_category::WIDTH>::KEY_TEXT);
	
		assert(cache_file_exists(KEY_SA, m_config));
		assert(cache_file_exists(KEY_TEXT, m_config));
		
		// Create the suffix array.
		typedef int_vector<t_csa_rao::alphabet_category::WIDTH> text_type;
		
		std::string const text_file(cache_file_name(KEY_TEXT, m_config));
		int_vector_buffer<t_csa_rao::alphabet_type::int_width> text_buf_tmp(text_file);
		int_vector_buffer<> sa_buf_tmp(cache_file_name(KEY_SA, m_config));
		
		m_text_buf = std::move(text_buf_tmp);
		
		// FIXME: this should be memory-mapped, too?
		m_sa_buf = std::move(sa_buf_tmp);
		m_sa_buf.buffersize(8 * 1024 * 1024);
	}
	
	
	template<class t_csa_rao>
	void csa_rao_builder<t_csa_rao>::check_parameters()
	{
		// Initialize the remaining instance variables;
		// set suitable values for t (m_level_count) and l (m_partition_count).
		
		auto const n(m_text_buf.size());
	
		if (0 == m_csa.m_level_count)
			m_csa.m_level_count = 1;
		
		if (0 == m_csa.m_partition_count)
		{
			m_csa.m_partition_count = util::find_divisor(
				n,
				std::ceil(std::pow(std::log2(n), 1.0 / (1 + m_csa.m_level_count)))
			);
		}
		
		uint64_t const l_t(util::ipow(m_csa.m_partition_count, m_csa.m_level_count));
		m_delegate.check_l_t(*this, n, l_t);
		
		// Assume that the size of the text is a multiple of partitions^levels.
		// (We can't easily make sure of this by e.g. adding ending characters in sdsl::construct.)
		assert(! (0 == m_csa.m_level_count || 0 == m_csa.m_partition_count));
		assert(0 == n % l_t);
	
		auto const d_item_bits(std::log2(m_csa.m_partition_count)); // Each item in d is l - (SA[i] mod l) ≤ l (3).
		m_d_size = std::max(static_cast<uint8_t>(8), static_cast<uint8_t>(util::upper_power_of_2(d_item_bits))); // Round up to 8, 16, 32, 64.
	
		typename t_csa_rao::template array<typename t_csa_rao::level> levels;
		levels.reserve(m_csa.m_level_count);
		m_csa.m_levels = std::move(levels);
	}
	
	
	template<class t_csa_rao>
	void csa_rao_builder<t_csa_rao>::compress_sa()
	{
		compress_level(0, m_sa_buf);
		for (uint8_t ll(1); ll < m_csa.m_level_count; ++ll)
			compress_level(ll, m_csa.m_sa);
	}
	
	
	template<class t_csa_rao>
	template <class t_builder, class t_sa_buf>
	int_vector<0>::size_type
	csa_rao_builder<t_csa_rao>::psi_k_support_builder_delegate<t_builder, t_sa_buf>
		::stored_count(t_builder &builder, uint32_t partition)
	{
		// Count (n / m_partition_count) from the beginning of p. 310; same as the number of items in the Ψ_k subsequences.
		return (builder.sa_buf().size() / builder.csa().m_partition_count);
	}
	
	
	template<class t_csa_rao>
	template <class t_builder, class t_sa_buf>
	uint8_t
	csa_rao_builder<t_csa_rao>::psi_k_support_builder_delegate<t_builder, t_sa_buf>
		::stored_width(t_builder &builder, uint32_t partition)
	{
		return util::upper_power_of_2(std::log2(1 + m_kc_max));
	}
	
	
	template<class t_csa_rao>
	template <class t_builder, class t_sa_buf>
	bool csa_rao_builder<t_csa_rao>::psi_k_support_builder_delegate<t_builder, t_sa_buf>
		::psi_k(t_builder &builder, uint32_t partition, uint64_t i, typename psi_k_index<t_sa_buf>::value_type &psi_k, uint64_t &j)
	{
		// Only consider the values that belong to the subsequences {Ψ_k(i) | d[i] = k} (3 (4), p. 310).
		if (m_d_values[i] != partition)
			return false;

		psi_k = builder.psi_k_fn()(partition, i);
		assert(psi_k); // psi_k shouldn't be zero here.

		// Calculate j for L^k_j in Lemma 3, i.e. the value in base-σ of the
		// k * m_partition_count^current_level symbols that appear before SA[Ψ_k(i)].
		auto const pos(builder.sa_buf()[psi_k - 1]); // psi_k_fn returns 1-based indices.
		auto decompressed_pos(builder.csa().decompress_sa(m_ll, pos));
		auto const nc(partition * util::ipow(builder.csa().m_partition_count, m_ll));
		// decompressed_pos won't be included.
		j = util::str_to_base_sigma(builder.text_buf(), builder.csa().m_alphabet, decompressed_pos, nc);
		
		assert(j <= m_kc_max);
		
		return true;
	}
	
	
	template<class t_csa_rao>
	template<class t_sa_buf_type>
	void csa_rao_builder<t_csa_rao>::compress_level(uint8_t const ll, t_sa_buf_type &sa_buf)
	{
		m_delegate.start_level(*this, ll, sa_buf);
		assert(0 < m_csa.m_partition_count);
		
		typename t_csa_rao::template array<typename t_csa_rao::psi_k_support_type> partitions;
		partitions.reserve(m_csa.m_partition_count);
	
		// Specialization of psi_k for the current level.
		psi_k_index<t_sa_buf_type> psi_k_fn(m_config, sa_buf);
		
		uint64_t const n(sa_buf.size());
		std::vector<typename t_sa_buf_type::value_type> sa_ll; // SA_n, calculated with 1-based indices.
	
		bit_vector b_values(n, 0);
		int_vector<0> d_values(n, 0, m_d_size);	// l − (SA[i] mod l), indices 0-based, SA[i] 1-based (3 (3)).
		
		// (3 (1)) Store values from SA divisible by m_partition_count.
		// (3 (2)) Create the B vector.
		// (3 (3)) Create the d array.
		{
			uint64_t j(0);
			for (auto it(sa_buf.cbegin()), end(sa_buf.cend()); it != end; ++it)
			{
				// Change to 1-based.
				auto const val(*it + 1);
				if (0 == (val % m_csa.m_partition_count))
				{
					auto div(val / m_csa.m_partition_count);
					assert(div);
					assert(div - 1 <= std::numeric_limits<typename decltype(sa_ll)::value_type>::max());
					sa_ll.push_back(div - 1);
					b_values[j] = 1;
				}
			
				d_values[j] = m_csa.m_partition_count - (val % m_csa.m_partition_count);
				++j;
			}
		}
		
		// (3 (4)) Create Ψ_k and compress it.
		// m_sa may still be used via sa_buf.
		{
			auto builder(construct_psi_k_support_builder(m_csa, m_text_buf, sa_buf, m_csa.m_alphabet, psi_k_fn));
			psi_k_support_builder_delegate<decltype(builder), decltype(sa_buf)> delegate(d_values, ll);
		
			for (uint64_t k(1); k < m_csa.m_partition_count; ++k)
			{
				// Maximum value for the k * m_partition_count^current_level symbols in base-σ.
				// (Similar to the maximum of “the value in binary of the k symbols”).
				delegate.set_kc_max(util::ipow(
					m_csa.m_alphabet.sigma,
					k * util::ipow(m_csa.m_partition_count, ll)
				) - 1);
			
				psi_k_support_type psi_k_support;
				builder.build(psi_k_support, k, delegate);
			
				// Debugging helper.
				m_delegate.finish_create_psi_k(*this, ll, k, d_values.size(), psi_k_fn);
			
				assert(partitions.size() == k - 1);
				partitions.emplace_back(std::move(psi_k_support));
			}
		}
		
		assert(partitions.size() == m_csa.m_partition_count - 1);
		typename t_csa_rao::level level(partitions, b_values, d_values);
	
		m_csa.m_sa.resize(sa_ll.size()); // Shrinks.
		std::copy(sa_ll.cbegin(), sa_ll.cend(), m_csa.m_sa.begin());
		assert(m_csa.m_levels.size() == ll);
		m_csa.m_levels.emplace_back(std::move(level));
	}
}

#endif
