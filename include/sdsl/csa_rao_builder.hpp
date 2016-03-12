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
#include <sdsl/int_vector.hpp>
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
		typedef psi_k_support<
			typename t_csa_rao::spec_type::bit_vector,
			typename t_csa_rao::spec_type::r_bit_vector,
			typename t_csa_rao::spec_type::s_bit_vector
		> psi_k_support_type;
		
	protected:
		friend typename t_csa_rao::spec_type::delegate_type;
		
		template <class t_builder, class t_sa_buf>
		class psi_k_support_builder_delegate
		{
		protected:
			int_vector<0> const &m_d_values;
			uint8_t const m_ll{0};
			
		public:
			psi_k_support_builder_delegate(int_vector<0> const &d_values, uint8_t ll):
				m_d_values(d_values),
				m_ll(ll)
			{
			}
			
			int_vector<0>::size_type stored_count(t_builder &builder, uint32_t partition);
			bool psi_k(
				t_builder &builder,
				uint32_t partition,
				uint64_t i,
				typename psi_k_index<t_sa_buf>::value_type &psi_k,
				typename t_builder::text_range &j
			);
		};
	
	protected:
		delegate_type m_delegate;
		
		// Builder lifetime is shorter than any of the two below.
		t_csa_rao &m_csa;
		cache_config &m_config;
		
		int_vector<t_csa_rao::alphabet_type::int_width> m_text_buf;
		int_vector<> m_sa_buf;
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
		
		decltype(m_text_buf) const &text_buf() const { return m_text_buf; }
		decltype(m_sa_buf) const &sa_buf() const { return m_sa_buf; }
	
		void build();
		
		// Called by build.
		void create_alphabet();
		void read_text_and_sa();
		void check_parameters(std::size_t const n);
		void compress_sa();
		
		template<class t_sa_buf_type>
		void compress_level(uint8_t const ll, t_sa_buf_type &sa_buf);
		
	protected:
		template <typename t_vec>
		static void open_and_read_vector_size(
			t_vec const &, isfstream &stream, std::string const &file_name, std::size_t &size, uint8_t &width
		);
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
			read_text_and_sa();
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
	template <typename t_vec>
	void csa_rao_builder<t_csa_rao>::open_and_read_vector_size(
		t_vec const &, isfstream &stream, std::string const &file_name, std::size_t &size, uint8_t &width
	)
	{
		auto status(open_file(stream, file_name));
		if (!status)
			throw std::runtime_error("Unable to open '" + file_name + "'");
		
		uint64_t file_size_bits(0);
		t_vec::read_header(file_size_bits, width, stream);
		auto const fixed_width(t_vec::fixed_int_width);
		if (fixed_width)
			width = fixed_width;
		size = file_size_bits / width;
	}
	
	
	template<class t_csa_rao>
	void csa_rao_builder<t_csa_rao>::read_text_and_sa()
	{
		auto const KEY_SA(conf::KEY_SA);
		auto const KEY_TEXT(key_text_trait<t_csa_rao::alphabet_category::WIDTH>::KEY_TEXT);
	
		assert(cache_file_exists(KEY_SA, m_config));
		assert(cache_file_exists(KEY_TEXT, m_config));
		
		// Check the text size.
		isfstream text_stream;
		isfstream sa_stream;
		std::size_t text_size(0);
		std::size_t sa_size(0);
		uint8_t text_width(0);
		uint8_t sa_width(0);
		
		open_and_read_vector_size(m_text_buf, text_stream, cache_file_name(KEY_TEXT, m_config), text_size, text_width);
		
		// Calculate the padding and load the text and SA.
		check_parameters(text_size);
		auto padding(m_csa.padding());
		
		{
			text_stream.seekg(0);
			m_text_buf.load(text_stream, text_width * padding);
			text_stream.close();
			
			if (padding)
				std::fill(m_text_buf.end() - padding, m_text_buf.end(), 0);
		}
		
		open_and_read_vector_size(m_sa_buf, sa_stream, cache_file_name(KEY_SA, m_config), sa_size, sa_width);

		{
			sa_stream.seekg(0);
			m_sa_buf.load(sa_stream, sa_width * padding);
			sa_stream.close();
			
			if (padding)
			{
				auto const sa_max(text_size - 1 + padding);
				if (sa_max <= m_sa_buf.max_value())
					std::move_backward(m_sa_buf.begin(), m_sa_buf.end() - padding, m_sa_buf.end());
				else
				{
					using std::swap;

					// Calculate the bits needed to represent sa_max.
					auto sa_bits(sa_max);
					
					{
						--sa_bits;
						sa_bits <<= 1;
						sa_bits = 1 + bits::hi(sa_bits);
						// Round up the number of bits to the next power of two.
						--sa_bits;
						sa_bits <<= 1;
						sa_bits = bits::hi(sa_bits);
						sa_bits = 1 << sa_bits;
					}

					int_vector<0> tmp(m_sa_buf.size(), 0, sa_bits);
					swap(tmp, m_sa_buf);
					
					std::copy(tmp.begin(), tmp.end() - padding, m_sa_buf.begin() + padding);
				}

				for (decltype(padding) i(0); i < padding; ++i)
					m_sa_buf[i] = text_size - 1 + (padding - i);
			}
		}
	}
	
	
	template<class t_csa_rao>
	void csa_rao_builder<t_csa_rao>::check_parameters(std::size_t const n)
	{
		// Initialize the remaining instance variables;
		// set suitable values for t (m_level_count) and l (m_partition_count).
		
		if (0 == m_csa.m_level_count)
			m_csa.m_level_count = 1;
		
		if (0 == m_csa.m_partition_count)
		{
			auto const start(std::ceil(std::pow(std::log2(n), 1.0 / (1 + m_csa.m_level_count))));
			m_csa.m_partition_count = util::find_divisor(n, (start < 1 ? 1 : start));
			assert(m_csa.m_partition_count);
		}

		uint64_t const partitions(m_csa.m_partition_count);
		uint64_t const l_t(util::ipow(partitions, m_csa.m_level_count));
		
		// Calculate padding.
		if (n <= l_t)
			m_csa.m_padding = l_t - n;
		else
		{
			auto const rem(n % l_t);
			if (rem)
				m_csa.m_padding = l_t - rem;
			else
				m_csa.m_padding = 0;
		}
		
		// Each item in d is l - (SA[i] mod l) ≤ l (3).
		uint64_t const d_item_bits(std::ceil(std::log2(partitions)));
		// Round up to 8, 16, 32, 64.
		m_d_size = std::max(static_cast<uint8_t>(8), static_cast<uint8_t>(util::upper_power_of_2(d_item_bits)));
	
		typename t_csa_rao::template array<typename t_csa_rao::level> levels;
		levels.reserve(m_csa.m_level_count);
		m_csa.m_levels = std::move(levels);
	}
	
	
	template<class t_csa_rao>
	void csa_rao_builder<t_csa_rao>::compress_sa()
	{
		compress_level(0, m_sa_buf);
		// Both buffers may still be needed for isa_lsw construction.
		for (uint8_t ll(1); ll < m_csa.m_level_count; ++ll)
			compress_level(ll, m_csa.m_sa);
	}
	
	
	template<class t_csa_rao>
	template <class t_builder, class t_sa_buf>
	int_vector<0>::size_type
	csa_rao_builder<t_csa_rao>::psi_k_support_builder_delegate<t_builder, t_sa_buf>::stored_count(
		t_builder &builder, uint32_t partition
	)
	{
		// Count (n / m_partition_count) from the beginning of p. 310; same as the number of items in the Ψ_k subsequences.
		return (builder.sa_buf().size() / builder.csa().m_partition_count);
	}
	
	
	template<class t_csa_rao>
	template <class t_builder, class t_sa_buf>
	bool csa_rao_builder<t_csa_rao>::psi_k_support_builder_delegate<t_builder, t_sa_buf>::psi_k(
		t_builder &builder, uint32_t partition, uint64_t i, typename psi_k_index<t_sa_buf>::value_type &psi_k, typename t_builder::text_range &j
	)
	{
		// Only consider the values that belong to the subsequences {Ψ_k(i) | d[i] = k} (3 (4), p. 310).
		if (m_d_values[i] != partition)
			return false;

		psi_k = builder.psi_k_fn()(partition, i);
		assert(psi_k); // psi_k shouldn't be zero here.

		// Calculate j for L^k_j in Lemma 3, i.e. the value in base-σ of the
		// k * m_partition_count^current_level symbols that appear before SA[Ψ_k(i)].
		uint64_t const partitions(builder.csa().m_partition_count);
		auto const pos(builder.sa_buf()[psi_k - 1]); // psi_k_fn returns 1-based indices.
		auto decompressed_pos(builder.csa().decompress_sa(m_ll, pos));
		auto const nc(partition * util::ipow(partitions, m_ll));
		// decompressed_pos won't be included.
		j = typename t_builder::text_range(builder.text_buf(), decompressed_pos, nc);
		
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
			
				auto const d_val(m_csa.m_partition_count - (val % m_csa.m_partition_count));
				assert(d_val <= d_values.max_value());
				d_values[j] = d_val;
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
