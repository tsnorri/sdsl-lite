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

#ifndef INCLUDED_SDSL_ISA_LSW
#define INCLUDED_SDSL_ISA_LSW

#include <sdsl/elias_inventory.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/iterators.hpp>
#include <sdsl/psi_k_index.hpp>
#include <sdsl/psi_k_support.hpp>
#include <sdsl/psi_k_support_builder.hpp>
#include <sdsl/rrr_vector.hpp>
#include <sdsl/select_support.hpp>


namespace sdsl
{
	template<class t_r_bit_vector, class t_s_bit_vector>
	class isa_lsw_base
	{
	public:
		typedef t_r_bit_vector								r_bit_vector;
		typedef t_s_bit_vector								s_bit_vector;
		typedef psi_k_support<r_bit_vector, s_bit_vector>	psi_k_support_type;
		
	protected:
		int_vector<0> m_isa;
		std::vector<psi_k_support_type> m_psi_k_support;
		
	public:
		isa_lsw_base() = default;
		isa_lsw_base(isa_lsw_base const &) = default;
		isa_lsw_base(isa_lsw_base &&) = default;
		isa_lsw_base &operator=(isa_lsw_base const &) & = default;
		isa_lsw_base &operator=(isa_lsw_base &&) & = default;
	};
	
	
	//! A class for the Inverse Suffix Array (ISA) proposed by T-W. Lam, W-K. Sung and S-S. Wong.
	/*! For use with CSAs that provide the Ψ_k function.
	 *  \tparam	t_csa	CSA class.
	 *  \sa sdsl::csa_rao
	 *  
	 *  \par Reference
	 *  Tak-Wah Lam, Wing-Kin Sung, Swee-Seong Wong:
	 *  Improved Approximate String Matching using Compressed Suffix Data Structures.
	 *  ISAAC 2005: 339–348
	 */
	// TODO: verify time and space complexity.
	template<class t_csa, class t_r_bit_vector = rrr_vector<>, class t_s_bit_vector = bit_vector>
	class isa_lsw : public isa_lsw_base<t_r_bit_vector, t_s_bit_vector>
	{
	public:
		typedef isa_lsw											isa_type;
		typedef random_access_const_iterator<isa_type>			const_iterator;
		typedef typename t_csa::value_type						value_type;
		typedef typename t_csa::size_type						size_type;
		typedef typename t_csa::difference_type					difference_type;
		
		typedef isa_lsw_base<t_r_bit_vector, t_s_bit_vector>	base_class;
		typedef typename base_class::r_bit_vector				r_bit_vector;
		typedef typename base_class::s_bit_vector				s_bit_vector;
		typedef typename base_class::psi_k_support_type			psi_k_support_type;
	
	protected:
		template <class t_builder, class t_sa_buf>
		class psi_k_support_builder_delegate
		{
		protected:
			uint64_t m_l{0};
			uint64_t m_kc_max{0};
		
		public:
			psi_k_support_builder_delegate(uint64_t const l): m_l(l) {}
			void set_kc_max(uint64_t kc_max) { m_kc_max = kc_max; }
			int_vector<0>::size_type stored_count(t_builder &builder, uint32_t partition);
			uint8_t stored_width(t_builder &builder, uint32_t partition);
			
			bool psi_k(
				t_builder &builder,
				uint32_t partition,
				uint64_t i,
				typename psi_k_index<t_sa_buf>::value_type &psi_k,
				uint64_t &j
			);
		};
		
	protected:
		t_csa const &m_csa;		// Pointer to the CSA that provides Ψ_k function.
		
	public:
		isa_lsw() = delete; // A CSA is needed.
	
		isa_lsw(t_csa const &csa, cache_config& config);
	
		isa_lsw(t_csa const &csa):
			base_class::isa_lsw_base(),
			m_csa(csa)
		{
		}
		
		isa_lsw(isa_type const &other):
			base_class::isa_lsw_base(other),
			m_csa(other.m_csa)
		{
		}
		
		isa_lsw(isa_type &&other):
			base_class::isa_lsw_base(std::move(other)),
			m_csa(other.m_csa)
		{
		}
		
		isa_type &operator=(isa_type const &other) &;
		isa_type &operator=(isa_type &&other) &;
		
		// i is the value from isa.
		uint64_t psi_k(uint64_t k, uint64_t i) const;
		
		value_type operator[](size_type i) const;
	
		//! Returns the size of the ISA.
		size_type size() const { return m_csa.size(); }
	
		//! Returns if the ISA is empty.
		size_type empty() const { return m_csa.empty(); }
	
		//! Returns a const_iterator to the first element.
		const_iterator begin() const { return const_iterator(this, 0); }
	
		//! Returns a const_iterator to the element after the last element.
		const_iterator end() const { return const_iterator(this, size()); }

		//! Returns a const_iterator to the first element.
		const_iterator cbegin() const { return const_iterator(this, 0); }
	
		//! Returns a const_iterator to the element after the last element.
		const_iterator cend() const { return const_iterator(this, size()); }
		
		auto serialize(std::ostream& out, structure_tree_node *v = nullptr, std::string name = "") const -> size_type;
		void load(std::istream& in);
	};
	
	
	template<class t_csa, class t_r_bit_vector, class t_s_bit_vector>
	template <class t_builder, class t_sa_buf>
	auto isa_lsw<t_csa, t_r_bit_vector, t_s_bit_vector>::psi_k_support_builder_delegate<t_builder, t_sa_buf>
		::stored_count(t_builder &builder, uint32_t partition) -> int_vector<0>::size_type
	{
		return builder.sa_buf().size() / m_l;
	}
	

	template<class t_csa, class t_r_bit_vector, class t_s_bit_vector>
	template <class t_builder, class t_sa_buf>
	auto isa_lsw<t_csa, t_r_bit_vector, t_s_bit_vector>::psi_k_support_builder_delegate<t_builder, t_sa_buf>
		::stored_width(t_builder &builder, uint32_t partition) -> uint8_t
	{
		return util::upper_power_of_2(std::log2(1 + m_kc_max));
	}
	
	
	template<class t_csa, class t_r_bit_vector, class t_s_bit_vector>
	template <class t_builder, class t_sa_buf>
	auto isa_lsw<t_csa, t_r_bit_vector, t_s_bit_vector>::psi_k_support_builder_delegate<t_builder, t_sa_buf>
		::psi_k(
			t_builder &builder,
			uint32_t partition,
			uint64_t i,
			typename psi_k_index<t_sa_buf>::value_type &psi_k,
			uint64_t &j
		) -> bool
	{
		auto val(builder.sa_buf()[i]);
		if (0 == val % m_l)
		{
			psi_k = builder.psi_k_fn().from_sa_val(partition, val);
			if (psi_k)
			{
				auto pos(builder.csa()[psi_k - 1]);
				assert(partition <= pos);
				j = util::str_to_base_sigma(builder.text_buf(), builder.alphabet(), pos, partition);
				
				assert(j <= m_kc_max);
				
				return true;
			}
		}
		return false;
	}
	
	
	template<class t_csa, class t_r_bit_vector, class t_s_bit_vector>
	isa_lsw<t_csa, t_r_bit_vector, t_s_bit_vector>::isa_lsw(t_csa const &csa, cache_config& config):
		isa_lsw(csa)
	{
		// Access the text and the suffix array. m_csa could be used instead.
		auto const KEY_SA(conf::KEY_SA);
		auto const KEY_TEXT(key_text_trait<t_csa::alphabet_category::WIDTH>::KEY_TEXT);
		
		assert(cache_file_exists(KEY_SA, config));
		assert(cache_file_exists(KEY_TEXT, config));
		std::string const text_file(cache_file_name(KEY_TEXT, config));

		int_vector_buffer<> sa_buf(cache_file_name(KEY_SA, config));
		int_vector_buffer<t_csa::alphabet_type::int_width> text_buf(text_file);
		
		uint64_t const l(m_csa.m_partition_count);
		
		{
			uint64_t const n(sa_buf.size());
			decltype(this->m_isa) isa_tmp(n / l, 0);
			for (uint64_t i(0); i < n; ++i)
			{
				auto val(sa_buf[i]);
				if (0 == val % l)
					isa_tmp[val / l] = i;
			}
			
			this->m_isa = std::move(isa_tmp);
		}

		{
			// Lemma 2.
			sdsl::psi_k_index<decltype(sa_buf)> psi_k_fn(config, sa_buf);
			
			this->m_psi_k_support.reserve(l - 1);
			auto builder(construct_psi_k_support_builder(m_csa, text_buf, sa_buf, m_csa.m_alphabet, psi_k_fn));
			
			psi_k_support_builder_delegate<decltype(builder), decltype(sa_buf)> delegate(l);
			for (uint64_t k(1); k < l; ++k)
			{
				delegate.set_kc_max(util::ipow(m_csa.m_alphabet.sigma, k) - 1);
				
				psi_k_support_type psi_k_support;
				builder.build(psi_k_support, k, delegate);
				this->m_psi_k_support.emplace_back(std::move(psi_k_support));
			}
		}
	}
	
	
	template<class t_csa, class t_r_bit_vector, class t_s_bit_vector>
	auto isa_lsw<t_csa, t_r_bit_vector, t_s_bit_vector>::operator=(isa_type const &other) & -> isa_type &
	{
		// m_csa must have already been set since it needs to be given in the constructor.
		base_class::operator=(other);
		return *this;
	}
	
	
	template<class t_csa, class t_r_bit_vector, class t_s_bit_vector>
	auto isa_lsw<t_csa, t_r_bit_vector, t_s_bit_vector>::operator=(isa_type &&other) & -> isa_type &
	{
		// m_csa must have already been set since it needs to be given in the constructor.
		base_class::operator=(std::move(other));
		return *this;
	}
	
	
	// i is the value from isa.
	template<class t_csa, class t_r_bit_vector, class t_s_bit_vector>
	uint64_t isa_lsw<t_csa, t_r_bit_vector, t_s_bit_vector>::psi_k(uint64_t k, uint64_t i) const
	{
		if (k == 0)
			return i;
		
		psi_k_support_type const &partition(this->m_psi_k_support[k - 1]);
		auto retval(partition[i]);
		return retval;
	}
	
	
	template<class t_csa, class t_r_bit_vector, class t_s_bit_vector>
	auto isa_lsw<t_csa, t_r_bit_vector, t_s_bit_vector>::operator[](size_type i) const -> value_type
	{
		uint64_t const l(m_csa.m_partition_count);
		size_type y(i / l);
		size_type yl(y * l);
		size_type k(i - yl);
		value_type z(this->m_isa[y]);
		value_type retval(psi_k(k, z));
		return retval;
	}
	
	
	template<class t_csa, class t_r_bit_vector, class t_s_bit_vector>
	auto isa_lsw<t_csa, t_r_bit_vector, t_s_bit_vector>::serialize(std::ostream& out, structure_tree_node *v, std::string name) const -> size_type
	{
		structure_tree_node *child(structure_tree::add_child(v, name, util::class_name(*this)));
		size_type written_bytes(0);
		auto psi_k_count(this->m_psi_k_support.size());

		written_bytes += this->m_isa.serialize(out, child, "m_isa");
		written_bytes += write_member(psi_k_count, out, child, "psi_k_count");
		
		typename decltype(this->m_psi_k_support)::size_type i(1);
		for (auto it(this->m_psi_k_support.cbegin()), end(this->m_psi_k_support.cend()); it != end; ++it)
		{
			written_bytes += it->serialize(out, child, "psi_k_" + std::to_string(i));
			++i;
		}
		
		structure_tree::add_size(child, written_bytes);
		return written_bytes;
	}


	template<class t_csa, class t_r_bit_vector, class t_s_bit_vector>
	void isa_lsw<t_csa, t_r_bit_vector, t_s_bit_vector>::load(std::istream& in)
	{
		this->m_isa.load(in);
		
		typename decltype(this->m_psi_k_support)::size_type psi_k_count(0);
		read_member(psi_k_count, in);
		
		this->m_psi_k_support.resize(psi_k_count);
		for (decltype(psi_k_count) i(0); i < psi_k_count; ++i)
			this->m_psi_k_support[i].load(in);
	}
}

#endif
