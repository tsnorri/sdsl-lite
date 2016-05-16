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
#include <sdsl/isa_simple.hpp>
#include <sdsl/psi_k_support.hpp>
#include <sdsl/psi_k_support_builder.hpp>
#include <sdsl/select_support.hpp>


namespace sdsl
{
	template<class t_bit_vector, class t_r_bit_vector, class t_s_bit_vector>
	class isa_lsw_base
	{
	public:
		typedef t_r_bit_vector								r_bit_vector;
		typedef t_s_bit_vector								s_bit_vector;
		typedef psi_k_support<
			t_bit_vector, r_bit_vector, s_bit_vector
		> psi_k_support_type;
		
	protected:
		int_vector<0> m_isa;
		std::vector<psi_k_support_type> m_psi_k_support;
		uint64_t m_l{0};
		
	public:
		isa_lsw_base() = default;
		isa_lsw_base(isa_lsw_base const &) = default;
		isa_lsw_base(isa_lsw_base &&) = default;
		isa_lsw_base &operator=(isa_lsw_base const &) & = default;
		isa_lsw_base &operator=(isa_lsw_base &&) & = default;
	};
	
	
	//! A class for the Inverse Suffix Array (ISA) proposed by T-W. Lam, W-K. Sung and S-S. Wong.
	/*! For use with CSAs that provide the Ψ^k function.
	 *  \tparam	t_csa	CSA class.
	 *  \sa sdsl::csa_rao
	 *  
	 *  \par Reference
	 *  Tak-Wah Lam, Wing-Kin Sung, Swee-Seong Wong:
	 *  Improved Approximate String Matching using Compressed Suffix Data Structures.
	 *  ISAAC 2005: 339–348
	 */
	// TODO: verify time and space complexity.
	template<class t_csa, class t_bit_vector, class t_r_bit_vector, class t_s_bit_vector>
	class isa_lsw : public isa_lsw_base<t_bit_vector, t_r_bit_vector, t_s_bit_vector>
	{
	public:
		typedef isa_lsw											isa_type;
		typedef random_access_const_iterator<isa_type>			const_iterator;
		typedef typename t_csa::value_type						value_type;
		typedef typename t_csa::size_type						size_type;
		typedef typename t_csa::difference_type					difference_type;
		
		typedef isa_lsw_base<
			t_bit_vector, t_r_bit_vector, t_s_bit_vector
		> base_class;
		typedef typename base_class::r_bit_vector				r_bit_vector;
		typedef typename base_class::s_bit_vector				s_bit_vector;
		typedef typename base_class::psi_k_support_type			psi_k_support_type;
	
	protected:
		template <class t_builder, class t_sa_buf>
		class psi_k_support_builder_delegate
		{
		protected:
			uint64_t m_l{0};
		
		public:
			psi_k_support_builder_delegate(uint64_t const l): m_l(l) {}
			int_vector<0>::size_type stored_count(t_builder &builder, uint32_t partition);
			
			bool psi_k(
				t_builder &builder,
				uint32_t partition,
				uint64_t i,
				typename isa_simple<t_sa_buf>::value_type &psi_k,
				typename t_builder::text_range &j
			);
		};
		
	protected:
		t_csa const &m_csa;		// Pointer to the CSA that provides Ψ^k function.
		
	public:
		isa_lsw() = delete; // A CSA is needed.
	
		isa_lsw(t_csa const &csa, csa_rao_builder<t_csa> const &builder, cache_config& config);
	
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
		
		value_type operator[](size_type i) const SDSL_HOT;
	
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
	
	
	template<class t_csa, class t_bit_vector, class t_r_bit_vector, class t_s_bit_vector>
	template <class t_builder, class t_sa_buf>
	auto isa_lsw<t_csa, t_bit_vector, t_r_bit_vector, t_s_bit_vector>::psi_k_support_builder_delegate<t_builder, t_sa_buf>::stored_count(
		t_builder &builder, uint32_t partition
	) -> int_vector<0>::size_type
	{
		uint64_t const n(builder.sa_buf().size());
		auto const isa_size(1 + uint64_t(std::floor(double(n) / m_l)));
		return isa_size;
	}
	

	template<class t_csa, class t_bit_vector, class t_r_bit_vector, class t_s_bit_vector>
	template <class t_builder, class t_sa_buf>
	bool isa_lsw<t_csa, t_bit_vector, t_r_bit_vector, t_s_bit_vector>::psi_k_support_builder_delegate<t_builder, t_sa_buf>::psi_k(
		t_builder &builder,
		uint32_t partition,
		uint64_t i,
		typename isa_simple<t_sa_buf>::value_type &psi_k,
		typename t_builder::text_range &j
	)
	{
		auto val(builder.sa_buf()[i]);
		if (0 == val % m_l)
		{
			psi_k = builder.isa().psi_k_from_sa_val(partition, val);
			if (psi_k)
			{
				auto pos(builder.csa().sa(psi_k - 1, 0));
				assert(partition <= pos);
				j = typename t_builder::text_range(builder.text_buf(), pos, partition);
				return true;
			}
		}
		return false;
	}
	
	
	template<class t_csa, class t_bit_vector, class t_r_bit_vector, class t_s_bit_vector>
	isa_lsw<t_csa, t_bit_vector, t_r_bit_vector, t_s_bit_vector>::isa_lsw(
		t_csa const &csa, csa_rao_builder<t_csa> const &builder, cache_config &config
	):
		isa_lsw(csa)
	{
		// Access the text and the suffix array. m_csa could be used instead.
		int_vector<> sa_buf(builder.sa_buf());
		int_vector<t_csa::alphabet_type::int_width> text_buf(builder.text_buf());
		uint64_t const n(sa_buf.size());
		
		{
			auto const log_n(util::log2_ceil(n));
			uint64_t const l(std::floor(std::sqrt(log_n)));
			this->m_l = std::max(uint64_t(1), l);
		}

		{
			// Find the maximum value to be stored. (Or just take log_2(n).)
			uint64_t max_val(0);
			for (uint64_t i(0); i < n; ++i)
			{
				auto const val(sa_buf[i]);
				if (0 == val % this->m_l)
					max_val = std::max(max_val, i);
			}
			
			auto const isa_size(1 + uint64_t(std::floor(double(n) / this->m_l)));
			decltype(this->m_isa) isa_tmp(isa_size, 0, util::log2_ceil(1 + max_val));
			
			// Create a sample of the inverse suffix array (Lemma 7).
			for (uint64_t i(0); i < n; ++i)
			{
				auto const val(sa_buf[i]);
				if (0 == val % this->m_l)
				{
					assert(i <= isa_tmp.max_value());
					isa_tmp[val / this->m_l] = i;
				}
			}
			
			this->m_isa = std::move(isa_tmp);
		}

		{
			// Compress the required Ψ^k values (Lemma 7 and 2).
			isa_simple<decltype(sa_buf)> isa(config, sa_buf);
			
			this->m_psi_k_support.reserve(this->m_l - 1);
			auto builder(construct_psi_k_support_builder(m_csa, text_buf, sa_buf, m_csa.m_alphabet, isa));
			
			psi_k_support_builder_delegate<decltype(builder), decltype(sa_buf)> delegate(this->m_l);
			for (uint64_t k(1); k < this->m_l; ++k)
			{
				psi_k_support_type psi_k_support;
				builder.build(psi_k_support, k, delegate);
				this->m_psi_k_support.emplace_back(std::move(psi_k_support));
			}
		}
	}
	
	
	template<class t_csa, class t_bit_vector, class t_r_bit_vector, class t_s_bit_vector>
	auto isa_lsw<t_csa, t_bit_vector, t_r_bit_vector, t_s_bit_vector>::operator=(isa_type const &other) & -> isa_type &
	{
		// m_csa must have already been set since it needs to be given in the constructor.
		base_class::operator=(other);
		return *this;
	}
	
	
	template<class t_csa, class t_bit_vector, class t_r_bit_vector, class t_s_bit_vector>
	auto isa_lsw<t_csa, t_bit_vector, t_r_bit_vector, t_s_bit_vector>::operator=(isa_type &&other) & -> isa_type &
	{
		// m_csa must have already been set since it needs to be given in the constructor.
		base_class::operator=(std::move(other));
		return *this;
	}
	
	
	// i is the value from isa.
	template<class t_csa, class t_bit_vector, class t_r_bit_vector, class t_s_bit_vector>
	uint64_t isa_lsw<t_csa, t_bit_vector, t_r_bit_vector, t_s_bit_vector>::psi_k(uint64_t k, uint64_t i) const
	{
		if (k == 0)
			return i;
		
		psi_k_support_type const &partition(this->m_psi_k_support[k - 1]);
		auto retval(partition[i]);
		return retval;
	}
	
	
	template<class t_csa, class t_bit_vector, class t_r_bit_vector, class t_s_bit_vector>
	auto isa_lsw<t_csa, t_bit_vector, t_r_bit_vector, t_s_bit_vector>::operator[](size_type i) const -> value_type
	{
		size_type const y(i / this->m_l);
		size_type const yl(y * this->m_l);
		size_type const k(i - yl);
		value_type const z(this->m_isa[y]);
		value_type const psi_val(psi_k(k, z));
		value_type const retval(psi_val - m_csa.padding());
		return retval;
	}
	
	
	template<class t_csa, class t_bit_vector, class t_r_bit_vector, class t_s_bit_vector>
	auto isa_lsw<t_csa, t_bit_vector, t_r_bit_vector, t_s_bit_vector>::serialize(std::ostream& out, structure_tree_node *v, std::string name) const -> size_type
	{
		structure_tree_node *child(structure_tree::add_child(v, name, util::class_name(*this)));
		size_type written_bytes(0);
		auto psi_k_count(this->m_psi_k_support.size());

		written_bytes += this->m_isa.serialize(out, child, "m_isa");
		written_bytes += write_member(this->m_l, out, child, "m_l");
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


	template<class t_csa, class t_bit_vector, class t_r_bit_vector, class t_s_bit_vector>
	void isa_lsw<t_csa, t_bit_vector, t_r_bit_vector, t_s_bit_vector>::load(std::istream& in)
	{
		this->m_isa.load(in);
		
		read_member(this->m_l, in);
		
		typename decltype(this->m_psi_k_support)::size_type psi_k_count(0);
		read_member(psi_k_count, in);
		
		this->m_psi_k_support.resize(psi_k_count);
		for (decltype(psi_k_count) i(0); i < psi_k_count; ++i)
			this->m_psi_k_support[i].load(in);
	}
}

#endif
