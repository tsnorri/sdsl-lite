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

#ifndef INCLUDED_SDSL_CSA_RAO
#define INCLUDED_SDSL_CSA_RAO

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <cmath>

#include <sdsl/csa_alphabet_strategy.hpp>
#include <sdsl/csa_rao_delegate.hpp>
#include <sdsl/elias_inventory.hpp>
#include <sdsl/iterators.hpp>
#include <sdsl/psi_k_support.hpp>
#include <sdsl/rrr_vector.hpp>
#include <sdsl/suffix_array_helper.hpp>


namespace sdsl
{
	template<class t_csa_rao> class csa_rao_builder;
	template<class t_csa, class t_r_bit_vector, class t_s_bit_vector> class isa_lsw;
	
	//! Template argument container for sdsl::csa_rao.
	/*! \tparam	t_levels			Number of levels. Zero indicates a value determined at run time.
	 *  \tparam	t_partitions		Number of partitions. The size of the text has to be a multiple of
									t_partitions^t_levels. Zero indicates a value determined at run time.
	 *  \tparam	t_alphabet_strat	Policy for alphabet representation.
	 *  \tparam t_r_bit_Vector		Type of bit vectors for which rank support is needed.
	 *  \tparam t_s_bit_vector		Type of bit vectors for which select support is needed.
	 *  \tparam t_rs_bit_vector		Type of bit vectors for which rank and select support are needed.
	 *  \tparam t_isa				Class that implements the inverse suffix array.
	 *  \tparam t_delegate			Helper class for debugging purposes.
	 */
	// TODO: verify time and space complexity.
	template<uint32_t t_levels							= 0,
			 uint32_t t_partitions						= 0,
			 class t_alphabet_strat						= byte_alphabet,
			 class t_r_bit_vector						= rrr_vector<>,
			 class t_s_bit_vector						= bit_vector,
			 class t_rs_bit_vector						= bit_vector,
			 template<class, class, class> class t_isa	= isa_lsw,
			 class t_delegate							= csa_rao_delegate
			>
	class csa_rao_spec
	{
	public:
		typedef t_alphabet_strat	alphabet_strat;
		typedef t_r_bit_vector		r_bit_vector;
		typedef t_s_bit_vector		s_bit_vector;
		typedef t_rs_bit_vector		rs_bit_vector;
		typedef t_delegate			delegate_type;
		
		template<class t_csa> using isa_type = t_isa<t_csa, r_bit_vector, s_bit_vector>;
		
		static uint32_t const s_levels{t_levels};
		static uint32_t const s_partitions{t_partitions};
	};
	
	
	//! A class for the Compressed Suffix Array (CSA) proposed by S. Srinivasa Rao.
	/*! Rao's CSA is a generalization of the CSA proposed by R. Grossi and J. S. Vitter.
	 *  \tparam	t_spec	Specification class.
	 *  \sa sdsl::csa_rao_spec, sdsl::isa_lsw, sdsl::csa_rao_builder
	 *  
	 *  \par Reference
	 *  S. Srinivasa Rao:
	 *  Time-space trade-offs for compressed suffix arrays
	 *  Information Processing Letters 82(6): 307–311 (2002)
	 *  
	 *  @ingroup csa
	 */
	template<class t_spec = csa_rao_spec<>>
	class csa_rao
	{
	protected:
		template<typename T> using array = std::vector<T>; // dynarray is not move-assignable, which we need.
		
	public:
		typedef csa_rao																		csa_type;
		typedef uint64_t																	value_type;
		typedef random_access_const_iterator<csa_type>										const_iterator;
		typedef const_iterator																iterator;
		typedef const value_type															const_reference;
		typedef const_reference																reference;
		typedef const_reference																*pointer;
		typedef const pointer																const_pointer;
		typedef int_vector<>::size_type														size_type;
		typedef size_type																	csa_size_type;
		typedef ptrdiff_t																	difference_type;
		typedef traverse_csa_saisa<csa_type, true>											psi_type;
		typedef traverse_csa_saisa<csa_type, false>											lf_type;
		typedef rank_bwt_of_csa_psi<csa_type>												rank_bwt_type;
		typedef select_bwt_of_csa_psi<csa_type>												select_bwt_type;
		typedef bwt_of_csa_psi<csa_type>													bwt_type;
		typedef text_of_csa<csa_type>														text_type;
		typedef first_row_of_csa<csa_type>													first_row_type;
		typedef typename t_spec::template isa_type<csa_type>								isa_type;
		typedef typename t_spec::alphabet_strat												alphabet_type;
		typedef typename alphabet_type::alphabet_category									alphabet_category;
		typedef typename alphabet_type::char_type											char_type; // Note: This is the char type of the CSA not the WT!
		typedef typename alphabet_type::comp_char_type										comp_char_type;
		typedef typename alphabet_type::string_type											string_type;
		
		typedef csa_rao_builder<csa_type>													builder_type;
		typedef std::remove_const_t<decltype(t_spec::s_levels)>								level_count_type;
		typedef std::remove_const_t<decltype(t_spec::s_partitions)>							partition_count_type;
		typedef psi_k_support<typename t_spec::r_bit_vector, typename t_spec::s_bit_vector>	psi_k_support_type;
		typedef t_spec																		spec_type;

		typedef csa_tag																		index_category;
		typedef psi_tag																		extract_category;

		// FIXME: Does this make sense? (Affects child() in CST, which calls get_char_pos, which opts
		// to use Psi instead of SA and ISA, if the densities are high enough. However, this class doesn't
		// provide a Psi function, so the values are calculated with SA and ISA anyway. Another option would
		// be to check for a suitable Ψ_k and use it.)
		enum {
			sa_sample_dens = 1,
			isa_sample_dens = 1
		};
	
		friend class csa_rao_builder<csa_type>;
		friend typename spec_type::delegate_type;
		friend isa_type;
		friend class traverse_csa_saisa<csa_type, true>;
		friend class traverse_csa_saisa<csa_type, false>;
		
	protected:
		class level_base;
		class level;
		
	protected:
		int_vector<0> m_sa; // The most recent (in terms of level) suffix array, 0-based indices. // FIXME: really <0>? Check the initializer that it actually chooses the smallest type possible.
		array<level> m_levels;
		alphabet_type m_alphabet;
		level_count_type m_level_count;			// FIXME: is this needed?
		partition_count_type m_partition_count;	// FIXME: is this needed?
		isa_type m_isa;
	
	protected:
		rank_bwt_type							const rank_bwt		= rank_bwt_type(*this);
		select_bwt_type							const select_bwt	= select_bwt_type(*this);
		
	public:
		typename alphabet_type::char2comp_type	const &char2comp	= m_alphabet.char2comp;
		typename alphabet_type::comp2char_type	const &comp2char	= m_alphabet.comp2char;
		typename alphabet_type::C_type			const &C			= m_alphabet.C;
		typename alphabet_type::sigma_type		const &sigma		= m_alphabet.sigma;
		psi_type								const psi			= psi_type(*this);
		lf_type									const lf			= lf_type(*this);
		bwt_type								const bwt			= bwt_type(*this);
		bwt_type								const L				= bwt_type(*this);
		isa_type								const &isa			= m_isa;
		first_row_type							const F				= first_row_type(*this);
		text_type								const text			= text_type(*this);

	protected:
		void copy(csa_type const &other);
		void move(csa_type &&other);
		
	public:
		csa_rao():
			m_levels(),
			m_level_count(spec_type::s_levels),
			m_partition_count(spec_type::s_partitions),
			m_isa(*this)
		{
		}
	
		csa_rao(cache_config& config);
		
		csa_rao(csa_type const &csa):
			csa_rao()
		{
			copy(csa);
		}
		
		csa_rao(csa_type &&csa):
			csa_rao()
		{
			move(std::move(csa));
		}

		csa_rao &operator=(csa_type const &csa) &;
		csa_rao &operator=(csa_type &&csa) &;
	
		size_type size() const { return m_levels[0].d_values().size(); }
		size_type max_size() const { return m_levels[0].d_values().max_size(); }
		bool empty()const { return 0 == size(); }
		const_iterator begin() const { return const_iterator(this, 0); }
		const_iterator end() const { return const_iterator(this, size()); }
		const_iterator cbegin() const { return const_iterator(this, 0); }
		const_iterator cend() const { return const_iterator(this, size()); }
		value_type operator[](size_type i) const;
		value_type sa(size_type i, level_count_type level) const;
		uint64_t psi_k(uint64_t k, uint64_t i) const;
		uint64_t psi_k(level_count_type lidx, uint64_t k, uint64_t i) const;
	
		size_type serialize(std::ostream &out, structure_tree_node *v = nullptr, std::string name = "") const;
	
		void load(std::istream &in);
		
		void swap(csa_type &other) { using std::swap; swap(other, *this); }
		void swap(csa_type &&other) { swap(other); }
		
		partition_count_type partition_count() const { return m_partition_count; }
		level_count_type level_count() const { return m_level_count; }
		uint64_t decompress_sa(uint8_t ll, uint64_t val) const;
	};
	
	
	template<class t_spec>
	class csa_rao<t_spec>::level_base
	{
	public:
		typedef typename spec_type::r_bit_vector r_bit_vector;
		
	protected:
		array<psi_k_support_type> m_partitions;
		r_bit_vector m_b_values;
		int_vector<0> m_d_values;	// l − (SA[i] mod l), SA[i] 1-based (3.3). // FIXME: really <0>?
		
	public:
		level_base():
			m_partitions(0)
		{
		}
		
		template<class t_p_bit_vector>
		level_base(
			array<psi_k_support_type> &partitions,	// Moved
			t_p_bit_vector const &b_values,			// Copied
			int_vector<0> &d_values					// Moved
		):
			m_partitions(std::move(partitions)),
			m_b_values(b_values),
			m_d_values(std::move(d_values))
		{
		}
		
		level_base(level_base const &) = default;
		level_base(level_base &&other) = default;
		level_base &operator=(level_base const &other) & = default;
		level_base &operator=(level_base &&other) & = default;
	};
	
	
	template<class t_spec>
	class csa_rao<t_spec>::level : public level_base
	{
	public:
		typedef typename level_base::r_bit_vector r_bit_vector;
		
	protected:
		typename level_base::r_bit_vector::rank_1_type m_b_r1_support;
	
	public:
		level(): level_base() {}
		
		template<class t_p_bit_vector>
		level(
			array<psi_k_support_type> &partitions,	// Moved
			t_p_bit_vector const &b_values,			// Copied
			int_vector<0> &d_values					// Moved
		):
			level_base(partitions, b_values, d_values),
			m_b_r1_support(&this->m_b_values)
		{
		}
		
		level(level const &other):
			level_base(other),
			m_b_r1_support(&this->m_b_values)
		{
		}
		
		level(level &&other):
			level_base(std::move(other)),
			m_b_r1_support(&this->m_b_values)
		{
		}
		
		level &operator=(level const &other) &;
		level &operator=(level &&other) &;
		psi_k_support_type const &partition(typename array<psi_k_support_type>::size_type i) const { return this->m_partitions[i]; }
		int_vector<0> const &d_values() const { return this->m_d_values; }
		auto b_rank_1(typename r_bit_vector::size_type i) const { return m_b_r1_support.rank(i); } // rank in [0, i-1].
		
		auto serialize(std::ostream& out, structure_tree_node *v = nullptr, std::string name = "") const -> size_type;
		void load(std::istream& in, partition_count_type partition_count);
	};
	
	
	template<class t_spec>
	auto csa_rao<t_spec>::level::operator=(level const &other) & -> level &
	{
		level_base::operator=(other);
		m_b_r1_support = other.m_b_r1_support;
		m_b_r1_support.set_vector(&this->m_b_values);
		return *this;
	}
	
	
	template<class t_spec>
	auto csa_rao<t_spec>::level::operator=(level &&other) & -> level &
	{
		level_base::operator=(std::move(other));
		m_b_r1_support = std::move(other.m_b_r1_support);
		m_b_r1_support.set_vector(&this->m_b_values);
		return *this;
	}
	
	
	template<class t_spec>
	auto csa_rao<t_spec>::level::serialize(std::ostream &out, structure_tree_node *v, std::string name) const -> size_type
	{
		structure_tree_node *child(structure_tree::add_child(v, name, util::class_name(*this)));
		size_type written_bytes(0);
		
		partition_count_type i(1);
		for (auto it(this->m_partitions.cbegin()), end(this->m_partitions.cend()); it != end; ++it)
		{
			written_bytes += it->serialize(out, child, "partition_" + std::to_string(i));
			++i;
		}
		
		written_bytes += this->m_b_values.serialize(out, child, "m_b_values");
		written_bytes += this->m_d_values.serialize(out, child, "m_d_values");
		written_bytes += m_b_r1_support(out, child, "m_b_r1_support");
		
		structure_tree::add_size(child, written_bytes);
		return written_bytes;
	}
	
	
	template<class t_spec>
	void csa_rao<t_spec>::level::load(std::istream &in, partition_count_type partition_count)
	{
		this->m_partitions.resize(partition_count);
		for (decltype(partition_count) i(0); i < partition_count; ++i)
			this->m_partitions[i].load(in);
		
		this->m_b_values.load(in);
		this->m_d_values.load(in);
		m_b_r1_support.load(in);
		m_b_r1_support.set_vector(&this->m_b_values);
	}
	
	
	template<class t_spec>
	csa_rao<t_spec>::csa_rao(cache_config& config):
		csa_rao()
	{
		builder_type builder(*this, config);
		builder.build();
		
		isa_type isa_tmp(*this, config);
		m_isa = std::move(isa_tmp);
	}
	
	
	template<class t_spec>
	void csa_rao<t_spec>::copy(csa_type const &other)
	{
		m_sa = other.m_sa;
		m_levels = other.m_levels;
		m_alphabet = other.m_alphabet;
		m_level_count = other.m_level_count;
		m_partition_count = other.m_partition_count;
		m_isa = other.m_isa;
	}
	
	
	template<class t_spec>
	void csa_rao<t_spec>::move(csa_type &&other)
	{
		m_sa = std::move(other.m_sa);
		m_levels = std::move(other.m_levels);
		m_alphabet = std::move(other.m_alphabet);
		m_level_count = std::move(other.m_level_count);
		m_partition_count = std::move(other.m_partition_count);
		m_isa = std::move(other.m_isa);
	}
	
	
	template<class t_spec>
	auto csa_rao<t_spec>::operator=(csa_type const &csa) & -> csa_type &
	{
		copy(csa);
		return *this;
	}
	
	
	template<class t_spec>
	auto csa_rao<t_spec>::operator=(csa_type &&csa) & -> csa_type &
	{
		move(std::move(csa));
		return *this;
	}
	
	
	template<class t_spec>
	auto csa_rao<t_spec>::operator[](size_type i) const -> value_type
	{
		return sa(i, 0);
	}
		
	
	template<class t_spec>
	auto csa_rao<t_spec>::sa(size_type i, level_count_type lidx) const -> value_type
	{
		assert(lidx <= m_level_count);
		if (lidx == m_level_count)
			return m_sa[i];
		
		level const &level(m_levels[lidx]);
		
		// Check if i is stored in SA_{i+1}. If not, use Ψ_k.
		auto const k(level.d_values()[i]);
		if (k == m_partition_count)
		{
			auto const r(level.b_rank_1(1 + i));
			assert(r);
			auto const val(sa(r - 1, 1 + lidx));
			return m_partition_count * (1 + val) - 1;
		}
		
		auto const psi_k_i(psi_k(lidx, k, i));
		
		// The value for psi_k_i should be stored in SA_{i+1}.
		assert(psi_k_i);
		assert(level.d_values()[psi_k_i - 1] == m_partition_count);
		
		auto const r(level.b_rank_1(psi_k_i)); // rank in [0, i-1], psi_k_i 1-based.
		assert(r);
		auto const retval(m_partition_count * (1 + sa(r - 1, 1 + lidx)) - k - 1);
		return retval;
	}

	
	template<class t_spec>
	uint64_t csa_rao<t_spec>::psi_k(uint64_t k, uint64_t i) const
	{
		return this->psi_k(0, k, i);
	}

		
	template<class t_spec>
	uint64_t csa_rao<t_spec>::psi_k(level_count_type lidx, uint64_t k, uint64_t i) const
	{
		assert(lidx < m_level_count);
		assert(k < m_partition_count);
		
		if (0 == k)
			return i;
		
		level const &level(m_levels[lidx]);
		psi_k_support_type const &partition(level.partition(k - 1));
		
		auto retval(1 + partition[i]);
		return retval;
	}
	
	
	template<class t_spec>
	auto csa_rao<t_spec>::serialize(std::ostream &out, structure_tree_node *v, std::string name) const -> size_type
	{
		structure_tree_node *child(structure_tree::add_child(v, name, util::class_name(*this)));
		size_type written_bytes(0);

		written_bytes += m_sa.serialize(out, child, "m_sa");
		written_bytes += m_alphabet.serialize(out, child, "m_alphabet");
		written_bytes += write_member(m_level_count, child, "m_level_count");
		written_bytes += write_member(m_partition_count, child, "m_partition_count");
		
		level_count_type i(1);
		for (auto it(m_levels.cbegin()), end(m_levels.cend()); it != end; ++it)
		{
			written_bytes += it->serialize(out, child, "level_" + std::to_string(i));
			++i;
		}
		
		written_bytes += m_isa.serialize(out, child, "isa");

		structure_tree::add_size(child, written_bytes);
		return written_bytes;
	}
	
	
	template<class t_spec>
	void csa_rao<t_spec>::load(std::istream &in)
	{
		m_sa.load(in);
		m_alphabet.load(in);
		read_member(m_level_count, in);
		read_member(m_partition_count, in);
		
		m_levels.resize(m_level_count);
		for (decltype(m_level_count) i(0); i < m_level_count; ++i)
			m_levels[i].load(in);
		
		m_isa.load(in);
	}
	
	
	template<class t_spec>
	uint64_t csa_rao<t_spec>::decompress_sa(uint8_t ll, uint64_t val) const
	{
		uint64_t const a(util::ipow(m_partition_count, ll));
		uint64_t const retval(a + a * val - 1);
		return retval;
	}
} // end namespace sdsl
#endif
