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

#ifndef INCLUDED_SDSL_PSI_K_SUPPORT
#define INCLUDED_SDSL_PSI_K_SUPPORT

#include <sdsl/elias_inventory.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/rrr_vector.hpp>


namespace sdsl
{
	template<class t_r_bit_vector, class t_s_bit_vector>
	class psi_k_support_base
	{
	public:
		typedef	t_r_bit_vector	r_bit_vector;
		typedef	t_s_bit_vector	s_bit_vector;
		
	protected:
		r_bit_vector m_v_values;
		// FIXME: check that these really aren't needed.
		//int_vector<0> m_l_k_values;
		//int_vector<64> m_c_k_values;
		elias_inventory<s_bit_vector> m_psi_k_values;
		
	public:
		psi_k_support_base() {}
		psi_k_support_base(psi_k_support_base const &) = default;
		psi_k_support_base(psi_k_support_base &&) = default;
		psi_k_support_base &operator=(psi_k_support_base const &) & = default;
		psi_k_support_base &operator=(psi_k_support_base &&) & = default;
		
		psi_k_support_base(
			bit_vector const &v_values,					// Copied
			int_vector<0> &l_k_values,					// Moved
			int_vector<64> &c_k_values,					// Moved
			elias_inventory<s_bit_vector> &psi_k_values	// Moved
		):
			m_v_values(v_values),
			//m_l_k_values(std::move(l_k_values)),
			//m_c_k_values(std::move(c_k_values)),
			m_psi_k_values(std::move(psi_k_values))
		{
		}
	};
	
	
	//! A class for providing compressed Ψ_k support as proposed by S. Srinivasa Rao.
	/*! Rao's CSA is a generalization of the CSA proposed by R. Grossi and J. S. Vitter.
	 *  \tparam t_r_bit_Vector		Type of bit vectors for which rank support is needed.
	 *  \tparam t_s_bit_vector		Type of bit vectors for which select support is needed.
	 *  \sa sdsl::psi_k_support_builder
	 *  
	 *  \par Reference
	 *  S. Srinivasa Rao:
	 *  Time-space trade-offs for compressed suffix arrays
	 *  Information Processing Letters 82(6): 307–311 (2002)
	 */
	// TODO: verify time and space complexity.
	template<
		class t_r_bit_vector = rrr_vector<>,
		class t_s_bit_vector = bit_vector
	>
	class psi_k_support : public psi_k_support_base<t_r_bit_vector, t_s_bit_vector>
	{
	public:
		typedef psi_k_support_base<t_r_bit_vector, t_s_bit_vector>	base_class;
		typedef typename base_class::r_bit_vector					r_bit_vector;
		typedef typename base_class::s_bit_vector					s_bit_vector;
		typedef int_vector<0>::size_type							size_type;
		
	protected:
		typename r_bit_vector::rank_1_type m_v_values_r1_support;
		
	public:
		psi_k_support(
			bit_vector const &v_values,					// Copied
			int_vector<0> &l_k_values,					// Moved
			int_vector<64> &c_k_values,					// Moved
			elias_inventory<s_bit_vector> &psi_k_values	// Moved
		):
			base_class::psi_k_support_base(v_values, l_k_values, c_k_values, psi_k_values),
			m_v_values_r1_support(&this->m_v_values)
		{
		}
		
		
		psi_k_support():
			base_class::psi_k_support_base(),
			m_v_values_r1_support(&this->m_v_values)
		{
		}
		
		
		psi_k_support(psi_k_support const &other):
			base_class::psi_k_support_base(other),
			m_v_values_r1_support(&this->m_v_values)
		{
		}
		
		
		psi_k_support(psi_k_support &&other):
			base_class::psi_k_support_base(std::move(other)),
			m_v_values_r1_support(&this->m_v_values)
		{
		}
		
		
		psi_k_support &operator=(psi_k_support const &other) &;
		psi_k_support &operator=(psi_k_support &&other) &;
		uint64_t operator[](size_type i) const;
		auto serialize(std::ostream& out, structure_tree_node *v = nullptr, std::string name = "") const -> size_type;
		void load(std::istream& in);
	};
	
	
	template<class t_r_bit_vector, class t_s_bit_vector>
	auto psi_k_support<t_r_bit_vector, t_s_bit_vector>::operator=(psi_k_support const &other) & -> psi_k_support &
	{
		base_class::operator=(other);
		m_v_values_r1_support = other.m_v_values_r1_support;
		m_v_values_r1_support.set_vector(&this->m_v_values);
		return *this;
	}
	
	
	template<class t_r_bit_vector, class t_s_bit_vector>
	auto psi_k_support<t_r_bit_vector, t_s_bit_vector>::operator=(psi_k_support &&other) & -> psi_k_support &
	{
		base_class::operator=(std::move(other));
		m_v_values_r1_support = std::move(other.m_v_values_r1_support);
		m_v_values_r1_support.set_vector(&this->m_v_values);
		return *this;
	}
	
	
	template<class t_r_bit_vector, class t_s_bit_vector>
	uint64_t psi_k_support<t_r_bit_vector, t_s_bit_vector>::operator[](size_type i) const
	{
		uint64_t const r(m_v_values_r1_support.rank(1 + i) - 1);
		auto retval(this->m_psi_k_values[r] - 1);
		return retval;
	}
	
	
	template<class t_r_bit_vector, class t_s_bit_vector>
	auto psi_k_support<t_r_bit_vector, t_s_bit_vector>::serialize(std::ostream& out, structure_tree_node *v, std::string name) const -> size_type
	{
		structure_tree_node *child(structure_tree::add_child(v, name, util::class_name(*this)));
		size_type written_bytes(0);
		
		written_bytes += this->m_v_values.serialize(out, child, "m_v_values");
		written_bytes += this->m_psi_k_values.serialize(out, child, "m_psi_k_values");
		written_bytes += m_v_values_r1_support.serialize(out, child, "m_v_values_r1_support");
	
		structure_tree::add_size(child, written_bytes);
		return written_bytes;
	}
	
	
	template<class t_r_bit_vector, class t_s_bit_vector>
	void psi_k_support<t_r_bit_vector, t_s_bit_vector>::load(std::istream& in)
	{
		this->m_v_values.load(in);
		this->m_psi_k_values.load(in);
		m_v_values_r1_support.load(in);
		m_v_values_r1_support.set_vector(&this->m_v_values);
	}
}

#endif
