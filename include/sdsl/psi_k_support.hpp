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


namespace sdsl { namespace detail {
		
	template<typename t_r_bit_vector>
	class psi_k_v
	{
	public:
		typedef t_r_bit_vector				r_bit_vector;
		typedef int_vector<0>::size_type	size_type;
		
	protected:
		r_bit_vector m_v_values;
		typename r_bit_vector::rank_1_type m_v_values_r1_support;
		
	protected:
		void copy(psi_k_v const &other)
		{
			m_v_values = other.m_v_values;
			m_v_values_r1_support = other.m_v_values_r1_support;
			m_v_values_r1_support.set_vector(&m_v_values);
		}
		
		void move(psi_k_v &&other)
		{
			m_v_values = std::move(other.m_v_values);
			m_v_values_r1_support = std::move(other.m_v_values_r1_support);
			m_v_values_r1_support.set_vector(&m_v_values);
		}
		
	public:
		psi_k_v(): m_v_values_r1_support(&m_v_values) {}
		psi_k_v(psi_k_v const &other) { copy(other); }
		psi_k_v(psi_k_v &&other) { move(std::move(other)); }
		psi_k_v &operator=(psi_k_v const &other) & { copy(other); return *this; }
		psi_k_v &operator=(psi_k_v &&other) & { move(std::move(other)); return *this; }

		template<typename t_bit_vector>
		psi_k_v(t_bit_vector const &v_values):
			m_v_values(v_values),
			m_v_values_r1_support(&m_v_values)
		{
		}

		r_bit_vector const &v_values() const { return m_v_values; }
		typename r_bit_vector::rank_1_type const &r1_support() const { return m_v_values_r1_support; }
		auto serialize(std::ostream& out, structure_tree_node *v = nullptr, std::string name = "") const -> size_type;
		void load(std::istream& in);
	};
	
	
	template<typename t_r_bit_vector>
	auto psi_k_v<t_r_bit_vector>::serialize(std::ostream& out, structure_tree_node *v, std::string name) const -> size_type
	{
		structure_tree_node *child(structure_tree::add_child(v, name, util::class_name(*this)));
		size_type written_bytes(0);
		
		written_bytes += m_v_values.serialize(out, child, "v_values");
		written_bytes += m_v_values_r1_support.serialize(out, child, "v_values_r1_support");
		
		structure_tree::add_size(child, written_bytes);
		return written_bytes;
	}
	
	
	template<class t_r_bit_vector>
	void psi_k_v<t_r_bit_vector>::load(std::istream& in)
	{
		m_v_values.load(in);
		m_v_values_r1_support.load(in);
		m_v_values_r1_support.set_vector(&m_v_values);
	}
}}
	
	
namespace sdsl {
	
	template<class t_s_bit_vector>
	class psi_k_support
	{
	public:
		typedef	t_s_bit_vector				s_bit_vector;
		typedef int_vector<0>::size_type	size_type;
		
	protected:
		elias_inventory<s_bit_vector> m_psi_k_values;
		
	public:
		psi_k_support() {}
		psi_k_support(psi_k_support const &) = default;
		psi_k_support(psi_k_support &&) = default;
		psi_k_support &operator=(psi_k_support const &) & = default;
		psi_k_support &operator=(psi_k_support &&) & = default;
		
		psi_k_support(
			bit_vector const &v_values,					// Unused
			int_vector<0> &l_k_values,					// Unused
			int_vector<0> &c_k_values,					// Unused
			elias_inventory<s_bit_vector> &psi_k_values	// Moved
		):
			m_psi_k_values(std::move(psi_k_values))
		{
		}

		elias_inventory<s_bit_vector> const &psi_k_values() const { return m_psi_k_values; }
		uint64_t operator[](size_type i) const SDSL_HOT;
		auto serialize(std::ostream& out, structure_tree_node *v = nullptr, std::string name = "") const -> size_type;
		void load(std::istream& in);
	};
	
	
	//! A class for providing compressed Ψ^k support as proposed by S. Srinivasa Rao.
	/*! Rao's CSA is a generalization of the CSA proposed by R. Grossi and J. S. Vitter.
	 *  \tparam t_r_bit_vector		Type of bit vectors for which rank support is needed.
	 *  \tparam t_s_bit_vector		Type of bit vectors for which select support is needed.
	 *  \sa sdsl::psi_k_support_builder
	 *  
	 *  \par Reference
	 *  S. Srinivasa Rao:
	 *  Time-space trade-offs for compressed suffix arrays
	 *  Information Processing Letters 82(6): 307–311 (2002)
	 */
	// TODO: verify time and space complexity.
	template<class t_bit_vector, class t_r_bit_vector, class t_s_bit_vector>
	class psi_k_support_v : public psi_k_support<t_s_bit_vector>
	{
	public:
		typedef psi_k_support<t_s_bit_vector>		base_class;
		typedef	t_r_bit_vector						r_bit_vector;
		typedef typename base_class::s_bit_vector	s_bit_vector;
		typedef typename base_class::size_type		size_type;
		
	protected:
		detail::psi_k_v<r_bit_vector> m_v;
		
	public:
		psi_k_support_v(
			bit_vector const &v_values,					// Copied
			int_vector<0> &l_k_values,					// Unused
			int_vector<0> &c_k_values,					// Unused
			elias_inventory<s_bit_vector> &psi_k_values	// Moved
		):
			base_class::psi_k_support(v_values, l_k_values, c_k_values, psi_k_values),
			m_v(v_values)
		{
		}
		
		psi_k_support_v() = default;
		psi_k_support_v(psi_k_support_v const &) = default;
		psi_k_support_v(psi_k_support_v &&other) = default;
		psi_k_support_v &operator=(psi_k_support_v const &other) & = default;
		psi_k_support_v &operator=(psi_k_support_v &&other) & = default;
		
		uint64_t operator[](size_type i) const SDSL_HOT;
		auto serialize(std::ostream& out, structure_tree_node *v = nullptr, std::string name = "") const -> size_type;
		void load(std::istream& in);
	};
	
	
	template<class t_s_bit_vector>
	uint64_t psi_k_support<t_s_bit_vector>::operator[](size_type i) const
	{
		auto retval(m_psi_k_values[i] - 1);
		return retval;
	}
	
	
	template<class t_bit_vector, class t_r_bit_vector, class t_s_bit_vector>
	uint64_t psi_k_support_v<t_bit_vector, t_r_bit_vector, t_s_bit_vector>::operator[](size_type i) const
	{
		uint64_t const r(m_v.r1_support().rank(1 + i) - 1);
		auto retval(base_class::operator[](r));
		return retval;
	}
	
	
	template<class t_s_bit_vector>
	auto psi_k_support<t_s_bit_vector>::serialize(std::ostream& out, structure_tree_node *v, std::string name) const -> size_type
	{
		structure_tree_node *child(structure_tree::add_child(v, name, util::class_name(*this)));
		size_type written_bytes(0);
		
		written_bytes += m_psi_k_values.serialize(out, child, "psi_k_values");
		
		structure_tree::add_size(child, written_bytes);
		return written_bytes;
	}

	
	template<class t_bit_vector, class t_r_bit_vector, class t_s_bit_vector>
	auto psi_k_support_v<t_bit_vector, t_r_bit_vector, t_s_bit_vector>::serialize(std::ostream& out, structure_tree_node *v, std::string name) const -> size_type
	{
		structure_tree_node *child(structure_tree::add_child(v, name, util::class_name(*this)));
		size_type written_bytes(0);
		
		written_bytes += base_class::serialize(out, child, "psi_k_support");
		written_bytes += m_v.serialize(out, child, "v");
	
		structure_tree::add_size(child, written_bytes);
		return written_bytes;
	}
	

	template<class t_s_bit_vector>
	void psi_k_support<t_s_bit_vector>::load(std::istream& in)
	{
		m_psi_k_values.load(in);
	}

	
	template<class t_bit_vector, class t_r_bit_vector, class t_s_bit_vector>
	void psi_k_support_v<t_bit_vector, t_r_bit_vector, t_s_bit_vector>::load(std::istream& in)
	{
		base_class::load(in);
		m_v.load(in);
	}
}

#endif
