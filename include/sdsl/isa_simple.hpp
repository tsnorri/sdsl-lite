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

#ifndef INCLUDED_SDSL_ISA_SIMPLE
#define INCLUDED_SDSL_ISA_SIMPLE

#include <sdsl/int_vector.hpp>


namespace sdsl
{
	//! A simple ISA.
	template<class t_sa_buf_type>
	class isa_simple
	{
	public:
		typedef uint64_t value_type;
	protected:
		t_sa_buf_type &m_sa_buf; // Check that isa_simple object lifetime does not exceed that of m_sa_buf.
		int_vector<0> m_isa;
	
	public:
		isa_simple() = delete;
		isa_simple(isa_simple const &) = delete;
		isa_simple(isa_simple &&) = delete;
		isa_simple &operator=(isa_simple const &) = delete;
		isa_simple &operator=(isa_simple &&) = delete;
		
		isa_simple(cache_config& config, t_sa_buf_type &sa_buf):
			m_sa_buf(sa_buf)
		{
			auto const count(m_sa_buf.size());
			int_vector<0> isa_tmp(count, 0, 1 + bits::hi(count));
		
			for (uint64_t i(0); i < count; ++i)
			{
				assert(i <= isa_tmp.max_value());
				isa_tmp[m_sa_buf[i]] = i;
			}
		
			m_isa = std::move(isa_tmp);
		}
		
		int_vector<0> &isa() { return m_isa; }
	
		// Accessors return values in range [0, n], i.e. 1-based indices.
		value_type psi_k_from_sa_val(uint64_t k, uint64_t v) const;
			
		value_type psi_k(uint64_t k, uint64_t i) const;
	};
	
	
	// Accessors return values in range [0, n], i.e. 1-based indices.
	template<class t_sa_buf_type>
	auto isa_simple<t_sa_buf_type>::psi_k_from_sa_val(uint64_t k, uint64_t v) const -> value_type
	{
		v += k;
		if (v < m_sa_buf.size())
			return 1 + m_isa[v];
		else
			return 0;
	}
	
	
	template<class t_sa_buf_type>
	auto isa_simple<t_sa_buf_type>::psi_k(uint64_t k, uint64_t i) const -> value_type
	{
		auto v(m_sa_buf[i]);
		return psi_k_from_sa_val(k, v);
	}
}

#endif
