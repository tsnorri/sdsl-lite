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

#ifndef INCLUDED_SDSL_PSI_K_INDEX
#define INCLUDED_SDSL_PSI_K_INDEX

#include <sdsl/int_vector_buffer.hpp>


namespace sdsl
{
	//! An index for calculating Ψ_k.
	template<class t_sa_buf_type>
	class psi_k_index
	{
	public:
		typedef uint64_t value_type;
	protected:
		t_sa_buf_type &m_sa_buf; // Check that psi_k object lifetime does not exceed that of sa_buf.
		int_vector_buffer<0> m_psi_idx;
	
	public:
		psi_k_index() = delete;
		psi_k_index(psi_k_index const &) = delete;
		psi_k_index(psi_k_index &&) = delete;
		psi_k_index &operator=(psi_k_index const &) = delete;
		psi_k_index &operator=(psi_k_index &&) = delete;
		
		psi_k_index(cache_config& config, t_sa_buf_type &sa_buf):
			m_sa_buf(sa_buf)
		{
			std::string tmp_key(util::to_string(util::pid()) + "_Psi_k_idx_" + util::to_string(util::id()));
			std::string tmp_file_name(cache_file_name(tmp_key, config));
			// FIXME: use int_vector_mapper instead.
			int_vector_buffer<0> psi_idx_tmp(tmp_file_name, std::ios::out, 16 * 1024 * 1024, m_sa_buf.width());
		
			for (uint64_t i(0), count(m_sa_buf.size()); i < count; ++i)
				psi_idx_tmp[m_sa_buf[i]] = i;
		
			m_psi_idx = std::move(psi_idx_tmp);
		}
		
		int_vector_buffer<0> &psi_idx() { return m_psi_idx; }
	
		// Accessors return values in range [0, n], i.e. 1-based indices.
		value_type from_sa_val(uint64_t k, uint64_t v);
			
		value_type operator()(uint64_t k, uint64_t i);
	};
	
	
	// Accessors return values in range [0, n], i.e. 1-based indices.
	template<class t_sa_buf_type>
	auto psi_k_index<t_sa_buf_type>::from_sa_val(uint64_t k, uint64_t v) -> value_type
	{
		v += k;
		if (v < m_sa_buf.size())
			return 1 + m_psi_idx[v];
		else
			return 0;
	}
	
	
	template<class t_sa_buf_type>
	auto psi_k_index<t_sa_buf_type>::operator()(uint64_t k, uint64_t i) -> value_type
	{
		auto v(m_sa_buf[i]);
		return from_sa_val(k, v);
	}
}

#endif
