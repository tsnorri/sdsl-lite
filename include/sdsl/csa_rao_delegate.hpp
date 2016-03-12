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

#ifndef INCLUDED_SDSL_CSA_RAO_DELEGATE
#define INCLUDED_SDSL_CSA_RAO_DELEGATE

#include <cstdint>

namespace sdsl
{
	template<class t_csa_rao>
	class csa_rao_builder;
	
	
	//! Debugging helper for sdsl::csa_rao.
	class csa_rao_delegate
	{
	public:
		template<class t_csa_rao>
		void setup(t_csa_rao &csa, csa_rao_builder<t_csa_rao> &builder) {}
		
		template<class t_csa_rao, class t_sa_buf_type>
		void start_level(csa_rao_builder<t_csa_rao> &builder, uint8_t ll, t_sa_buf_type &sa_buf) {}
		
		template<class t_csa_rao, class t_psi_k_fn>
		void finish_create_psi_k(csa_rao_builder<t_csa_rao> &builder, uint8_t ll, uint64_t k, uint64_t count, t_psi_k_fn &psi_k_fn) {}
	};
}

#endif
