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

#ifndef INCLUDED_SDSL_ELIAS_INVENTORY
#define INCLUDED_SDSL_ELIAS_INVENTORY

#include <sdsl/int_vector.hpp>
#include <sdsl/select_support.hpp>


namespace sdsl {
	template<class t_s_bit_vector>
	class elias_inventory_base
	{
	public:
		typedef	t_s_bit_vector					s_bit_vector;
		typedef int_vector<0>::size_type		size_type;
		typedef uint64_t						value_type;
		
	protected:
		s_bit_vector m_values_high;
		int_vector<0> m_values_low;
		value_type m_max{0};
		uint8_t m_low_bits{0};
	
	public:
		elias_inventory_base() = default;
		elias_inventory_base(elias_inventory_base const &) = default;
		elias_inventory_base(elias_inventory_base &&) = default;
		elias_inventory_base &operator=(elias_inventory_base const &) & = default;
		elias_inventory_base &operator=(elias_inventory_base &&) & = default;
	};
	
	
	//! A class for storing ordered sequences of integers in compressed form.
	/*! The implementation follows Grossi and Vitter's and Rao's papers.
	 *  \tparam t_s_bit_vector		Type of bit vectors for which select support is needed.
	 *  
	 *  \par References
	 *  – Roberto Grossi, Jeffrey Scott Vitter:
	 *    Compressed Suffix Arrays and Suffix Trees With Applications to Text Indexing and String Matching.
	 *    SIAM Journal on Computing 35(2): 378–407 (2005)
	 *  – S. Srinivasa Rao:
	 *    Time-space trade-offs for compressed suffix arrays
	 *    Information Processing Letters 82(6): 307–311 (2002)
	 */
	// TODO: compare Grossi and Vitter's and Rao's variants to the one proposed by Elias in Efficient Storage and Retrieval by Content and Address of Static Files, Journal of the ACM, 21(2): 246–260 (1974).
	// TODO: verify time and space complexity.
	template<class t_s_bit_vector>
	class elias_inventory : public elias_inventory_base<t_s_bit_vector>
	{
	public:
		typedef elias_inventory_base<t_s_bit_vector>	base_class;
		typedef typename base_class::s_bit_vector		s_bit_vector;
		typedef typename base_class::size_type			size_type;
		typedef typename base_class::value_type			value_type;
		typedef	typename s_bit_vector::select_1_type	select_1_support_type;
		
	protected:
		select_1_support_type m_values_high_s1_support;
		
		
	public:
		elias_inventory():
			base_class::elias_inventory_base(),
			m_values_high_s1_support(&this->m_values_high)
		{
		}
		
		
		elias_inventory(elias_inventory const &other):
			base_class::elias_inventory_base(other),
			m_values_high_s1_support(&this->m_values_high)
		{
		}
		
		
		elias_inventory(elias_inventory &&other):
			base_class::elias_inventory_base(std::move(other)),
			m_values_high_s1_support(&this->m_values_high)
		{
		}
		
		
		// vec 0-based, ck 1-based.
		template<class t_vec, class t_ck>
		elias_inventory(t_vec const &vec, t_ck const &ck);
		
		elias_inventory &operator=(elias_inventory const &other) &;
		elias_inventory &operator=(elias_inventory &&other) &;
		value_type operator[](size_type i) const SDSL_HOT;
		
		auto serialize(std::ostream& out, structure_tree_node *v = nullptr, std::string name = "") const -> size_type;
		void load(std::istream& in);
	};
	
	
	template<class t_s_bit_vector>
	auto elias_inventory<t_s_bit_vector>::operator=(elias_inventory const &other) & -> elias_inventory &
	{
		base_class::operator=(other);
		m_values_high_s1_support = other.m_values_high_s1_support;
		m_values_high_s1_support.set_vector(&this->m_values_high);
		return *this;
	}
	
	
	template<class t_s_bit_vector>
	auto elias_inventory<t_s_bit_vector>::operator=(elias_inventory &&other) & -> elias_inventory &
	{
		base_class::operator=(std::move(other));
		m_values_high_s1_support = std::move(other.m_values_high_s1_support);
		m_values_high_s1_support.set_vector(&this->m_values_high);
		return *this;
	}
	
	
	template<class t_s_bit_vector>
	auto elias_inventory<t_s_bit_vector>::operator[](size_type i) const -> value_type
	{
		value_type low_val(this->m_values_low[i]);
		value_type high_val(m_values_high_s1_support.select(1 + i) - i);
		high_val <<= this->m_low_bits;
		high_val |= low_val;
		return high_val % this->m_max;
	}
	
	
	// vec 0-based, ck 1-based.
	template<class t_s_bit_vector>
	template<class t_vec, class t_ck>
	elias_inventory<t_s_bit_vector>::elias_inventory(t_vec const &vec, t_ck const &ck):
		elias_inventory()
	{
		decltype(this->m_max) max(0);
		size_type const total_count(ck[ck.size() - 1]);
		
		// Find a value greater than any of the stored.
		for (size_type j(0), set_count(vec.size()); j < set_count; ++j)
		{
			for (size_type i(0), count(vec[j].size()); i < count; ++i)
				max = std::max(vec[j][i], max);
		}
		++max;
		
		// Combine the values in vec into one set.
		int_vector<64> combined(total_count, 0);
		value_type max_sum(0);
		for (size_type j(0), set_count(vec.size()); j < set_count; ++j)
		{
			for (size_type i(0), count(vec[j].size()); i < count; ++i)
			{
				// ck[j] is the cumulative sum of the counts of the items in the previous lists.
				auto current_val(j * max + vec[j][i]);
				assert(current_val <= combined.max_value());
				combined[ck[j] + i] = current_val;
				max_sum = std::max(current_val, max_sum);
			}
		}
		
		// Rao's lemma 1.
		// Calculate the minimum number of bits.
		size_type const bits(std::ceil(std::log2(max_sum)));
		// Grossi and Vitter use floor (lemma 2), Rao uses ceil (lemma 1).
		size_type const high_bits(std::floor(std::log2(total_count)));
		decltype(this->m_low_bits) low_bits = bits - high_bits;
		bit_vector high_values(2 * total_count, 0); // Size from lemma 1.
		int_vector<0> low_values(total_count, 0, util::upper_power_of_2(low_bits));
		size_type ptr(0);
		value_type prev_val(0);
		
		for (size_type i(0), count(combined.size()); i < count; ++i)
		{
			auto const val(combined[i]);
			assert(prev_val <= val);
			auto const high_val((val >> low_bits) - (prev_val >> low_bits));
			auto const low_val((val << (64 - low_bits)) >> (64 - low_bits));
			
			ptr += high_val;
			high_values[ptr] = 1;
			assert(low_val <= low_values.max_value());
			low_values[i] = low_val;
			
			prev_val = val;
			++ptr;
		}
		
		this->m_values_high = high_values;
		m_values_high_s1_support = std::move(select_1_support_type(&this->m_values_high));
		this->m_values_low = std::move(low_values);
		this->m_max = max;
		this->m_low_bits = low_bits;
	}
	
	
	template<class t_s_bit_vector>
	auto elias_inventory<t_s_bit_vector>::serialize(std::ostream& out, structure_tree_node *v, std::string name) const -> size_type
	{
		structure_tree_node *child(structure_tree::add_child(v, name, util::class_name(*this)));
		size_type written_bytes(0);
		
		written_bytes += this->m_values_high.serialize(out, child, "m_values_high");
		written_bytes += this->m_values_low.serialize(out, child, "m_values_low");
		written_bytes += write_member(this->m_max, out, child, "m_max");
		written_bytes += write_member(this->m_low_bits, out, child, "m_low_bits");
		written_bytes += m_values_high_s1_support.serialize(out, child, "m_values_high_s1_support");
	
		structure_tree::add_size(child, written_bytes);
		return written_bytes;
	}
	
	
	template<class t_s_bit_vector>
	void elias_inventory<t_s_bit_vector>::load(std::istream& in)
	{
		this->m_values_high.load(in);
		this->m_values_low.load(in);
		read_member(this->m_max, in);
		read_member(this->m_low_bits, in);
		m_values_high_s1_support.load(in);
		m_values_high_s1_support.set_vector(&this->m_values_high);
	}
}

#endif
