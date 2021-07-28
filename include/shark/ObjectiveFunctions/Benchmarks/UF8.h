//===========================================================================
/*!
 * 
 *
 * \brief       Objective function UF8
 * 
 * 
 *
 * \author      T.Voss, T. Glasmachers, O.Krause
 * \date        2010-2011
 *
 *
 * \par Copyright 1995-2017 Shark Development Team
 * 
 * <BR><HR>
 * This file is part of Shark.
 * <http://shark-ml.org/>
 * 
 * Shark is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published 
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Shark is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with Shark.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
//===========================================================================
#ifndef SHARK_OBJECTIVEFUNCTIONS_BENCHMARK_UF8_H
#define SHARK_OBJECTIVEFUNCTIONS_BENCHMARK_UF8_H

#include <shark/ObjectiveFunctions/AbstractObjectiveFunction.h>
#include <shark/ObjectiveFunctions/BoxConstraintHandler.h>
#include <shark/ObjectiveFunctions/Benchmarks/cec09.h>

namespace shark {namespace benchmarks{

/**
* \brief Implements the benchmark function UF8.
*
* See: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.18.7531&rep=rep1&type=pdf
* The benchmark function exposes the following features:
*	- Scalable w.r.t. the searchspace and w.r.t. the objective space.
*	- Highly multi-modal.
* \ingroup benchmarks
*/
struct UF8 : public MultiObjectiveFunction
{
	UF8(std::size_t numVariables = 0) : m_objectives(3), m_handler(numVariables,0,1 ){
		announceConstraintHandler(&m_handler);
	}

	/// \brief From INameable: return the class name.
	std::string name() const
	{ return "UF8"; }

	std::size_t numberOfObjectives()const{
		return m_objectives;
	}
	bool hasScalableObjectives()const{
		return true;
	}
	
	void setNumberOfObjectives( std::size_t numberOfObjectives ){
		m_objectives = numberOfObjectives;
	}
	
	
	std::size_t numberOfVariables()const{
		return m_handler.dimensions();
	}
	
	bool hasScalableDimensionality()const{
		return true;
	}

	void setNumberOfVariables( std::size_t numberOfVariables ){
		SearchPointType lb(numberOfVariables,-2);
		SearchPointType ub(numberOfVariables, 2);
		lb(0)=0;
		ub(0)=1;
		lb(1)=0;
		ub(1)=1;

		m_handler.setBounds(lb, ub);
	}
	ResultType eval( const SearchPointType & x ) const {
		m_evaluationCounter++;

		ResultType value( numberOfObjectives() );
		CEC09::UF8(&(*(x.begin())), &(*(value.begin())), numberOfVariables());
		return value;
	}
private:
	std::size_t m_objectives;
	BoxConstraintHandler<SearchPointType> m_handler;
	
};

}}

#endif
