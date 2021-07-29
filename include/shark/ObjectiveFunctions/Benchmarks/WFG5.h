//===========================================================================
/*!
 * 
 *
 * \brief       Objective function WFG5
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
#ifndef SHARK_OBJECTIVEFUNCTIONS_BENCHMARK_WFG5_H
#define SHARK_OBJECTIVEFUNCTIONS_BENCHMARK_WFG5_H

#include <shark/ObjectiveFunctions/AbstractObjectiveFunction.h>
#include <shark/ObjectiveFunctions/BoxConstraintHandler.h>
#include <shark/ObjectiveFunctions/Benchmarks/cec09.h>
#include <shark/ObjectiveFunctions/Benchmarks/Toolkit/ExampleProblems.h>
#include <shark/ObjectiveFunctions/Benchmarks/Toolkit/TransFunctions.h>
using namespace WFG::Toolkit;
using namespace WFG::Toolkit::Examples;
namespace shark {namespace benchmarks{

/**
* \brief Implements the benchmark function WFG5.
*
* See: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.18.7531&rep=rep1&type=pdf
* The benchmark function exposes the following features:
*	- Scalable w.r.t. the searchspace and w.r.t. the objective space.
*	- Highly multi-modal.
* \ingroup benchmarks
*/
struct WFG5 : public MultiObjectiveFunction
{
	WFG5(std::size_t numVariables = 0) : m_objectives(2), m_handler(numVariables,0,1){
		announceConstraintHandler(&m_handler);
	}

	/// \brief From INameable: return the class name.
	std::string name() const
	{ return "WFG5"; }

	std::size_t numberOfObjectives()const{
		return m_objectives;
	}
	bool hasScalableObjectives()const{
		return true;
	}
	
	void setNumberOfObjectives( std::size_t numberOfObjectives ){
		m_objectives = numberOfObjectives;
	}
	void setPositionParameters(std::size_t position_k)
	{
	       m_position_k=position_k;
	}
	std::size_t numberOfVariables()const{
		return m_handler.dimensions();
	}
	
	bool hasScalableDimensionality()const{
		return true;
	}

	void setNumberOfVariables( std::size_t numberOfVariables ){
		SearchPointType lb(numberOfVariables,0);
		SearchPointType ub(numberOfVariables, 1);
		for(std::size_t i=0; i < numberOfVariables; i++) ub(i) = 2.0*(i+1);
		m_handler.setBounds(lb, ub);
	}
	ResultType eval( const SearchPointType & x ) const {
		m_evaluationCounter++;
		std::vector<double> f= Problems::WFG5(std::vector<double> ( &(*(x.begin())), &(*(x.begin()))+numberOfVariables()), m_position_k, numberOfObjectives());
		ResultType value(&(*(f.begin())), &(*(f.begin())) + numberOfObjectives());
		return value;
	}
private:
	std::size_t m_objectives, m_position_k;
	BoxConstraintHandler<SearchPointType> m_handler;
	
};

}}

#endif
