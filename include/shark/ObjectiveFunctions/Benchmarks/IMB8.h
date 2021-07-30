//===========================================================================
/*!
 * 
 *
 * \brief       Objective function IMB8
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
#ifndef SHARK_OBJECTIVEFUNCTIONS_BENCHMARK_IMB8_H
#define SHARK_OBJECTIVEFUNCTIONS_BENCHMARK_IMB8_H

#include <shark/ObjectiveFunctions/AbstractObjectiveFunction.h>
#include <shark/ObjectiveFunctions/BoxConstraintHandler.h>

namespace shark {namespace benchmarks{

/**
* \brief Implements the benchmark function IMB8.
*
* See: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.18.7531&rep=rep1&type=pdf
* The benchmark function exposes the following features:
*	- Scalable w.r.t. the searchspace and w.r.t. the objective space.
*	- Highly multi-modal.
* \ingroup benchmarks
*/
struct IMB8 : public MultiObjectiveFunction
{
	IMB8(std::size_t numVariables = 0) : m_objectives(2), m_handler(numVariables,0,1 ){
		announceConstraintHandler(&m_handler);
	}

	/// \brief From INameable: return the class name.
	std::string name() const
	{ return "IMB8"; }

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
		SearchPointType lb(numberOfVariables,0);
		SearchPointType ub(numberOfVariables, 1);
		m_handler.setBounds(lb, ub);
	}
	ResultType eval( const SearchPointType & x ) const {
		m_evaluationCounter++;
		ResultType value( numberOfObjectives() );
		int n = numberOfVariables();
                double gx = 0.0;
                for(int i = 1; i < n; i++)
                {
              	  double ti = x[i] - 0.5, si=x[i]-sin(0.5*M_PI*x[0]);
              	  if(x[0]>=0.5 && x[0]<=0.8)
              	    gx += (-0.9*si*si + pow(fabs(si), 0.6)); 
              	  else 
              	    gx += pow(fabs(ti), 0.6);
                }
                value[0] = (1.0+gx)*x[0];
                value[1] = (1.0+gx)*(1.0-x[0]);
                return value;
	}
private:
	std::size_t m_objectives;
	BoxConstraintHandler<SearchPointType> m_handler;
	
};

}}

#endif
