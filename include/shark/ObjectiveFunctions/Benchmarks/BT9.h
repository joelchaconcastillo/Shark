//===========================================================================
/*!
 * 
 *
 * \brief       Objective function BT9
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
#ifndef SHARK_OBJECTIVEFUNCTIONS_BENCHMARK_BT9_H
#define SHARK_OBJECTIVEFUNCTIONS_BENCHMARK_BT9_H

#include <shark/ObjectiveFunctions/AbstractObjectiveFunction.h>
#include <shark/ObjectiveFunctions/BoxConstraintHandler.h>

namespace shark {namespace benchmarks{

/**
* \brief Implements the benchmark function BT9.
*
* See: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.18.7531&rep=rep1&type=pdf
* The benchmark function exposes the following features:
*	- Scalable w.r.t. the searchspace and w.r.t. the objective space.
*	- Highly multi-modal.
* \ingroup benchmarks
*/
struct BT9 : public MultiObjectiveFunction
{
	BT9(std::size_t numVariables = 0) : m_objectives(3), m_handler(numVariables,0,1 ){
		announceConstraintHandler(&m_handler);
	}

	/// \brief From INameable: return the class name.
	std::string name() const
	{ return "BT9"; }

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
	double D1(double g, double theta)const
        {
            return (g*g) + (1.0-std::exp(-(g*g)/theta))/5.0;
        }
	ResultType eval( const SearchPointType & x ) const {
		m_evaluationCounter++;
		ResultType value( numberOfObjectives() );
		int n = numberOfVariables();
                double sum1 = 0.0, sum2 = 0.0, sum3=0.0, theta = 1.0e-9;
                for(int j = 3; j <= n; j++)
                {
                   double yj = x[j-1] - sin((j*M_PI)/(2.0*n));
                   if((j%3)==0) sum1 += D1(yj, theta);
                   if((j%3)==1) sum2 += D1(yj, theta);
                   if((j%3)==2) sum3 += D1(yj, theta);
                } 
                value[0] = cos(0.5*x[0]*M_PI)*cos(0.5*x[1]*M_PI) + 10.0*sum1;
                value[1] = cos(0.5*x[0]*M_PI)*sin(0.5*x[1]*M_PI) + 10.0*sum2;
                value[2] = sin(0.5*x[0]*M_PI) + 10.0*sum3;

		return value;
	}
private:
	std::size_t m_objectives;
	BoxConstraintHandler<SearchPointType> m_handler;
	
};

}}

#endif
