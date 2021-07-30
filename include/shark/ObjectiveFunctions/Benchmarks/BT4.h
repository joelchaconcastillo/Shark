//===========================================================================
/*!
 * 
 *
 * \brief       Objective function BT4
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
#ifndef SHARK_OBJECTIVEFUNCTIONS_BENCHMARK_BT4_H
#define SHARK_OBJECTIVEFUNCTIONS_BENCHMARK_BT4_H

#include <shark/ObjectiveFunctions/AbstractObjectiveFunction.h>
#include <shark/ObjectiveFunctions/BoxConstraintHandler.h>

namespace shark {namespace benchmarks{

/**
* \brief Implements the benchmark function BT4.
*
* See: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.18.7531&rep=rep1&type=pdf
* The benchmark function exposes the following features:
*	- Scalable w.r.t. the searchspace and w.r.t. the objective space.
*	- Highly multi-modal.
* \ingroup benchmarks
*/
struct BT4 : public MultiObjectiveFunction
{
	BT4(std::size_t numVariables = 0) : m_objectives(2), m_handler(numVariables,0,1 ){
		announceConstraintHandler(&m_handler);
	}

	/// \brief From INameable: return the class name.
	std::string name() const
	{ return "BT4"; }

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
        double S2(double x, double gamma)const
        {
           if( x>=0.0 && x<0.25) return (1.0 - pow(1.0-4.0*x, gamma))/4.0;
           else if( x>=0.25 && x<0.5) return (1.0 + pow(4.0*x-1.0, gamma))/4.0;
           else if( x>=0.5 && x<0.75) return (3.0 - pow(3.0-4.0*x, gamma))/4.0;
           else if( x>=0.75 && x<=1.0) return (3.0 + pow(4.0*x-3.0, gamma))/4.0;
           return 0.0;
        }
	ResultType eval( const SearchPointType & x ) const {
		m_evaluationCounter++;

		ResultType value( numberOfObjectives() );
		int n = numberOfVariables();
                double sum1 = 0.0, sum2 = 0.0, theta = 1.0e-8;
                for(int j = 2; j <= n; j++)
                {
                   double yj = x[j-1] - sin((j*M_PI)/(2.0*n));
                   if(!(j%2)) sum1 += D1(yj, theta);
                   if(j%2) sum2 += D1(yj, theta);
               } 
                value[0] = S2(x[0], 0.06) + sum1;
                value[1] = 1.0 - sqrt(S2(x[0], 0.06)) + sum2;

		return value;
	}
private:
	std::size_t m_objectives;
	BoxConstraintHandler<SearchPointType> m_handler;
	
};

}}

#endif
