/*!
 * 
 *
 * \brief       Example for running MO-CMA-ES on an exemplary benchmark function.

 * 
 *
 * \author      tvoss
 * \date        -
 *
 *
 * \par Copyright 1995-2015 Shark Development Team
 * 
 * <BR><HR>
 * This file is part of Shark.
 * <http://image.diku.dk/shark/>
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
// Implementation of the SMS-EMOA
#include <shark/Algorithms/DirectSearch/SMS-EMOA.h>
// Access to benchmark functions
#include <shark/ObjectiveFunctions/Benchmarks/Benchmarks.h>
using namespace shark;		
using namespace shark::benchmarks;		
using namespace std;
int main( int argc, char ** argv ) {
       omp_set_num_threads(1);
        std::size_t mu=100;
        std::size_t iterations=250000;
//       	RealVector reference(3);//	reference(0) = 100;
	UF1 f;
        f.setNumberOfVariables(30);
        f.setNumberOfObjectives(3);

	SMSEMOA smsemoa;
	smsemoa.mu() = mu;
	smsemoa.crossoverProbability()=0.9;
        smsemoa.nc()=2;
        smsemoa.nm()=50;
 
//	smsemoa.indicator().setReference(reference);
	f.init();
	smsemoa.init(f);
	
	for(std::size_t i = 0; i != iterations; ++i){
		smsemoa.step(f);
		//~ if(i%(iterations/100)==0)
			//~ std::clog<<"\r"<<i<<" "<<std::flush;
	}
	for( std::size_t i = 0; i < smsemoa.solution().size(); i++ ) {
		for( std::size_t j = 0; j < f.numberOfObjectives(); j++ ) {
			std::cout<< smsemoa.solution()[ i ].value[j]<<" ";
		}
		std::cout << std::endl;
	}

//	// Iterate the optimizer
//	while( dtlz2.evaluationCounter() < 25000 ) {
//		mocma.step( dtlz2 );
//	}
//
}
