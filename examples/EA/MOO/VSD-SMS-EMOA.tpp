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
#include <shark/Algorithms/DirectSearch/VSD-SMS-EMOA.h>
// Access to benchmark functions
#include <shark/ObjectiveFunctions/Benchmarks/Benchmarks.h>
using namespace shark;		
using namespace shark::benchmarks;		
using namespace std;
void save_front(char saveFilename[1024], VSDSMSEMOA &vsdsmsemoa, bool overwrite=true)
{
        std::fstream fout;
	if(overwrite)
	  fout.open(saveFilename);
	else
	  fout.open(saveFilename,fstream::app|fstream::out );
	for( std::size_t i = 0; i < vsdsmsemoa.solution().size(); i++ ) {
		for( std::size_t j = 0; j < vsdsmsemoa.solution()[i].value.size(); j++ ) {
			fout<< vsdsmsemoa.solution()[ i ].value[j]<<" ";
		}
		fout << std::endl;
	}
	fout.close();
}
int main( int argc, char ** argv ) {
       omp_set_num_threads(1);
        std::size_t mu=100;
        std::size_t iterations=250000;
//       	RealVector reference(3);//	reference(0) = 100;
	IMB1 f;
        f.setNumberOfVariables(10);
        f.setNumberOfObjectives(2);
//	f.setPositionParameters(4);
	VSDSMSEMOA vsdsmsemoa;
	vsdsmsemoa.mu() = mu;
	vsdsmsemoa.crossoverProbability()=0.4;
        vsdsmsemoa.nc()=2;
        vsdsmsemoa.nm()=50;
	vsdsmsemoa.d0()=0.8;
	vsdsmsemoa.df()=0.5;
	vsdsmsemoa.nIte()=iterations;
//	vsdsmsemoa.indicator().setReference(reference);
	f.init();
	vsdsmsemoa.init(f);
	long long acum=0;	
        save_front("out.txt", vsdsmsemoa);
	for(std::size_t i = 0; i != iterations; ++i){
		vsdsmsemoa.currentIte()=i;
		vsdsmsemoa.step(f);
		if(acum > 0.01*iterations)
		{
		  acum -=0.01*iterations;
		  save_front("out.txt", vsdsmsemoa, false);
		}
		acum++;
	}
//	// Iterate the optimizer
//	while( dtlz2.evaluationCounter() < 25000 ) {
//		mocma.step( dtlz2 );
//	}
//
}
