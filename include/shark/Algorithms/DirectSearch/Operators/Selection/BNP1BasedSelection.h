/*!
 *
 *
 * \brief       Indicator-based selection strategy for multi-objective selection.
 *
 *
 *
 * \author      T.Voss
 * \date        2010
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
#ifndef SHARK_ALGORITHMS_DIRECT_SEARCH_OPERATORS_SELECTION_BNP_BASED_SELECTION_H
#define SHARK_ALGORITHMS_DIRECT_SEARCH_OPERATORS_SELECTION_BNP_BASED_SELECTION_H

#include <shark/Algorithms/DirectSearch/Operators/Domination/NonDominatedSort.h>

#include <map>
#include <vector>

namespace shark {

/**
* \brief Implements the well-known indicator-based selection strategy.
*
* See
* Kalyanmoy Deb and Amrit Pratap and Sameer Agarwal and T. Meyarivan,
* A Fast Elitist Multi-Objective Genetic Algorithm: NSGA-II,
* IEEE Transactions on Evolutionary Computation
* Year 2000, Volume 6, p. 182-197
*
* \tparam Indicator The second-level sorting criterion.
*/
template<typename Indicator>
struct BNPBasedSelection {

       typedef shark::Individual<RealVector,RealVector> IndividualType;
       double norm2(RealVector &A, RealVector &B, RealVector & lowerBounds, RealVector & upperBounds)
	{
	   double sum=0.0;
	   int nvar=A.size();
	  for(int i = 0; i < nvar; i++)
	  {
		double low=upperBounds[i]-lowerBounds[i], na=A[i]/low, nb=B[i]/low;
		sum +=  (na-nb)*(na-nb);
	  }
	 return sqrt(sum/(double)nvar);//
	}
	/**
	* \brief Executes the algorithm and assigns each member of the population
	* its level non-dominance (rank) and its individual contribution to the front
	* it belongs to (share).
	*
	* \param [in,out] population The population to be ranked.
	* \param [in,out] mu the number of individuals to select
	*/
////This is a variant of the BNP algorithm
////	template<typename PopulationType>
////	void operator()( PopulationType & population, std::size_t mu, double dt, RealVector & lowerBounds, RealVector & upperBounds){
////		if(population.empty()) return;
////		typedef std::vector< view_reference<typename PopulationType::value_type > > View;
////		while(population.size() > mu)
////		{
////			//compute fronts and ranks
////		        nonDominatedSort(penalizedFitness(population), ranks(population));
////			//detect the nearest pair of individuals..
////			std::pair<double, std::pair<int, int>> near(DBL_MAX, std::make_pair(-1, -1));
////			double adi=0.0, cont=0.0;
////		  	for(int i = 0; i < population.size(); i++)	
////			{
////			  for(int j = i+1; j< population.size(); j++)
////			  {
////				near = std::min(near, std::make_pair(norm2(population[i].searchPoint(), population[j].searchPoint(), lowerBounds, upperBounds),std::make_pair(i,j)));
////			  adi+=norm2(population[i].searchPoint(), population[j].searchPoint(), lowerBounds, upperBounds);
////				cont++;
////			  }
////			}
//////			std::cout << near.first<<" "<<dt<<" " <<adi/cont<<std::endl;
////			if(near.first>=dt)break;
////			std::pair<int,int> cp = near.second;
////			// remove the element with biggest rank  ties are broken by HV contributions in that front
////			if(population[cp.first].rank()>population[cp.second].rank())
////			{
////			   std::iter_swap(population.begin()+cp.first, population.end()-1);
////			   population.pop_back();
////			}
////			else if(population[cp.first].rank()<population[cp.second].rank())
////			{
////			   std::iter_swap(population.begin()+cp.second, population.end()-1);
////			   population.pop_back();
////			}
////			else//same rank
////			{
////			  View front;
////			  std::vector<int>idxs;
////			  for( unsigned int i = 0; i < population.size(); i++ ) 
////				if( population[cp.first].rank() == population[i].rank()) 
////				   front.push_back(population[i]), idxs.push_back(i);
////			  std::size_t frontSize = front.size();
////			  std::vector<KeyValuePair<double,std::size_t> > contrFront = m_indicator.leastContributors(penalizedFitness(front), frontSize);
////			  for(auto idx:contrFront)
////			  {
////			       if(idxs[idx.value]==cp.first)
////				{
////			   	   std::iter_swap(population.begin()+cp.first, population.end()-1);
////			      	   population.pop_back();
////				   break;
////				}
////				if(idxs[idx.value]==cp.second)
////				{
////			   	   std::iter_swap(population.begin()+cp.second, population.end()-1);
////			      	   population.pop_back();
////				   break;
////				}
////			  }
////			}
////		}
////
////		
////		//This is still necessary since that the binary torunament is performed
////		nonDominatedSort(penalizedFitness(population), ranks(population));
////		unsigned int maxRank = 0;
////		std::map< unsigned int, View > fronts;
////
////		for( unsigned int i = 0; i < population.size(); i++ ) {
////			maxRank = std::max( maxRank, population[i].rank());
////			fronts[population[i].rank()].push_back( population[i] );
////			population[i].selected() = true;
////		}
////
////		//deselect the highest rank fronts until we would end up with less or equal mu elements
////		unsigned int rank = maxRank;
////		std::size_t popSize = population.size();
////		
////		while(popSize-fronts[rank].size() >= mu){
////			//deselect all elements in this front
////			View & front = fronts[rank];
////			for(std::size_t i = 0; i != front.size(); ++i){
////				front[i].selected() = false;
////			}
////			popSize -= front.size();
////			--rank;
////		}
////		//now use the indicator to deselect the worst approximating elements of the last selected front
////		View& front = fronts[rank];
////		
////		//create an archive of points which are surely selected because of their domination rank
////		View archive;
////		archive.reserve(popSize - front.size());
////		for(unsigned int r = 1; r != rank; ++r){
////			archive.insert(archive.end(),fronts[r].begin(), fronts[r].end());
////		}
////		//deselect 
////		std::vector<std::size_t> deselected = m_indicator.leastContributors(penalizedFitness(front),penalizedFitness(archive), popSize - mu);
////		for(auto lc:deselected){
////			front[lc].selected() = false;
////		}
////	}

	template<typename PopulationType>
	void operator()( PopulationType & population, std::size_t mu, double dt, RealVector & lowerBounds, RealVector & upperBounds){
		if(population.empty()) return;
		typedef std::vector< view_reference<typename PopulationType::value_type > > View;
		while(population.size() > mu)
		{
			//compute fronts and ranks
		        nonDominatedSort(penalizedFitness(population), ranks(population));
			//detect the nearest pair of individuals..
			std::pair<double, std::pair<int, int>> near(DBL_MAX, std::make_pair(-1, -1));
			double adi=0.0, cont=0.0;
		  	for(int i = 0; i < population.size(); i++)	
			{
			  for(int j = i+1; j< population.size(); j++)
			  {
				near = std::min(near, std::make_pair(norm2(population[i].searchPoint(), population[j].searchPoint(), lowerBounds, upperBounds),std::make_pair(i,j)));
			  adi+=norm2(population[i].searchPoint(), population[j].searchPoint(), lowerBounds, upperBounds);
				cont++;
			  }
			}
//			std::cout << near.first<<" "<<dt<<" " <<adi/cont<<std::endl;
			if(near.first>=dt)break;
			std::pair<int,int> cp = near.second;
			// remove the element with biggest rank  ties are broken by HV contributions in that front
			if(population[cp.first].rank()>population[cp.second].rank())
			{
			   std::iter_swap(population.begin()+cp.first, population.end()-1);
			   population.pop_back();
			}
			else if(population[cp.first].rank()<population[cp.second].rank())
			{
			   std::iter_swap(population.begin()+cp.second, population.end()-1);
			   population.pop_back();
			}
			else//same rank
			{
			  View front;
			  std::vector<int>idxs;
			  for( unsigned int i = 0; i < population.size(); i++ ) 
				if( population[cp.first].rank() == population[i].rank()) 
				   front.push_back(population[i]), idxs.push_back(i);
			  std::size_t frontSize = front.size();
			  std::vector<KeyValuePair<double,std::size_t> > contrFront = m_indicator.leastContributors(penalizedFitness(front), frontSize);
			  for(auto idx:contrFront)
			  {
			       if(idxs[idx.value]==cp.first)
				{
			   	   std::iter_swap(population.begin()+cp.first, population.end()-1);
			      	   population.pop_back();
				   break;
				}
				if(idxs[idx.value]==cp.second)
				{
			   	   std::iter_swap(population.begin()+cp.second, population.end()-1);
			      	   population.pop_back();
				   break;
				}
			  }
			}
		}

		
		//This is still necessary since that the binary torunament is performed
		nonDominatedSort(penalizedFitness(population), ranks(population));
		unsigned int maxRank = 0;
		std::map< unsigned int, View > fronts;

		for( unsigned int i = 0; i < population.size(); i++ ) {
			maxRank = std::max( maxRank, population[i].rank());
			fronts[population[i].rank()].push_back( population[i] );
			population[i].selected() = true;
		}

		//deselect the highest rank fronts until we would end up with less or equal mu elements
		unsigned int rank = maxRank;
		std::size_t popSize = population.size();
		
		while(popSize-fronts[rank].size() >= mu){
			//deselect all elements in this front
			View & front = fronts[rank];
			for(std::size_t i = 0; i != front.size(); ++i){
				front[i].selected() = false;
			}
			popSize -= front.size();
			--rank;
		}
		//now use the indicator to deselect the worst approximating elements of the last selected front
		View& front = fronts[rank];
		
		//create an archive of points which are surely selected because of their domination rank
		View archive;
		archive.reserve(popSize - front.size());
		for(unsigned int r = 1; r != rank; ++r){
			archive.insert(archive.end(),fronts[r].begin(), fronts[r].end());
		}
		//deselect 
		std::vector<std::size_t> deselected = m_indicator.leastContributors(penalizedFitness(front),penalizedFitness(archive), popSize - mu);
		for(auto lc:deselected){
			front[lc].selected() = false;
		}
	}


	/**
	 * \brief Stores/restores the serializer's state.
	 * \tparam Archive Type of the archive.
	 * \param [in,out] archive The archive to serialize to.
	 * \param [in] version number, currently unused.
	 */
	template<typename Archive>
	void serialize( Archive & archive, const unsigned int version )
	{
		archive & BOOST_SERIALIZATION_NVP( m_indicator );
	}
	Indicator& indicator(){
		return m_indicator;
	}
	Indicator const& indicator()const{
		return m_indicator;
	}
private:
	Indicator m_indicator; ///< Instance of the second level sorting criterion.

	/** \cond */
	template<typename T>
	struct view_reference {
		T * mep_value;
	public:
		typedef RealVector FitnessType;
	
		view_reference() : mep_value( NULL ){}
		view_reference(T & value) : mep_value( &value ){}

		operator T & ()
		{
			return( *mep_value );
		}

		operator const T & () const
		{
			return( *mep_value );
		}

		view_reference<T> operator=( const T & rhs )
		{
			*mep_value = rhs;
			return *this;
		}

		RealVector& penalizedFitness() {
			return mep_value->penalizedFitness();
		}
		
		RealVector& unpenalizedFitness(){
			return mep_value->unpenalizedFitness();
		}
		
		RealVector const& penalizedFitness() const{
			return mep_value->penalizedFitness();
		}
		
		RealVector const& unpenalizedFitness() const{
			return mep_value->unpenalizedFitness();
		}
		
		bool& selected(){
			return mep_value->selected();
		}
	};
	/** \endcond */
};
}

#endif
