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
	template<typename PopulationType>
	void operator()( PopulationType & population, std::size_t mu, double dt, RealVector & lowerBounds, RealVector & upperBounds){
		if(population.size() == mu)return;
		if(population.empty()) return; // I don't know why this is here!...
		typedef std::vector< view_reference<typename PopulationType::value_type > > View;
		int nobj = population[0].unpenalizedFitness().size(), npops=population.size();
	      if(dt>0)
	      {
		PopulationType penalized, selected, candidates=population;
		while(selected.size() < mu && !candidates.empty())
  	        {
		   PopulationType all;
		   for(auto sind:selected)all.push_back(sind);
		   for(auto cind:candidates)all.push_back(cind);
		   nonDominatedSort(unpenalizedFitness(all), ranks(all));
	
		   for(int i = 0; i < selected.size(); i++)selected[i].rank() = all[i].rank();
		   for(int i = 0; i < candidates.size(); i++)candidates[i].rank() = all[selected.size()+i].rank();
		  
		   int lowestfront=npops;
		   RealVector ref( nobj, -DBL_MAX);
		   for(auto cind:candidates)
		   {
			 if(lowestfront > (int)cind.rank())
			 {
			    lowestfront = (int) cind.rank();
			    for(int m = 0; m < nobj; m++)
			      ref[m]= cind.unpenalizedFitness()[m];
		 	 }
			 else if(lowestfront == (int)cind.rank())
			    for(int m = 0; m < nobj; m++) ref[m]=std::max(ref[m], cind.unpenalizedFitness()[m]);
		   }
		
		   for(int m = 0; m < nobj; m++) ref[m] +=1.0;
		   m_indicator.setReference(ref);
		   std::pair<double, int> maxHV(-DBL_MAX, -1);
		   PopulationType frontSelected;
		   for(auto sind:selected) if(sind.rank()==lowestfront) frontSelected.push_back(sind); 

		   for(int i = 0; i < candidates.size(); i++)
		   {
		       if(candidates[i].rank()!=lowestfront) continue;
		        auto cind=candidates[i];
			frontSelected.push_back(cind);
		        std::vector<KeyValuePair<double,std::size_t> > contribSelec = m_indicator.leastContributors(unpenalizedFitness(frontSelected), frontSelected.size());
		       for(auto cs:contribSelec){
			 if(cs.value==frontSelected.size()-1) maxHV=max(maxHV, std::make_pair(cs.key, i));
			}
		       frontSelected.pop_back();
		   }
		   selected.push_back(candidates[maxHV.second]);
		   std::iter_swap(candidates.begin()+maxHV.second, candidates.end()-1);
		   candidates.pop_back();
		   for(int i = candidates.size()-1; i>=0; i--)
		   {
		        if(norm2(selected.back().searchPoint(), candidates[i].searchPoint(), lowerBounds, upperBounds) <= dt)
		        {
		           penalized.push_back(candidates[i]);
		           std::iter_swap(candidates.begin()+i, candidates.end()-1);
		           candidates.pop_back();
		        }
		   }
		} 
		std::vector<double> dists(penalized.size(), DBL_MAX);
		for(int i = 0; i < penalized.size(); i++)
		  for(int j = 0; j < selected.size(); j++)
			dists[i]=std::min(dists[i], norm2(selected[j].searchPoint(), penalized[i].searchPoint(), lowerBounds, upperBounds));
		while(selected.size() < mu)
		{
		   std::pair<double, int> maxdcn(-1, -1);
		   for(int i = 0; i < penalized.size(); i++)
			maxdcn=max(maxdcn, std::make_pair(dists[i], i));
		  selected.push_back(penalized[maxdcn.second]);
		  std::iter_swap(penalized.begin()+maxdcn.second, penalized.end()-1);
		  std::iter_swap(dists.begin()+maxdcn.second, dists.end()-1);
		  penalized.pop_back();
		  dists.pop_back();
		}
		  // static long ite=0;
		  // ite++;
		  // if(ite%1000==0)std::cout<<ite<<std::endl;
		population=selected;

		   nonDominatedSort(unpenalizedFitness(population), ranks(population));
		return;
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
	       for(int i = population.size()-1; i>=0; i--)
		if(! population[i].selected())
		{
		   iter_swap(population.begin()+i, population.end()-1);
		   population.pop_back();
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
