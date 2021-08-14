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
// Implementation of the VSD-SMS-EMOA
#include <shark/Algorithms/DirectSearch/VSD-SMS-EMOA.h>
// Access to benchmark functions
#include <shark/ObjectiveFunctions/Benchmarks/Benchmarks.h>
using namespace shark;		
using namespace shark::benchmarks;		
using namespace std;
void save_front(string saveFilename, VSDSMSEMOA &vsdsmsemoa, bool overwrite=true)
{
        std::fstream fout;
//	if(overwrite)
//	  fout.open(saveFilename.c_str());
//	else
	  fout.open(saveFilename.c_str(),fstream::app|fstream::out );
	for( std::size_t i = 0; i < vsdsmsemoa.solution().size(); i++ ) {
		for( std::size_t j = 0; j < vsdsmsemoa.solution()[i].value.size(); j++ ) {
			fout<< vsdsmsemoa.solution()[ i ].value[j]<<" ";
		}
		fout << std::endl;
	}
	fout.close();
}
void save_pos(string saveFilename, VSDSMSEMOA &vsdsmsemoa, bool overwrite=true)
{
        std::fstream fout;
//	if(overwrite)
//	fout.open(saveFilename.c_str());
//	else
	fout.open(saveFilename.c_str(), fstream::app|fstream::out);
	for( std::size_t i = 0; i < vsdsmsemoa.solution().size(); i++ ) {
		for( std::size_t j = 0; j < vsdsmsemoa.solution()[i].point.size(); j++ ) {
			fout<< vsdsmsemoa.solution()[ i ].point[j]<<" ";
		}
		fout << std::endl;
	}

	fout.close();
}
MultiObjectiveFunction *obj;
VSDSMSEMOA vsdsmsemoa;
int run, iterations;
map<string, MultiObjectiveFunction*> mp;

string strTestInstance, currentPATH;
void setInstance()
{
   mp["DTLZ1"] = new DTLZ1;
   mp["DTLZ2"] = new DTLZ2;
   mp["DTLZ3"] = new DTLZ3;
   mp["DTLZ4"] = new DTLZ4;
   mp["DTLZ5"] = new DTLZ5;
   mp["DTLZ6"] = new DTLZ6;
   mp["DTLZ7"] = new DTLZ7;
   mp["UF1"] = new UF1;
   mp["UF2"] = new UF2;
   mp["UF3"] = new UF3;
   mp["UF4"] = new UF4;
   mp["UF5"] = new UF5;
   mp["UF6"] = new UF6;
   mp["UF7"] = new UF7;
   mp["UF8"] = new UF8;
   mp["UF9"] = new UF9;
   mp["UF10"] = new UF10;
   mp["BT1"] = new BT1;
   mp["BT2"] = new BT2;
   mp["BT3"] = new BT3;
   mp["BT4"] = new BT4;
   mp["BT5"] = new BT5;
   mp["BT6"] = new BT6;
   mp["BT7"] = new BT7;
   mp["BT8"] = new BT8;
   mp["BT9"] = new BT9;
   mp["IMB1"] = new IMB1;
   mp["IMB2"] = new IMB2;
   mp["IMB3"] = new IMB3;
   mp["IMB4"] = new IMB4;
   mp["IMB5"] = new IMB5;
   mp["IMB6"] = new IMB6;
   mp["IMB7"] = new IMB7;
   mp["IMB8"] = new IMB8;
   mp["IMB9"] = new IMB9;
   mp["IMB10"] = new IMB10;
   mp["WFG1"] = new WFG1;
   mp["WFG2"] = new WFG2;
   mp["WFG3"] = new WFG3;
   mp["WFG4"] = new WFG4;
   mp["WFG5"] = new WFG5;
   mp["WFG6"] = new WFG6;
   mp["WFG7"] = new WFG7;
   mp["WFG8"] = new WFG8;
   mp["WFG9"] = new WFG9;
}
void SetConfiguration(int argc, char*argv[])
{
        setInstance();
	currentPATH=".";
	for(int i = 1; i < argc ; i++)
    	{
		string Terminal(argv[i]);
		if( Terminal == "--Instance")
		{
			strTestInstance = string(argv[++i]);
			obj=mp[strTestInstance];
		}
		else if(Terminal == "--Seed")
			run = atoi(argv[++i]);
		else if(Terminal == "--Px")
			vsdsmsemoa.crossoverProbability()=atof(argv[++i]);
		else if(Terminal == "--Di")
			vsdsmsemoa.d0()=atof(argv[++i]);
		else if(Terminal == "--Df")
			vsdsmsemoa.df()=atof(argv[++i]);
		else if(Terminal == "--Path")
			currentPATH =string(argv[++i]);
		else if(Terminal =="--n")
			vsdsmsemoa.mu() = (size_t) atoi(argv[++i]);
		else if(Terminal =="--nobj")
			obj->setNumberOfObjectives(atoi(argv[++i]));
		else if(Terminal == "--nfes")
			iterations= atoi(argv[++i]);
		else if(Terminal == "--nvar")
		{
			obj->setNumberOfVariables(atoi(argv[++i]));
		}
		else if(Terminal == "--param_k")
			obj->setPositionParameters(atoi(argv[++i]));
		else
		{
			cout << Terminal<<endl;
			cout << "Unknown Argument...";
			exit(0);
		}
	    }
        vsdsmsemoa.nc()=2;//eta sbx
        vsdsmsemoa.nm()=50;//eta polynomial
}
int main( int argc, char * argv[] ) {


	srand(run);
        omp_set_num_threads(1);
	SetConfiguration(argc, argv);
	
	obj->init();
	vsdsmsemoa.init(*obj);

  	long long acum = 0;
	string posfix="_"+strTestInstance+"_nobj_"+to_string(obj->numberOfObjectives())+"_nvar_"+to_string(obj->numberOfVariables())+"_seed_"+to_string(run)+"_px_"+to_string(vsdsmsemoa.crossoverProbability())+"_nfes_"+to_string(iterations);
	string filepof=currentPATH+"/POF"+posfix, filepos=currentPATH+"/POS"+posfix;
        save_front(filepof, vsdsmsemoa);
        save_pos(filepos, vsdsmsemoa);
	int i = 0;
	vsdsmsemoa.nIte()=iterations;
	while(i < iterations)
	{
		vsdsmsemoa.currentIte()=i;
		vsdsmsemoa.step(*obj);
		i++;
		acum++;
		if(acum > 0.1*iterations)
		{
		  acum -=0.1*iterations;
		  save_front(filepof, vsdsmsemoa, false);
		  save_pos(filepos, vsdsmsemoa, false);
		}
	}
        save_front(filepof, vsdsmsemoa, false);
	save_pos(filepos, vsdsmsemoa, false);


}
