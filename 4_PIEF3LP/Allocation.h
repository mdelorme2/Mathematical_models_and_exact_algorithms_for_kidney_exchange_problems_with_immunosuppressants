#ifndef ALLOCATION_H
	#define ALLOCATION_H
	
	using namespace std;
	#include <iostream> 
	#include <iomanip> 
	#include <fstream>
	#include <sstream>
	#include <vector>
	#include <string>
	#include <set>
	#include <iostream> 
	#include <math.h> 
	#include <cstdlib>
	#include <algorithm>

	bool sortpair(const pair<int,int>& p1, const pair<int,int>& p2); 
		
	class Info;
	class Allocation;

/*	*************************************************************************************
	************************************* INFO ******************************************
	************************************************************************************* */

	class Info{
	public:
		bool opt;
		vector<double> timeCPU;
		double contUB;
		int LB;
		int UB;
		int nbCons;
		int nbVar;
		int nbNZ;
	};

/*	*************************************************************************************
	********************************** ALLOCATION ***************************************
	************************************************************************************* */

	class Allocation{
	public:
		
		// Data read from the file
		string name;
		int nbNodes;
		vector<vector<int> > edges;
		vector<vector<int> > adj;
		vector<vector<int> > adm;
		vector<int> itx;
		vector<int> xti;

		// Data to be decided by the user
		int K;
		int B;
		
		// Given by the ILP model
		Info infos;

		void load(const string& path, const string& filein);
		void printProb();
		void printInfo(const string& pathAndFileout);
	};
	
	
#endif 