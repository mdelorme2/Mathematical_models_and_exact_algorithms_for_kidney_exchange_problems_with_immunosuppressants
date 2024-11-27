#include "main.h"

/*	*************************************************************************************
	*************************************  MAIN *****************************************
	************************************************************************************* */

double initTimeModelCPU;

int main(int argc, char **argv){
	
	initTimeModelCPU = getCPUTime();
	
	// local variables
	Allocation allo ;
	string filein = argv[2];
	string path = argv[1];
	string pathAndFileout = argv[3];
	int K = atoi(argv[4]);
	int B = atoi(argv[5]);
	
	// functions
	allo.load(path,filein);
	allo.K = K;
	allo.B = B;
	allo.printProb();
	allo.infos.timeCPU.push_back(0);
	
	EEF(allo);
	allo.printInfo(pathAndFileout);
}

int EEF(Allocation& allo){

	// Local variables
	vector<vector<vector<int> > > anedges(allo.nbNodes);
	vector<vector<int> > arcs;
	
	// Create cycles 
	for(int i = 0; i < allo.nbNodes;i++){
		vector<vector<bool> > ee (allo.nbNodes,vector<bool>(allo.nbNodes,false));
				
		// Back
		vector<vector<int> > heads (allo.K,vector<int>(allo.nbNodes,1));
		for(int j = allo.K-1; j >= 1; j--){
			heads[j][i] = 0;
			for(int k = 0; k < allo.edges.size();k++){
				if(heads[j][allo.itx[allo.edges[k][1]]] == 0 && allo.itx[allo.edges[k][0]] >= i){
					heads[j-1][allo.itx[allo.edges[k][0]]] = 0;
				}
			}
		}
		heads[0][i] = 0;
		
		// Front 
		vector<vector<int> > tails (allo.K,vector<int>(allo.nbNodes,allo.K+1));
		tails[0][i] = 0;
		for(int j = 0; j < allo.K; j++){
			for(int k = 0; k < allo.edges.size();k++){
				if(tails[j][allo.itx[allo.edges[k][0]]] + heads[j][allo.itx[allo.edges[k][1]]] == 0){
					if(j < allo.K - 1) tails[j+1][allo.itx[allo.edges[k][1]]] = 0;
					if(ee[allo.itx[allo.edges[k][0]]][allo.itx[allo.edges[k][1]]] == false){
						ee[allo.itx[allo.edges[k][0]]][allo.itx[allo.edges[k][1]]] = true;
					//	cout << "ADD " << allo.xti[i] << " " << allo.edges[k][0] << " " << allo.edges[k][1] << " pos " << j << endl;
						anedges[i].push_back({allo.itx[allo.edges[k][0]],allo.itx[allo.edges[k][1]]});
					}
				}
			}
		}		
	}	

	// Create chains	
	if (allo.B >= 1){
		vector<bool> tails(allo.nbNodes,true);
		// Add initial closing arcs
		for(int k = 0; k < allo.nbNodes;k++){
			arcs.push_back({k,0,-1,-1});
			//cout << "ADD " << allo.xti[k] << " " << 0 << " " << "-1" << " " << "-1" << endl;
		}
		for(int j = 0; j < allo.K - 1; j++){
			vector<bool> tailsN(allo.nbNodes,false);
			// Add arcs in position j
			for(int k = 0; k < allo.edges.size();k++){
				if(tails[allo.itx[allo.edges[k][0]]]){
					tailsN[allo.itx[allo.edges[k][1]]] = true;
					arcs.push_back({allo.itx[allo.edges[k][0]],j,allo.itx[allo.edges[k][1]],j+1});
					//cout << "ADD " << allo.edges[k][0] << " " << j << " " << allo.edges[k][1] << " " << j + 1 << endl;
				}
			}
			// Add closing arcs
			for(int k = 0; k < allo.nbNodes;k++){
				if(tailsN[k]){
					arcs.push_back({k,j+1,-1,-1});
					//cout << "ADD " << allo.xti[k] << " " << j+1 << " " << "-1" << " " << "-1" << endl;
				}
			}
			// Update the reachable nodes
			tails = tailsN;
		}
	}
	
	cout << arcs.size() << endl;

	allo.infos.timeCPU.push_back(getCPUTime() - initTimeModelCPU);
	
	GRBEnv env = GRBemptyenv;
	env.set(GRB_DoubleParam_MemLimit, 60);
	env.start();

	// Model
	try{
		// Local variables
		GRBModel model = GRBModel(env);
		GRBLinExpr objFun = 0;
		GRBLinExpr costA = 0;
		
		vector<vector<GRBVar> > isEdgeUsed (allo.nbNodes);
		vector<vector<GRBLinExpr> > fiE(allo.nbNodes, vector<GRBLinExpr>(allo.nbNodes,0));
		vector<vector<GRBLinExpr> > foE(allo.nbNodes, vector<GRBLinExpr>(allo.nbNodes,0));
		vector<vector<bool> > fbE(allo.nbNodes, vector<bool>(allo.nbNodes,false));
		vector<GRBLinExpr> isNodeUsed(allo.nbNodes,0);

		vector<GRBVar> isArcUsed (arcs.size());
		vector<vector<GRBLinExpr> > fi(allo.nbNodes, vector<GRBLinExpr>(allo.K,0));
		vector<vector<GRBLinExpr> > fo(allo.nbNodes, vector<GRBLinExpr>(allo.K,0));
		vector<vector<bool> > fb(allo.nbNodes, vector<bool>(allo.K,false));
		
		// Initialization
		for (int i = 0; i < allo.nbNodes; i++){
			isEdgeUsed[i].resize(anedges[i].size());
			for (int j = 0; j < anedges[i].size();j++){
				isEdgeUsed[i][j] = model.addVar(0, 1, 0, GRB_BINARY);
			}
		}

		for (int i = 0; i < arcs.size(); i++){
			isArcUsed[i] = model.addVar(0, 1, 0, GRB_BINARY);
		}
		
		model.update();

		// Perform values
		for (int i = 0; i < allo.nbNodes; i++){
			for (int j = 0; j < anedges[i].size(); j++){
				isNodeUsed[anedges[i][j][1]] += isEdgeUsed[i][j];
				fiE[i][anedges[i][j][1]] += isEdgeUsed[i][j];
				foE[i][anedges[i][j][0]] += isEdgeUsed[i][j];
				fbE[i][anedges[i][j][1]] = true;
				objFun += isEdgeUsed[i][j];
			}
		}

		for (int j = 0; j < arcs.size(); j++){
			isNodeUsed[arcs[j][0]] += isArcUsed[j];
			if(arcs[j][3] != -1){
				fi[arcs[j][2]][arcs[j][3]] += isArcUsed[j];
				fb[arcs[j][2]][arcs[j][3]] = true;
			}
			else costA += isArcUsed[j];
			fo[arcs[j][0]][arcs[j][1]] += isArcUsed[j];
			fb[arcs[j][0]][arcs[j][1]] = true;
			objFun += isArcUsed[j];
		}

		// Flow conservation in cycle structure and use node once
		for (int i = 0; i < allo.nbNodes; i++){
			GRBLinExpr sum = 0;
			for (int j = 0; j < allo.nbNodes; j++){
				if(fbE[i][j]){
					sum += fiE[i][j];
					model.addConstr(fiE[i][j] == foE[i][j]);
					// model.addConstr(fiE[i][j] <= foE[i][i]); // Version 1 -> without, Version 2 -> with
				}
			}
			model.addConstr(isNodeUsed[i] <= 1);
			model.addConstr(sum <= allo.K * foE[i][i]);
		}
		
		// Flow conservation in chain structure and compute budget
		GRBLinExpr sum = 0;
		for (int j = 0; j < allo.nbNodes; j++){
			for (int k = 1; k < allo.K; k++){
				if(fb[j][k]){
					model.addConstr(fi[j][k] == fo[j][k]);
				}
			}
			sum += fo[j][0];
		}
		model.addConstr(sum == costA);
		
		// Budget constraint	
		model.addConstr(costA <= allo.B);
		
		// Objective function
		model.setObjective(objFun, GRB_MAXIMIZE);
				
		// Setting of Gurobi
		model.getEnv().set(GRB_DoubleParam_TimeLimit, 3600);
		model.getEnv().set(GRB_IntParam_Method, 2);
		model.getEnv().set(GRB_IntParam_Threads, 1);
		model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
		model.optimize();
		
		// Filling Info
		allo.infos.timeCPU[0] = getCPUTime() - initTimeModelCPU;
		allo.infos.UB = ceil(model.get(GRB_DoubleAttr_ObjBound) - EPSILON);
		allo.infos.opt = false;

		// Get Info pre preprocessing
		allo.infos.nbVar =  model.get(GRB_IntAttr_NumVars);
		allo.infos.nbCons = model.get(GRB_IntAttr_NumConstrs);
		allo.infos.nbNZ = model.get(GRB_IntAttr_NumNZs);

		// If no solution found
		if (model.get(GRB_IntAttr_SolCount) < 1){
			cout << "Failed to optimize ILP. " << endl;
			allo.infos.LB  = 0;
			return -1;
		}

		// If solution found
		allo.infos.LB = ceil(model.get(GRB_DoubleAttr_ObjVal) - EPSILON);	
		if(allo.infos.LB == allo.infos.UB) allo.infos.opt = true;

		// Filling Solution
		for (int i = 0; i < allo.nbNodes; i++){
			for (int j = 0; j < anedges[i].size(); j++){			
				if(ceil(isEdgeUsed[i][j].get(GRB_DoubleAttr_X) - EPSILON) == 1)
					cout << allo.xti[i] << " " << allo.xti[anedges[i][j][0]] << " " << allo.xti[anedges[i][j][1]] << endl;
			}
		}
		cout << "-----" << endl;		
		for (int j = 0; j < arcs.size(); j++){
			if(ceil(isArcUsed[j].get(GRB_DoubleAttr_X) - EPSILON) == 1){
				if(arcs[j][3] != -1)
					cout << allo.xti[arcs[j][0]] << " " << arcs[j][1] << " " << allo.xti[arcs[j][2]] << " " << arcs[j][3] << endl;
				else
					cout << allo.xti[arcs[j][0]] << " " << arcs[j][1] << " " << arcs[j][2] << " " << arcs[j][3] << endl;
			}
		}
	}

	// Exceptions
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Exception during optimization" << endl;
	}


	// End
	return 0;
}
