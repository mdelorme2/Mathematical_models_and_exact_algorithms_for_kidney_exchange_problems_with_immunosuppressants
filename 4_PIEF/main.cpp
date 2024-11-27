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
	
	PIEF(allo);
	allo.printInfo(pathAndFileout);
}

int PIEF(Allocation& allo){

	// Local variables
	vector<vector<vector<int> > > edges(allo.nbNodes);
		
	for(int i = 0; i < allo.nbNodes;i++){
		vector<vector<vector<bool> > > ee (allo.nbNodes,vector<vector<bool> >(allo.nbNodes,vector<bool>(allo.K,false)));
				
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
			for(int k1 = i; k1 < allo.nbNodes;k1++){
				for(int k2 = i; k2 < allo.nbNodes;k2++){
					if(j == 0 && k1 != i) continue;
					if(j == allo.K - 1 && k2 != i) continue;
					if(j > 0 && k1 == i) continue;
					if(tails[j][k1] + (1 - allo.adm[k1][k2]) + heads[j][k2] <= allo.B){
						if(j < allo.K - 1) tails[j+1][k2] = min(tails[j+1][k2],tails[j][k1] + (1 - allo.adm[k1][k2]));
						if(ee[k1][k2][j] == false){
							ee[k1][k2][j] = true;
						//	cout << "ADD " << allo.xti[i] << " " << allo.xti[k1] << " " << allo.xti[k2]  << " pos " << j << endl;
						//	cout << "BECAUSE " << tails[j][k1] << " + " << (1 - allo.adm[k1][k2]) << " + " << heads[j][k2] << " <= " << allo.B << endl;
							if(k2 == i) edges[i].push_back({k1,j,k2,0});
							else edges[i].push_back({k1,j,k2,j+1});
						}
					}
				}
			}
		}		
	}	

	allo.infos.timeCPU.push_back(getCPUTime() - initTimeModelCPU);
	
	GRBEnv env = GRBemptyenv;
	env.set(GRB_DoubleParam_MemLimit, 60);
	env.start();

	// Model
	try{
		// Local variables
		GRBModel model = GRBModel(env);
		GRBLinExpr objFun = 0;
		GRBLinExpr costE = 0;
		vector<vector<GRBVar> > isEdgeUsed (allo.nbNodes);
		vector<vector<vector<GRBLinExpr> > > fi(allo.nbNodes, vector<vector<GRBLinExpr> > (allo.nbNodes, vector<GRBLinExpr> (allo.K,0)));
		vector<vector<vector<GRBLinExpr> > > fo(allo.nbNodes, vector<vector<GRBLinExpr> > (allo.nbNodes, vector<GRBLinExpr> (allo.K,0)));
		vector<vector<vector<bool> > > fb(allo.nbNodes, vector<vector<bool> > (allo.nbNodes, vector<bool> (allo.K,false)));
		vector<GRBLinExpr> isNodeUsed(allo.nbNodes,0);

		// Initialization
		for (int i = 0; i < allo.nbNodes; i++){
			isEdgeUsed[i].resize(edges[i].size());
			for (int j = 0; j < edges[i].size();j++){
				isEdgeUsed[i][j] = model.addVar(0, 1, 0, GRB_BINARY);
			}
		}

		model.update();

		// Perform values
		for (int i = 0; i < allo.nbNodes; i++){
			for (int j = 0; j < edges[i].size(); j++){
				isNodeUsed[edges[i][j][2]] += isEdgeUsed[i][j];
				fi[i][edges[i][j][2]][edges[i][j][3]] += isEdgeUsed[i][j];
				fo[i][edges[i][j][0]][edges[i][j][1]] += isEdgeUsed[i][j];
				fb[i][edges[i][j][0]][edges[i][j][1]] = true;
				fb[i][edges[i][j][2]][edges[i][j][3]] = true;
				objFun += isEdgeUsed[i][j];
				costE += (1 - allo.adm[edges[i][j][0]][edges[i][j][2]]) * isEdgeUsed[i][j];
			}
		}

		// Flow conservation and use node once
		for (int i = 0; i < allo.nbNodes; i++){
			for (int k = 0; k <= allo.K-1; k++){
				for (int j = 0; j < allo.nbNodes; j++){
					if(fb[i][j][k]){
						model.addConstr(fi[i][j][k] == fo[i][j][k]);
					}
				}
			}
			model.addConstr(isNodeUsed[i] <= 1);
		}
				
		// Budget constraint	
		model.addConstr(costE <= allo.B);
		
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
			for (int j = 0; j < edges[i].size(); j++){			
				if(ceil(isEdgeUsed[i][j].get(GRB_DoubleAttr_X) - EPSILON) == 1)
					cout << allo.xti[i] << " " << allo.xti[edges[i][j][0]] << " " << edges[i][j][1] << " " << allo.xti[edges[i][j][2]] << " " << edges[i][j][3] << endl;
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
