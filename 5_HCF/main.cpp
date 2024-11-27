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

	HCF(allo);
	allo.printInfo(pathAndFileout);
}

int HCF(Allocation& allo){

	// Local variables
	vector<vector<int> > hcycles;
	vector<int> costs;
	vector<vector<vector<vector<vector<int> > > > > khcycles(allo.nbNodes);
	vector<vector<vector<vector<int> > > > khcyclesC(allo.nbNodes);
	vector<vector<vector<int> > > khcyclesMinC(allo.nbNodes);
	for(int i = 0; i < allo.nbNodes;i++){
		khcycles[i].resize(allo.nbNodes);
		khcyclesC[i].resize(allo.nbNodes);
		khcyclesMinC[i].resize(allo.nbNodes);
		for(int j = 0; j < allo.nbNodes;j++){
			khcycles[i][j].resize(allo.K);
			khcyclesC[i][j].resize(allo.K);
			khcyclesMinC[i][j].resize(allo.K,9999);
		}
	}

	// Create Floyd's matrix
	vector<vector<int> > flop (allo.nbNodes, vector<int> (allo.nbNodes, allo.nbNodes));
	for(int i = 0; i < allo.adj.size();i++){
		for(int j = 0; j < allo.adj[i].size();j++)
			flop[i][allo.adj[i][j]] = 1; 	
	}
	
	for(int i = 0; i < allo.nbNodes;i++){
		for(int j = 0; j < allo.nbNodes;j++){
			for(int k = 0; k < allo.nbNodes;k++){
				flop[j][k] = min(flop[j][k],flop[j][i] + flop[i][k]);
			}
		}
	}

	// Create half cycles
	for(int i = 0; i < allo.nbNodes;i++){
		vector<vector<int> > curC; 
		vector<int> curCost; 
		curC.push_back({i});
		curCost.push_back(0);
		for(int j = 1; j <= (allo.K+1)/2; j++){
			vector<vector<int> > newC;
			vector<int> newCost;
			for(int k = 0; k < curC.size();k++){
				for(int l = 0; l< allo.nbNodes;l++){		
					int potNC = 1 - allo.adm[curC[k].back()][l];
					if(l != i && (curCost[k] + potNC < allo.B || (curCost[k] + potNC == allo.B && flop[l][i] <= allo.K - j))){
						bool add = true;
						for(int m = 1; m <= j -1; m++){
							if(curC[k][m] == l){
								add = false;
								break;
							}
						}
						if(add){
							newC.push_back(curC[k]);
							newCost.push_back(curCost[k] + potNC);
							newC.back().push_back(l);
						}
					}
				}
			}
			for(int k = 0; k < newC.size();k++){
				bool add = true;
				for(int l = 1 ; l < newC[k].size()-1; l++){
					if(newC[k][l] < newC[k][0] && newC[k][l] < newC[k].back())
						add = false;
				}
				if(add && (allo.K%2 == 0 || j < (allo.K+1)/2 || newC[k][0] < newC[k].back())){
					khcycles[newC[k][0]][newC[k].back()][newC[k].size()-1].push_back(newC[k]);
					khcyclesC[newC[k][0]][newC[k].back()][newC[k].size()-1].push_back(newCost[k]);
					khcyclesMinC[newC[k][0]][newC[k].back()][newC[k].size()-1] = min(khcyclesMinC[newC[k][0]][newC[k].back()][newC[k].size()-1], newCost[k]);
					/*cout << "Add ";
					for(int l = 0 ; l < newC[k].size(); l++){
						cout << allo.xti[newC[k][l]] << " ";
					}
					cout << "with cost " << newCost[k] << endl;*/
				}
			}
			curC = newC;
			curCost = newCost;
		}
	}		

	// Only add the HC that can be completed
	for(int i = 0; i < allo.nbNodes;i++){
		for(int j = i+1; j < allo.nbNodes;j++){
			for(int k = 1; k < allo.K; k++){
				if(khcycles[i][j][k].size() > 0 && khcycles[j][i][k].size() + khcycles[j][i][k-1].size() > 0 && khcyclesMinC[i][j][k] + min(khcyclesMinC[j][i][k], khcyclesMinC[j][i][k-1]) <= allo.B){
					for(int l = 0; l < khcycles[i][j][k].size();l++){
						hcycles.push_back(khcycles[i][j][k][l]);
						costs.push_back(khcyclesC[i][j][k][l]);
					}
					for(int l = 0; l < khcycles[j][i][k].size();l++){
						hcycles.push_back(khcycles[j][i][k][l]);
						costs.push_back(khcyclesC[j][i][k][l]);
					}
					if(khcycles[i][j][k-1].size() == 0){
						for(int l = 0; l < khcycles[j][i][k-1].size();l++){
							hcycles.push_back(khcycles[j][i][k-1][l]);
							costs.push_back(khcyclesC[j][i][k-1][l]);
						}
					}
				}
			}
		}
	}

	cout << hcycles.size() << endl;

	/*for(int i = 0; i < hcycles.size();i++){
		for(int j = 0; j < hcycles[i].size();j++)
			cout << allo.xti[hcycles[i][j]] << " ";
		cout << endl;
	}*/
	
	allo.infos.timeCPU.push_back(getCPUTime() - initTimeModelCPU);
	
	GRBEnv env = GRBemptyenv;
	env.set(GRB_DoubleParam_MemLimit, 60);
	env.start();
		
	// Model
	try{
		// Local variables
		GRBModel model = GRBModel(env);
		GRBLinExpr objFun = 0;
		GRBLinExpr costC = 0;
		vector<GRBVar> isSelfLoop (allo.nbNodes);
		vector<GRBVar> isHCycleUsed (hcycles.size());
		vector<GRBLinExpr> isNodeUsedS(allo.nbNodes,0);
		vector<GRBLinExpr> isNodeUsedM(allo.nbNodes,0);
		vector<GRBLinExpr> isNodeUsedE(allo.nbNodes,0);
		vector<vector<GRBLinExpr> > isPNodeUsed(allo.nbNodes,vector<GRBLinExpr>(allo.nbNodes,0));
		vector<vector<bool> > isPNodeUsedB(allo.nbNodes,vector<bool>(allo.nbNodes,false));

		// Initialization
		for (int i = 0; i < hcycles.size(); i++){
			isHCycleUsed[i] = model.addVar(0, 1, 0, GRB_BINARY);
		}

		for (int i = 0; i < allo.nbNodes; i++){
			isSelfLoop[i] = model.addVar(0, 1, 0, GRB_BINARY);
		}
		
		model.update();

		// Perform values
		for (int i = 0; i < hcycles.size(); i++){		
			isNodeUsedS[hcycles[i][0]] += isHCycleUsed[i];
			for(int j=1;j<hcycles[i].size()-1;j++){
				isNodeUsedM[hcycles[i][j]] += isHCycleUsed[i];
			}
			isNodeUsedE[hcycles[i].back()] += isHCycleUsed[i];
			isPNodeUsed[hcycles[i][0]][hcycles[i].back()] += isHCycleUsed[i];
			isPNodeUsedB[hcycles[i][0]][hcycles[i].back()] = true;
			objFun += (hcycles[i].size()-1) * isHCycleUsed[i];
			costC += costs[i] * isHCycleUsed[i];
		}
		
		for (int i = 0; i < allo.nbNodes; i++){
			costC += (1 - allo.adm[i][i]) * isSelfLoop[i];
		}			

		// Unique assignment for patients
		for (int i = 0; i < allo.nbNodes; i++){
			model.addConstr(0.5*isNodeUsedS[i] + isNodeUsedM[i] + 0.5*isNodeUsedE[i] + isSelfLoop[i] <= 1);
		}

		// Pair correspondance 
		for (int i = 0; i < allo.nbNodes; i++){
			for (int j = i+1; j < allo.nbNodes; j++){
				if(isPNodeUsedB[i][j] || isPNodeUsedB[j][i])
					model.addConstr(isPNodeUsed[i][j] == isPNodeUsed[j][i]);
			}
		}
		
		// Budget constraint
		model.addConstr(costC <= allo.B);
			
		// Objective function
		model.setObjective(objFun, GRB_MAXIMIZE);
		
		// Setting of Gurobi
		model.getEnv().set(GRB_DoubleParam_TimeLimit,  3600);
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
		for (int i = 0; i < hcycles.size(); i++){
			if(ceil(isHCycleUsed[i].get(GRB_DoubleAttr_X) - EPSILON) == 1){
				for(int j = 0; j < hcycles[i].size();j++)
					cout << allo.xti[hcycles[i][j]] << " ";
				cout << endl;
			}
		}
		cout << " --- and selfloop --- " << endl;
		for (int i = 0; i < allo.nbNodes; i++){
			if(ceil(isSelfLoop[i].get(GRB_DoubleAttr_X) - EPSILON) == 1){
				cout << allo.xti[i] << endl;
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
