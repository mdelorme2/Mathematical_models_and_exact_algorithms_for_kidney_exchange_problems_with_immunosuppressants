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

	cycle(allo);
	allo.printInfo(pathAndFileout);
}

int cycle(Allocation& allo){

	// Local variables
	vector<vector<int> > cycles;
	vector<int> costs;
		
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
	
	// Create all cycles 
	for(int i = 0; i < allo.nbNodes;i++){
		vector<vector<int> > curC; 
		vector<int> curCost; 
		curC.push_back({i});
		curCost.push_back(0);
		for(int j = 1; j <= allo.K; j++){
			vector<vector<int> > newC;
			vector<int> newCost;
			for(int k = 0; k < curC.size();k++){
				for(int l = i; l< allo.nbNodes;l++){		
					int potNC = 1 - allo.adm[curC[k].back()][l];
					if(curCost[k] + potNC > allo.B) continue;
					if(l == i){
						cycles.push_back(curC[k]);
						costs.push_back(curCost[k]+potNC);
					}
					else{
						bool add = true;
						for(int m = 1; m <= j -1; m++){
							if(curC[k][m] == l){
								add = false;
								break;
							}
						}
						if(add && (curCost[k] + potNC < allo.B || flop[l][i] <= allo.K - j)){
							newC.push_back(curC[k]);
							newCost.push_back(curCost[k] + potNC);
							newC.back().push_back(l);
						}
					}
				}
			}
			curC = newC;
			curCost = newCost;
		}
	}		
	
	cout << cycles.size() << endl;
	
/*	for(int i = 0; i < cycles.size();i++){
		for(int j = 0; j < cycles[i].size();j++)
			cout << allo.xti[cycles[i][j]] << " ";
		cout << "with cost " << costs[i] << endl;
	} */

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
		vector<GRBVar> isCycleUsed (cycles.size());
		vector<GRBLinExpr> isNodeUsed(allo.nbNodes,0);

		// Initialization
		for (int i = 0; i < cycles.size(); i++){
			isCycleUsed[i] = model.addVar(0, 1, 0, GRB_BINARY);
		}

		model.update();

		// Perform values
		for (int i = 0; i < cycles.size(); i++){
			for(int j=0;j<cycles[i].size();j++){
				isNodeUsed[cycles[i][j]] += isCycleUsed[i];
			}
			objFun += cycles[i].size() * isCycleUsed[i];
			costC += costs[i] * isCycleUsed[i];
		}

		// Unique assignment for patients
		for (int i = 0; i < allo.nbNodes; i++){
			model.addConstr(isNodeUsed[i] <= 1);
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
		for (int i = 0; i < cycles.size(); i++){
			if(ceil(isCycleUsed[i].get(GRB_DoubleAttr_X) - EPSILON) == 1){
				for(int j = 0; j < cycles[i].size();j++)
					cout << allo.xti[cycles[i][j]] << " ";
				cout << "with cost " << costs[i] << endl;
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
