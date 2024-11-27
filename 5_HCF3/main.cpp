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
	vector<vector<int> > arcs;
	
	// Create cycles 		
	vector<vector<vector<vector<vector<int> > > > > khcycles(allo.nbNodes);
	for(int i = 0; i < allo.nbNodes;i++){
		khcycles[i].resize(allo.nbNodes);
		for(int j = 0; j < allo.nbNodes;j++){
			khcycles[i][j].resize(allo.K);
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
		curC.push_back({i});
		for(int j = 1; j <= (allo.K+1)/2; j++){
			vector<vector<int> > newC;
			for(int k = 0; k < curC.size();k++){
				for(int l = 0; l< allo.adj[curC[k].back()].size();l++){
					if(allo.adj[curC[k].back()][l] != i && flop[allo.adj[curC[k].back()][l]][i] <= allo.K - j){
						bool add = true;
						for(int m = 1; m <= j -1; m++){
							if(curC[k][m] == allo.adj[curC[k].back()][l]){
								add = false;
								break;
							}
						}
						if(add){
							newC.push_back({curC[k]});
							newC.back().push_back(allo.adj[curC[k].back()][l]);
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
				}
			}
			curC = newC;
		}
	}		
	
	// Only add the HC that can be completed	
	for(int i = 0; i < allo.nbNodes;i++){
		for(int j = i+1; j < allo.nbNodes;j++){
			for(int k = 1; k < allo.K; k++){
				if(khcycles[i][j][k].size() > 0 && khcycles[j][i][k].size() + khcycles[j][i][k-1].size() > 0){
					for(int l = 0; l < khcycles[i][j][k].size();l++){
						hcycles.push_back(khcycles[i][j][k][l]);
					}
					for(int l = 0; l < khcycles[j][i][k].size();l++){
						hcycles.push_back(khcycles[j][i][k][l]);
					}
					if(khcycles[i][j][k-1].size() == 0){
						for(int l = 0; l < khcycles[j][i][k-1].size();l++){
							hcycles.push_back(khcycles[j][i][k-1][l]);
						}
					}
				}
			}
		}
	}
	
	cout << hcycles.size() << endl;
	
/*	for(int i = 0; i < hcycles.size();i++){
		for(int j = 0; j < hcycles[i].size();j++)
			cout << allo.xti[hcycles[i][j]] << " ";
		cout << endl;
	}*/

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

		vector<GRBVar> isHCycleUsed (hcycles.size());
		vector<GRBLinExpr> isNodeUsedS(allo.nbNodes,0);
		vector<GRBLinExpr> isNodeUsedM(allo.nbNodes,0);
		vector<GRBLinExpr> isNodeUsedE(allo.nbNodes,0);
		vector<vector<GRBLinExpr> > isPNodeUsed(allo.nbNodes,vector<GRBLinExpr>(allo.nbNodes,0));
		vector<vector<bool> > isPNodeUsedB(allo.nbNodes,vector<bool>(allo.nbNodes,false));

		vector<GRBVar> isArcUsed (arcs.size());
		vector<vector<GRBLinExpr> > fi(allo.nbNodes, vector<GRBLinExpr>(allo.K,0));
		vector<vector<GRBLinExpr> > fo(allo.nbNodes, vector<GRBLinExpr>(allo.K,0));
		vector<vector<bool> > fb(allo.nbNodes, vector<bool>(allo.K,false));

		// Initialization
		for (int i = 0; i < hcycles.size(); i++){
			isHCycleUsed[i] = model.addVar(0, 1, 0, GRB_BINARY);
		}

		for (int i = 0; i < arcs.size(); i++){
			isArcUsed[i] = model.addVar(0, 1, 0, GRB_BINARY);
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
		}
	
		for (int j = 0; j < arcs.size(); j++){
			isNodeUsedM[arcs[j][0]] += isArcUsed[j];
			if(arcs[j][3] != -1){
				fi[arcs[j][2]][arcs[j][3]] += isArcUsed[j];
				fb[arcs[j][2]][arcs[j][3]] = true;
			}
			else costA += isArcUsed[j];
			fo[arcs[j][0]][arcs[j][1]] += isArcUsed[j];
			fb[arcs[j][0]][arcs[j][1]] = true;
			objFun += isArcUsed[j];
		}		

		// Unique assignment for patients
		for (int i = 0; i < allo.nbNodes; i++){
			model.addConstr(0.5*isNodeUsedS[i] + isNodeUsedM[i] + 0.5*isNodeUsedE[i] <= 1);
		}

		// Pair correspondance 
		for (int i = 0; i < allo.nbNodes; i++){
			for (int j = i+1; j < allo.nbNodes; j++){
				if(isPNodeUsedB[i][j] || isPNodeUsedB[j][i])
					model.addConstr(isPNodeUsed[i][j] == isPNodeUsed[j][i]);
			}
		}

		// Flow conservation
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
