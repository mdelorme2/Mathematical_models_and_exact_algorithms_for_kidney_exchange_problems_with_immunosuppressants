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
	vector<vector<int> > arcs;
	vector<double> RCC;
	vector<double> RCA;
	
	// Create Floyd's matrix
	vector<vector<int> > flop (allo.nbNodes, vector<int> (allo.nbNodes, allo.nbNodes));
	for(int i = 0; i < allo.adj.size();i++){
		for(int j = 0; j < allo.adj[i].size();j++)
			flop[i][allo.adj[i][j]] = 1; 	
	}
	
	for(int i = 0; i < allo.nbNodes;i++){
		for(int j = 0; j < allo.nbNodes;j++){
			for(int k = 0; k < allo.nbNodes;k++)
				flop[j][k] = min(flop[j][k],flop[j][i] + flop[i][k]);
		}
	}

	/*for(int i = 0; i < flop.size();i++){
		for(int j = 0; j < flop[i].size();j++)
			cout << flop[i][j] << "\t";
		cout << endl;
	}*/
	
	// Create cycles
	for(int i = 0; i < allo.nbNodes;i++){
		vector<vector<int> > curC;
		curC.push_back({i});
		for(int j = 1; j <= allo.K; j++){
			vector<vector<int> > newC;
			for(int k = 0; k < curC.size();k++){
				for(int l = 0; l< allo.adj[curC[k].back()].size();l++){
					if(allo.adj[curC[k].back()][l] >= i){
						if(allo.adj[curC[k].back()][l] == i){
							cycles.push_back({curC[k]});
						}
						else{
							if(flop[allo.adj[curC[k].back()][l]][i] <= allo.K - j){
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
				}
			}
			curC = newC;
		}
	}		
	
	cout << cycles.size() << endl;
	RCC.resize(cycles.size(),0.0);
	
	/*for(int i = 0; i < cycles.size();i++){
		for(int j = 0; j < cycles[i].size();j++)
			cout << allo.xti[cycles[i][j]] << " ";
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
	RCA.resize(arcs.size(),0.0);
	
	allo.infos.timeCPU.push_back(getCPUTime() - initTimeModelCPU);
		
	// Model
	try{
		GRBEnv env = GRBemptyenv;
		env.set(GRB_DoubleParam_MemLimit, 30);
		env.start();
	
		// Local variables
		GRBModel model = GRBModel(env);
		GRBLinExpr objFun = 0;
		GRBLinExpr costA = 0;
		vector<GRBLinExpr> isNodeUsed(allo.nbNodes,0);
		
		vector<GRBVar> isCycleUsed (cycles.size());		
		vector<GRBVar> isArcUsed (arcs.size());
		vector<vector<GRBLinExpr> > fi(allo.nbNodes, vector<GRBLinExpr>(allo.K,0));
		vector<vector<GRBLinExpr> > fo(allo.nbNodes, vector<GRBLinExpr>(allo.K,0));
		vector<vector<bool> > fb(allo.nbNodes, vector<bool>(allo.K,false));

		// Initialization
		for (int i = 0; i < cycles.size(); i++){
			isCycleUsed[i] = model.addVar(0, 1, 0, GRB_CONTINUOUS);
		}

		for (int i = 0; i < arcs.size(); i++){
			isArcUsed[i] = model.addVar(0, 1, 0, GRB_CONTINUOUS);
		}
		
		model.update();

		// Perform values
		for (int i = 0; i < cycles.size(); i++){
			for(int j=0;j<cycles[i].size();j++){
				isNodeUsed[cycles[i][j]] += isCycleUsed[i];
			}
			objFun += cycles[i].size() * isCycleUsed[i];
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

		// Unique assignment for patients
		for (int i = 0; i < allo.nbNodes; i++){
			model.addConstr(isNodeUsed[i] <= 1);
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
		model.getEnv().set(GRB_IntParam_Crossover, 0); 
		model.getEnv().set(GRB_IntParam_Threads, 1);
		model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
		model.optimize();

		// Filling Solution
		for (int i = 0; i < cycles.size(); i++){
			if(isCycleUsed[i].get(GRB_DoubleAttr_X) < EPSILON){
			//	cout << "Cycle " << i << " " << isCycleUsed[i].get(GRB_DoubleAttr_RC) << endl;
				RCC[i] = isCycleUsed[i].get(GRB_DoubleAttr_RC);
			}
			else
				RCC[i] = 0.0;
		}
		// cout << "-----" << endl;		
		for (int j = 0; j < arcs.size(); j++){
			if(isArcUsed[j].get(GRB_DoubleAttr_X) < EPSILON){
				// cout << "Arc " << j << " " << isArcUsed[j].get(GRB_DoubleAttr_RC) << endl;
				RCA[j] = isArcUsed[j].get(GRB_DoubleAttr_RC);
			}
			else
				RCA[j] = 0.0;				
		}

		allo.infos.contUB = model.get(GRB_DoubleAttr_ObjVal);
		allo.infos.UB = floor(model.get(GRB_DoubleAttr_ObjVal) + EPSILON);	
	}

	// Exceptions
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
		
		// Filling Info
		allo.infos.timeCPU[0] = getCPUTime() - initTimeModelCPU;
		allo.infos.opt = false;
		allo.infos.LB = 0;
		return -1;
	}
	catch (...) {
		cout << "Exception during optimization" << endl;
	}
	
	while(1){	
		// Local variables
		try{
			GRBEnv env = GRBemptyenv;
			env.set(GRB_DoubleParam_MemLimit, 60);
			env.start();
		
			// Local variables
			GRBModel model = GRBModel(env);
			GRBLinExpr objFun = 0;
			GRBLinExpr costA = 0;
			vector<GRBLinExpr> isNodeUsed(allo.nbNodes,0);
			
			vector<GRBVar> isCycleUsed (cycles.size());	
			vector<bool> isCycleActivated (cycles.size(),false);			
			vector<GRBVar> isArcUsed (arcs.size());
			vector<bool> isArcActivated (arcs.size(),false);
			vector<vector<GRBLinExpr> > fi(allo.nbNodes, vector<GRBLinExpr>(allo.K,0));
			vector<vector<GRBLinExpr> > fo(allo.nbNodes, vector<GRBLinExpr>(allo.K,0));
			vector<vector<bool> > fb(allo.nbNodes, vector<bool>(allo.K,false));

			// Initialization
			for (int i = 0; i < cycles.size(); i++){
				if(allo.infos.contUB + RCC[i] + EPSILON >= allo.infos.UB){
					isCycleUsed[i] = model.addVar(0, 1, 0, GRB_BINARY);
					isCycleActivated[i] = true;
				}
			}

			for (int j = 0; j < arcs.size(); j++){
				if(allo.infos.contUB + RCA[j] + EPSILON >= allo.infos.UB){
					isArcUsed[j] = model.addVar(0, 1, 0, GRB_BINARY);
					isArcActivated[j] = true;
				}
			}
			
			model.update();

			// Perform values
			for (int i = 0; i < cycles.size(); i++){
				if(isCycleActivated[i]){ 
					for(int j=0;j<cycles[i].size();j++){
						isNodeUsed[cycles[i][j]] += isCycleUsed[i];
					}
					objFun += cycles[i].size() * isCycleUsed[i];
				}
			}

			for (int j = 0; j < arcs.size(); j++){
				if(isArcActivated[j]){
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
			}

			// Unique assignment for patients
			for (int i = 0; i < allo.nbNodes; i++){
				model.addConstr(isNodeUsed[i] <= 1);
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
			model.getEnv().set(GRB_DoubleParam_TimeLimit,  3600 - (getCPUTime() - initTimeModelCPU));
			model.getEnv().set(GRB_IntParam_Method, 2);
			model.getEnv().set(GRB_IntParam_Threads, 1);
			model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
			model.getEnv().set(GRB_IntParam_MIPFocus, 1);
			model.optimize();
			
			if(ceil(model.get(GRB_DoubleAttr_ObjVal) - EPSILON) == allo.infos.UB){
				// Filling Info
				allo.infos.timeCPU[0] = getCPUTime() - initTimeModelCPU;
				allo.infos.nbVar =  model.get(GRB_IntAttr_NumVars);
				allo.infos.nbCons = model.get(GRB_IntAttr_NumConstrs);
				allo.infos.nbNZ = model.get(GRB_IntAttr_NumNZs);
				allo.infos.LB = allo.infos.UB;	
				allo.infos.opt = true;

				// Filling Solution
				for (int i = 0; i < cycles.size(); i++){
					if(isCycleActivated[i]){
						if(ceil(isCycleUsed[i].get(GRB_DoubleAttr_X) - EPSILON) == 1){
							for(int j = 0; j < cycles[i].size();j++)
								cout << allo.xti[cycles[i][j]] << " ";
							cout << endl;
						}
					}
				}
				cout << "-----" << endl;		
				for (int j = 0; j < arcs.size(); j++){
					if(isArcActivated[j]){
						if(ceil(isArcUsed[j].get(GRB_DoubleAttr_X) - EPSILON) == 1){
							if(arcs[j][3] != -1)
								cout << allo.xti[arcs[j][0]] << " " << arcs[j][1] << " " << allo.xti[arcs[j][2]] << " " << arcs[j][3] << endl;
							else
								cout << allo.xti[arcs[j][0]] << " " << arcs[j][1] << " " << arcs[j][2] << " " << arcs[j][3] << endl;
						}
					}
				}
				break;
			}
			else{
				if(model.get(GRB_IntAttr_Status) == 9){
					// Filling Info
					allo.infos.timeCPU[0] = getCPUTime() - initTimeModelCPU;
					allo.infos.opt = false;

					// Get Info
					allo.infos.nbVar =  model.get(GRB_IntAttr_NumVars);
					allo.infos.nbCons = model.get(GRB_IntAttr_NumConstrs);
					allo.infos.nbNZ = model.get(GRB_IntAttr_NumNZs);
					allo.infos.LB = 0;
					return -1;
				}
				else allo.infos.UB--;
			}
		}
		
		// Exceptions
		catch (GRBException e) {
			cout << "Error code = " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
		
			// Filling Info
			allo.infos.timeCPU[0] = getCPUTime() - initTimeModelCPU;
			allo.infos.opt = false;
			allo.infos.LB = 0;
			return -1;
		}
		catch (...) {
			cout << "Exception during optimization" << endl;
		}
	}

	// End
	return 0;
}