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
	double p = atof(argv[6]);
	int seed = atoi(argv[7]);
	
	// functions
	allo.load(path,filein,p,seed);
	allo.K = K;
	allo.B = B;
	allo.p = p;
	allo.printProb();
	allo.infos.timeCPU.push_back(0);
	
	cycle(allo);
	allo.printInfo(pathAndFileout);
}

int cycle(Allocation& allo){

	// Local variables
	vector<vector<int> > cycles;
	vector<vector<int> > cyclesB; vector<int> costs; set<vector<int> > cyclesBS;
	vector<vector<int> > arcs;
	vector<vector<int> > cuts;
	
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

	allo.infos.timeCPU.push_back(getCPUTime() - initTimeModelCPU);
				
	GRBEnv env = GRBemptyenv;
	env.set(GRB_DoubleParam_MemLimit, 60);
	env.start();
	allo.infos.nbIter = 0;
	// Model
	while(1){
		allo.infos.nbIter++;
		try{
			// Local variables
			GRBModel model = GRBModel(env);
			GRBLinExpr objFun = 0;
			GRBLinExpr costA = 0;
			vector<GRBLinExpr> isNodeUsed(allo.nbNodes,0);
			
			vector<GRBVar> isCycleUsed (cycles.size());		
			vector<GRBVar> isCycleBUsed (cyclesB.size());
			vector<GRBVar> isArcUsed (arcs.size());
			vector<vector<GRBLinExpr> > fi(allo.nbNodes, vector<GRBLinExpr>(allo.K,0));
			vector<vector<GRBLinExpr> > fo(allo.nbNodes, vector<GRBLinExpr>(allo.K,0));
			vector<vector<bool> > fb(allo.nbNodes, vector<bool>(allo.K,false));

			// Initialization
			for (int i = 0; i < cycles.size(); i++){
				isCycleUsed[i] = model.addVar(0, 1, 0, GRB_BINARY);
			}
			for (int i = 0; i < cyclesB.size(); i++){
				isCycleBUsed[i] = model.addVar(0, 1, 0, GRB_BINARY);
			}
			for (int i = 0; i < arcs.size(); i++){
				isArcUsed[i] = model.addVar(0, 1, 0, GRB_BINARY);
			}
			
			model.update();

			// Perform values
			for (int i = 0; i < cycles.size(); i++){
				for(int j=0;j<cycles[i].size();j++){
					isNodeUsed[cycles[i][j]] += isCycleUsed[i];
				}
				objFun += cycles[i].size() * isCycleUsed[i];
			}

			for (int i = 0; i < cyclesB.size(); i++){
				for(int j=0;j<cyclesB[i].size();j++){
					isNodeUsed[cyclesB[i][j]] += isCycleBUsed[i];
				}
				objFun += cyclesB[i].size() * isCycleBUsed[i];
				costA += costs[i] * isCycleBUsed[i];
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
			for (int j = 0; j < allo.nbNodes; j++){
				for (int k = 1; k < allo.K; k++){
					if(fb[j][k]){
						model.addConstr(fi[j][k] == fo[j][k]);
					}
				}
			}
					
			// Budget constraint
			model.addConstr(costA <= allo.B);
			
			// No-good cuts
			for (int i = 0; i < cuts.size(); i++){	
				GRBLinExpr LHS = 0;
				for (int j = 0; j < cuts[i].size(); j++){
					LHS += isArcUsed[cuts[i][j]];				
				}
				model.addConstr(LHS <= cuts[i].size() - 1);
			}			
					
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
			vector<int> init;
			vector<int> succ (allo.nbNodes,-1);
			vector<int> succAI (allo.nbNodes,-1);
			for (int j = 0; j < arcs.size(); j++){
				if(ceil(isArcUsed[j].get(GRB_DoubleAttr_X) - EPSILON) == 1){
					succ[arcs[j][0]] = arcs[j][2]; succAI[arcs[j][0]] = j;
					if(arcs[j][1] == 0) init.push_back(arcs[j][0]); 
				}
			}
			
			// prepare the cuts
			cout << "CHAINS: " << endl;
			bool isV = true;
			for (int i = 0; i < init.size(); i++){	
				vector<int> cut;			
				int first = init[i]; int last = first;
				while(succ[last] != -1){
					cout <<  allo.xti[last] << " ";
					cut.push_back(succAI[last]);
					last = succ[last];
				}
				cout <<  allo.xti[last] << " ";
				cut.push_back(succAI[last]);
				if(allo.adm[last][first] == -1){
					cout << "INFEASIBLE, ADD CUT ";
					vector<int> cycleB;
					vector<vector<int> > curC;
					vector<int> curCost; 
					for(int j = 0; j < cut.size();j++){
						cout << cut[j] << "+";
						cycleB.push_back(arcs[cut[j]][0]);
					}
					cout << "<=" << cut.size()-1 << endl;
					isV = false;
					cuts.push_back(cut);
					// Add the right cyclesB, only if B >= 2
					if(allo.B >= 2){
						// Add the first RA
						for(int j = 0; j < allo.nbNodes;j++){
							if(allo.adm[cycleB.back()][j] == 0 && find(cycleB.begin(), cycleB.end(), j) == cycleB.end()){
								curC.push_back(cycleB); curC.back().push_back(j);
								curCost.push_back(1);
							}
						}
						for(int j = cycleB.size(); j <= allo.K-1; j++){
							vector<vector<int> > newC;
							vector<int> newCost;
							for(int k = 0; k < curC.size();k++){
								for(int l = 0; l< allo.nbNodes;l++){		
									int potNC = 1 - allo.adm[curC[k].back()][l];
									if(allo.adm[curC[k].back()][l] == -1 || curCost[k] + potNC > allo.B) continue;
									if(l == cycleB[0] && allo.adm[curC[k].back()][l] == 0){
										// Set in normal form with minimum first
										vector<int> copy = curC[k];
										rotate(copy.begin(), min_element(copy.begin(), copy.end()), copy.end());
										if(cyclesBS.find(copy) == cyclesBS.end()){							
											cyclesB.push_back(copy); cyclesBS.insert(copy);
											costs.push_back(curCost[k]+potNC);
											cout << "EXTRA CYCLEB ->";
											for(int m = 0; m < copy.size();m++){
												cout << allo.xti[copy[m]] << " ";
											}
											cout << "(cost " << costs.back() << ")" << endl;
										}
									}
									else{
										if(find(curC[k].begin(), curC[k].end(), l) == curC[k].end() && curCost[k] + potNC < allo.B){
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
				}
				else{
					cout << "FEASIBLE" << endl;
				}
			}
			if(isV){
				allo.infos.nbAdd = cyclesB.size();
				allo.infos.nbCuts = cuts.size();
				for (int i = 0; i < cycles.size(); i++){
					if(ceil(isCycleUsed[i].get(GRB_DoubleAttr_X) - EPSILON) == 1){
						for(int j = 0; j < cycles[i].size();j++)
							cout << allo.xti[cycles[i][j]] << " ";
						cout << endl;
					}
				}
				cout << "-----" << endl;	
				for (int i = 0; i < cyclesB.size(); i++){
					if(ceil(isCycleBUsed[i].get(GRB_DoubleAttr_X) - EPSILON) == 1){
						for(int j = 0; j < cyclesB[i].size();j++)
							cout << allo.xti[cyclesB[i][j]] << " ";
						cout << "with cost " << costs[i] << endl;
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
				break;
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
	}

	// End
	return 0;
}
