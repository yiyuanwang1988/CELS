#include <iostream>
#include "Graph.h"
#include "Solution.h"
#include "AGENS.h"
#include "file.h"
#include "iostream"
#include <fstream>

using namespace std;
int main(int argc, char** argv) {
    Param Parameters(argc, argv);
   
    Graph* G = Graph::GetGraph();
    G->GraphInit(Parameters.File_In);
	G->fname = Parameters.File_Name;

	/*
	vector<int> s;
	for (int node = 1; node <= G->Node_Num; node++) {

		bool flag = false;
		for (int i = 0; i < G->Adj_List[node].size(); i++) {//判断当前点是否与他的邻居的非共享顶点都<=x
			int a = G->Adj_List[node][i];

			if (G->contrast_Matrix[a][node] <= G->average) {
				flag = true;
				break;
			}
		}
		if (flag) {
			s.push_back(node);
		}
	}
	cout << G->fname << " " << s.size() << " " << G->Node_Num << endl;
	
	*/
	srand(Parameters.rseed);
    Solution Sol(Parameters, Parameters.rseed);
	file = Parameters.File_Name;
    if (Parameters.InitMethod == "HCluster") Sol.Timer->StartClock();
    AGENS A(20);
		
	
	//A.HierarchicalCluster_Min();
	/*for (int i = 0; i < 10000; i++) {
		A.findMaximalCliqueRandom();
		A.findMaximalCliqueRandom();
		for (int i = 0; i <= G->Node_Num; i++) {
			clique_flag[i] = 0;
		}
	}*/

	if (Parameters.InitMethod != "HCluster")Sol.Timer->StartClock();
	if (Parameters.InitMethod == "Random")Sol.GenerateInitSol_Random();
	else if (Parameters.InitMethod == "Greedy")Sol.GenerateInitSol_Greedy();
	else if (Parameters.InitMethod == "C1")Sol.GenerateInitSol_C1(Parameters.alpha);
	else if (Parameters.InitMethod == "C2")Sol.GenerateInitSol_C2(Parameters.alpha);
	else if (Parameters.InitMethod == "HCluster") {
		Sol.GenerateInitSol_HCluster();
	}
	else if (Parameters.InitMethod == "Check") {
		if (Sol.CheckBestSolution()) {
			std::cout << "Right Solution" << std::endl;
			std::cout << "---------------------------------------" << std::endl;
			return 0;
		}
		else {
			std::cout << "Wrong Solution" << std::endl;
			std::cout << "---------------------------------------" << std::endl;
			return -1;
		}
	}
	else {
		std::cout << "Error In InitMethod!" << std::endl;
		exit(-1);
	}
    

    Sol.Timer->RunTime();
    Sol.InitSol = Sol.CurSol;
    Sol.RLS();

	//std::ofstream file("example.txt", std::ios::app);

	/*
		if (file.is_open()) {
		// 向文件写入内容
		file << Parameters.File_Name << ' ' << Sol.BestSol  << " " << Sol.CurEpoch << endl;

		// 关闭文件
		file.close();

	}
	else {
		std::cout << "Unable to open the file." << std::endl;
	}
	*/

	//for (int i = 1; i <= G->Node_Num; i++) {
	//	cout << Sol.freqq[i] << " ";
	//}
	//cout << endl;
	//cout << Parameters.File_Name << ' ' << Sol.BestSol << " " << ave1 / sol_vec1.size() << endl;


//	cout << Parameters.File_Name << ' ' << Sol.BestSol << " " << aver_a / sol_vec_a.size() << " " << aver_b / sol_vec_b.size() << endl;

	

	cout << "File name:" << Parameters.File_Name << endl;
	cout << "Best sol: " << Sol.BestSol << endl;
	cout<< "The best time is: " << Sol.Timer->GetTime() << " " << endl;

	cout << "B is as follows :";
	Graph *G1 = Graph::GetGraph();
	int node, neighbor;
	int testSol = 0;
	std::vector<int>::iterator iter;
	for (node = 1; node <= G1->Node_Num; ++node) {
		if (Sol.BestSolutions[node] == ACTIVE)
			cout << node << " ";
		
	}
	cout << endl;
	cout << "B' is as follows :";

	for (node = 1; node <= G1->Node_Num; ++node) {
		if (Sol.BestSolutions[node] == PASSIVE)
			cout << node << " ";

	}
	cout << endl;


    return 0;
}