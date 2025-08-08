#pragma once
#include <fstream>
#include "Solution.h"
#include <cstring>
#include <iostream>
#include <ctime>

#include "Graph.h"
using namespace std;
void preFileName(char *FileName);
void getFilePath(char *File_Path,const char *Path,char *FileName);

void preFileName(char *FileName) {
	int i, j = 0;
	int len = strlen(FileName);
	for (i = len - 1; i >= 0; --i) {
		if (FileName[i] == '/')break;
	}
	if (i == 0)return;
	for (++i; i < len; ++i) {
		FileName[j++] = FileName[i];
	}
	FileName[j] = '\0';
}

void getFilePath(char *File_Path, const char *Path, char *FileName) {
	strcpy(File_Path, Path);
	strcat(File_Path, FileName);
}

//std::ofstream Monitor(char *FileName) {
//	std::ofstream fout;
//	boost::filesystem::path filepath = boost::filesystem::current_path() / "Data_Out/";
//	if (!boost::filesystem::exists(filepath)) {
//		boost::filesystem::create_directory(filepath);
//	}
//	//    filepath/=FileName;
//	//    if(!boost::filesystem::exists(filepath)){
//	//        boost::filesystem::create_directory(filepath);
//	//    }
//	//    time_t t=time(nullptr);
//		//char tmp[100];
//		//strftime(tmp, sizeof(tmp),"%Y_%m_%d_%H.%M.%S",localtime(&t));
//		//filepath/=tmp;
//	filepath /= "PerTurbCompare.txt";
//	std::string Fileout = filepath.string();
//	fout.open(Fileout, std::ios::app);
//	if (fout.fail()) {
//		std::cout << "Error In Monitor Path" << std::endl;
//		exit(-1);
//	}
//	return fout;
//}

//bool CheckSolution(Solution &Sol) {
//	Graph *G = Graph::GetGraph();
//	int node, neighbor;
//	int testSol = 0;
//	std::vector<int>::iterator iter;
//	std::cout << "------------------Check Solution-------------------" << std::endl;
//	for (node = 1; node <= G->Node_Num; ++node) {
//		if (Sol.BestSolutions[node] == PASSIVE)continue;
//		iter = G->Adj_List[node].begin();
//		G->Node_State[node] = 0;
//		while (iter != G->Adj_List[node].end()) {
//			neighbor = *iter++;
//			if (Sol.BestSolutions[neighbor] == PASSIVE) {
//				++G->Node_State[node];
//			}
//		}
//		if (G->Node_State[node] != 0) {
//			std::cout << node << "(" << G->Node_State[node] << "|" << G->Gain_Out[node] << "),";
//			++testSol;
//		}
//		G->Node_State[node] = 0;
//	}
//	std::cout << std::endl;
//	std::cout << "-------------------Check Finish--------------------" << std::endl;
//	if (testSol == Sol.BestSol) {
//		return true;
//	}
//	else {
//		return false;
//	}
//}