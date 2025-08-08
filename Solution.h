#pragma once

#include <fstream>
#include "Graph.h"
#include "Bucket.h"
#include "tool.h"
#include <algorithm>
#include <iostream>
#include <string>
#include <cstring>
#include "file.h"
#include "AGENS.h"

class Param {
public:
	int argc;
	char **argv;
	char *File_Name;
	char *File_In;
	const char* InPath = "Instancias/";
	int cutoffTime;
	int MaxBatch;
	int unbetterBound;
	int test_time;
	double alpha;
	std::string InitMethod;
	std::string SelectMethod;
	std::string PerturbMethod;
	unsigned int rseed = 0;
	std::string outPutDir;

	Param(int paramNum, char **params);
	static unsigned int getRandSeed();
	~Param();
};

Param::Param(int paramNum, char **params) {
	/**
	 InitMethod: Random; Greedy; C1; C2;
	 SelectMethod: First; Random; Prob;
	 PerturbMethod: Random; MAB; Prob;
	**/
	argc = paramNum;
	argv = params;
	if (argc == 10) {
		File_Name = argv[1];
		cutoffTime = atoi(argv[2]);
		MaxBatch = atoi(argv[3]);
		unbetterBound = atoi(argv[4]);
		test_time = atoi(argv[5]);
		InitMethod = argv[6];
		alpha = atof(argv[7]);
		SelectMethod = argv[8];
		PerturbMethod = argv[9];
	}
	else if (argc == 6) {
		File_Name = argv[1];
		cutoffTime = atoi(argv[5]);
		MaxBatch = 10000;
		unbetterBound = 150;
		test_time = 0;
		InitMethod = argv[4];
		alpha = 0.5;
		SelectMethod = "First";
		PerturbMethod = argv[3];
		rseed = atoi(argv[2]);
	}
	else if (paramNum == 4) {
		File_Name = argv[1];
		cutoffTime = atoi(argv[3]);
		MaxBatch = 10000;
		unbetterBound = 250;
		test_time = atoi(argv[3]);
		InitMethod = "HCluster";
		alpha = 0.5;
		// some instances are run with Random
		SelectMethod = "First";
		PerturbMethod = "HCluster";
		outPutDir = "1";
		// use the random seed
		rseed = atoi(argv[2]);
	}
	else if (argc == 2) {
		File_Name = argv[1];
		cutoffTime = 100;
		MaxBatch = 10000;
		unbetterBound = 250;
		test_time = 0;
		InitMethod = "HCluster";
		alpha = 0.5;
		SelectMethod = "First";
		PerturbMethod = "HCluster";
		rseed = getRandSeed();;
	}
	else if (argc == 3) {
		File_Name = argv[1];
		cutoffTime = 100;
		MaxBatch = 10000;
		unbetterBound = atoi(argv[2]);
		test_time = 0;
		InitMethod = "HCluster";
		alpha = 0.5;
		SelectMethod = "First";
		PerturbMethod = "HCluster";
		rseed = getRandSeed();;
	}
	else {
		std::cout << "Error in parameters!" << std::endl;
		exit(-1);
	}
	preFileName(File_Name);
	File_In = (char*)malloc((strlen(InPath) + strlen(File_Name) + 1) * sizeof(char));
	getFilePath(File_In, InPath, File_Name);
}

unsigned int Param::getRandSeed() {
	FILE *fs_p = NULL;
	unsigned int seed = 0;

	fs_p = fopen("/dev/urandom", "r");
	if (NULL == fs_p)
	{
		printf("Can not open /dev/urandom\n");
	}
	else {
		fread(&seed, sizeof(int), 1, fs_p);  //obtain one unsigned int data
		fclose(fs_p);
	}
	return seed;
}

Param::~Param() {
	free(File_In);
}

class Solution{
public:
    std::vector<int> Solutions;//标记点是在左边还是右边，B是ACTIVE B'是PASSIVE
    std::vector<int> BestSolutions;//最优解的状态
    std::vector<int> LocalBestSolutions;//当前最优
    long long TotalIter;//当前是第几循环
    int CurSol,BestSol,LocalBestSol,InitSol;
    int CurIter,BestIter;
    int CurBatch,BestBatch;
    int CurEpoch,BestEpoch;

    int SolutionNum,CurSolutionNum;
    int Max_Batch;
    int PerTurbNum;
    int UnBetterBound;
    int AMoveNodeNum;
    int BMoveNodeNum;

	vector <int> freqq;
	long long total_tabu = 0;
	long long total_cur_sol = 0;
	long long	 stepp = 0;
    std::string SelectMethod,InitMethod;
    timeTool *Timer;
    randTool *Rand;
    Bucket *Bk_Head;
    Bucket *BK_Tail_A;
    Bucket *BK_Tail_B;
    Graph *G;
	AGENS *A;
    std::ofstream fMonitor;
    std::string PerturbMethod;
    int testNum;
    explicit Solution(Param &param);
    Solution(Param &param,unsigned int rseed);
    ~Solution();
    bool CheckCurSol();
    bool CheckBucket();
    bool CheckGainOut();

	bool checkAMoveNodeNum(int a,int b);
	int SelectMoveNodePair(int TabuStep);
	int CalPairScore(int v1, int v2);//v1在B,v2在B
    bool CheckBestSolution();
    void StoreInitSolutions();
    void StoreSolutions();
    void StoreLocalSolutions();
    void InitBucket();
    void GenerateInitSol_Random();
    void GenerateInitSol_Greedy();
    void GenerateInitSol_C1(double a);
    void GenerateInitSol_C2(double a);
    void GenerateInitSol_HCluster();
//    void GenerateInitSol_HClusterSize();
//    void GenerateInitSol_HClusterGreedy();

    void MoveBucketNode(Bucket *point,int node,int direction);
    void AdjustBucket(int set,int direction,int movecount);
    void AdjustVertexList(int node);
    int SelectMoveNode(int direct);
	void SelectMoveNode3(int direct);
	int SelectMoveNode1(int direct);
	int SelectMoveNode2(int direct);
    void DropNode(int dropnode);
    void AddNode(int addnode);
    void TabuSearch();
    void PerTurbMultiDA();
    void PerTurbHCluster();
    void RLS();
	void testbug();
};

void Solution::testbug() {
	for (int i = 1; i < G->Clustering.size(); ++i) {
		if (G->ClusterState[i] == ACTIVE)continue;
		if (G->ClusterNumInASet[i] != G->Clustering[i].size() && G->ClusterNumInASet[i] != 0 && perturb_flag[G->Clustering[i][0]] == 0) {
			std::cout << file << " G->ClusterNumInASet[i] != G->Clustering[i].size() && G->ClusterNumInASet[i] != 0 "<< G->ClusterNumInASet[i] <<" " << G->Clustering[i].size() << std::endl;
		}
	}
}

bool cmp_degree(int a, int b) {
	Graph *G = Graph::GetGraph();
	return G->Node_Degree[a] < G->Node_Degree[b];
}

bool cmp_MoveGain(int a, int b) {
	Graph *G = Graph::GetGraph();
	return G->MoveGain[a] < G->MoveGain[b];
}

bool cmp_ClusterNum(int a, int b) {
	Graph *G = Graph::GetGraph();
	return G->Clustering[a].size() < G->Clustering[b].size();
}

bool cmp_GainIn(int a, int b) {
	Graph *G = Graph::GetGraph();
	return G->Gain_In[a] > G->Gain_In[b];
}

bool cmp_GainOut(int a, int b) {
	Graph *G = Graph::GetGraph();
	return G->Gain_Out[a] < G->Gain_Out[b];
}

Solution::~Solution() {
	delete Timer;
	delete Rand;
	while (Bk_Head->pre != NULL) {
		delete Bk_Head->pre;
	}
	while (Bk_Head->next != NULL) {
		delete Bk_Head->next;
	}
	delete Bk_Head;
}

Solution::Solution(Param &param) {
	G = Graph::GetGraph();
	Timer = new timeTool(param.cutoffTime);
	Rand = new randTool();
	Bk_Head = new Bucket(None);
	BK_Tail_A = new Bucket(G->Node_Num);
	BK_Tail_B = new Bucket(-G->Node_Num);
	A = new AGENS(100);
	Bk_Head->next = BK_Tail_A;
	BK_Tail_A->pre = Bk_Head;
	Bk_Head->pre = BK_Tail_B;
	BK_Tail_B->next = Bk_Head;
	Max_Batch = param.MaxBatch;
	UnBetterBound = param.unbetterBound;
	testNum = param.test_time;
	Solutions.resize(G->Node_Num + 1, -1);
	BestSolutions.resize(G->Node_Num + 1);
	LocalBestSolutions.resize(G->Node_Num + 1);
	TotalIter = 0;
	PerTurbNum = 0;
	SolutionNum = G->Node_Num / 2;
	CurSolutionNum = 0;
	CurIter = 0;
	CurBatch = 0;
	CurEpoch = 0;
	InitMethod = param.InitMethod;
	SelectMethod = param.SelectMethod;
	PerturbMethod = param.PerturbMethod;
}

Solution::Solution(Param &param, unsigned int rseed) {
	G = Graph::GetGraph();
	Timer = new timeTool(param.cutoffTime);
	Rand = new randTool(rseed);
	Bk_Head = new Bucket(None);
	BK_Tail_A = new Bucket(G->Node_Num);
	BK_Tail_B = new Bucket(-G->Node_Num);
	Bk_Head->next = BK_Tail_A;
	BK_Tail_A->pre = Bk_Head;
	Bk_Head->pre = BK_Tail_B;
	BK_Tail_B->next = Bk_Head;
	Max_Batch = param.MaxBatch;
	testNum = param.test_time;
	Solutions.resize(G->Node_Num + 1, -1);
	BestSolutions.resize(G->Node_Num + 1, 0);
	LocalBestSolutions.resize(G->Node_Num + 1);
	TotalIter = 0;
	PerTurbNum = 0;
	SolutionNum = G->Node_Num / 2;
	//    UnBetterBound=param.unbetterBound;
	UnBetterBound = SolutionNum;
	CurSolutionNum = 0;
	CurIter = 0;
	CurBatch = 0;
	CurEpoch = 0;
	InitMethod = param.InitMethod;
	SelectMethod = param.SelectMethod;
	PerturbMethod = param.PerturbMethod;
}

void Solution::StoreInitSolutions() {
	int i;
	for (i = 1; i <= G->Node_Num; ++i) {
		BestSolutions[i] = Solutions[i];
	}
	BestSol = CurSol;
	BestEpoch = CurEpoch;
	BestBatch = CurBatch;
	BestIter = CurIter;
	Timer->FinishClock();
}

void Solution::StoreSolutions() {
	int i;
	for (i = 1; i <= G->Node_Num; ++i) {
		BestSolutions[i] = Solutions[i];
	}
	BestSol = CurSol;
	BestEpoch = CurEpoch;
	BestBatch = CurBatch;
	BestIter = CurIter;
	Timer->FinishClock();
}

void Solution::StoreLocalSolutions() {
	int i;
	for (i = 1; i <= G->Node_Num; ++i) {
		LocalBestSolutions[i] = Solutions[i];
	}
	LocalBestSol = CurSol;
	if (LocalBestSol < BestSol) {
		BestIter = CurIter;
		Timer->FinishClock();
	}
}

void Solution::InitBucket() {
	int i, node, neighbor;
	std::vector<int>::iterator iter;
	Bucket *point = Bk_Head;
	CurSol = 0;
	for (i = 0; i < SolutionNum; ++i) {
		node = G->VertexList[i];//// VertexList[i - 1] = i;  Node_Index[i] = i - 1;
		if (Solutions[node] != ACTIVE) {
			std::cout << "Error IN Solutions" << "InitBucket" << std::endl;
			exit(-1);
		}
		if (G->Node_Index[node] != i) {
			std::cout << "Error IN NodeIndex" << std::endl;
			exit(-1);
		}
		G->Node_State[i] = node;//
		if (G->Gain_Out[node] != 0) {//node与B'由连接,则node移动到B'CurSol--
			--G->MoveGain[node];
			++CurSol;
		}
		if (G->Gain_Out[node] == 0) {//若node不与B'中任意顶点相邻，则其邻居移除过去会导致CurSol++
			for (iter = G->Adj_List[node].begin(); iter != G->Adj_List[node].end(); ++iter) {
				neighbor = *iter;
				++G->MoveGain[neighbor];
			}
		}
		else if (G->Gain_Out[node] == 1) {//node仅仅与B'中的一个顶点v'相邻,则将v'从B'移动到B,B中的点不属于C,CurSol--
			for (iter = G->Adj_List[node].begin(); iter != G->Adj_List[node].end(); ++iter) {
				neighbor = *iter;
				if (Solutions[neighbor] == PASSIVE) {//B是ACTIVE B'是PASSIVE
					--G->MoveGain[neighbor];
					break;
				}
			}
		}
		else {
			continue;
		}
	}

	for (; i < G->Node_Num; i++) {
		node = G->VertexList[i];
		if (Solutions[node] != PASSIVE) {//B'中的点一定是PASSIVE
			std::cout << "Error IN Solutions" << std::endl;
			exit(-1);
		}
		if (G->Node_Index[node] != i) {
			std::cout << "Error IN NodeIndex" << std::endl;
			exit(-1);
		}
		G->Node_State[i] = node;//下标对应的点
		if (G->Gain_Out[node] != 0)++G->MoveGain[node];//node为B'中的顶点，且跟B'中其他顶点有连接，则移过去node为C中的一个顶点，CurSol++;
	}
	std::sort(G->Node_State.begin(), G->Node_State.begin() + SolutionNum, cmp_MoveGain);//B与B'分别排序
	std::sort(G->Node_State.begin() + SolutionNum, G->Node_State.begin() + G->Node_Num, cmp_MoveGain);//B与B'分别排序
	for (i = 0; i < SolutionNum; ++i) {
		node = G->Node_State[i];//可以看作Vertexlist的copy数组
		G->Node_State[i] = None;
		if (G->MoveGain[node] + 2 > point->BucketId) {
			point = new Bucket(point, G->MoveGain[node] + 2);//     当前Bucket newBucket
			point->BucketList.push_back(node);
		}
		else if (G->MoveGain[node] + 2 == point->BucketId) {
			point->BucketList.push_back(node);
		}
		else {
			std::cout << "Error In InitBucket!" << std::endl;
			exit(-1);
		}
	}
	point = Bk_Head;
	for (i = G->Node_Num - 1; i >= SolutionNum; --i) {
		node = G->Node_State[i];
		G->Node_State[i] = None;
		if (G->MoveGain[node] - 2 < point->BucketId) {
			point = new Bucket(G->MoveGain[node] - 2, point);
			point->BucketList.push_back(node);
		}
		else if (G->MoveGain[node] - 2 == point->BucketId) {
			point->BucketList.push_back(node);
		}
		else {
			std::cout << "Error In InitBucket!" << std::endl;
			exit(-1);
		}
	}
}

void Solution::GenerateInitSol_Random() {
	int i, node, neighbor;
	while (CurSolutionNum < SolutionNum) {
		i = Rand->RandInt(CurSolutionNum, G->Node_Num - 1);
		node = G->VertexList[i];
		G->VertexList[i] = G->VertexList[CurSolutionNum];
		G->Node_Index[G->VertexList[i]] = i;
		G->VertexList[CurSolutionNum] = node;
		G->Node_Index[node] = CurSolutionNum++;
		Solutions[node] = ACTIVE;
		G->AdjustGain(node, 1);
	}
	InitBucket();
	StoreInitSolutions();
}


void Solution::GenerateInitSol_Greedy() {
	int i = 0, j, node, neighbor;
	std::vector<int>::iterator iter;
	std::sort(G->VertexList.begin(), G->VertexList.begin() + G->Node_Num, cmp_degree);
	while (CurSolutionNum < SolutionNum) {
		node = G->VertexList[i++];
		if (Solutions[node] == PASSIVE) {
			Solutions[node] = ACTIVE;
			++CurSolutionNum;
			if (CurSolutionNum == SolutionNum) {
				for (iter = G->Adj_List[node].begin(); iter != G->Adj_List[node].end(); ++iter) {
					neighbor = *iter;
					++G->Gain_In[neighbor];
					--G->Gain_Out[neighbor];
				}
				break;
			}
			for (iter = G->Adj_List[node].begin(); iter != G->Adj_List[node].end(); ++iter) {
				neighbor = *iter;
				++G->Gain_In[neighbor];
				--G->Gain_Out[neighbor];
				if (Solutions[neighbor] == ACTIVE || CurSolutionNum == SolutionNum)continue;
				Solutions[neighbor] = ACTIVE;
				G->AdjustGain(neighbor, 1);
				++CurSolutionNum;
			}
		}
		else {
			for (iter = G->Adj_List[node].begin(); iter != G->Adj_List[node].end(); ++iter) {
				neighbor = *iter;
				if (Solutions[neighbor] == ACTIVE)continue;
				Solutions[neighbor] = ACTIVE;
				G->AdjustGain(neighbor, 1);
				++CurSolutionNum;
				if (CurSolutionNum == SolutionNum)break;
			}
		}
	}
	i = 0;
	j = G->Node_Num - 1;
	for (node = 1; node <= G->Node_Num; ++node) {
		if (Solutions[node] == ACTIVE) {
			G->VertexList[i] = node;
			G->Node_Index[node] = i++;
		}
		else {
			G->VertexList[j] = node;
			G->Node_Index[node] = j--;
		}
	}
	InitBucket();
	StoreInitSolutions();
}

void Solution::GenerateInitSol_C1(double a) {
	int i, node, index;
	int maxGainIn, minGainIn, th, RCLnum;
	i = Rand->RandInt(0, G->Node_Num - 1);
	node = G->VertexList[i];
	G->VertexList[i] = G->VertexList[CurSolutionNum];
	G->Node_Index[G->VertexList[i]] = i;
	G->VertexList[CurSolutionNum] = node;
	G->Node_Index[node] = CurSolutionNum++;
	Solutions[node] = ACTIVE;
	G->AdjustGain(node, 1);
	while (CurSolutionNum < SolutionNum) {
		maxGainIn = 0;
		minGainIn = G->Node_Num;
		for (i = CurSolutionNum; i < G->Node_Num; ++i) {
			node = G->VertexList[i];
			if (G->Gain_In[node] > maxGainIn)maxGainIn = G->Gain_In[node];
			if (G->Gain_In[node] < minGainIn)minGainIn = G->Gain_In[node];
		}
		th = minGainIn + static_cast<int>(a*(maxGainIn - minGainIn));
		RCLnum = 0;
		for (i = CurSolutionNum; i < G->Node_Num; ++i) {
			node = G->VertexList[i];
			if (G->Gain_In[node] >= th)G->Node_State[RCLnum++] = node;
		}
		i = Rand->RandInt(0, RCLnum - 1);
		node = G->Node_State[i];
		index = G->Node_Index[node];
		G->VertexList[index] = G->VertexList[CurSolutionNum];
		G->Node_Index[G->VertexList[index]] = index;
		G->VertexList[CurSolutionNum] = node;
		G->Node_Index[node] = CurSolutionNum++;
		Solutions[node] = ACTIVE;
		G->AdjustGain(node, 1);
	}
	/*for(i=0;i<G->Node_Num;i++){
		G->Node_State[i]=0;
	}*/
	InitBucket();
	StoreInitSolutions();
}

void Solution::GenerateInitSol_C2(double a) {
	int i, node, num, tmp;
	int RCLnum, index, minGainIn, canNode;
	i = Rand->RandInt(0, G->Node_Num - 1);
	node = G->VertexList[i];
	G->VertexList[i] = G->VertexList[CurSolutionNum];
	G->Node_Index[G->VertexList[i]] = i;
	G->VertexList[CurSolutionNum] = node;
	G->Node_Index[node] = CurSolutionNum++;
	Solutions[node] = ACTIVE;
	G->AdjustGain(node, 1);
	while (CurSolutionNum < SolutionNum) {
		tmp = static_cast<int>(a*(G->Node_Num - CurSolutionNum));
		num = tmp > 1 ? tmp : 1;
		minGainIn = G->Node_Num;
		canNode = 0;
		for (i = 0; i < num; ++i) {
			index = Rand->RandInt(CurSolutionNum, G->Node_Num - 1);
			if (minGainIn > G->Gain_In[G->VertexList[index]]) {
				minGainIn = G->Gain_In[G->VertexList[index]];
				canNode = G->VertexList[index];
			}
		}
		index = G->Node_Index[canNode];
		G->VertexList[index] = G->VertexList[CurSolutionNum];
		G->Node_Index[G->VertexList[index]] = index;
		G->VertexList[CurSolutionNum] = canNode;
		G->Node_Index[canNode] = CurSolutionNum++;
		Solutions[canNode] = ACTIVE;
		G->AdjustGain(canNode, 1);
	}
	InitBucket();
	StoreInitSolutions();
}


//产生初始解，
void Solution::GenerateInitSol_HCluster() {
	Graph* G = Graph::GetGraph();
	int i, j, node, neighbor, num;
	int clusterNo, index, curSolNum;
	std::vector<int> ClusterOrder;
	std::vector<int>::iterator iter;

	for (i = 1; i < G->Clustering.size(); ++i) {//wang
		std::sort(G->Clustering[i].begin(), G->Clustering[i].end(), cmp_degree);//每个Cluster按照度排序
		G->ClusterNumInASet[i] = 0;
	}

	for (i = 1; i <= 3; ++i) {//wang
		if (G->Clustering[i].size() == 0)
			continue;
		ClusterOrder.push_back(G->Clustering[i][0]);//找到每个组里面度数最小的点
	}

	std::sort(ClusterOrder.begin(), ClusterOrder.end(), cmp_degree);
	for (iter = ClusterOrder.begin(); iter != ClusterOrder.end(); ++iter) {
		node = *iter;
		clusterNo = G->Node_Cluster[node];

		num = G->Clustering[clusterNo].size() - G->ClusterNumInASet[clusterNo];////非固定的cluster在A的顶点个数 G->ClusterNumInASet[clusterNo]此时为0 wang
		
	//	cout << clusterNo << " " << G->Clustering[clusterNo].size() << " " << G->ClusterNumInASet[clusterNo] << " " << num << " " << CurSolutionNum << " " << SolutionNum << endl;
		if (num + CurSolutionNum > SolutionNum)continue;
		else {
			G->ClusterNumInASet[clusterNo] = G->Clustering[clusterNo].size();//index为clusterNo的Cluster在B中的点的个数
			for (i = 0; i < G->Clustering[clusterNo].size(); ++i) {
				node = G->Clustering[clusterNo][i];
				if (Solutions[node] == ACTIVE)continue;
				else {
					Solutions[node] = ACTIVE;//若顶点在B,则为ACTIVE，若顶点在B',则为PASSIVE
					G->AdjustGain(node, 1);//调整邻居的Gain_in Gain_out
					index = G->Node_Index[node];// VertexList[i - 1] = i;  Node_Index[i] = i - 1;
					G->VertexList[index] = G->VertexList[CurSolutionNum];//CurSolutionNum初始化为0
					G->Node_Index[G->VertexList[index]] = index;
					G->VertexList[CurSolutionNum] = node;
					G->Node_Index[node] = CurSolutionNum++;
				}
			}
			if (CurSolutionNum == SolutionNum)break;
			for (i = 0; i < G->Clustering[clusterNo].size(); ++i) {
				node = G->Clustering[clusterNo][i];
				if (G->Gain_In[node] == G->Node_Degree[node])continue;
				else {
					for (j = 0; j < G->Adj_List[node].size(); ++j) {
						neighbor = G->Adj_List[node][j];
						if (Solutions[neighbor] == ACTIVE)continue;
						++G->ClusterNumInASet[G->Node_Cluster[neighbor]];
						Solutions[neighbor] = ACTIVE;
						G->AdjustGain(neighbor, 1);
						index = G->Node_Index[neighbor];
						G->VertexList[index] = G->VertexList[CurSolutionNum];
						G->Node_Index[G->VertexList[index]] = index;
						G->VertexList[CurSolutionNum] = neighbor;
						G->Node_Index[neighbor] = CurSolutionNum++;
						if (CurSolutionNum == SolutionNum)break;
					}
				}
				if (CurSolutionNum == SolutionNum)break;
			}
		}
	}

	if (CurSolutionNum < SolutionNum) {
		std::sort(G->VertexList.begin() + CurSolutionNum, G->VertexList.end(), cmp_GainIn);
		for (i = CurSolutionNum; i < SolutionNum; ++i) {
			node = G->VertexList[i];
			++G->ClusterNumInASet[G->Node_Cluster[node]];
			Solutions[node] = ACTIVE;
			G->AdjustGain(node, 1);
			G->Node_Index[node] = i;
		}
		for (; i < G->Node_Num; ++i) {
			node = G->VertexList[i];
			G->Node_Index[node] = i;
		}
		CurSolutionNum = SolutionNum;
	}

	if (PerturbMethod == "HCluster") {
		AMoveNodeNum = G->ClusterNumInASet[0];
		BMoveNodeNum = G->Clustering[0].size() - AMoveNodeNum;
		for (i = 1; i <= 3;  ++i) {
			//如果这个组在左边wang

			if (G->ClusterNumInASet[i] == 0 || G->Clustering[i].size() == G->ClusterNumInASet[i])G->ClusterState[i] = PASSIVE;
			else {
				//std::cout << G->fname << " G->ClusterNumInASet[i] != G->Clustering[i].size() && G->ClusterNumInASet[i] != 0 " << std::endl;
				AMoveNodeNum += G->ClusterNumInASet[i];
				BMoveNodeNum += G->Clustering[i].size() - G->ClusterNumInASet[i];
			}			
		}
	
	}

	int ta = SolutionNum;
	int tb = G->Node_Num - SolutionNum;
	for (int node = 1; node <= G->Node_Num; node++) {
		if (G->ClusterState[G->Node_Cluster[node]] == PASSIVE && perturb_flag[node] == 0) {
			if (Solutions[node] == ACTIVE)ta--;
			else tb--;
		}
	}
	AMoveNodeNum = ta;
	BMoveNodeNum = tb;

	//if (checkAMoveNodeNum(AMoveNodeNum, BMoveNodeNum)) {
	//	//cout << "init is ok" << endl;
	//}

	InitBucket();
	StoreInitSolutions();
	//cout << 444 << endl;
	//testbug();

}

void Solution::AdjustVertexList(int node) {
	int tmp;
	if (Solutions[node] == ACTIVE) {
		tmp = G->VertexList[CurSolutionNum];
		G->VertexList[CurSolutionNum] = node;
		G->VertexList[G->Node_Index[node]] = tmp;
		G->Node_Index[tmp] = G->Node_Index[node];
		G->Node_Index[node] = CurSolutionNum++;
	}
	else {
		tmp = G->VertexList[CurSolutionNum - 1];
		G->VertexList[CurSolutionNum - 1] = node;
		G->VertexList[G->Node_Index[node]] = tmp;
		G->Node_Index[tmp] = G->Node_Index[node];
		G->Node_Index[node] = --CurSolutionNum;
	}
}

void Solution::MoveBucketNode(Bucket *point, int node, int direction) {
	int MG_index = G->MoveGain[node] + Solutions[node] * 2;
	if (direction == 1) {
		while (point->BucketId < MG_index) {
			point = point->next;
		}
		if (point->BucketId == MG_index)point->BucketList.push_back(node);
		else {
			point = new Bucket(MG_index, point);
			point->BucketList.push_back(node);
		}
	}
	else {
		while (point->BucketId > MG_index) {
			point = point->pre;
		}
		if (point->BucketId == MG_index)point->BucketList.push_back(node);
		else {
			point = new Bucket(point, MG_index);
			point->BucketList.push_back(node);
		}
	}
}

void Solution::AdjustBucket(int set, int direction, int movecount) {
	Bucket *point;
	int node, offset, end;
	std::vector<int>::iterator iter;
	if (set == ACTIVE) {
		offset = 2;
		if (direction == 1) {
			point = Bk_Head->next;
			end = G->Node_Num;
		}
		else {
			point = BK_Tail_A->pre;
			end = None;
		}
	}
	else {
		offset = -2;
		if (direction == 1) {
			point = BK_Tail_B->next;
			end = None;
		}
		else {
			point = Bk_Head->pre;
			end = -G->Node_Num;
		}
	}
	if (direction == 1) {
		while (point->BucketId != end) {
			iter = point->BucketList.begin();
			while (iter != point->BucketList.end()) {
				node = *iter;
				if (G->MoveGain[node] + offset != point->BucketId) {//检查所有点的move_gain
					iter = point->BucketList.erase(iter);
					MoveBucketNode(point->next, node, direction);
					movecount -= G->MoveGain[node] + offset - point->BucketId;
					if (movecount == 0) {
						if (point->BucketList.empty()) {
							point = point->next;
							delete point->pre;
						}
						return;
					}
				}
				else {
					++iter;
				}
			}
			if (point->BucketList.empty()) {
				point = point->next;
				delete point->pre;
			}
			else {
				point = point->next;
			}
		}
	}
	else {
		while (point->BucketId != end) {
			iter = point->BucketList.begin();
			while (iter != point->BucketList.end()) {
				node = *iter;
				if (G->MoveGain[node] + offset != point->BucketId) {
					iter = point->BucketList.erase(iter);
					MoveBucketNode(point->pre, node, direction);
					movecount -= point->BucketId - G->MoveGain[node] - offset;
					if (movecount == 0) {
						if (point->BucketList.empty()) {
							point = point->pre;
							delete point->next;
						}
						return;
					}
				}
				else {
					++iter;
				}
			}
			if (point->BucketList.empty()) {
				point = point->pre;
				delete point->next;
			}
			else {
				point = point->pre;
			}
		}
	}
}

void Solution::PerTurbMultiDA() {
	int i, node, addnode, dropnode;
	std::vector<int> MoveNodes;
	for (i = 0; i < G->Node_Num; ++i) {
		node = G->VertexList[i];
		G->Node_Tabu[node] = 0;
	}
	for (i = 0; i < PerTurbNum; ++i) {
		dropnode = SelectMoveNode(1);
		G->Node_Tabu[dropnode] = CurIter + 1;
		MoveNodes.push_back(dropnode);
		DropNode(dropnode);
		AdjustVertexList(dropnode);
	}
	for (i = 0; i < PerTurbNum; ++i) {
		addnode = SelectMoveNode(0);
		AddNode(addnode);
		AdjustVertexList(addnode);
	}
	for (i = 0; i < PerTurbNum; ++i) {
		node = MoveNodes[i];
		G->Node_Tabu[node] = 0;
	}
}

void Solution::PerTurbHCluster() {
	Graph* G = Graph::GetGraph();
	std::vector<int> ASetClusters;//A集合中固定的类
	std::vector<int> BSetClusters;//B集合中固定的类
	std::vector<int> ActiveClusters;//非固定的类

	//AMoveNodeNum表示左边可移动的点
	int i, node;
	int clusterNo;
	int NodesInAset, NodesInBset;
	int  tag, count = 0;
	int m;

	for (node = 1; node <= G->Node_Num; ++node) {
		G->Node_Tabu[node] = 0;
		G->Node_Cluster[node] = 0;
	}

//	checkAMoveNodeNum(AMoveNodeNum, BMoveNodeNum);
	for (i = 1; i < G->Clustering.size(); ++i) {
	//	if (G->ClusterState[i] != PASSIVE) cout << i << "error!!!!!!!!!!!!!!!!!!!" << endl;//???????
		G->ClusterNumInASet[i] = 0;
		for(int j = 0; j < G->Clustering[i].size(); j++) {
			node = G->Clustering[i][j];
			G->Node_Cluster[node] = 0;//每个点不指向任何一个组

			if (G->ClusterState[i] == PASSIVE) {
				if (Solutions[node] == ACTIVE && perturb_flag[node] == 0) {
					AMoveNodeNum++;
				}
				else if (Solutions[node] == PASSIVE && perturb_flag[node] == 0) {
					BMoveNodeNum++;
				}
			}	
		}
		G->ClusterState[i] = ACTIVE;//将固定的组全部Unlock
	}


	for (int i = 0; i < record_large.size(); i++) {
		perturb_flag[record_large[i]] = 0;
	}
	record_large.clear();
	//for (int i = 1; i <= G->Node_Num; i++) {
	//	if (perturb_flag[i] == 1)
	//	{
	//		cout << file << 11111111111 << endl;
	//		exit(0);
	//	}
	//}
	//checkAMoveNodeNum(AMoveNodeNum, BMoveNodeNum);
	//int generate_times = rand() % 3 + 1;?
	int generate_times = 3;//wang
	
	/*
	int seed = rand() % 100;
	if (seed < 25)generate_times = 2;
	else if (seed < 50)generate_times = 4;
	else generate_times = 3;
	*/

	A->set_0_clique_flag();

	lockvec1.clear();
	lockvec2.clear();
	lockvec3.clear();
	
	//然后随机生成三个团进行lock，然后讨论放在左边还是右边
	for (m = 1; m <= generate_times; m++) {

		int ccc = rand() % 10;
		int iniver;
		if (ccc < 5) {
			iniver = SelectMoveNode1(1);
			if(iniver = -1)
				iniver = SelectMoveNode1(0);
		}
		else {
			iniver = SelectMoveNode1(0);
			if (iniver = -1)
				iniver = SelectMoveNode1(1);
		}
		if (iniver == -1)
			return;

		if (rand() % 10 < 5 || iniver == -1) {
			int selectini = rand() % G->Node_Num + 1;
			do {
				selectini = rand() % G->Node_Num + 1;
			} while (clique_flag[selectini] != 0);
			iniver = selectini;
		}
		//cout << iniver << " ";
		freqq[iniver]++;

		if (G->flag_23 < 0){
			if (A->findcontrast_kRandom(m, iniver)) {
				//	if (m == 1)cout << lock_clique.size() << endl;
				break;//1.第一个顶点选不到 2.第二个顶点选不到,则不进行集合枚举
			}
		}
		else {
			if (A->findcontrast_kRandom1(m, iniver)) {
				//	if (m == 1)cout << lock_clique.size() << endl;
				break;//1.第一个顶点选不到 2.第二个顶点选不到,则不进行集合枚举
			}
		}

		count = 0;

		for (int j = 0; j < G->Clustering[m].size(); j++) {
			node = G->Clustering[m][j];
			G->Node_Cluster[node] = m;//标记每个点所在的组号
		}
	
		NodesInAset = 0;
		NodesInBset = 0;

		for (int n = 0; n < G->Clustering[m].size(); n++) {
			if (Solutions[G->Clustering[m][n]] == ACTIVE && perturb_flag[G->Clustering[m][n]] == 0) {//
				NodesInAset++;
			}
			if (Solutions[G->Clustering[m][n]] == PASSIVE && perturb_flag[G->Clustering[m][n]] == 0) {// 
				NodesInBset++;
			}
		}

		if (G->Clustering[m].size() <= BMoveNodeNum && G->Clustering[m].size() <= AMoveNodeNum) {
			tag = Rand->RandInt(0, 1);
		}
		else if (G->Clustering[m].size() <= BMoveNodeNum) {
			tag = 0;
		}
		else  if (G->Clustering[m].size() <= AMoveNodeNum) {
			tag = 1;
		}
		else {
			return;
		}
		//cout << AMoveNodeNum << " " << BMoveNodeNum << endl;
		AMoveNodeNum -= NodesInAset;
		BMoveNodeNum -= NodesInBset;
		G->ClusterNumInASet[m] = NodesInAset;
		
		G->ClusterState[m] = PASSIVE;//表示被固定了


		if (flag2 == 1 && m == 2) {//跟第一个块保持一致
			if (Solutions[G->Clustering[1][0]] == PASSIVE && G->Clustering[m].size() <= BMoveNodeNum) {
				tag = 0;
			}
			else if (Solutions[G->Clustering[1][0]] == ACTIVE && G->Clustering[m].size() <= AMoveNodeNum)
			{
				tag = 1;
			}
		}

		if (flag3 == 1 && m == 3) {//跟第一个块保持一致
			if (Solutions[G->Clustering[1][0]] == PASSIVE && G->Clustering[m].size() <= BMoveNodeNum) {
				tag = 0;
			}
			else if (Solutions[G->Clustering[1][0]] == ACTIVE && G->Clustering[m].size() <= AMoveNodeNum)
			{
				tag = 1;
			}
		}

		if (flag3 == 2 && m == 3) {//跟第二个块保持一致
			if (Solutions[G->Clustering[2][0]] == PASSIVE && G->Clustering[m].size() <= BMoveNodeNum) {
				tag = 0;
			}
			else if (Solutions[G->Clustering[2][0]] == ACTIVE && G->Clustering[m].size() <= AMoveNodeNum)
			{
				tag = 1;
			}
		}

		switch (tag) {
			//将类固定至B集合
		case 0:
			for (i = 0; i < G->Clustering[m].size(); ++i) {//Clusters的A到B ACTIVE是A，PASSIVE是B
				node = G->Clustering[m][i];
				if (Solutions[node] == ACTIVE) {
					++count;
					Bk_Head->DeleteNode(node, G->MoveGain[node] + 2, 1);
					DropNode(node);
					AdjustVertexList(node);
					--G->ClusterNumInASet[m];//
					if (perturb_flag[node] == 1) {
						--AMoveNodeNum;
						++BMoveNodeNum;
					}

				}
			}
			for (i = 0; i < count; ++i) {//B到A
				node = SelectMoveNode(0);
				if (Solutions[node] == ACTIVE) {
					cout << "error" << endl;
				}
				AddNode(node);
				AdjustVertexList(node);
				++G->ClusterNumInASet[G->Node_Cluster[node]];
				++AMoveNodeNum;
				--BMoveNodeNum;
			}
			break;
			//将类固定至A集合
		case 1:
			for (i = 0; i < G->Clustering[m].size(); ++i) {
				node = G->Clustering[m][i];
				if (Solutions[node] == PASSIVE) {
					++count;
					Bk_Head->DeleteNode(node, G->MoveGain[node] - 2, -1);
					AddNode(node);
					AdjustVertexList(node);
					++G->ClusterNumInASet[m];
					if (perturb_flag[node] == 1) {
						++AMoveNodeNum;
						--BMoveNodeNum;
					}
				}
			}
			for (i = 0; i < count; ++i) {
				node = SelectMoveNode(1);
				if (Solutions[node] == PASSIVE) {
					cout << "error" << endl;
				}
				DropNode(node);
				AdjustVertexList(node);
				--G->ClusterNumInASet[G->Node_Cluster[node]];
				--AMoveNodeNum;
				++BMoveNodeNum;
			}
			break;
		default:
			std::cout << "Error In Switch!" << std::endl;
			break;
		}
	}
}

void Solution::DropNode(int dropnode) {
	int neighbor, neighbor2;
	int movecountA = 0, movecountB = 0;
	std::vector<int>::iterator iter, iter2;
	CurSol += G->MoveGain[dropnode];//当前解的objective function CurSol越小越好,MoveGain小于0,越小越好，每次选择选择最小的MoveGain
	G->MoveGain[dropnode] = -G->MoveGain[dropnode];//顶点由B移到B'
	Solutions[dropnode] = PASSIVE;//若顶点在B,则为ACTIVE，若顶点在B',则为PASSIVE
	MoveBucketNode(Bk_Head->pre, dropnode, -1);
	for (iter = G->Adj_List[dropnode].begin(); iter != G->Adj_List[dropnode].end(); ++iter) {
		neighbor = *iter;
		--G->Gain_In[neighbor];//neighbor中在B中的顶点个数
		++G->Gain_Out[neighbor];//neighbor中在B'中的顶点个数
		if (Solutions[neighbor] == ACTIVE) {//由B移到B',若其邻居在B
			if (G->Gain_Out[dropnode] == 0) {//若dropnode的所有邻居全在B; 若其邻居没有B’的连接倒是可以理解
				--G->MoveGain[neighbor];
				++movecountA;
			}
			if (G->Gain_Out[neighbor] == 1) {//若neighbor的邻居中原先没有任何顶点在B,现在仅有dropnode在B
				--G->MoveGain[neighbor];//neighbor只有一个Gain_Out,则
				++movecountA;
				for (iter2 = G->Adj_List[neighbor].begin(); iter2 != G->Adj_List[neighbor].end(); ++iter2) {//neighbor的Gainout由0变为了1
					neighbor2 = *iter2;
					if (neighbor2 != dropnode) {
						--G->MoveGain[neighbor2];
						++movecountA;
					}
				}
			}
			else if (G->Gain_Out[neighbor] == 2) {
				for (iter2 = G->Adj_List[neighbor].begin(); iter2 != G->Adj_List[neighbor].end(); ++iter2) {
					neighbor2 = *iter2;
					if (Solutions[neighbor2] == PASSIVE && neighbor2 != dropnode) {
						++G->MoveGain[neighbor2];
						++movecountB;
						break;
					}
				}
			}
		}
		else {
			if (G->Gain_Out[dropnode] == 1) {
				++G->MoveGain[neighbor];
				++movecountB;
			}
			if (G->Gain_Out[neighbor] == 1) {
				++G->MoveGain[neighbor];
				++movecountB;
			}
		}
	}
	if (movecountA != 0) {
		AdjustBucket(1, -1, movecountA);
	}
	if (movecountB != 0) {
		AdjustBucket(-1, 1, movecountB);
	}
}

void Solution::AddNode(int addnode) {
	int neighbor, neighbor2;
	int movecountA = 0, movecountB = 0;
	std::vector<int>::iterator iter, iter2;
	CurSol += G->MoveGain[addnode];
	G->MoveGain[addnode] = -G->MoveGain[addnode];
	Solutions[addnode] = ACTIVE;
	MoveBucketNode(BK_Tail_A->pre, addnode, -1);
	for (iter = G->Adj_List[addnode].begin(); iter != G->Adj_List[addnode].end(); ++iter) {
		neighbor = *iter;
		++G->Gain_In[neighbor];
		--G->Gain_Out[neighbor];
		if (Solutions[neighbor] == ACTIVE) {
			if (G->Gain_Out[addnode] == 0) {
				++G->MoveGain[neighbor];
				++movecountA;
			}
			if (G->Gain_Out[neighbor] == 0) {
				++G->MoveGain[neighbor];
				++movecountA;
				for (iter2 = G->Adj_List[neighbor].begin(); iter2 != G->Adj_List[neighbor].end(); ++iter2) {
					neighbor2 = *iter2;
					if (neighbor2 != addnode) {
						++G->MoveGain[neighbor2];
						++movecountA;
					}
				}
			}
			else if (G->Gain_Out[neighbor] == 1) {
				for (iter2 = G->Adj_List[neighbor].begin(); iter2 != G->Adj_List[neighbor].end(); ++iter2) {
					neighbor2 = *iter2;
					if (Solutions[neighbor2] == PASSIVE && neighbor2 != addnode) {
						--G->MoveGain[neighbor2];
						++movecountB;
						break;
					}
				}
			}
		}
		else {
			if (G->Gain_Out[addnode] == 1) {
				--G->MoveGain[neighbor];
				++movecountB;
			}
			if (G->Gain_Out[neighbor] == 0) {
				--G->MoveGain[neighbor];
				++movecountB;
			}
		}
	}
	if (movecountA != 0) {
		AdjustBucket(1, 1, movecountA);
	}
	if (movecountB != 0) {
		AdjustBucket(-1, -1, movecountB);
	}
}

int Solution::SelectMoveNode2(int direct) {
	Bucket* Point;
	int finish;
	std::vector<int>::iterator iter, tmp_iter;
	std::vector<int> CandNode;
	std::vector<int>::iterator bestiter;
	double Maxprob = -1;
	int bestnode = -1;
	int movenode = -1, i;
	bool tag = false;

	if (direct == 1) {
		//DropNode
		Point = Bk_Head->next;
		finish = G->Node_Num;
	}
	else {
		//AddNode
		Point = BK_Tail_B->next;
		finish = 0;
	}

	int cnt1 = 0;
	int cnt2 = 0;
	bool st = false;
	int localbestnode = 0;
	int localNode_Tabu = 0x3f3f3f3f;
	//cout << SelectMethod << endl;
	int outercount = 0;
	int bbbbbbbbbb = 0;
	while (Point->BucketId < finish) {

		for (iter = Point->BucketList.begin(); iter != Point->BucketList.end(); ++iter) {
			movenode = *iter;
			if(Point->BucketList.size() != 0)
				bbbbbbbbbb++;
			if (clique_flag[movenode] == 1) {
				continue;
			}
			if (1) {
				tag = true;
				CandNode.push_back(movenode);
			}
		}
		if (tag) {
			outercount++;
		}
		else {
			Point = Point->next;
		}
		if (outercount == 2)
			break;
	}

	//cout << 1111111111 << endl;
	if (1) {
		if (CandNode.size() == 0) {
			return -1;
		}
		bestnode = CandNode[0];
		for (int i = 1; i < CandNode.size(); i++) {
			if (freqq[bestnode] > freqq[CandNode[i]]) {
				bestnode = CandNode[i];
			}
		}
	}
	return bbbbbbbbbb;
}

int Solution::SelectMoveNode1(int direct) {
	Bucket* Point;
	int finish;
	std::vector<int>::iterator iter, tmp_iter;
	std::vector<int> CandNode;
	std::vector<int>::iterator bestiter;
	double Maxprob = -1;
	int bestnode = -1;
	int movenode = -1, i;
	bool tag = false;

	if (direct == 1) {
		//DropNode
		Point = Bk_Head->next;
		finish = G->Node_Num;
	}
	else {
		//AddNode
		Point = BK_Tail_B->next;
		finish = 0;
	}

	int cnt1 = 0;
	int cnt2 = 0;
	bool st = false;
	int localbestnode = 0;
	int localNode_Tabu = 0x3f3f3f3f;
	//cout << SelectMethod << endl;
	int outercount = 0;

	int counttttt = 0;
	while (Point->BucketId < finish) {
		
		for (iter = Point->BucketList.begin(); iter != Point->BucketList.end(); ++iter) {
			movenode = *iter;
			counttttt++;
			if (clique_flag[movenode] == 1) {
				continue;
			}
			if (1) {
				tag = true;
				CandNode.push_back(movenode);
			}

		}

		if (tag) {
			outercount++;
		}
		else {
			Point = Point->next;
		}
		if (outercount == 2)
			break;
	}
	//if (CurEpoch == 10042) {
	//	SelectMoveNode3(0);
	//	SelectMoveNode3(1);
	//}
		

	if (1) {
		if (CandNode.size() == 0) {
			return -1;
		}
		bestnode = CandNode[0];
		for (int i = 1; i < CandNode.size(); i++) {
			if (freqq[bestnode] > freqq[CandNode[i]]) {
				bestnode = CandNode[i];
			}
		}
	}
	return bestnode;
}

int Solution::SelectMoveNode(int direct) {
	Bucket* Point;
	int finish;
	std::vector<int>::iterator iter, tmp_iter;
	std::vector<int> CandNode;
	std::vector<int>::iterator bestiter;
	double Maxprob = -1;
	int bestnode = -1;
	int movenode = -1, i;
	bool tag = false;

	if (direct == 1) {
		//DropNode
		Point = Bk_Head->next;
		finish = G->Node_Num;
	}
	else {
		//AddNode
		Point = BK_Tail_B->next;
		finish = 0;
	}

	int cnt1 = 0;
	int cnt2 = 0;
	bool st = false;
	int localbestnode = 0;
	int localNode_Tabu = 0x3f3f3f3f;
	
	while (Point->BucketId < finish) {
		for (iter = Point->BucketList.begin(); iter != Point->BucketList.end(); ++iter) {
			movenode = *iter;	
			if (direct == 1 && Solutions[movenode] != ACTIVE)
				cout << "SSSSSS1" << endl;
			if (direct == 0 && Solutions[movenode] != PASSIVE)
				cout << "SSSSSS2" << endl;
		/*
				if (!st) {
				if (G->ClusterState[G->Node_Cluster[movenode]] != PASSIVE) {
					if (G->Node_Tabu[movenode] < localNode_Tabu) {
						localNode_Tabu = G->Node_Tabu[movenode];
						localbestnode = movenode;
						tmp_iter = iter;
					}
				}
				else cnt1++;
			}
		*/
			
			if (G->Node_Tabu[movenode] > CurIter || (G->ClusterState[G->Node_Cluster[movenode]] == PASSIVE && perturb_flag[movenode] == 0)) {
				continue;
			}
			if (SelectMethod == "First") {
				Point->BucketList.erase(iter);
				return movenode;
			}
			else {
				tag = true;
				CandNode.push_back(movenode);
			}
		}
	//	if (localbestnode != 0)st = true;


		if (tag) {
			break;
		}
		else {
			Point = Point->next;
		}
	}


	if (SelectMethod == "Random") {
		if (CandNode.size() == 0) {
			std::cout << "Error IN SelectMoveNode!" << std::endl;
			exit(-1);
		}
		i = Rand->RandInt(0, CandNode.size() - 1);
		bestnode = CandNode[i];
		for (bestiter = Point->BucketList.begin(); bestiter != Point->BucketList.end(); ++bestiter) {
			if (bestnode == *bestiter)break;
		}
		Point->BucketList.erase(bestiter);
	}

	else {
		if (Point->BucketId == finish) {		
			/*
			if (st) {
				Point->BucketList.erase(tmp_iter);
				return localbestnode;
			}
			*/
			std::cout << SolutionNum << std::endl;
			cout << G->fname << " cnt: "  << " " << cnt1 << " " << cnt2 << endl;
			std::cout << "Error IN SelectMoveNode!" << std::endl;

			for (int i = 1; i < G->Node_Num; i++) {
				cout << Solutions[i] << G->Node_Cluster[i] << endl;
			}
				
			cout << "ssssss" << endl;

			exit(-1);
		}
	}
	return bestnode;
}

int Solution::CalPairScore(int v1, int v2) {//v1在B ACTIVE, v2在B'PASSIVE v1对应u，v2对应v gain_in是B，gain_out是B'
	Graph *G = Graph::GetGraph();
	std::vector<int>::iterator iter;
	if (Solutions[v1] != ACTIVE)
		cout << "NMSL" << endl;
	if (Solutions[v2] != PASSIVE)
		cout << "NMSL" << endl;
	int additional_minus = 0;
	if (G->Adj_Matrix[v1][v2] == 1 && G->Gain_Out[v2] == 0) {
		additional_minus++;
	}
	for (iter = G->Adj_List[v2].begin(); iter != G->Adj_List[v2].end(); ++iter) {
		int neighbor = *iter;
		if ((G->Adj_Matrix[v1][neighbor] == 1 || neighbor == v1) && G->Gain_Out[neighbor] == 1 && Solutions[neighbor] == ACTIVE) {
			additional_minus++;
		}
	}

	return  additional_minus;//??????????

}

int Solution::SelectMoveNodePair(int TabuStep) {
	Bucket* Point;
	Bucket* Point1;
	int finish = G->Node_Num;
	int finish1 = 0;;
	std::vector<int> CandNode;
	std::vector<int> CandNode1;
	Point = Bk_Head->next;
	Point1 = BK_Tail_B->next;

	std::vector<int>::iterator iter, tmp_iter;
	std::vector<int>::iterator bestiter;
	double Maxprob = -1;
	int bestnode = -1;
	int bestnode1 = -1;
	int movenode = -1, i;
	bool tag = false;
	vector<pair<int, int> > best_cand_ver;
	int PreCurSol = CurSol;
	while (Point->BucketId < finish) {
		for (iter = Point->BucketList.begin(); iter != Point->BucketList.end(); ++iter) {
			movenode = *iter;
			if (Solutions[movenode] != ACTIVE)
				cout << " if(Solutions[movenode] != ACTIVE) " << endl;

			if (G->Node_Tabu[movenode] > CurIter || (G->ClusterState[G->Node_Cluster[movenode]] == PASSIVE && perturb_flag[movenode] == 0)) {
				continue;
			}
			tag = true;
			CandNode.push_back(movenode);
		}

		if (tag) {
			break;
		}
		else {
			Point = Point->next;
		}
	}
	int kkkkkk = 0;
	tag = false;
	while (Point1->BucketId < finish1) {
		for (iter = Point1->BucketList.begin(); iter != Point1->BucketList.end(); ++iter) {
			movenode = *iter;
			kkkkkk++;
			if (Solutions[movenode] != PASSIVE)
				cout << " if(Solutions[movenode] != PASSIVE) " << endl;
			if ((G->ClusterState[G->Node_Cluster[movenode]] == PASSIVE && perturb_flag[movenode] == 0)) {//G->Node_Tabu[movenode] > CurIter || 
				continue;
			}
			tag = true;
			CandNode1.push_back(movenode);
		}

		if (tag) {
			break;
		}
		else {
			Point1 = Point1->next;
		}
	}
	if(CandNode1.size() == 0)
		cout << kkkkkk << endl;
	int best_result = G->Node_Num;

	if (CandNode.size() > 4) {
		CandNode.resize(4);
	}
	if (CandNode1.size() > 4) {
		CandNode1.resize(4);
	}

	if (CandNode.size() == 1 && CandNode1.size() == 1) {
		bestnode = CandNode[0];
		bestnode1 = CandNode1[0];
	}
	else {
		int can;
		int can1;
		best_result = G->Node_Num * 2;
		for (int i = 0; i < CandNode.size(); i++) {
			can = CandNode[i];
			for (int j = 0; j < CandNode1.size(); j++) {
				can1 = CandNode1[j];

				int curtwo_result = G->MoveGain[can] + G->MoveGain[can1] + CalPairScore(can, can1);
				if (curtwo_result < best_result) {
					best_result = curtwo_result;
					best_cand_ver.clear();
					best_cand_ver.push_back(make_pair(can, can1));
				}
				else if (curtwo_result == best_result) {
					best_cand_ver.push_back(make_pair(can, can1));
				}

			}
		}

		pair<int, int> vvvv = best_cand_ver[rand() % best_cand_ver.size()];
		bestnode = vvvv.first;
		bestnode1 = vvvv.second;
	}
	

	if (CandNode.size() == 0 || CandNode1.size() == 0) {
		std::cout << "Error IN SelectMoveNode! " << CandNode.size() << " " << CandNode1.size() << std::endl;
		exit(-1);
	}
	for (bestiter = Point->BucketList.begin(); bestiter != Point->BucketList.end(); ++bestiter) {
		if (bestnode == *bestiter)break;
	}
	Point->BucketList.erase(bestiter);

	for (bestiter = Point1->BucketList.begin(); bestiter != Point1->BucketList.end(); ++bestiter) {
		if (bestnode1 == *bestiter)break;
	}
	Point1->BucketList.erase(bestiter);


	//dropnode = SelectMoveNode(1);
	--G->ClusterNumInASet[G->Node_Cluster[bestnode]];
	DropNode(bestnode);
	AdjustVertexList(bestnode);
	if (TabuStep > 2) {
		G->Node_Tabu[bestnode] = CurIter + 1 + Rand->RandInt(1, TabuStep - 1);
	}
	else {
		G->Node_Tabu[bestnode] = CurIter + TabuStep;
	}


	++G->ClusterNumInASet[G->Node_Cluster[bestnode1]];
	AddNode(bestnode1);
	AdjustVertexList(bestnode1);
	if (TabuStep > 2) {
		G->Node_Tabu[bestnode1] = CurIter + 1 + Rand->RandInt(1, TabuStep - 1);
	}
	else {
		G->Node_Tabu[bestnode1] = CurIter + TabuStep;
	}
	if (best_result != G->Node_Num && PreCurSol + best_result != CurSol)
		cout << file << "score wrong" << endl;
	//CheckBucket();
	//cout << bestnode << " " << bestnode1 << endl;
	return bestnode;
}

void Solution::SelectMoveNode3(int direct) {
	Bucket* Point;
	int finish;
	std::vector<int>::iterator iter, tmp_iter;
	std::vector<int> CandNode;
	std::vector<int>::iterator bestiter;
	double Maxprob = -1;
	int bestnode = -1;
	int movenode = -1, i;
	bool tag = false;

	if (direct == 1) {
		//DropNode
		Point = Bk_Head->next;
		finish = G->Node_Num;
	}
	else {
		//AddNode
		Point = BK_Tail_B->next;
		finish = 0;
	}

	int cnt1 = 0;
	int cnt2 = 0;
	bool st = false;
	int localbestnode = 0;
	int localNode_Tabu = 0x3f3f3f3f;

	while (Point->BucketId < finish) {
		for (iter = Point->BucketList.begin(); iter != Point->BucketList.end(); ++iter) {
			movenode = *iter;

			if(Point->BucketList.size() > 0)
				CandNode.push_back(movenode);
			if (G->Node_Tabu[movenode] > CurIter || (G->ClusterState[G->Node_Cluster[movenode]] == PASSIVE && perturb_flag[movenode] == 0)) {
				continue;
			}
			if (CurEpoch == 10042) {
				cout << clique_flag[movenode] << " ";
			}
		}
		//	if (localbestnode != 0)st = true;



		Point = Point->next;

	}
	cout << endl;
	cout << CandNode.size() << " " << endl;

	//return bestnode;
}

bool Solution::checkAMoveNodeNum(int a, int b) {
	Graph* G = Graph::GetGraph();

	int ta = SolutionNum;
	int tb = G->Node_Num - SolutionNum;
	for (int node = 1; node <= G->Node_Num; node++) {
		if (G->ClusterState[G->Node_Cluster[node]] == PASSIVE && perturb_flag[node] == 0) {
			if (Solutions[node] == ACTIVE)ta--;
			else tb--;
		}
	}
	
	if (a != ta || b != tb) {
		cout << G->fname << " " << CurSolutionNum << " " << a << " " << ta << " " << b << " " << tb << endl;
		exit(-1);
		return false;
	}
	else return true;
}

bool Solution::CheckCurSol() {
	int testSol = 0;
	int node, i, neighbor;
	std::vector<int>::iterator iter;
	for (node = 1; node <= G->Node_Num; ++node) {
		iter = G->Adj_List[node].begin();
		G->Node_State[node] = 0;
		while (iter != G->Adj_List[node].end()) {
			neighbor = *iter++;
			if (Solutions[neighbor] == PASSIVE) {
				++G->Node_State[node];
			}
		}
		if (Solutions[node] == ACTIVE && G->Node_State[node] != 0)++testSol;
		if (G->Node_State[node] != G->Gain_Out[node]) {
			std::cout << "Error Gain_Out:" << node << "(" << G->Node_State[node] << "," << G->Gain_Out[node] << ")" << std::endl;
			return false;
		}
		G->Node_State[node] = 0;
	}
	if (testSol != CurSol) {
		return false;
	}
	else {
		return true;
	}
}

bool Solution::CheckGainOut() {
	int i, node, neighbor, count;
	std::vector<int>::iterator iter;
	for (node = 1; node < G->Node_Num; ++node) {
		count = 0;
		for (iter = G->Adj_List[node].begin(); iter != G->Adj_List[node].end(); ++iter) {
			neighbor = *iter;
			if (Solutions[neighbor] == PASSIVE) {
				++count;
			}
		}
		if (count != G->Gain_Out[node]) {
			std::cout << "Error In GainOut:" << node << std::endl;
			return false;
		}
	}
	return true;
}

bool Solution::CheckBucket() {
	Bucket *Point;
	std::vector<int>::iterator iter;
	int node;
	for (node = 1; node <= G->Node_Num; ++node) {
		G->Node_State[node] = 0;
	}
	Point = Bk_Head->next;
	int aaaaa = 0;
	int bbbbb = 0;
	while (Point != BK_Tail_A) {
		for (iter = Point->BucketList.begin(); iter != Point->BucketList.end(); ++iter) {
			node = *iter;
			aaaaa++;
			if (Solutions[node] == PASSIVE || Solutions[node] != ACTIVE) {
				std::cout << "Error In BK" << std::endl;
				return false;
			}
			if (G->MoveGain[node] + 2 != Point->BucketId) {
				std::cout << "Error In BKid" << std::endl;
				return false;
			}
		}
		Point = Point->next;
	}
	Point = Bk_Head->pre;
	while (Point != BK_Tail_B) {
		for (iter = Point->BucketList.begin(); iter != Point->BucketList.end(); ++iter) {
			node = *iter;
			bbbbb++;
			if (Solutions[node] == ACTIVE || Solutions[node] != PASSIVE) {
				std::cout << "Error In BK" << std::endl;
				return false;
			}
			if (G->MoveGain[node] - 2 != Point->BucketId) {
				std::cout << "Error In BKid" << std::endl;
				return false;
			}
		}
		Point = Point->pre;
	}
	//cout << " (" << aaaaa << " " << bbbbb << ") ";
	return true;
}

bool Solution::CheckBestSolution() {
	int i, node, neighbor;
	int star, finish;
	int solNum, CheckSol;
	int StoreBestSolution = 2;
	std::string Sol = "3,9,11,14,16,17,18,20,21,22,24,26,29,30,34,35,36,39,40,41,42,43,45,46,47,49,54,55,56,";
	std::string words;
	for (node = 1; node <= G->Node_Num; ++node) {
		if (Solutions[node] == ACTIVE) {
			std::cout << "Error in Solutions" << std::endl;
			return false;
		}
	}
	solNum = 0;
	star = 0;
	while (star < Sol.size()) {
		for (finish = star + 1; finish < Sol.size(); ++finish) {
			if (Sol[finish] == ',')break;
		}
		words = Sol.substr(star, finish - star);
		node = std::stoi(words, 0);
		Solutions[node] = ACTIVE;
		++solNum;
		star = finish + 1;
	}
	if (solNum != SolutionNum) {
		std::cout << "Error in SolutionNum" << std::endl;
		return false;
	}
	CheckSol = 0;
	for (node = 1; node <= G->Node_Num; ++node) {
		if (Solutions[node] == PASSIVE)continue;
		for (i = 0; i < G->Adj_List[node].size(); ++i) {
			neighbor = G->Adj_List[node][i];
			if (Solutions[neighbor] == PASSIVE) {
				++CheckSol;
				break;
			}
		}
	}
	std::cout << "BestSol:" << CheckSol << std::endl;
	if (CheckSol != StoreBestSolution) {
		return false;
	}
	return true;
}

void Solution::TabuSearch() {
	int unbetter = 0, i;
	int addnode = 0, dropnode = 0;
	int TabuStep;
	std::vector<int>::iterator iter;
	CurIter = 0;
	LocalBestSol = G->Node_Num;
	if (PerturbMethod == "HCluster") {
		TabuStep = AMoveNodeNum < BMoveNodeNum ? AMoveNodeNum : BMoveNodeNum;
		if (TabuStep == 0)return;
		total_cur_sol += CurSol;
		total_tabu += TabuStep;
		stepp++;
	}
	else {
		TabuStep = CurSol;
	}
	UnBetterBound = AMoveNodeNum + BMoveNodeNum;
	while (unbetter <= UnBetterBound) {
		++TotalIter;
		++CurIter;

		if (CurSol <= BestSol ) {
			SelectMoveNodePair(TabuStep);
		}
		else {
			dropnode = SelectMoveNode(1);
			--G->ClusterNumInASet[G->Node_Cluster[dropnode]];
			DropNode(dropnode);
			AdjustVertexList(dropnode);
			//        G->Node_Tabu[dropnode]=CurIter+TabuStep;
			if (TabuStep > 2) {
				G->Node_Tabu[dropnode] = CurIter + 1 + Rand->RandInt(1, TabuStep - 1);
			}
			else {
				G->Node_Tabu[dropnode] = CurIter + TabuStep;
			}
			addnode = SelectMoveNode(0);
			++G->ClusterNumInASet[G->Node_Cluster[addnode]];
			AddNode(addnode);
			AdjustVertexList(addnode);
			//        G->Node_Tabu[dropnode]=CurIter+TabuStep;
			if (TabuStep > 2) {
				G->Node_Tabu[dropnode] = CurIter + 1 + Rand->RandInt(1, TabuStep - 1);
			}
			else {
				G->Node_Tabu[dropnode] = CurIter + TabuStep;
			}


		}
		if (CurSol < LocalBestSol) {
			LocalBestSol = CurSol;
			unbetter = 0;
		}
		else {
			++unbetter;
		}

		if (CurSol < BestSol) {
			StoreSolutions();
		}
		
		//testbug();
		//if (checkAMoveNodeNum(AMoveNodeNum, BMoveNodeNum)) {
		//	//cout << "init is ok" << endl;
		//}
	}
}

void Solution::RLS() {
	//生成三个团
	//扰动
	freqq.resize(G->Node_Num + 1);
	while (true) {
		++CurEpoch;
		total_count++;

		if(!CheckBucket()){
		    std::cout<<"Error IN BUcket"<<std::endl;
		}
		
		if (Timer->OverCutoff()) {
			return;
		}

		TabuSearch();

	//	std::cout<<"Epoch "<<CurEpoch<<":"<<LocalBestSol<<std::endl;
		//if (LocalBestSol < BestSol) {
		//	StoreSolutions();
		//}
		//        if(!CheckBucket()){
		//            std::cout<<"Error IN BUcket"<<std::endl;
		//        }
		//		  fMonitor<<LocalBestSol<<std::endl;
		if (PerturbMethod == "Random") {
			PerTurbNum = Rand->RandInt(1, G->Node_Num / 2);
			PerTurbMultiDA();
		}
		//进行扰动
		else if (PerturbMethod == "HCluster") {
			//cout << " | ";
			PerTurbHCluster();
			//cout << " | ";
			//删除已有团，生成新的三个团
			//扰动
		}
		else {
			std::cout << "Error PurTurbMethod" << std::endl;
			exit(-1);
		}
		if (CurSolutionNum != SolutionNum) cout << "error" << endl;
	
	}
}
