#pragma once

#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <numeric>
#include <cmath>
#include <tuple>

using namespace std;
#define None 0
#define ACTIVE 1
#define PASSIVE -1

std::vector<int> sol_vec_a;
std::vector<int> sol_vec_b;
vector<int> perturb_flag;
int minn_LocalSol = 0x3f3f3f3f;
int local_cnt = 0;
int tag1 = 1;//如果当前tag为1，则表示选择阈值为a， 为2， 则表示选择阈值为 b
double aver_a = 0, aver_b = 0;

int last_con;
int current_con;
double best_value;

double kkkk = 0;
std::vector<int> sol_vec_c;
int lasta= 0, lastb = 0;
int threshold;
int total_count = 0;

int thre3;

class Graph {
public:
	/*图的信息*/
	string fileaa;
	int Node_Num;
	int Edge_Num;
	int MaxDegree;
	int MinDegree;

	int minn = 0x3f3f3f3f;
	int average = 0;
	int tmp_k = 0;
	double sum = 0;
	double flag_23;

	std::vector<double> Ft;
	std::vector<double> copy_Ft;
	std::vector<double> Pre_Pt;

	std::vector<int> potens;;//预处理初始点的候选解

	std::vector<int> VertexList;//VertexList[i - 1] = i;
	std::vector<int> Node_Index;//Node_Index[i] = i - 1;
	std::vector<int> Node_State;
	std::vector<int> Node_Degree; //点的度数
	std::vector<int> Node_Tabu;
	std::vector<int> Node_Cluster;
	std::vector<std::vector<int> > Adj_List;//邻接表
	std::vector<std::vector<int> > Adj_Matrix;//邻接表

	std::vector<int> Gain_In;
	std::vector<int> Gain_Out;
	std::vector<int> MoveGain;
	std::vector<std::vector<int>> Clustering;//存每个分组里面的点
	std::vector<int> ClusterNumInASet;//组在左边的点的数量
	std::vector<int> ClusterState;//组的状态
	std::vector<std::vector<int> > Matrix;//相似度矩阵
	std::vector<std::vector<int> > contrast_Matrix;//非共享的顶点

	//std::vector<std::vector<int> > Edge_Degree;//边的权重
	//std::vector<std::vector<int> > Edge_Weight;//边的权重

	std::string fname;
	static Graph* GetGraph() {
		static Graph GInstance;
		return &GInstance;
	}
	void GraphInit(char *FileIn);
	void AdjustGain(int node, int toSet);
	~Graph();

private:
	Graph() = default;
	Graph(const Graph&);
	Graph &operator=(const Graph&);
};

double calculateDegreeAssortativity(const vector<vector<int>>& adjList, int Node_Num) {
	vector<int> degrees(Node_Num + 1, 0); // Store degree of each node (1-based index)

	// Calculate the degree of each node
	for (size_t i = 1; i <= Node_Num; ++i) {
		degrees[i] = adjList[i].size();
	}

	// Store joint-degree pairs for all edges
	vector<pair<int, int>> degreePairs;
	for (size_t u = 1; u <= Node_Num; ++u) {
		for (int v : adjList[u]) {
			// Avoid double-counting edges
			if (u < v) {
				degreePairs.emplace_back(degrees[u], degrees[v]);
			}
		}
	}

	// Calculate sums for the numerator and denominator of assortativity formula
	double sum_jk = 0.0, sum_j = 0.0, sum_k = 0.0, sum_j2 = 0.0, sum_k2 = 0.0;
	for (const auto&[j, k] : degreePairs) {
		sum_jk += j * k;
		sum_j += j;
		sum_k += k;
		sum_j2 += j * j;
		sum_k2 += k * k;
	}

	size_t m = degreePairs.size(); // Number of edges

	// Assortativity coefficient formula
	double numerator = (sum_jk / m) - (sum_j / m) * (sum_k / m);
	double denominator = sqrt((sum_j2 / m) - pow(sum_j / m, 2)) * sqrt((sum_k2 / m) - pow(sum_k / m, 2));

	return (denominator == 0.0) ? 0.0 : numerator / denominator;
}

void Graph::GraphInit(char* FileIn) {

	std::ifstream fin(FileIn);
	char Buff[256];
	int l_node, r_node, offset = 0, node = 0, edge = 0;
	int i, j;
	if (fin.fail()) {
		std::cout << "Error In reading File!" << FileIn << std::endl;
		exit(-1);
	}
	fin.getline(Buff, 256);
	fin >> l_node >> r_node >> Edge_Num;
	Node_Num = l_node > r_node ? l_node : r_node;
	Node_Degree.resize(Node_Num + 2);
	while (!fin.eof()) {
		fin >> l_node >> r_node;
		if (l_node == r_node)
			continue;
		if (l_node > node)node = l_node;
		if (r_node > node)node = r_node;
		if (l_node == 0 || r_node == 0)offset = 1;
		Node_Degree[l_node]++;
		Node_Degree[r_node]++;
		edge++;
	}
	/*if((Node_Num!=(node+offset))||(edge!=Edge_Num)){
		std::cout<<"Error In reading Node Edge"<<std::endl;
		exit(-1);
	}*/
	Node_State.resize(Node_Num + 1);
	Node_Index.resize(Node_Num + 1, -1);
	Node_Tabu.resize(Node_Num + 1);
	VertexList.resize(Node_Num);
	Gain_In.resize(Node_Num + 1);
	Gain_Out.resize(Node_Num + 1);
	MoveGain.resize(Node_Num + 1);
	Adj_List.resize(Node_Num + 1);
	Node_Cluster.resize(Node_Num + 1);
	MinDegree = Node_Num;
	MaxDegree = 0;
	for (i = 1; i <= Node_Num; i++) {
		Adj_List[i].resize(Node_Degree[i - offset]);
		Node_Degree[i - offset] = 0;
	}
	fin.clear();
	fin.seekg(0);
	edge = 0;
	fin.getline(Buff, 256);
	fin.getline(Buff, 256);


	std::vector<bool> st;
	st.resize(Node_Num + 1);
	//Edge_Degree.resize(Node_Num + 1, std::vector<int>(Node_Num + 1));
	//Edge_Weight.resize(Node_Num + 1, std::vector<int>(Node_Num + 1));
	Adj_Matrix.resize(Node_Num + 1);
	for (int i = 0; i <= Node_Num; i++) {
		Adj_Matrix[i].resize(Node_Num + 1);
	}

	while (!fin.eof()) {
		fin >> l_node >> r_node;
		Adj_Matrix[l_node][r_node] = 1;
		Adj_Matrix[r_node][l_node] = 1;
		if (l_node != r_node) {
			l_node += offset;
			r_node += offset;

			//处理点集
			if (!st[l_node]) {
				st[l_node] = true;
			}
			if (!st[r_node]) {
				st[r_node] = true;
			}

			for (i = 0; i < Node_Degree[l_node]; i++) {
				if (Adj_List[l_node][i] == r_node)break;
			}
			if (i == Node_Degree[l_node]) {
				Adj_List[l_node][Node_Degree[l_node]++] = r_node;
				Adj_List[r_node][Node_Degree[r_node]++] = l_node;
				edge++;
				//Edge_Degree[l_node][r_node]++;
				//Edge_Degree[r_node][l_node]++;
			}
		}
	}
	Edge_Num = edge;
	for (i = 1; i <= Node_Num; i++) {
		Gain_Out[i] = Node_Degree[i];
		VertexList[i - 1] = i;
		Node_Index[i] = i - 1;
		if (Node_Degree[i] > MaxDegree)MaxDegree = Node_Degree[i];
		if (Node_Degree[i] < MinDegree)MinDegree = Node_Degree[i];
		Adj_List[i].resize(Node_Degree[i]);
	}

	//Edge_Weight = Edge_Degree;

	best_value = Node_Num;
	int n, m;
	Matrix.resize(Node_Num + 1);
	contrast_Matrix.resize(Node_Num + 1);
	for (i = 0; i <= Node_Num; ++i) {
		Matrix[i].resize(Node_Num + 1, 0);
		contrast_Matrix[i].resize(Node_Num + 1, 0);
	}
	//minn一邻居 minnn二邻居 二邻居平均 所有顶点平均


	int two_min = Node_Num;
	int two_max = 0;
	int two_avg = 0;
	long long add_sum = 0;
	int cntt = 0;
	int cnt = 0;

	long long add_sum2 = 0;
	for (m = 1; m <= Node_Num; ++m) {
		Matrix[m][m] = 1;
		Node_State[m] = ACTIVE;
		for (i = 0; i < Adj_List[m].size(); ++i) {
			node = Adj_List[m][i];
			Node_State[node] = ACTIVE;//标记邻居
		}
		for (n = m + 1; n <= Node_Num; ++n) {
			if (Node_State[n] == ACTIVE)++Matrix[m][n];
			for (i = 0; i < Adj_List[n].size(); ++i) {
				node = Adj_List[n][i];
				if (Node_State[node] == ACTIVE) {
					++Matrix[m][n];
				}
			}

			contrast_Matrix[m][n] = (Node_Degree[m] + 1) + (Node_Degree[n] + 1) - 2 * Matrix[m][n];
			contrast_Matrix[n][m] = contrast_Matrix[m][n];

			if (Node_State[n] == ACTIVE && contrast_Matrix[n][m] < best_value) {
				best_value = contrast_Matrix[n][m];
			}

			Matrix[n][m] = Matrix[m][n];

			cnt++;
			add_sum += contrast_Matrix[m][n];
			add_sum2 += Matrix[m][n];
			if (Node_State[n] == ACTIVE) {
				minn = std::min(minn, contrast_Matrix[m][n]);
			}
			if (Matrix[m][n] > 0) {
				if (contrast_Matrix[m][n] > two_max) {
					two_max = contrast_Matrix[m][n];
				}

				if (contrast_Matrix[m][n] < two_min) {
					two_min = contrast_Matrix[m][n];
				}

				two_avg += contrast_Matrix[m][n];
				cntt++;
			}
		}
		Node_State[m] = None;
		for (i = 0; i < Adj_List[m].size(); ++i) {
			node = Adj_List[m][i];
			Node_State[node] = None;
		}
	}
	for (i = 0; i <= Node_Num; ++i) {
		Node_State[i] = None;
	}
	two_avg = two_avg / cntt;

	//cout << FileIn << " " << 1.0 * 2 * Edge_Num / (Node_Num*(Node_Num - 1)) << endl;
	//double kkkkkkk = 0;
	//for (m = 1; m <= Node_Num; ++m) {
	//	for (n = m + 1; n <= Node_Num; ++n) {
	//		if (Matrix[m][n] > 0) {
	//			kkkkkkk += (contrast_Matrix[m][n] - two_avg)*(contrast_Matrix[m][n] - two_avg);
	//		}
	//	}
	//}

	//cout << FileIn << " " <<  two_min << " " << kkkkkkk / cntt << " " << two_max << " " << endl;
	//exit(0);
 	/*
	if (cnt != (Node_Num - 1) * Node_Num / 2) {
		std::cout << "error in cnt" << std::endl;
		exit(-1);
	}
	*/
	flag_23 = calculateDegreeAssortativity(Adj_List, Node_Num);
	if (Edge_Num / Node_Num > 70) {
		flag_23 = 1;
	}

	average = add_sum / cnt;
	threshold = std::max(average * 2 / 5, minn + 1);//另外一个阈值为 * 3 / 5 today, 影响potens， maxn , Ft[]

	int average2 = add_sum2 / cnt;

	thre3 = average2 * 3 / 5;
	//cout << FileIn << " " << threshold << endl;
	//exit(0);
	//cout << threshold << endl;
	Ft.resize(Node_Num + 1, 0);
	copy_Ft.resize(Node_Num + 1, 0);
	perturb_flag.resize(Node_Num + 1, 0);
	
	for (i = 1; i <= Node_Num; ++i) {
		for (j = 1; j <= Node_Num ; ++j) {
			if (i == j)
				continue;
			copy_Ft[i] += contrast_Matrix[i][j];//越小越好
		}
	}
	//cout << flag_23 << endl;

	double maxn = 0;
	for (int node = 1; node <= Node_Num; node++) {

		bool flag = false;
		for (int i = 0; i < Adj_List[node].size(); i++) {//判断当前点是否存在 与他的邻居的非共享顶点都<=x
			int a = Adj_List[node][i];

			if (contrast_Matrix[a][node] <= threshold && a != node) {
				flag = true;
				break;
			}
			if (flag)break;
			for (int j = 0; j < Adj_List[a].size(); j++) {
				int b = Adj_List[a][j];
				if (contrast_Matrix[b][node] <= threshold && b != node) {//判断当前点是否与他的两层邻居的非共享顶点都<=x
					flag = true;
					break;
				}
			}
			if (flag)break;
		}
		if (flag) {
			potens.push_back(node);
			maxn = std::max(maxn, copy_Ft[node]);
		}
	}


	for (i = 1; i <= Node_Num; ++i) {
		int best_valu11 = Node_Num;
		for (j = 0; j < Adj_List[i].size(); ++j) {
			int one_adj = Adj_List[i][j];
			if (contrast_Matrix[i][one_adj] < best_valu11) {
				best_valu11 = contrast_Matrix[i][one_adj];
			}
		}

		Ft[i] = 1.0 * (maxn + 1 - copy_Ft[i]);//越小越好

		//cout << maxn << " " << copy_Ft[i] << " " <<  1.0 * (maxn + 1 - copy_Ft[i]) << " " << 1 + best_valu11 - best_value << endl;
		//cout << endl;
	}

	Pre_Pt.resize(potens.size() + 1, 0);
	for (int i = 0; i < potens.size(); i++) {
		sum += Ft[potens[i]];
		Pre_Pt[i + 1] = Pre_Pt[i] + Ft[potens[i]];
		if (Ft[potens[i]] < 0)
			cout << FileIn << " (Ft[potens[i]] < 0) " << endl;
	}
	//cout << potens.size() << " " << Node_Num << endl;
	if (sum != Pre_Pt[potens.size()]) {
		std::cout << FileIn << " " << "error in sum" << std::endl;
		exit(-1);
	}
	fin.close();
}

void Graph::AdjustGain(int node, int toSet) {
	int neighbor;
	std::vector<int>::iterator iter;
	for (iter = Adj_List[node].begin(); iter != Adj_List[node].end(); ++iter) {
		neighbor = *iter;
		Gain_Out[neighbor] -= toSet;
		Gain_In[neighbor] += toSet;
	}
}

Graph::~Graph() {
	// std::cout << "Graph Destroy!" << std::endl;
}
