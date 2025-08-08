#pragma once

#include <vector>
#include <list>
#include <string>
#include "Graph.h"
#include <iostream>
#include <assert.h>
#include <set>
using namespace std;
const double eps = 1e-10;

vector<int> nonlockvec;

vector<int> lockvec1;
vector<int> lockvec2;
vector<int> lockvec3;

int flag2;
int flag3;

int iter_time;
int** Neighbor_Table1;
int* Neighbor_Table1_Degree;
int** Weight_Table1;
int* group_w_inside_degree;//组内权重度
int* grouping_Map;//分组时候每个超节点对应一个组
int* ori_compact_map;
int super_node_num;//超节点个数
int* ver_degree;//超节点内外权重度
int* group_size;//每个组的元素个数
int* mark;
int** Edge_Weight1;
vector<int> lock_clique;
vector<int> record_thre_cri;
vector<int> record_large;
int last_record_thre;
string file;

int* clique_flag;

struct three_ele {
	int v1;
	int v2;
	int weight;
};

class AGENS {
public:
	AGENS();
	AGENS(int n);
	~AGENS() = default;
	int ClusterNum;
	double MaxSim;
	double ThSim;
	int cnt;
	std::vector<std::vector<float> > SimilarMatrix;//相似度矩阵
	std::vector<std::vector<int> > Clusters;


	int CC = 0;
	int sum_m;
	bool mod_inc;
	

	void check_sum();
	void check_EDGE_degree();
	void check1();
	void check2();
	void check3(int group, int inside_num);
	void check4(int kin, int ps_ver, int group);
	void check5(int tmp_pos, int result_inside);
	void check6(int tmp_pos, int kj);
	void InitCluster();
	void HierarchicalCluster_Min();
	void HierarchicalCluster_Max();
	void HierarchicalCluster_Avg();
	void StoreCluster();

	void findMaximalCliqueRandom(int x);
	void findMaximal2PlexRandom(int x);

	bool findcontrast_kRandom(int x, int iniver);

	bool findcontrast_kRandom1(int x, int iniver);

	int get_initvertex();

	int getRandomElementFromSet(const set<int>& s);
	int getMinmaxElementFromSet(const set<int>& s, int best_ttt, int x);

	int getMax_minElementFromSet(const set<int>& s, int best_ttt, int x);
	int getdouble_ElementFromSet(const set<int>& s);


	void StoreCluster_louvain();
	double cal_modularity();

	void init();
	void louvain();//
	void Structural_Diagram();
	void Point_Grouping();

	void set_0_clique_flag();
	void set_1_clique_flag();
};


AGENS::AGENS() {
	ClusterNum = 5;
	InitCluster();
}

AGENS::AGENS(int n) :ClusterNum(n) {
	InitCluster();
}

void AGENS::set_0_clique_flag(){
	Graph *G = Graph::GetGraph();
	for (int i = 0; i <= G->Node_Num; i++) {
		clique_flag[i] = 0;
	}
}

void AGENS::set_1_clique_flag() {
	Graph *G = Graph::GetGraph();
	for (int i = 0; i <= G->Node_Num; i++) {
		clique_flag[i] = 1;
	}
}



//初始化相似度矩阵和类
void AGENS::InitCluster() {
	Graph *G = Graph::GetGraph();
	clique_flag = new int[G->Node_Num + 1];
	for (int i = 0; i <= G->Node_Num; i++) {
		clique_flag[i] = 0;
		G->Node_Cluster[i] = 0;
	}

	int ClusterCount = 4;//组的数量

/*
	int seed = rand() % 100;
	if (seed < 25)ClusterCount = 2;
	else if (seed < 50)ClusterCount = 4;
	else ClusterCount = 3;
*/

	int i, j, node;
	G->ClusterNumInASet.resize(ClusterCount + 1);
	G->ClusterState.resize(ClusterCount + 1);//组是否被lock掉，ACTIVE表示没有
	G->Clustering.resize(ClusterCount + 1);
	Clusters.resize(ClusterCount + 1);

	for (int i = 0; i <= ClusterCount; i++) {
		G->ClusterState[i] = ACTIVE;
	}

	for (i = 1; i <= 3; i++) {
		int selectini = rand() % G->Node_Num + 1;
		do {
			selectini = rand() % G->Node_Num + 1;
		} while (clique_flag[selectini] != 0);

		if(G->flag_23 < 0)
			findcontrast_kRandom(i, selectini);
		else {
			findcontrast_kRandom1(i, selectini);
		}
		
		//G->ClusterState[i] = PASSIVE;
		G->ClusterNumInASet[i] = 0;

		for (j = 0; j < G->Clustering[i].size(); j++) {
			node = G->Clustering[i][j];
			G->Node_Cluster[node] = i;//标记每个点所在的组号
		}

		/*
			if (CurSolutionNum + Clusters[i - 1].size() <= SolutionNum && B0 + Clusters[i - 1].size() <= SolutionNum) {//如果放在左右都不会超过 G->Node_Num / 2,则随机选择一边固定
			int x = rand() % 2;
			if (x == 0) {
				CurSolutionNum += G->Clustering[i - 1].size();
				for (j = 0; j < G->Clustering[i - 1].size(); j++) {
					node = G->Clustering[i - 1][j];
					Solutions[node] = ACTIVE;//标记表示在左边
				}
			}
			else if (x == 1) {
				B0 += G->Clustering[i - 1].size();
			}
		}
		*/
	}
	
	
	//cout << " || ";
	for (int i = 1; i <= G->Node_Num; i++) {
		if (G->Node_Cluster[i] == 0)G->Clustering[0].push_back(i);
	}
	//cout << G->Clustering[0].size() << endl;

	//int i, m, n, node;
	//int total_simi = 0;
	//int count_simi = 0;
	//int count_simi2 = 0;
	//SimilarMatrix.resize(G->Node_Num);
	//for (i = 0; i < G->Node_Num; ++i) {
	//	SimilarMatrix[i].resize(G->Node_Num);
	//}
	//MaxSim = G->Node_Num * 10;
	//Clusters.resize(G->Node_Num);
	//for (m = 1; m <= G->Node_Num; ++m) {
	//	SimilarMatrix[m - 1][m - 1] = 1;
	//	Clusters[m - 1].push_back(m);
	//	G->Node_State[m] = ACTIVE;
	//	for (i = 0; i < G->Adj_List[m].size(); ++i) {
	//		node = G->Adj_List[m][i];
	//		G->Node_State[node] = ACTIVE;
	//	}
	//	for (n = m + 1; n <= G->Node_Num; ++n) {
	//		if (G->Node_State[n] == ACTIVE)++SimilarMatrix[m - 1][n - 1];
	//		for (i = 0; i < G->Adj_List[n].size(); ++i) {
	//			node = G->Adj_List[n][i];
	//			if (G->Node_State[node] == ACTIVE) {
	//				++SimilarMatrix[m - 1][n - 1];
	//			}
	//		}
	//		int temp = SimilarMatrix[m - 1][n - 1];
	//		SimilarMatrix[m - 1][n - 1] = (G->Node_Degree[m] + 1) + (G->Node_Degree[n] + 1) - 2 * temp;
	//		total_simi += SimilarMatrix[m - 1][n - 1];
	//		count_simi++;
	//		if (SimilarMatrix[m - 1][n - 1] == 0) {
	//			count_simi2++;
	//		}
	//		SimilarMatrix[n - 1][m - 1] = SimilarMatrix[m - 1][n - 1];
	//		if (SimilarMatrix[n - 1][m - 1] < MaxSim) MaxSim = SimilarMatrix[n - 1][m - 1];
	//	}
	//	G->Node_State[m] = None;
	//	for (i = 0; i < G->Adj_List[m].size(); ++i) {
	//		node = G->Adj_List[m][i];
	//		G->Node_State[node] = None;
	//	}
	//}
	//double avg_simi = total_simi / count_simi;
	////	cout << avg_simi << " " <<  count_simi << " " <<  avg_simi/ count_simi << " " << count_simi2  << endl;
	// //   ThSim= max(3, avg_simi * 2);
	//ThSim = (avg_simi * 3) / 5;
	//if (ThSim < 3) {
	//	ThSim = 3;
	//}
	//else if (ThSim > 20) {
	//	ThSim = 20;
	//}
	/*
	if(ThSim<2){
		std::cout<<"ThSim <2"<<std::endl;
		ThSim=2;
	}
	if(ThSim>5){
		std::cout<<"ThSim >5"<<std::endl;
		ThSim=5;
	}
	std::cout<<"MaxSim:"<<MaxSim<<std::endl;
	std::cout<<"ThSim:"<<ThSim<<std::endl;
	std::cout<<"----------------"<<std::endl;
	*/
}

//层次聚类，求出Clusters，然后进行存储
void AGENS::HierarchicalCluster_Min() {
/*
		int i, j;
	// int MaxDis;
	int MinDis;
	int mi, mj;
	Graph *G = Graph::GetGraph();
	int ClusterCount = G->Node_Num;
	while (Clusters.size() > ClusterNum) {
		// MaxDis=0;
		MinDis = 0x3f3f3f3f;
		mi = 0;
		mj = 0;
		for (i = 0; i < Clusters.size(); ++i) {  //找到最小的SimilarMatrix
			for (j = i + 1; j < Clusters.size(); ++j) {
				if (SimilarMatrix[i][j] < MinDis) {
					MinDis = SimilarMatrix[i][j];
					mi = i;
					mj = j;
				}
			}
		}
		//        std::cout<<MaxDis<<"|"<<Clusters.size()<<std::endl;
		if (MinDis >= ThSim)break;
		for (i = 0; i < Clusters[mj].size(); ++i) {
			Clusters[mi].push_back(Clusters[mj][i]);
		}
		Clusters.erase(Clusters.begin() + mj);
		SimilarMatrix[mi][mi] += SimilarMatrix[mj][mj];
		//合并
		for (i = 0; i < SimilarMatrix.size(); ++i) {
			if (i == mi || i == mj)continue;
			//最大合并
			if (SimilarMatrix[i][mi] < SimilarMatrix[i][mj]) {
				SimilarMatrix[i][mi] = SimilarMatrix[i][mj];
				SimilarMatrix[mi][i] = SimilarMatrix[i][mi];
			}
		}
		SimilarMatrix.erase(SimilarMatrix.begin() + mj);
		for (i = 0; i < SimilarMatrix.size(); ++i) {
			SimilarMatrix[i].erase(SimilarMatrix[i].begin() + mj);
		}
	}
*/
	StoreCluster();
}

void AGENS::HierarchicalCluster_Max() {
	//int i, j;
	//int MaxDis;
	//int mi, mj;
	//Graph *G = Graph::GetGraph();
	//int ClusterCount = G->Node_Num;
	//while (Clusters.size() > ClusterNum) {
	//	MaxDis = 0;
	//	mi = 0;
	//	mj = 0;
	//	for (i = 0; i < Clusters.size(); ++i) {
	//		for (j = i + 1; j < Clusters.size(); ++j) {
	//			if (SimilarMatrix[i][j] > MaxDis) {
	//				MaxDis = SimilarMatrix[i][j];
	//				mi = i;
	//				mj = j;
	//			}
	//		}
	//	}
	//	// std::cout<<MaxDis<<"|"<<Clusters.size()<<std::endl;
	//	for (i = 0; i < Clusters[mj].size(); ++i) {
	//		Clusters[mi].push_back(Clusters[mj][i]);
	//	}
	//	Clusters.erase(Clusters.begin() + mj);
	//	SimilarMatrix[mi][mi] += SimilarMatrix[mj][mj];
	//	//合并
	//	for (i = 0; i < SimilarMatrix.size(); ++i) {
	//		if (i == mi || i == mj)continue;
	//		//最大合并
	//		if (SimilarMatrix[i][mi] < SimilarMatrix[i][mj]) {
	//			SimilarMatrix[i][mi] = SimilarMatrix[i][mj];
	//			SimilarMatrix[mi][i] = SimilarMatrix[i][mi];
	//		}
	//	}
	//	SimilarMatrix.erase(SimilarMatrix.begin() + mj);
	//	for (i = 0; i < SimilarMatrix.size(); ++i) {
	//		SimilarMatrix[i].erase(SimilarMatrix[i].begin() + mj);
	//	}
	//}
	//StoreCluster();
}

void AGENS::HierarchicalCluster_Avg() {
	//int i, j;
	//int MaxDis;
	//int mi, mj;
	//Graph *G = Graph::GetGraph();
	//int ClusterCount = G->Node_Num;
	//while (Clusters.size() > ClusterNum) {
	//	MaxDis = 0;
	//	mi = i;
	//	mj = j;
	//	for (i = 0; i < Clusters.size(); ++i) {
	//		for (j = i + 1; j < Clusters.size(); ++j) {
	//			if (SimilarMatrix[i][j] > MaxDis) {
	//				MaxDis = SimilarMatrix[i][j];
	//				mi = i;
	//				mj = j;
	//			}
	//			if (SimilarMatrix[j][i] > MaxDis) {
	//				MaxDis = SimilarMatrix[j][i];
	//				mi = i;
	//				mj = j;
	//			}
	//		}
	//	}
	//	for (i = 0; i < Clusters[mj].size(); ++i) {
	//		Clusters[mi].push_back(Clusters[mj][i]);
	//	}
	//	Clusters.erase(Clusters.begin() + mj);
	//	SimilarMatrix[mi][mi] += SimilarMatrix[mj][mj];
	//	//合并
	//	for (i = 0; i < SimilarMatrix.size(); ++i) {
	//		if (i == mi || i == mj)continue;
	//		//平均合并
	//		SimilarMatrix[mi][i] = (SimilarMatrix[i][i] * SimilarMatrix[i][mi] + SimilarMatrix[i][i] * SimilarMatrix[i][mj]) / SimilarMatrix[mi][mi];
	//		SimilarMatrix[i][mi] = SimilarMatrix[i][mi] + SimilarMatrix[i][mj];
	//	}
	//	SimilarMatrix.erase(SimilarMatrix.begin() + mj);
	//	for (i = 0; i < SimilarMatrix.size(); ++i) {
	//		SimilarMatrix[i].erase(SimilarMatrix[i].begin() + mj);
	//	}
	//}
	//StoreCluster();
}

int AGENS::getRandomElementFromSet(const set<int>& s) {
	if (s.empty()) {
		return -1;
	}

	// 生成一个随机索引，范围在 [0, s.size() - 1] 之间
	int randomIndex = rand() % s.size();

	// 使用迭代器定位到随机索引的位置
	auto it = s.begin();
	advance(it, randomIndex);  // 移动迭代器到随机索引处

	return *it;
}

int AGENS::getdouble_ElementFromSet(const set<int>& s) {
	//Graph* G = Graph::GetGraph();
	//int max = 0;
	//int max_node = 0;
	//for (int a : s) {//遍历候选解
	//	if (a >= 1 && a <= G->Node_Num) {

	//	}
	//	else {
	//		cout << "error int a s" << endl;
	//		exit(-1);
	//	}

	//	int minn = 0x3f3f3f3f;
	//	for (int i = 0; i < lock_clique.size(); i++) {
	//		int b = lock_clique[i];//遍历当前解

	//		if (a >= 1 && a <= G->Node_Num && b >= 1 && b <= G->Node_Num) {
	//			if (G->Matrix[a][b] < minn) {//算出共享邻居最小的值??
	//				minn = G->Matrix[a][b];
	//			}
	//		}
	//		else {
	//			cout << "error int a and b" << endl;
	//			exit(-1);
	//		}
	//	}

	//	if (minn > max) {//找出最小值最大的点
	//		max = minn;
	//		max_node = a;
	//	}
	//}
	//return max_node;
	Graph* G = Graph::GetGraph();

	int maxn = 0;
	vector<int> candidate;
	for (int a : s) {//遍历候选解
		if (a >= 1 && a <= G->Node_Num) {

		}
		else {
			cout << "error int a s" << endl;
			exit(-1);
		}

		int minn = 0x3f3f3f3f;
		for (int i = 0; i < lock_clique.size(); i++) {
			int b = lock_clique[i];//遍历当前解

			if (a >= 1 && a <= G->Node_Num && b >= 1 && b <= G->Node_Num) {
				if (G->Matrix[a][b] < minn) {//算出共享邻居最小的值
					minn = G->Matrix[a][b];
				}
			}
			else {
				cout << "error int a and b" << endl;
				exit(-1);
			}
		}
		if (minn > maxn) {//找出最大值
			candidate.clear();
			candidate.push_back(a);
			maxn = minn;
		}
		else if (minn == maxn) {//找出最小值最大的点
			candidate.push_back(a);
		}
	}

	if (candidate.size() == 1)return candidate[0];

	vector<int> res;
	int minn = 0x3f3f3f3f;
	for (int i = 0; i < candidate.size(); i++) {
		int a = candidate[i];

		if (a >= 1 && a <= G->Node_Num) {

		}
		else {
			cout << "error int a s" << endl;
			exit(-1);
		}

		int maxn = 0;
		for (int j = 0; j < lock_clique.size(); j++) {
			int b = lock_clique[j];//遍历当前解

			if (a >= 1 && a <= G->Node_Num && b >= 1 && b <= G->Node_Num) {
				if (G->contrast_Matrix[a][b] > maxn) {//算出非共享邻居最大的值
					maxn = G->contrast_Matrix[a][b];
				}
			}
			else {
				cout << "error int a and b" << endl;
				exit(-1);
			}
		}


		if (maxn < minn) {//找出最小值
			res.clear();
			res.push_back(a);
			minn = maxn;
		}
		else if (maxn == minn) {
			res.push_back(a);
		}
	}

	if (res.size() == 1)return res[0];
	else {
		int seed = rand() % res.size();
		return res[seed];
	}
}

int AGENS::get_initvertex() {
	Graph* G = Graph::GetGraph();

	long long add_sum = 0;
	
	for (int i = 0; i < G->potens.size(); i++) {
		int node = G->potens[i];
		if (clique_flag[node] == 0) {//如果没有被选择过
			add_sum += G->Ft[node];
		}
	}

	if (add_sum == 0)return 0;//如果没有可满足的点，则返回

	bool flag = false;
	while (!flag) {
		long long seed = rand() % add_sum + 1;
		long long l = 1, r = 0;

		for (int i = 0; i < G->potens.size(); i++) {
			int node = G->potens[i];
			if (clique_flag[node] != 0)continue;//当前点没有被选择过
			
			r += G->Ft[node];

			if (seed >= l && seed <= r) {
				return node;
				flag = true;
				break;
			}
			l += G->Ft[node];
		}
	}

	return 0;
}

int AGENS::getMinmaxElementFromSet(const set<int>& s, int best_ttt, int x) {
	Graph* G = Graph::GetGraph();
	int max = 0;
	int max_node = 0;
	vector<int> res;
	for (int a : s) {//遍历候选解
		if (a >= 1 && a <= G->Node_Num) {

		}
		else {
			cout << "error int a s" << endl;
			exit(-1);
		}

		int minn = 0x3f3f3f3f;
		for (int i = 0; i < lock_clique.size(); i++) {
			int b = lock_clique[i];//遍历当前解

			if (a >= 1 && a <= G->Node_Num && b >= 1 && b <= G->Node_Num) {
				if (G->Matrix[a][b] < minn) {//算出共享邻居最小的值??
					minn = G->Matrix[a][b];
				}
			}
			else {
				cout << "error int a and b" << endl;
				exit(-1);
			}
		}



		if (minn > max) {//找出最大值最小的点
			max = minn;
			res.clear();
			res.push_back(a);
		}
		else if (max == minn) {
			res.push_back(a);
		}
	}

	if (res.size() == 0)return -1;
	if (res.size() == 1) {
		if (best_ttt != max) {
			nonlockvec.push_back(res[0]);
		}

		if (x == 1 && best_ttt == max) {
			lockvec1.push_back(res[0]);
		}
		if (x == 2 && best_ttt == max) {
			lockvec2.push_back(res[0]);
		}
		if (x == 3 && best_ttt == max) {
			lockvec3.push_back(res[0]);
		}

		return res[0];
	}
	else {
		int seed = rand() % res.size();

		if (best_ttt != max) {
			nonlockvec.push_back(res[seed]);
		}

		if (x == 1 && best_ttt == max) {
			lockvec1.push_back(res[seed]);
		}
		if (x == 2 && best_ttt == max) {
			lockvec2.push_back(res[seed]);
		}
		if (x == 3 && best_ttt == max) {
			lockvec3.push_back(res[seed]);
		}

		return res[seed];
	}
	return max_node;
}





int AGENS::getMax_minElementFromSet(const set<int>& s, int best_ttt, int x) {
	Graph* G = Graph::GetGraph();
	int minn = 0x3f3f3f3f;
	int min_node = 0;
	current_con = 10;
	vector<int> res;
	int max = G->Node_Num;
	for (int a : s) {//遍历候选解
		if (a >= 1 && a <= G->Node_Num) {

		}
		else {
			cout << "error int a s" << endl;
			cout << a << endl;
			exit(-1);
		}

		if (clique_flag[a] != 0)continue;//当前候选解不能被选择过
		int maxn = 0;
		for (int i = 0; i < lock_clique.size(); i++) {
			int b = lock_clique[i];//遍历当前解

			if (a >= 1 && a <= G->Node_Num && b >= 1 && b <= G->Node_Num) {
				if (G->contrast_Matrix[a][b] > maxn) {//算出共享邻居最大的值
					maxn = G->contrast_Matrix[a][b];
				}
			}
			else {
				cout << "error int a and b" << endl;
				exit(-1);
			}
		}

		if (maxn < minn) {//找出最大值最小的点
			minn = maxn;
			res.clear();
			res.push_back(a);
		}
		else if (maxn == minn) {
			res.push_back(a);
		}
	}

	if (minn != last_record_thre) {
		record_thre_cri.push_back(lock_clique.size() + 1);
		last_record_thre = minn;
	}
	
	max = minn;

	if (res.size() == 0)return -1;
	if (res.size() == 1) {
		if (best_ttt != max) {
			nonlockvec.push_back(res[0]);
		}

		if (x == 1 && best_ttt == max) {
			lockvec1.push_back(res[0]);
		}
		if (x == 2 && best_ttt == max) {
			lockvec2.push_back(res[0]);
		}
		if (x == 3 && best_ttt == max) {
			lockvec3.push_back(res[0]);
		}

		return res[0];
	}
	else {
		int seed = rand() % res.size();

		if (best_ttt != max) {
			nonlockvec.push_back(res[seed]);
		}

		if (x == 1 && best_ttt == max) {
			lockvec1.push_back(res[seed]);
		}
		if (x == 2 && best_ttt == max) {
			lockvec2.push_back(res[seed]);
		}
		if (x == 3 && best_ttt == max) {
			lockvec3.push_back(res[seed]);
		}

		return res[seed];
	}
}

bool AGENS::findcontrast_kRandom(int x, int iniver) {
	Graph* G = Graph::GetGraph();

	lock_clique.clear();//lock_clique 
	G->Clustering[x].clear();
	record_thre_cri.clear();
	

	int initialVertex = iniver;

	lock_clique.push_back(initialVertex);
	G->Clustering[x].push_back(initialVertex);
	clique_flag[initialVertex] = 1;

	
	nonlockvec.clear();


	int best_ttt = G->Node_Num;
	set<int> potentialVertices;//满足条件的候选解的集合
	for (int i = 0; i < G->Adj_List[initialVertex].size(); i++) {//更新候选解，遍历初始点的邻居
		int a = G->Adj_List[initialVertex][i];
		if (clique_flag[a] == 0 ) {//没有被选过
			potentialVertices.insert(a);
			if (G->contrast_Matrix[initialVertex][a] < best_ttt)
				best_ttt = G->contrast_Matrix[initialVertex][a];

		}
		for (int j = 0; j < G->Adj_List[a].size(); j++) {//更新候选解，遍历初始点的邻居
			int bb = G->Adj_List[a][j];
			if (clique_flag[bb] == 0) {//没有被选过
				potentialVertices.insert(bb);
				if (G->contrast_Matrix[initialVertex][bb] < best_ttt)
					best_ttt = G->contrast_Matrix[initialVertex][bb];
			}
		}
	}
	last_record_thre = best_ttt;
	record_thre_cri.push_back(lock_clique.size());
	last_con = best_ttt;

	int thre1 = max(threshold, best_ttt + 1);

	if (potentialVertices.empty()) {
		return true;//如果不存在这样的点，或者非共享顶点的数量超过threshold，这提前退出
	}

	while (!potentialVertices.empty()) {
		int locak_vertex = getMax_minElementFromSet(potentialVertices, best_ttt, x);//选择最大最小的点，如果有多个，则随机进行选择
		if (locak_vertex == -1) {
			cout << G->fname << "error in vetex" << endl;
			break;
		}

		potentialVertices.erase(locak_vertex);//候选加入的顶点移除vertex

		lock_clique.push_back(locak_vertex);
		G->Clustering[x].push_back(locak_vertex);//wang
		clique_flag[locak_vertex] = 1;


		auto it = potentialVertices.begin();//加入一个点之后，先判断候选解中不满足average的点去掉

		while (it != potentialVertices.end()) {
			int v = *it;

			if (G->contrast_Matrix[v][locak_vertex] > thre1 || G->Matrix[v][locak_vertex] == 0) {// 删除元素，并更新迭代器
				it = potentialVertices.erase(it);  // erase 会返回指向下一个元素的迭代器
			}
			else {
				++it;  // 如果没有删除，则正常移动迭代器
			}
		}

		//for (int m = 0; m < G->Adj_List[locak_vertex].size(); m++) {
		//	int a = G->Adj_List[locak_vertex][m];
		//	if (clique_flag[a] != 0)continue; //将当前点的邻居存入其中，同时保证没有被选过，以及满足min_k
		//	bool flag = true;
		//	for (int j = 0; j < lock_clique.size(); j++) {//
		//		int b = lock_clique[j];
		//		if (G->contrast_Matrix[a][b] > thre1 || G->Matrix[a][b] == 0) {
		//			flag = false;
		//			break;
		//		}
		//	}
		//	if (flag)potentialVertices.insert(a);
		//}	
	}

	for (int i = 0; i < nonlockvec.size(); i++) {
		record_large.push_back(nonlockvec[i]);
		perturb_flag[nonlockvec[i]] = 1;
		//clique_flag[lock_clique[i]] = 0;
	}

	for (int i = 0; i < lock_clique.size(); i++) {
		int a = lock_clique[i];

		if (G->Clustering[x][i] != lock_clique[i]) {
			cout << "error in lock_clique" << endl;
			exit(-1);
		}
	}


	for (int i = 0; i < lock_clique.size(); i++) {
		int a = lock_clique[i];

		for (int j = i; j < lock_clique.size(); j++) {
			bool flag = false;
			if (i == j)continue;
			int b = lock_clique[j];

			if (a == b) {
				cout << "error in common" << endl;
				exit(-1);
			}
		}
	}

	if (x == 2) {
		flag2 = -1;
		for (int i = 0; i < lockvec1.size(); i++) {
			for (int j = 0; j < lockvec2.size(); j++) {
				if (G->contrast_Matrix[lockvec1[i]][lockvec2[j]] == best_ttt) {
					flag2 = 1;
				}
			}
		}
	}

	if (x == 3) {
		flag3 = -1;
		int a = -1;
		int b = -1;
		for (int i = 0; i < lockvec1.size(); i++) {
			for (int j = 0; j < lockvec3.size(); j++) {
				if (G->contrast_Matrix[lockvec1[i]][lockvec3[j]] == best_ttt) {
					a = 1;
				}
			}
		}
		for (int i = 0; i < lockvec2.size(); i++) {
			for (int j = 0; j < lockvec3.size(); j++) {
				if (G->contrast_Matrix[lockvec2[i]][lockvec3[j]] == best_ttt) {
					b = 1;
				}
			}
		}
		if (a == 1 && b == -1) {
			flag3 = 1;
		}
		if (a == -1 && b == 1) {
			flag3 = 2;
		}

	}
	
	//for (int i = 0; i < lock_clique.size(); i++) {
	//	cout << lock_clique[i] << " ";
	//}
	//cout << endl;
	if (lock_clique.size() == 1) {
		clique_flag[lock_clique[0]] = 0;
		return true;
	}
	else
		return false;
}

bool AGENS::findcontrast_kRandom1(int x, int iniver) {
	Graph* G = Graph::GetGraph();

	lock_clique.clear();//lock_clique 
	G->Clustering[x].clear();
	record_thre_cri.clear();

	//根据轮盘法选择初始点
	//int initialVertex;
	////根据轮盘法选择初始点
	//bool flag = false;
	//while (!flag) {
	//	double l = 0, r = 0;

	//	double value = static_cast<double>(rand()) / RAND_MAX *  G->sum ;
	//	//cout << value << " ";
	//	for (int i = 1; i <= G->potens.size(); ++i) {
	//		l = G->Pre_Pt[i - 1];
	//		r = G->Pre_Pt[i];

	//		if (value <= r && value >= l && clique_flag[G->potens[i - 1]] == 0) {
	//			initialVertex = G->potens[i - 1];
	//			flag = true;
	//			break;
	//		}
	//	}
	//}

	//if (initialVertex == 0) {
	//	
	//	return true;//表示找不到初始点，直接返回
	//}
	int initialVertex = iniver;

	lock_clique.push_back(initialVertex);
	G->Clustering[x].push_back(initialVertex);
	clique_flag[initialVertex] = 1;


	nonlockvec.clear();


	int best_ttt = 0;
	int top_unadj = 0;
	set<int> potentialVertices;//满足条件的候选解的集合

	for (int i = 0; i < G->Adj_List[initialVertex].size(); i++) {//更新候选解，遍历初始点的邻居
		int a = G->Adj_List[initialVertex][i];
		if (clique_flag[a] == 0) {//没有被选过
			if (G->Matrix[initialVertex][a] > best_ttt)
				best_ttt = G->Matrix[initialVertex][a];
			if (G->contrast_Matrix[initialVertex][a] > top_unadj)
				top_unadj = G->contrast_Matrix[initialVertex][a];

		}
		for (int j = 0; j < G->Adj_List[a].size(); j++) {//更新候选解，遍历初始点的邻居
			int bb = G->Adj_List[a][j];
			if (clique_flag[bb] == 0) {//没有被选过
				if (G->Matrix[initialVertex][bb] > best_ttt)
					best_ttt = G->Matrix[initialVertex][bb];
				if (G->contrast_Matrix[initialVertex][bb] > top_unadj)
					top_unadj = G->contrast_Matrix[initialVertex][bb];
			}
		}
	}

	for (int i = 0; i < G->Adj_List[initialVertex].size(); i++) {//更新候选解，遍历初始点的邻居
		int a = G->Adj_List[initialVertex][i];
		if (clique_flag[a] == 0) {//没有被选过
			potentialVertices.insert(a);
			if (G->Matrix[initialVertex][a] > best_ttt && G->contrast_Matrix[initialVertex][a] < top_unadj)
				best_ttt = G->Matrix[initialVertex][a];

		}
		for (int j = 0; j < G->Adj_List[a].size(); j++) {//更新候选解，遍历初始点的邻居
			int bb = G->Adj_List[a][j];
			if (clique_flag[bb] == 0) {//没有被选过
				potentialVertices.insert(bb);
				if (G->Matrix[initialVertex][bb] > best_ttt && G->contrast_Matrix[initialVertex][bb] < top_unadj)
					best_ttt = G->Matrix[initialVertex][bb];
			}
		}
	}

	last_record_thre = best_ttt;
	record_thre_cri.push_back(lock_clique.size());
	last_con = best_ttt;




	int thre1;//thre越小越容易，越小越宽松，如果thre3太大了，则用best_ttt约束
	if (best_ttt < thre3) {
		thre1 = best_ttt;
	}
	else {
		thre1 = thre3;
	}

	if (thre1 < 2)
		thre1 = 2;



	if (potentialVertices.empty()) {
		return true;//如果不存在这样的点，或者非共享顶点的数量超过threshold，这提前退出
	}

	while (!potentialVertices.empty()) {
		int locak_vertex = getMinmaxElementFromSet(potentialVertices, best_ttt, x);//选择最大最小的点，如果有多个，则随机进行选择
		if (locak_vertex == -1) {
			cout << G->fname << "error in vetex" << endl;
			break;
		}

		potentialVertices.erase(locak_vertex);//候选加入的顶点移除vertex

		lock_clique.push_back(locak_vertex);
		G->Clustering[x].push_back(locak_vertex);//wang
		clique_flag[locak_vertex] = 1;


		auto it = potentialVertices.begin();//加入一个点之后，先判断候选解中不满足average的点去掉

		while (it != potentialVertices.end()) {
			int v = *it;

			if (G->Matrix[v][locak_vertex] < thre1 || G->contrast_Matrix[v][locak_vertex] >= top_unadj) {// 删除元素，并更新迭代器
				it = potentialVertices.erase(it);  // erase 会返回指向下一个元素的迭代器
			}
			else {
				++it;  // 如果没有删除，则正常移动迭代器
			}
		}

		/*if (current_con == -1 && lock_clique.size() > 10) {
			cout << lock_clique.size() << " ";
			return false;
		}*/


		//for (int m = 0; m < G->Adj_List[locak_vertex].size(); m++) {
		//	int a = G->Adj_List[locak_vertex][m];
		//	if (clique_flag[a] != 0)continue; //将当前点的邻居存入其中，同时保证没有被选过，以及满足min_k
		//	bool flag = true;
		//	for (int j = 0; j < lock_clique.size(); j++) {//
		//		int b = lock_clique[j];
		//		if (G->Matrix[a][b] < thre1 ||  G->contrast_Matrix[a][b] >= top_unadj) {
		//			flag = false;
		//			break;
		//		}
		//	}
		//	if (flag)potentialVertices.insert(a);
		//}
	}

	for (int i = 0; i < nonlockvec.size(); i++) {
		record_large.push_back(nonlockvec[i]);
		perturb_flag[nonlockvec[i]] = 1;
		//clique_flag[lock_clique[i]] = 0;
	}

	for (int i = 0; i < lock_clique.size(); i++) {
		int a = lock_clique[i];

		if (G->Clustering[x][i] != lock_clique[i]) {
			cout << "error in lock_clique" << endl;
			exit(-1);
		}
	}

	for (int i = 0; i < lock_clique.size(); i++) {
		int a = lock_clique[i];

		for (int j = i; j < lock_clique.size(); j++) {
			bool flag = false;
			if (i == j)continue;
			int b = lock_clique[j];

			if (a == b) {
				cout << "error in common" << endl;
				exit(-1);
			}
		}
	}

	if (x == 2) {
		flag2 = -1;
		for (int i = 0; i < lockvec1.size(); i++) {
			for (int j = 0; j < lockvec2.size(); j++) {
				if (G->Matrix[lockvec1[i]][lockvec2[j]] == best_ttt) {
					flag2 = 1;
				}
			}
		}
	}
	if (x == 3) {
		flag3 = -1;
		int a = -1;
		int b = -1;
		for (int i = 0; i < lockvec1.size(); i++) {
			for (int j = 0; j < lockvec3.size(); j++) {
				if (G->Matrix[lockvec1[i]][lockvec3[j]] == best_ttt) {
					a = 1;
				}
			}
		}
		for (int i = 0; i < lockvec2.size(); i++) {
			for (int j = 0; j < lockvec3.size(); j++) {
				if (G->Matrix[lockvec2[i]][lockvec3[j]] == best_ttt) {
					b = 1;
				}
			}
		}
		if (a == 1 && b == -1) {
			flag3 = 1;
		}
		if (a == -1 && b == 1) {
			flag3 = 2;
		}

	}
	/*
		for (int i = 0; i < lock_clique.size(); i++) {
		int a = lock_clique[i];

		bool flag = false;
		for (int j = 0; j < lock_clique.size(); j++) {
			if (i == j)continue;
			int b = lock_clique[j];

			for (int k = 0; k < G->Adj_List[b].size(); k++) {
				int c = G->Adj_List[b][k];

				if (a == c) {
					flag = true;
					break;
				}
			}
		}
		if (!flag) {
			cout << "do not connection " << lock_clique.size() << endl;
			exit(-1);
		}
	}
	*/
	if (lock_clique.size() == 1) {
		clique_flag[lock_clique[0]] = 0;
		return true;
	}
	else
		return false;
}

void AGENS::StoreCluster() {
	/*
	int ClusterCount = 0;

	for (i = 0; i < Clusters.size(); ++i) {
		if (Clusters[i].size() > 1) {
			++ClusterCount;
		}
	}
	G->ClusterNumInASet.resize(ClusterCount + 1);
	G->ClusterState.resize(ClusterCount + 1);
	G->Clustering.resize(ClusterCount + 1);
	ClusterCount = 0;
	G->ClusterState[0] = ACTIVE;
	for (i = 0; i < Clusters.size(); ++i) {
		if (Clusters[i].size() > 1) {
			++ClusterCount;
			G->ClusterState[ClusterCount] = ACTIVE;
			for (j = 0; j < Clusters[i].size(); ++j) {
				node = Clusters[i][j];
				G->Clustering[ClusterCount].push_back(node);
				G->Node_Cluster[node] = ClusterCount;
			}
			cout << "wrong" << endl;
		}
		else {
			node = Clusters[i][0];
			G->Clustering[0].push_back(node);
			G->Node_Cluster[node] = 0;
		}
	}
	*/
}

//AGENS::AGENS() {
//	ClusterNum = 5;
//	//InitCluster();
//}
//
//AGENS::AGENS(int n) :ClusterNum(n) {
//	//InitCluster();
//}

