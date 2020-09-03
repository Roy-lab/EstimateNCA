#ifndef _EDGELIST
#define _EDGELIST

#include "CommonTypes.H"
#include <vector>
#include <map>
#include <string>
using namespace std;

class EdgeList
{
	public:
		EdgeList();
		~EdgeList();
		int addEdge(string& tname, string& gname, double v);
		NET_T* getNet();
		int setNet(NET_T*);
		int countEdges();
		int setConsensus(vector<EdgeList*>*);
		EdgeList* compyMe();
		int getOverlap(EdgeList* onet);
		int addSuffix(const char* suff);
		int setUnion(vector<EdgeList*>* outnets);
		int writeNet(char* oname);
		EdgeList* removeSuffix(const char* suff);
		EdgeList* copyMe();
	private:
		int deleteNet();
		map<string,map<string,double>* >* net;
};

#endif //_EDGELIST
