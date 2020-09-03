#include "EdgeList.H"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>

EdgeList::EdgeList()
{
	net = NULL;
}

EdgeList::~EdgeList()
{
	deleteNet();
}

int
EdgeList::deleteNet()
{
	if (net != NULL)
	{
		for (map<string,map<string,double>* >::iterator itr=net->begin(); itr!=net->end(); itr++)
		{
			map<string,double>* elist = itr->second;
			elist->clear();
			delete elist;
		}
		net->clear();
		delete net;
		net == NULL;
	}
	return 0;
}

NET_T*
EdgeList::getNet()
{
	return net;
}

int
EdgeList::setNet(NET_T* n)
{
	deleteNet();
	net = n;
	return 0;
}

int
EdgeList::addEdge(string& tname, string& gname, double v)
{
	if (net == NULL)
	{
		net = new NET_T;
	}
	map<string,double>* elist;
	if (net->find(tname) == net->end())
	{
		elist = new map<string,double>;
		(*net)[tname] = elist;
	}
	else
	{
		elist = (*net)[tname];
	}
	(*elist)[gname] = v;
	return 0;
}

int
EdgeList::countEdges()
{
	int ecount = 0;
	if (net == NULL)
	{
		return 0;
	}
	for (map<string,map<string,double>* >::iterator itr=net->begin(); itr!=net->end(); itr++)
	{
		ecount += itr->second->size();
	}
	return ecount;
}

EdgeList*
EdgeList::copyMe()
{
	NET_T* cnet = new NET_T;
	cnet->clear();
	for (map<string,map<string,double>* >::iterator itr=net->begin(); itr!=net->end(); itr++)
	{
		string name = itr->first;
		map<string,double>* elist = itr->second;
		map<string,double>* clist = new map<string,double>;
		clist->clear();
		for (map<string,double>::iterator eitr=elist->begin(); eitr!=elist->end(); eitr++)
		{
			(*clist)[eitr->first] = eitr->second;
		}
		(*cnet)[name] = clist;
	}
	EdgeList* e = new EdgeList;
	e->net = cnet;
	return e;
}

int
EdgeList::addSuffix(const char* suff)
{
	vector<string> names;
	for (map<string,map<string,double>* >::iterator itr=net->begin(); itr!=net->end(); itr++)
	{
		names.push_back(itr->first);
	}
	for (int i=0;i<names.size();i++)
	{
		string n = names[i];
		map<string,double>* elist = (*net)[n];
		n.append(suff);
		map<string,double>* nelist = new map<string,double>;
		nelist->clear();
		for (map<string,double>::iterator itr=elist->begin(); itr!=elist->end(); itr++)
		{
			(*nelist)[itr->first] = itr->second;
		}
		(*net)[n] = nelist;
	}
	return 0;
}

EdgeList* 
EdgeList::removeSuffix(const char* suff)
{
	NET_T* cnet = new NET_T;
	EdgeList* e = new EdgeList;
	e->net = cnet;
	cnet->clear();
	for (map<string,map<string,double>* >::iterator itr=net->begin(); itr!=net->end(); itr++)
	{
		string name = itr->first;
		int p = name.find(suff);
		string tf;
		if (p==-1)
		{
			tf = name;
		}
		else
		{
			tf = name.substr(0,p);
		}
		map<string,double>* elist = itr->second;
		map<string,double>* clist = NULL;
		if (cnet->find(tf) == cnet->end())
		{
			clist = new map<string,double>;
			clist->clear();
			(*cnet)[tf] = clist;
		}
		else
		{
			clist = (*cnet)[tf];
		}
		for (map<string,double>::iterator eitr=elist->begin(); eitr!=elist->end(); eitr++)
		{
			string gname = eitr->first;
			double val = eitr->second;
			if (clist->find(gname) == clist->end())
			{
				(*clist)[gname] = val;
			}
			else
			{
				if ((*clist)[gname] < val)
				{
					(*clist)[gname] = val;
				}
			}
		}
	}
	return e;
}

int
EdgeList::writeNet(char* oname)
{
	ofstream oFile(oname);
	for (map<string,map<string,double>* >::iterator itr=net->begin(); itr!=net->end(); itr++)
	{
		string name = itr->first;
		map<string,double>* elist = itr->second;
		for (map<string,double>::iterator eitr=elist->begin(); eitr!=elist->end(); eitr++)
		{
			oFile << name << "\t" << eitr->first << "\t" << eitr->second << endl;
		}
	}
	oFile.close();
	return 0;
}

int
EdgeList::setUnion(vector<EdgeList*>* outnets)
{
	deleteNet();
	net = new NET_T;
	net->clear();
	int osize = outnets->size();
	for (int i=0;i<osize;i++)
	{
		EdgeList* e = outnets->at(i);
		NET_T* onet = e->net;
		for (map<string,map<string,double>* >::iterator itr=onet->begin();itr!=onet->end();itr++)
		{
			string name = itr->first;
			map<string,double>* elist = itr->second;

			map<string,double>* rlist;
			if (net->find(name) == net->end())
			{
				rlist = new map<string,double>;
				rlist->clear();
				(*net)[name] = rlist;
			}
			else
			{
				rlist = (*net)[name];
			}
			for (map<string,double>::iterator gitr=elist->begin(); gitr!=elist->end(); gitr++)
			{
				string gname = gitr->first;
				if (rlist->find(gname) == rlist->end())
				{
					(*rlist)[gname] = fabs(gitr->second);
					//(*rlist)[gname] = 1;
				}
				else
				{
					cerr << "ERROR! in merging networks, edge " << name << "," << gname << " was repeated!" << endl;
				}
			}
		}
	}
	return 0;
}

int
EdgeList::getOverlap(EdgeList* onet)
{
	int oldcount = onet->countEdges();
	int newcount = countEdges();
	int overlap  = 0;
	NET_T* oldNet = onet->net;
	if (net == NULL || oldNet == NULL)
	{
		return 0;
	}
	for (map<string,map<string,double>* >::iterator itr=oldNet->begin(); itr!=oldNet->end(); itr++)
	{
		string name = itr->first;
		if (net->find(name) == net->end())
		{
			continue;
		}
		map<string,double>* oldelist = itr->second;
		map<string,double>* newelist = (*net)[name];
		for (map<string,double>::iterator eitr=oldelist->begin(); eitr!=oldelist->end(); eitr++)
		{
			string gname = eitr->first;
			if (newelist->find(gname)==newelist->end())
			{
				continue;
			}
			overlap++;
		}
	}
	return overlap;
}

int
EdgeList::setConsensus(vector<EdgeList*>* outnets)
{
	deleteNet();
	net = new NET_T;
	net->clear();
	int osize = outnets->size();
	for (int i=0;i<osize;i++)
	{
		EdgeList* e = outnets->at(i);
		NET_T* onet = e->net;
		for (map<string,map<string,double>* >::iterator itr=onet->begin();itr!=onet->end();itr++)
		{
			string name = itr->first;
			map<string,double>* elist = itr->second;

			map<string,double>* rlist;
			if (net->find(name) == net->end())
			{
				rlist = new map<string,double>;
				rlist->clear();
				(*net)[name] = rlist;
			}
			else
			{
				rlist = (*net)[name];
			}
			for (map<string,double>::iterator gitr=elist->begin(); gitr!=elist->end(); gitr++)
			{
				string gname = gitr->first;
				if (rlist->find(gname) == rlist->end())
				{
					(*rlist)[gname] = 0;
				}
				(*rlist)[gname] = (*rlist)[gname]+(1.0/((double)osize));
			}
		}
	}
	return 0;
}
