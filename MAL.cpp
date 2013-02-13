// MAL.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <bitset>
#include <map>
#include <algorithm>
using namespace std;

const int MAX=16;

class ReservationTable
{
	int stages,time;
	multimap<int,int> rtable;
protected:
	list<int> forbidden;
public:
	void readtable();
};

struct edge;

struct header
{
	bitset<MAX> bset;
	header *next;
	edge *pedge;
};

struct edge
{
	list<int> weights;
	edge *next;
	header *to;
};

class StateDiagram :public ReservationTable
{
	int cvsize;
	bitset<MAX> initialcv;
	void addstate(header *state);
	void insertedge(header *state);
protected:
	header *start;
public:
	void buildsd();
};

class Cycle :public StateDiagram
{
	double mal;
	list<list<int>> simplecycles;
	list<list<int>> greedycycles;
	void callsimple(header *node, list<bitset<MAX>> trav, list<int> cycle);
	void callgreedy(header *node, list<bitset<MAX>> trav, list<int> cycle);
public:
	void process();
};

void ReservationTable::readtable()
{
	ifstream fin;
	fin.open("ReservationTable.txt",ifstream::in);
	if(!fin) 
	{
		cerr<<"Cannot open INVENTORY file!"<<endl;
		system("pause");
		exit(1);	
	}
	string str;
	getline(fin,str);
	stages=stoi(str);
	getline(fin,str);
	time=stoi(str);
	int i=0,j;
	cout<<"Reservation Table:\n  ";
	for(j=1; j<=time; ++j)
		cout<<"  "<<j;
	cout<<endl;
	while(!fin.eof())
	{
		getline(fin,str);
		++i;
		j=0;
		cout<<"S"<<i<<" ";
		for (string::iterator it=str.begin(); it!=str.end(); ++it)
		{
			cout<<" ";
			if(*it=='1')
			{ cout<<"X"; rtable.insert(pair<int,int>(i,++j)); }
			else if(*it=='0')
			{ cout<<" "; ++j; }
		}
		cout<<endl;
	}
	fin.close();
	vector<int> row;
	for(i=1; i<=stages; i++)
    {
		multimap<int,int>::iterator it;
		for(it=rtable.equal_range(i).first; it!=rtable.equal_range(i).second; ++it)
		row.push_back((*it).second);
		for(size_t p=0; p<row.size(); ++p)
		for(size_t q=p+1; q<row.size(); ++q)
		{
			list<int>::iterator f=find(forbidden.begin(),forbidden.end(),row[q]-row[p]);
			if(f==forbidden.end())
			forbidden.push_back(row[q]-row[p]);
		}
		row.clear();
	}
	stable_sort(forbidden.begin(),forbidden.end());
	cout<<"\nForbidden latencies: ";
	for (list<int>::iterator it=forbidden.begin(); it!=forbidden.end(); ++it)
    cout<<" "<<*it;
	cout<<endl;
}

void StateDiagram::addstate(header *state)
{
	bitset<MAX> temp;
	header *ptr;
	for(int i=0;i<cvsize;i++)
	{
		if(!state->bset[i])
		{
			temp=state->bset;
			temp>>=i+1;
			temp|=start->bset;
			header *node=new header;
			node->bset=temp;
			node->next=nullptr;
			node->pedge=nullptr;
			for(ptr=start; ptr->bset!=node->bset && ptr->next; ptr=ptr->next);
			if(ptr->bset!=node->bset)
			ptr->next=node;
		}
	}
}

void StateDiagram::insertedge(header *state)
{
	bitset<MAX> temp;
	header *ptr;
	edge *head=nullptr;
	edge *e;
	for(int i=0;i<cvsize;i++)
	{
		if(!state->bset[i])
		{
			temp=state->bset;
			temp>>=i+1;
			temp|=start->bset;
			for(ptr=start; ptr->bset!=temp; ptr=ptr->next);
			if(!head)
			{
				head=new edge;
				head->weights.push_back(i+1);
				head->to=ptr;
				head->next=nullptr;
			}
			else
			{
				for(e=head; e->to!=ptr,e->next; e=e->next);
				if(e->to!=ptr)
				{
					edge *node=new edge;
					node->weights.push_back(i+1);
					node->next=nullptr;
					node->to=ptr;
					e->next=node;
				}
				else
				{
					e->weights.push_back(i+1);
				}
			}
		}
	}
	if(!head)
	{
		head=new edge;
		head->weights.push_back(cvsize+1);
		head->next=nullptr;
		head->to=start;
	}
	else
	{
		for(e=head;e->to!=start && e->next;e=e->next);
		if(e->to==start)
			e->weights.push_back(cvsize+1);
		else
		{
			edge *node=new edge;
			node->weights.push_back(cvsize+1);
			node->next=nullptr;
			node->to=start;
			e->next=node;
		}
	}
	state->pedge=head;
	cout<<"State: "<<state->bset<<endl;
	for(e=head;e;e=e->next)
	{
		cout<<"\t";
		for(list<int>::iterator it=e->weights.begin(); it!=e->weights.end();)
		{
			cout<<*it;
			if(++it!=e->weights.end()) cout<<",";
		}
		cout<<"\t=> "<<e->to->bset;
		cout<<endl;
	}
}

void StateDiagram::buildsd()
{
	for (list<int>::iterator it=forbidden.begin(); it!=forbidden.end(); ++it)
	initialcv.set(*it-1,1);
	cout<<"\nInitial Collision Vector: "<<initialcv<<endl;
	start=new header;
	start->bset=initialcv;
	start->next=nullptr;
	start->pedge=nullptr;
	cvsize=forbidden.back();
	header *ptr;
	cout<<"\nState Diagram:"<<endl;
	for(ptr=start;ptr;ptr=ptr->next)
	{
		addstate(ptr);
		insertedge(ptr);
	}
	cout<<endl;
}

void Cycle::callsimple(header *node, list<bitset<MAX>> trav,list<int> cycle)
{
	list<bitset<MAX>>::iterator it;
	trav.push_back(node->bset);
	for(edge *e=node->pedge; e; e=e->next)
	{
		cycle.push_back(e->weights.front());
		it=find(trav.begin(),trav.end(),e->to->bset);
		if(it==trav.end())
			callsimple(e->to,trav,cycle);
		else if(it==trav.begin())
			simplecycles.push_back(cycle);
		cycle.pop_back();
	}
}
void Cycle::callgreedy(header *node, list<bitset<MAX>> trav,list<int> cycle)
{
	list<bitset<MAX>>::iterator it;
	trav.push_back(node->bset);
	edge *e=node->pedge;
	cycle.push_back(e->weights.front());
	it=find(trav.begin(),trav.end(),e->to->bset);
	if(it==trav.end())
		callgreedy(e->to,trav,cycle);
	else if(it==trav.begin())
		greedycycles.push_back(cycle);
	cycle.pop_back();
}
void Cycle::process()
{
	list<bitset<MAX>> traversed;
	list<int> cycle;
	for(header *ptr=start; ptr; ptr=ptr->next)
	{
		cycle.clear();
		traversed.clear();
		traversed.push_back(ptr->bset);
		for(header *p=start;p!=ptr;p=p->next)
			traversed.push_back(p->bset);
		callsimple(ptr,traversed,cycle);
		callgreedy(ptr,traversed,cycle);
	}
	list<list<int>>::iterator it1;
	list<int>::iterator it2;
	cout<<"\nSimple Cycles:"<<endl;
	for(it1=simplecycles.begin(); it1!=simplecycles.end(); ++it1)
	{
		cout<<"(";
		for(it2=(*it1).begin(); it2!=(*it1).end();)
		{
			cout<<*it2;
			if(++it2!=(*it1).end()) cout<<",";
		}
		cout<<")"<<endl;
	}
	double min;
	mal=INT_MAX;
	cout<<"\nGreedy Cycles:"<<endl;
	for(it1=greedycycles.begin(); it1!=greedycycles.end(); ++it1)
	{
		min=0;
		cout<<"(";
		for(it2=(*it1).begin(); it2!=(*it1).end();)
		{
			cout<<*it2;
			min+=(*it2);
			if(++it2!=(*it1).end()) cout<<",";
		}
		min/=(*it1).size();
		if(min<mal) mal=min;
		cout<<")"<<endl;
	}
	cout<<"\nMinimum Average Latency: "<<mal<<endl;
}

int main()
{
	Cycle obj;
	obj.readtable();
	obj.buildsd();
	obj.process();
	system("pause");
	return 0;
}

