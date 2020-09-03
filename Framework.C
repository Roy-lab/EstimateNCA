#include "Framework.H"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <unistd.h>
#include <algorithm>
#include <cstdlib>
#include <ctime>

#include <pthread.h>
#include <semaphore.h>

#include "ThreadManager.H"
#include "MemoryCheck.H"

Framework::Framework()
{
	regMngr = NULL;
	tgtMngr = NULL;
	evidMngr = NULL;
	priorNet = NULL;
	lambda = 0;
	CV=0;
	numNcaRep=10;
}

Framework::~Framework()
{
	if (regMngr!=NULL)
	{
		delete regMngr;
	}
	if (tgtMngr!=NULL)
	{
		delete tgtMngr;
	}
	if (evidMngr!=NULL)
	{
		delete evidMngr;
	}
	if (priorNet!=NULL)
	{
		delete priorNet;
	}
}

int
Framework::inferNCA(EdgeList* initNet)
{
	//In case of duplicate (both expression and TFA regulator),
	//I am taking Max. Maybe I should take mean?
	EdgeList* initNet2 = initNet->removeSuffix("_nca");
	char ioutdir[1025];
	sprintf(ioutdir,"%s/",outdir);
	Graph* updatedPrior = new Graph;
	updatedPrior->setVariableManagers(regMngr,tgtMngr);
	updatedPrior->setNet(initNet2);
	//NCALearner* ncalrnr = new NCALearner(regMngr, tgtMngr, evidMngr, priorNet, lambda, 10, numThread, ioutdir);
	NCALearner* ncalrnr = new NCALearner(regMngr, tgtMngr, evidMngr, updatedPrior, lambda, CV, numNcaRep, ioutdir);
	ncalrnr->start();
	VariableManager* nregMngr;
	EvidenceManager* ntfaMngr;
	ncalrnr->getData(nregMngr,ntfaMngr,"_nca");
	delete ncalrnr;

	regMngr->mergeVarSets(nregMngr);

	delete ntfaMngr;
	delete initNet2;

	return 0;
}

int
Framework::start()
{
	MemoryCheck mc;
	//Start with MotifTFA
	mc.begin();
	cout << "Starting NCA:" << endl;
	inferNCA(pnet);
	mc.end();
	mc.print("inferNCA");

	return 0;
}

int
Framework::init(int argc, char** argv)
{
	bool inputG=false;
	bool inputR=false;
	bool inputE=false;
	bool inputN=false;
	bool inputL=false;
	bool inputC=false;
	bool inputO=false;

	char gname[1025];
	char rname[1025];
	char pname[1025];
	char dname[1025];

	int optret='-';
	opterr=1;
	int oldoptind=optind;
	int condCnt=1;
	while(optret=getopt(argc,argv,"r:g:p:o:d:l:c:h")!=-1)
	{
		if(optret=='?')
		{
			cout <<"Option error " << optopt << endl;
			return -1;
		}
		char c;
		char* my_optarg=NULL;
		c=*(argv[oldoptind]+1);
		if(optind-oldoptind == 2)        //for -v 1
		{
			my_optarg=argv[oldoptind+1];	
		}
		else                             //for -v1
		{
			my_optarg=argv[oldoptind]+2;
		}
		switch(c)
		{
			case 'h': //Help
			{
				printHelp(argv[0]);
				return -1;
			}
			case 'o': //Output directory
			{
				strcpy(outdir,my_optarg);
				inputO=true;
				break;
			}
			case 'l'://Lambda
			{
				lambda = atof(my_optarg);
				inputL=true;
				break;
			}
			case 'c'://CV
			{
				CV = atoi(my_optarg);
				inputC=true;
				break;
			}
			case 'd'://Expression
			{
				strcpy(dname,my_optarg);
				inputE=true;
				break;
			}
			case 'p'://Prior net
			{
				strcpy(pname,my_optarg);
				inputN=true;
				break;
			}
			case 'g'://Target genes
			{
				strcpy(gname,my_optarg);
				inputG=true;
				break;
			}
			case 'r'://Regulators
			{
				strcpy(rname,my_optarg);
				inputR=true;
				break;
			}
			default:
			{
				cerr <<"Unhandled option " << c  << endl;
				printHelp(argv[0]);
				return -1;
			}
		}
		oldoptind=optind;
	}
	if (!inputE)
	{
		cerr << "Expression matrix was not provided!" << endl;
		printHelp(argv[0]);
		return -1;
	}
	if (!inputN)
	{
		cerr << "Prior network was not provided!" << endl;
		printHelp(argv[0]);
		return -1;
	}
	if (!inputG)
	{
		cerr << "List of target genes was not provided!" << endl;
		printHelp(argv[0]);
		return -1;
	}
	if (!inputR)
	{
		cerr << "List of regulators was not provided!" << endl;
		printHelp(argv[0]);
		return -1;
	}
	if (!inputO)
	{
		cerr << "Output directory was not provided!" << endl;
		printHelp(argv[0]);
		return -1;
	}
	numNcaRep = 1;
	if (!inputL && !inputC)
	{
		cerr << "Lambda was not provided" << endl;
		cerr << "Setting lambda to 0" << endl;
		lambda = 0;
	}
	if (inputL && inputC)
	{
		cerr << "Both Lambda and CV were provided" << endl;
		cerr << "Ignoring CV" << endl;
		CV=0;
	}

	regMngr = new VariableManager;
	regMngr->readVariables(rname);

	tgtMngr = new VariableManager;
	tgtMngr->readVariables(gname);

	evidMngr = new EvidenceManager;
	evidMngr->setVariableManager(tgtMngr);
	evidMngr->loadEvidenceFromFile_Continuous(dname);

	priorNet = new Graph;
	priorNet->setVariableManagers(regMngr,tgtMngr);
	priorNet->readNet(pname);

	pnet = priorNet->getNetMap();
	pnet->addSuffix("_nca");

	return 0;
}

int
Framework::printHelp(char* name)
{
	cerr << "Usage for " << name << endl;
	cerr << "-h\t\tprint this message" << endl;
	cerr << "-d\t\tInput expression, tab delimited, [GeneName][\\t][val][\\t][val]..." << endl;
	cerr << "-r\t\tList of regulators." << endl;
	cerr << "-g\t\tList of target genes." << endl;
	cerr << "-p\t\tPrior network, [TF][\\t][TG][\\t][Confidence]" << endl;
	cerr << "-l\t\tLambda penalty, for LASSO, higher means more sparse (default 0)" << endl;
	cerr << "-c\t\tNumber of cross validation folds (instead of fixed Lambda)" << endl;
	cerr << "-o\t\tOutput directory, it will over write existing files" << endl;
	return 0;
}

/*
int
Framework::makeRandomSubs(int maxCnt)
{
	Matrix* d = allEvidMngr->getDataMat();
	int maxDim = d->getColCnt();
	for (int r=0;r<maxCnt;r++)
	{
		vector<int> allids;
		for (int i=0;i<maxDim;i++)
		{
			allids.push_back(i);
		}
		random_shuffle(allids.begin(),allids.end());
		vector<int> remids;
		for (int i=0;i<maxDim/2;i++)
		{
			remids.push_back(allids[i]);
		}
		allRMIDs.push_back(remids);
	}
	return 0;
}
*/
