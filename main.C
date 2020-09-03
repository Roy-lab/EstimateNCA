#include "Framework.H"
#include <iostream>

int
main (int argc, char **argv)
{
	srand ( unsigned ( time(0) ) );
	Framework* fw = new Framework;
	if(fw->init(argc,argv)==0)
	{
		fw->start();
	}
	delete fw;
	return 0;
}
