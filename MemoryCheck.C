#include "MemoryCheck.H"

MemoryCheck::MemoryCheck()
{
	b = false;
	e = false;
}

MemoryCheck::~MemoryCheck()
{
}

int
MemoryCheck::begin()
{
	b = true;
	getrusage(RUSAGE_SELF,&b_usage);
	return 0;
}

int
MemoryCheck::end()
{
	e = true;
	getrusage(RUSAGE_SELF,&e_usage);
	return 0;
}

int
MemoryCheck::print(const char* pref)
{
	if (b && e)
	{
		unsigned long int t1;
		unsigned long int t2;
		t1 = b_usage.ru_maxrss;
		t2 = e_usage.ru_maxrss;
		printf("Memory usage (%s) = b: %lu, e: %lu, %lu\n",pref, t1, t2, t2-t1);
	}
	b = false;
	e = false;
	return 0;
}
