#ifndef _THREADMANAGER
#define _THREADMANAGER

#include <pthread.h>
#include <semaphore.h>
#include <vector>

using namespace std;

class ThreadManager
{
	public:
		ThreadManager(int n, void *(*r)(void *));
		~ThreadManager();
		int addInput(void *);
		int run();
		sem_t& getSem();
		pthread_mutex_t& getLock();
	private:
		int maxpcnt;
		pthread_mutex_t qlock;
		sem_t qsem;
		void *(*start_routine) (void *);
		vector<void*> inputs;
};

#endif //_THREADMANAGER
