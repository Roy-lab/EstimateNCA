#include "ThreadManager.H"

ThreadManager::ThreadManager(int n, void *(*r)(void *))
{
	start_routine = r;
	maxpcnt = n;
	sem_init(&qsem, 0, 0);
	pthread_mutex_init(&qlock, NULL);
	inputs.clear();
}

ThreadManager::~ThreadManager()
{
	sem_destroy(&qsem);
	pthread_mutex_destroy(&qlock);
}

int 
ThreadManager::addInput(void *v)
{
	//void** ptr = new void*[2];
	//ptr[0] = (void*)this;
	//ptr[1] = v;
	//inputs.push_back((void*)ptr);
	inputs.push_back(v);
	return 0;
}

int
ThreadManager::run()
{
	int pcnt = maxpcnt;
	int netcnt = inputs.size();
	while (netcnt>0)
	{
		if (pcnt > 0)
		{
			void * v = inputs[netcnt-1];
			pthread_t one_t;
			pthread_create(&one_t, NULL, start_routine, v);
			pcnt--;
			netcnt--;
		}
		else
		{
			sem_wait(&qsem);
			pcnt++;
		}
	}
	while(pcnt < maxpcnt)
	{
		sem_wait(&qsem);
		pcnt++;
	}
	inputs.clear();
	return 0;
}

pthread_mutex_t&
ThreadManager::getLock()
{
	return qlock;
}

sem_t&
ThreadManager::getSem()
{
	return qsem;
}
