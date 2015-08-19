#include <cstdlib>
#include <sys/time.h>
#include <stopwatch.h>

struct timeval myStartTime;

void stopwatchStart()
{
	gettimeofday(&myStartTime, NULL);
}

double stopwatchReadSeconds()
{
	struct timeval endTime;
	gettimeofday(&endTime, 0);
    
	long ds = endTime.tv_sec - myStartTime.tv_sec;
	long dus = endTime.tv_usec - myStartTime.tv_usec;
	return ds + 0.000001*dus;
}


