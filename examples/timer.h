
#pragma once
#include <sys/time.h>
#include <unistd.h>

unsigned long status_ms(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    unsigned long long msSinceEpoch = (unsigned long long)(tv.tv_sec) * 1000
       + (unsigned long long)(tv.tv_usec) / 1000;
    return msSinceEpoch;
}

