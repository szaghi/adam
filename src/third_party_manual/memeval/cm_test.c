#define _GNU_SOURCE
#include <stdio.h>
#include <sys/sysinfo.h>
#include <unistd.h>

void getmemory(long int* totmem, long int* avmem);

int main(void) {
    long int a, b;
    getmemory(&a, &b);
    printf("memory: %ld/%ld ", b,a);
}

