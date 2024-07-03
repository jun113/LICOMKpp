#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <crts.h>
//#include <swblas_common.h>

#include <signal.h>

// #include <mpi.h>

#include "spawn.h"

extern void slave_spawn_proxy_stub(void *);
__uncached struct spawn_args g_spawn_args[SPAWN_CGS];
int spawn_proxy_n_cg = 1;

int spawn_proxy_initialized = 0;

static void sighandler(int _) {
    printf("Sig 55\n");
    while (1);
}

void spawn_proxy_auto_init();

void spawn_proxy_init() {
    athread_init();
    for (int i = 0; i < 1; ++i)
        g_spawn_args[i].flag = 0;
    athread_spawn(spawn_proxy_stub, NULL);
}
void spawn_proxy_init_cgs() {
    athread_init_cgs();
    for (int i = 0; i < SPAWN_CGS; ++i)
        g_spawn_args[i].flag = 0;
    athread_spawn_cgs(spawn_proxy_stub, NULL);
}
void __real_spawn_proxy_run(void *fnptr, void *arg) {
//    spawn_proxy_auto_init();
    g_spawn_args[0].fnptr = fnptr;
    g_spawn_args[0].arg = arg;
    asm volatile("memb");
    g_spawn_args[0].flag = 1;
}
void __real_spawn_proxy_run_cgs(void *fnptr, void *arg) {
    spawn_proxy_auto_init();
    for (int i = 0; i < SPAWN_CGS; ++i) {
        g_spawn_args[i].fnptr = fnptr;
        g_spawn_args[i].arg = arg;
    }
    asm volatile("memb");
    for (int i = 0; i < SPAWN_CGS; ++i)
        g_spawn_args[i].flag = 1;
}
void __real_spawn_proxy_run_at(void *fnptr, void *arg, int cpeg) {
    spawn_proxy_auto_init();
    g_spawn_args[cpeg].fnptr = fnptr;
    g_spawn_args[cpeg].arg = arg;
    asm volatile("memb");
    g_spawn_args[cpeg].flag = 1;
}
void __real_spawn_proxy_run_mpecg(void *fnptr, void *arg, int mpecg) {
    spawn_proxy_auto_init();
    for (int i = 0; i < mpecg; ++i) {
        g_spawn_args[i].fnptr = fnptr;
        g_spawn_args[i].arg = arg;
    }
    asm volatile("memb");
    for (int i = 0; i < mpecg; ++i)
        g_spawn_args[i].flag = 1;
}    
void spawn_proxy_join() {
    for (int i = 0; i < 1; ++i)
        while (g_spawn_args[i].flag)
            ;
}
void spawn_proxy_join_at(int cpeg) {
    while (g_spawn_args[cpeg].flag)
        ;
}
void spawn_proxy_join_cgs() {
    for (int i = 0; i < SPAWN_CGS; ++i)
        while (g_spawn_args[i].flag)
            ;
}
void spawn_proxy_join_mpecg(int mpecg) {
    for (int i = 0; i < mpecg; ++i)
        while (g_spawn_args[i].flag)
            ;
}    
void spawn_proxy_finalize() {
    for (int i = 0; i < spawn_proxy_n_cg; ++i)
        g_spawn_args[i].flag = -1;
    asm volatile("memb");
    if (spawn_proxy_n_cg == 1)
        athread_join();
    else if (spawn_proxy_n_cg == 6)
        athread_join_cgs();
    else if (spawn_proxy_n_cg != 0)
        assert(0);
        //swblas_unimplemented();
}

__attribute__((constructor)) void spawn_proxy_prepare_env() {
    spawn_proxy_n_cg = atoi(getenv("RMS_MPE_CGS"));
}

void spawn_proxy_auto_init() {
    if (spawn_proxy_initialized)
        return;
    spawn_proxy_initialized = 1;

    if (spawn_proxy_n_cg == 1) 
        spawn_proxy_init();
    else if (spawn_proxy_n_cg == 6)
        spawn_proxy_init_cgs();
    else 
        assert(0);
}
