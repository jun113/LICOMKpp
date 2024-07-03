#pragma once
#ifdef __cplusplus
extern "C" {
#endif

#define SPAWN_CGS 6

struct spawn_args {
    void *fnptr;
    void *arg;
    volatile int flag;
};

extern __uncached struct spawn_args g_spawn_args[SPAWN_CGS];
extern __uncached int spawn_proxy_n_cg;

void spawn_proxy_init();
void spawn_proxy_init_cgs();
void spawn_proxy_run(void *fnptr, void *arg);
void spawn_proxy_run_cgs(void *fnptr, void *arg);
void spawn_proxy_run_mpecg(void *fnptr, void *arg, int mpecg);
void spawn_proxy_run_at(void *fnptr, void *arg, int mpecg);
void spawn_proxy_join();
void spawn_proxy_join_at(int cpeg);
void spawn_proxy_join_cgs();
void spawn_proxy_join_mpecg(int mpecg);
void spawn_proxy_finalize();

#define spawn_proxy_run(a, b) __real_spawn_proxy_run(slave_##a, b)
#define spawn_proxy_run_cgs(a, b) __real_spawn_proxy_run_cgs(slave_##a, b)
#define spawn_proxy_run_mpecg(a, b, c) __real_spawn_proxy_run_mpecg(slave_##a, b, c)
#define spawn_proxy_run_at(a, b, c) __real_spawn_proxy_run_at(slave_##a, b, c)

void __real_spawn_proxy_run_cgs(void *fnptr, void *arg);
void __real_spawn_proxy_run(void *fnptr, void *arg);
void __real_spawn_proxy_run_mpecg(void *fnptr, void *arg, int mpecg);
void __real_spawn_proxy_run_at(void *fnptr, void *arg, int cpeg);

#ifdef __cplusplus
}
#endif
