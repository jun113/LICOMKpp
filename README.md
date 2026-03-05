# LICOMK++: A performance-portable ocean general circulation model (LICOM) using Kokkos, to facilitate global Kilometer-scale ocean simulations

**ACM Gordon Bell Prize Finalist for Climate Modelling'2024**

A Performance-Portable Kilometer-Scale Global Ocean Model on ORISE and New Sunway Heterogeneous Supercomputers

## 📊 Performance Highlights

| System | Resolution | Cores/GPUs | Performance (SYPD) | Parallel Efficiency |
|--------|------------|------------|-------------------|-------------------|
| Sunway OceanLight | 1 km | 38,366,250 cores | 1.05 | 54.8% |
| ORISE | 1 km | 16,000 HIP GPUs | 1.70 | 55.6% |

## 🏛️ Architecture Support

| Architecture | Backend | Status |
|-------------|------------------|---------|
| NVIDIA GPUs | CUDA | ✅ Supported |
| AMD GPUs | HIP | ✅ Supported |
| Hygon DCUs | HIP | ✅ Supported |
| Intel CPUs | OpenMP | ✅ Supported |
| ARM CPUs | OpenMP | ✅ Supported |
| Sunway Many-cores | Athread | ✅ Supported |

## 📄 Citation

If you use LICOMK++ in your research, please cite:
```bibtex
@inproceedings{10.1109/SC41406.2024.00009,
author = {Wei, Junlin and Han, Xiang and Yu, Jiangfeng and Jiang, Jinrong and Liu, Hailong and Lin, Pengfei and Yu, Maoxue and Xu, Kai and Zhao, Lian and Wang, Pengfei and Zheng, Weipeng and Xie, Jingwei and Zhou, Yanzhi and Zhang, Tao and Zhang, Feng and Zhang, Yehong and Yu, Yue and Wang, Yuzhu and Bai, Yidi and Li, Chen and Yu, Zipeng and Deng, Haoyu and Li, Yaxin and Chi, Xuebin},
title = {A Performance-Portable Kilometer-Scale Global Ocean Model on ORISE and New Sunway Heterogeneous Supercomputers},
year = {2024},
isbn = {9798350352917},
publisher = {IEEE Press},
url = {https://doi.org/10.1109/SC41406.2024.00009},
doi = {10.1109/SC41406.2024.00009},
abstract = {Ocean general circulation models (OGCMs) are indispensable for studying the multi-scale oceanic processes and climate change. High-resolution ocean simulations require immense computational power and thus become a challenge in climate science. We present LICOMK++, a performance-portable OGCM using Kokkos, to facilitate global kilometer-scale ocean simulations. The breakthroughs include: (1) we enhance cutting-edge Kokkos with the Sunway architecture, enabling LICOMK++ to become the first performance-portable OGCM on diversified architectures, i.e., Sunway processors, CUDA/HIP-based GPUs, and ARM CPUs. (2) LICOMK++ overcomes the one simulated-years-per-day (SYPD) performance challenge for global realistic OGCM at 1-km resolution. It records 1.05 and 1.70 SYPD with a parallel efficiency of 54.8\% and 55.6\% scaling on almost the entire new Sunway supercomputer and two-thirds of the ORISE supercomputer. (3) LICOMK++ is the first global 1-km-resolution realistic OGCM to generate scientific results. It successfully reproduces mesoscale and submesoscale structures that have considerable climate effects.},
booktitle = {Proceedings of the International Conference for High Performance Computing, Networking, Storage, and Analysis},
articleno = {3},
numpages = {12},
keywords = {High-Performance Computing, Ocean General Circulation Model, Performance-Portable, Sunway Architecture},
location = {Atlanta, GA, USA},
series = {SC '24}
}
```
```bibtex
@article{WEI2024901,
title = {Accelerating LASG/IAP climate system ocean model version 3 for performance portability using Kokkos},
journal = {Future Generation Computer Systems},
volume = {160},
pages = {901-917},
year = {2024},
issn = {0167-739X},
doi = {https://doi.org/10.1016/j.future.2024.06.029},
url = {https://www.sciencedirect.com/science/article/pii/S0167739X24003285},
author = {Junlin Wei and Pengfei Lin and Jinrong Jiang and Hailong Liu and Lian Zhao and Yehong Zhang and Xiang Han and Feng Zhang and Jian Huang and Yuzhu Wang and Youyun Li and Yue Yu and Xuebin Chi},
keywords = {High-performance computing, Ocean general circulation model, Performance portability},
abstract = {In this paper, the performance portability of the LASG/IAP Climate System Ocean Model version 3 (LICOM3) is demonstrated based on the C++ library Kokkos. Kokkos enables application execution in various High-Performance Computing (HPC) architectures for on-node parallelism. This study employs Kokkos to expose on-node parallelism and reuses pre-existing Message-Passing Interface (MPI) for internode parallelism. By porting to Kokkos, the single-source code LICOM3 is successfully executed on ARM CPUs, Tesla V100, and HIP-based GPUs. To this end, the characteristics and mechanisms of LICOM3 and Kokkos are considered, and the model is then optimized comprehensively in terms of data management, computation, and memory transmission. The proposed Kokkos optimization code at a 1∘ resolution accelerates operation by factors of 1.9, 1.2, and 1.1 compared to the raw Compute Unified Device Architecture (CUDA), Heterogeneous Interface for Portable (HIP) and OpenMP codes, respectively. Further, it exhibits 3.4 Simulated Years Per Day (SYPD) at a resolution of 0.05∘ when executed on 4096 HIP-based GPUs for large-scale simulations.}
}
```
