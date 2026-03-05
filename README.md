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
@article{10.1016/j.future.2024.06.029,
author = {Wei, Junlin and Lin, Pengfei and Jiang, Jinrong and Liu, Hailong and Zhao, Lian and Zhang, Yehong and Han, Xiang and Zhang, Feng and Huang, Jian and Wang, Yuzhu and Li, Youyun and Yu, Yue and Chi, Xuebin},
title = {Accelerating LASG/IAP climate system ocean model version 3 for performance portability using Kokkos},
year = {2024},
issue_date = {Nov 2024},
publisher = {Elsevier Science Publishers B. V.},
address = {NLD},
volume = {160},
number = {C},
issn = {0167-739X},
url = {https://doi.org/10.1016/j.future.2024.06.029},
doi = {10.1016/j.future.2024.06.029},
journal = {Future Gener. Comput. Syst.},
month = nov,
pages = {901–917},
numpages = {17},
keywords = {High-performance computing, Ocean general circulation model, Performance portability}
}
```
