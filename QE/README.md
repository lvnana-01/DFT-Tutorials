# Quantum ESPRESSO 中文教程：从入门到进阶（Quantum ESPRESSO Tutorial in Chinese）

系统讲解开源平面波赝势程序 Quantum ESPRESSO（QE）的第一性原理计算：conda-forge/源码/GPU 安装，pw.x 输入文件结构（&CONTROL/&SYSTEM/&ELECTRONS、ATOMIC_SPECIES、K_POINTS），SSSP/PseudoDojo/pslibrary 赝势选择，ecutwfc/ecutrho 与 k 点收敛测试，relax/vc-relax 与状态方程，bands.x/dos.x/projwfc.x 能带与 PDOS，磁性、非共线与 SOC，DFT+U/U+V（HUBBARD 卡）与 hp.x，范德华修正与功函数，差分电荷/ELF/Bader，neb.x 过渡态，ph.x/q2r.x/matdyn.x 声子（DFPT），介电/红外/拉曼，杂化泛函与 meta-GGA，Wannier90，EPW 电声耦合，pw.x/cp.x 分子动力学，以及并行参数（-nk/-nd/-nt）、自动化与报错排查。

**适用版本**：QE 7.6（2026-07 发布的最新正式版），兼顾 7.3.1 – 7.5；7.1 起 DFT+U 改用 `HUBBARD` 卡，7.6 起 CMake GPU 选项改为 `QE_GPU`/`QE_GPU_ARCHS`，相关章节均标注最低版本。7.6.1 与 8.0 截至 2026-11 尚未发布，文中只作预告。

> [!TIP]
> 建议先阅读 [通用基础](../Common/README.md) 了解 Linux 与集群作业提交，再按章节顺序学习；每章末尾都有参考资料与常见报错。QE 是 GPL 开源软件，没有 VASP 许可证时可以用它完整学习平面波 DFT 的全部流程。

> [!NOTE]
> 本板块只提供教程、输入模板与脚本，**不附带任何计算结果**；文中的输出片段都标注为“示意输出，数值仅供参考”。赝势文件（UPF）不随仓库分发，下载方法见 [第3章](tutorials/03-pseudopotentials.md)。

## 教程目录

| 章节 | 配套示例 |
|---|---|
| [第1章 QE 概述、版本与安装（Installing Quantum ESPRESSO: conda-forge, Source Build, GPU & WSL）](tutorials/01-install.md) | [01-install](examples/01-install/README.md) |
| [第2章 pw.x 输入文件结构与第一个 SCF（pw.x Input Anatomy & Your First SCF Calculation）](tutorials/02-input-anatomy.md) | [02-input-anatomy](examples/02-input-anatomy/README.md) |
| [第3章 赝势选择：SSSP、PseudoDojo 与 pslibrary（Pseudopotentials in QE: SSSP, PseudoDojo, pslibrary）](tutorials/03-pseudopotentials.md) | [03-pseudopotentials](examples/03-pseudopotentials/README.md) |
| [第4章 收敛测试：ecutwfc/ecutrho、k 点与 smearing（Convergence Tests in QE: Cutoffs, k-points & Smearing）](tutorials/04-convergence.md) | [04-convergence](examples/04-convergence/README.md) |
| [第5章 结构优化与状态方程（Structure Relaxation & Equation of State: relax, vc-relax, ev.x）](tutorials/05-relax-eos.md) | [05-relax-eos](examples/05-relax-eos/README.md) |
| [第6章 能带结构、态密度与投影态密度（Band Structure, DOS & PDOS: bands.x, dos.x, projwfc.x）](tutorials/06-band-dos-pdos.md) | [06-band-dos-pdos](examples/06-band-dos-pdos/README.md) |
| [第7章 磁性、非共线与自旋轨道耦合（Spin Polarization, Noncollinear Magnetism & SOC）](tutorials/07-magnetism-soc.md) | [07-magnetism-soc](examples/07-magnetism-soc/README.md) |
| [第8章 DFT+U/U+V 与 hp.x 线性响应（DFT+U, DFT+U+V & Linear-Response Hubbard Parameters with hp.x）](tutorials/08-dftu-hp.md) | [08-dftu-hp](examples/08-dftu-hp/README.md) |
| [第9章 范德华修正、二维材料、表面与功函数（vdW Corrections, 2D Materials, Surfaces & Work Function）](tutorials/09-vdw-surface-workfunction.md) | [09-vdw-surface-workfunction](examples/09-vdw-surface-workfunction/README.md) |
| [第10章 电荷密度、差分电荷、ELF、STM 与 Bader 分析（Charge Density, Difference, ELF, STM & Bader）](tutorials/10-charge-elf-bader.md) | [10-charge-elf-bader](examples/10-charge-elf-bader/README.md) |
| [第11章 过渡态搜索：neb.x（Nudged Elastic Band with neb.x）](tutorials/11-neb.md) | [11-neb](examples/11-neb/README.md) |
| [第12章 声子：ph.x、q2r.x 与 matdyn.x（Phonons with DFPT: ph.x, q2r.x, matdyn.x）](tutorials/12-phonon-dfpt.md) | [12-phonon-dfpt](examples/12-phonon-dfpt/README.md) |
| [第13章 介电常数、Born 有效电荷、红外与拉曼（Dielectric Tensor, Born Charges, IR & Raman）](tutorials/13-dielectric-ir-raman.md) | [13-dielectric-ir-raman](examples/13-dielectric-ir-raman/README.md) |
| [第14章 杂化泛函与 meta-GGA（Hybrid Functionals & meta-GGA: HSE, PBE0, ACE, SCAN）](tutorials/14-hybrid-metagga.md) | [14-hybrid-metagga](examples/14-hybrid-metagga/README.md) |
| [第15章 Wannier90 接口与插值能带（Maximally Localized Wannier Functions with Wannier90）](tutorials/15-wannier90.md) | [15-wannier90](examples/15-wannier90/README.md) |
| [第16章 电声耦合与超导：ph.x 与 EPW 入门（Electron-Phonon Coupling & Superconductivity: ph.x and EPW）](tutorials/16-elph-epw.md) | [16-elph-epw](examples/16-elph-epw/README.md) |
| [第17章 分子动力学：pw.x MD 与 Car-Parrinello cp.x（Ab Initio MD: Born-Oppenheimer with pw.x and Car-Parrinello with cp.x）](tutorials/17-md-cp.md) | [17-md-cp](examples/17-md-cp/README.md) |
| [第18章 自动化、并行性能与报错排查（Automation, Parallel Performance & Troubleshooting）](tutorials/18-automation-parallel-troubleshooting.md) | [18-automation-parallel](examples/18-automation-parallel/README.md) |

## 全部示例

| 示例目录 | 说明 |
|---|---|
| [01-install](examples/01-install/README.md) | 示例 01：QE 安装脚本与 Si 冒烟测试 |
| [02-input-anatomy](examples/02-input-anatomy/README.md) | 示例 02：pw.x 输入结构与第一个 SCF（Si、Al、H2O） |
| [03-pseudopotentials](examples/03-pseudopotentials/README.md) | 示例 03：赝势下载、检查与 MgO 多元素截断（SSSP / PseudoDojo） |
| [04-convergence](examples/04-convergence/README.md) | 示例 04：Si 与 Al 的截断能、k 点与 smearing 收敛测试 |
| [05-relax-eos](examples/05-relax-eos/README.md) | 示例 05：Si 变胞优化、CO 分子优化与 Cu 状态方程 |
| [06-band-dos-pdos](examples/06-band-dos-pdos/README.md) | 示例 06：Si 的能带、态密度与投影态密度 |
| [07-magnetism-soc](examples/07-magnetism-soc/README.md) | 示例 07：自旋极化、反铁磁、非共线与自旋轨道耦合（bcc Fe / fcc Ni / NiO / GaAs） |
| [08-dftu-hp](examples/08-dftu-hp/README.md) | 示例 08：DFT+U、DFT+U+V 与 hp.x 线性响应（NiO / LiCoO2） |
| [09-vdw-surface-workfunction](examples/09-vdw-surface-workfunction/README.md) | 示例 09：范德华修正、二维材料、表面 slab 与功函数（graphite / MoS2 / Al(111)） |
| [10-charge-elf-bader](examples/10-charge-elf-bader/README.md) | 示例 10：电荷密度、差分电荷、ELF、STM 与 Bader 分析 |
| [11-neb](examples/11-neb/README.md) | 示例 11：H 原子在 Al(100) 表面相邻四重空位之间扩散的 CI-NEB |
| [12-phonon-dfpt](examples/12-phonon-dfpt/README.md) | 示例 12：Si 与 Al 的 DFPT 声子（ph.x → q2r.x → matdyn.x） |
| [13-dielectric-ir-raman](examples/13-dielectric-ir-raman/README.md) | 示例 13：GaAs 介电张量、Born 有效电荷、LO-TO、拉曼与 H2O 红外/拉曼 |
| [14-hybrid-metagga](examples/14-hybrid-metagga/README.md) | 示例 14：Si/MgO 的 HSE、PBE0 与 r2SCAN，以及 HSE 能带的 Wannier 插值 |
| [15-wannier90](examples/15-wannier90/README.md) | 示例 15：Si 的最大局域 Wannier 函数与插值能带 |
| [16-elph-epw](examples/16-elph-epw/README.md) | 示例 16：电声耦合与超导（ph.x + lambda.x 路线 / EPW 路线） |
| [17-md-cp](examples/17-md-cp/README.md) | 示例 17：pw.x 分子动力学与 cp.x Car-Parrinello 动力学 |
| [18-automation-parallel](examples/18-automation-parallel/README.md) | 示例 18：并行参数扫描、Slurm 模板、ASE 驱动与 XML 解析 |

## 辅助脚本

完整说明见 [scripts/README.md](scripts/README.md)。所有脚本都有 `-h` 帮助；绘图与解析脚本只读取已有输出，不会启动计算。

- [qe_inputgen.py](scripts/qe_inputgen.py)：从 CIF/POSCAR 生成 pw.x 输入（ibrav=0、显式单位、SSSP 截断、k 间距、crystal_b 路径）
- [qe_check_input.py](scripts/qe_check_input.py)：静态检查 pw.x 输入（namelist、CRLF、单位、7.1 起删除的 DFT+U 旧关键字、赝势文件）
- [qe_conv_test.sh](scripts/qe_conv_test.sh)：批量生成 ecut / k 点 / smearing 扫描输入并汇总总能
- [qe_parse_out.py](scripts/qe_parse_out.py)：解析 pw.x 输出（总能、费米能、力、应力、磁矩、最终结构），导出 JSON/CSV
- [qe_plot_bands.py](scripts/qe_plot_bands.py)：能带作图与带隙
- [qe_plot_dos.py](scripts/qe_plot_dos.py)：DOS/PDOS 作图
- [qe_phonon_plot.py](scripts/qe_phonon_plot.py)：声子色散与声子 DOS 作图，标记虚频
- [qe_eos_fit.py](scripts/qe_eos_fit.py)：Birch-Murnaghan / Murnaghan 状态方程拟合
- [qe_workfunction.py](scripts/qe_workfunction.py)：平面平均势与功函数
- [qe_cube_diff.py](scripts/qe_cube_diff.py)：cube 差分电荷与平面平均
- [job_qe_slurm.sh](scripts/job_qe_slurm.sh)：Slurm 作业模板
- [job_qe_pbs.sh](scripts/job_qe_pbs.sh)：PBS/Torque 作业模板

## 推荐的外部课程

下列课程只给出链接，不复制其讲义或输入文件；使用前请留意各自的版本与许可证。

- [Cornell PARADIM Summer School 2016: Introduction to DFT for Experimentalists](https://www.paradim.org/cu_ss_2016)（Giustino）：pw.x/ph.x 上机，覆盖平衡结构、弹性常数、声子与能带。基于 QE 5.x/6.x 与 LDA，原理讲解好，下载/编译/赝势部分已过时。
- [Cornell PARADIM Summer School 2018](https://www.paradim.org/cu_ss_2018)（Giustino、Pamuk）：Si、Cl2、金刚石、GaAs 的声子、LO-TO、介电与光谱，每讲有视频。
- [Cornell PARADIM Computational Tutorials](https://www.paradim.org/toolbox/computation/tutorials)：QE 配置、输入与收敛参数、平衡结构、弹性常数、声子、能带，另有 Yambo GW/BSE 系列。
- [Cornell PARADIM Summer School 2023](https://www.paradim.org/2023_DFT_SS)：理论讲座值得看；注意上机用的是 JDFTx，不是 QE。
- [Imperial College MSE404-MM: Modelling Materials with DFT](https://imperialcollegelondon.github.io/MSE404-MM/)：9 个 QE 实验（输入输出、收敛、能带、金属/DOS、自旋、优化、声子），持续维护，强烈推荐。
- [Computational Materials Physics（Ghent，Cottenier）](https://www.compmatphys.org)：零基础 DFT 概念课，练习用 QE。
- [MaX e-School 2021（ICTP / MaX / QEF）](https://indico.ictp.it/event/9616/)：两周 QE 课程与视频。
- [Advanced QE School 2023（EPFL THEOS / Materials Cloud）](https://github.com/materialscloud-org/QuantumESPRESSO-school-2023)：杂化泛函、Wannier90、DFPT、DFT+U+V 与 hp.x、Koopmans、AiiDA。
- [ENCCS: Efficient materials modelling on HPC](https://enccs.github.io/efficient-materials-modelling-on-hpc/)：-nk/-nd/-nt 并行、GPU 与 Yambo。
- [EPW School / Tutorials](https://docs.epw-code.org/doc/Tutorial.html)：电声耦合、迁移率与 Eliashberg 超导。
- [pranabdas/espresso](https://pranabdas.github.io/espresso/)：从安装到 EPW 的个人教程站，已使用新 HUBBARD 语法。
- [QE 官方 Learn 页](https://www.quantum-espresso.org/learn/)：官方推荐的课程与学校列表。

## 相关资源

- [GitHub 相关仓库精选（含 QE 专节）](../Resources/github-repos.md)
- [计算案例：酶催化 QM 簇模型等完整研究流程](../Cases/README.md)
- [返回项目首页](../README.md)
