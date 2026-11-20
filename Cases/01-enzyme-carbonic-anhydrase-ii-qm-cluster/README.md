# 酶催化计算实例：Human carbonic anhydrase II (hCA II) 反应机理的 QM 簇模型研究（Enzyme Catalysis Case Study）

本案例以人碳酸酐酶 II（human carbonic anhydrase II, hCA II）为例，完整走一遍酶反应机理的 **QM 簇模型（quantum chemical cluster approach）** 研究流程。流程从 PDB 晶体结构（2CBA，以及 CO2 结合态 3D92）出发，依次是质子化态判断（PROPKA/PDB2PQR）、截取约 50–60 原子的 Zn(II) 活性位点簇模型（4-甲基咪唑截断、H 封端、冻结 Cβ），GFN2-xTB 预优化，B3LYP-D3(BJ)/def2-SVP + CPCM(ε=4) 几何优化、过渡态搜索、IRC 与频率，最后用 def2-TZVP 单点、能量剖面、模型收敛性检验、ONIOM/QM/MM 替代方案和 NBO/IGMH 分析收尾。关键词：enzyme mechanism、cluster model、Zn-hydroxide、CO2 hydration、transition state、Gaussian 16、ORCA 6、xTB。

> [!IMPORTANT]
> **本案例没有运行任何量子化学计算。** 目录中只有工作流、输入模板和脚本；模板中的坐标用 `@COORDS@` 等占位符表示，需由 [`00-structure/cut_cluster.py`](00-structure/cut_cluster.py) 和你自己的预优化结果生成。文中凡出现的数值，要么是晶体结构中可直接量出的距离（已注明 PDB 编号），要么是带 DOI 的**文献值**，要么明确标注为“示意输出，数值仅供参考”。请不要把任何示意数值当作结果引用。

---

## 目录

- [1. 本章目标与适用场景](#1-本章目标与适用场景)
- [2. 背景与反应机理](#2-背景与反应机理)
- [3. 为什么用 QM 簇模型，而不是 QM/MM](#3-为什么用-qm-簇模型而不是-qmmm)
- [4. 计算路线图](#4-计算路线图)
- [5. 文件说明](#5-文件说明)
- [Step 1 获取与处理 PDB 结构](#step-1-获取与处理-pdb-结构)
- [Step 2 截取活性位点簇模型](#step-2-截取活性位点簇模型)
- [Step 3 xTB 预优化](#step-3-xtb-预优化)
- [Step 4 DFT 几何优化](#step-4-dft-几何优化)
- [Step 5 过渡态搜索](#step-5-过渡态搜索)
- [Step 6 IRC 验证](#step-6-irc-验证)
- [Step 7 频率与热校正](#step-7-频率与热校正)
- [Step 8 高精度单点](#step-8-高精度单点)
- [Step 9 能垒与能量剖面](#step-9-能垒与能量剖面)
- [Step 10 模型收敛性检验](#step-10-模型收敛性检验)
- [Step 11 进阶：ONIOM(B3LYP:Amber) 与 CP2K QM/MM 替代方案](#step-11-进阶oniomb3lypamber-与-cp2k-qmmm-替代方案)
- [Step 12 结果分析](#step-12-结果分析)
- [参数选择与收敛建议](#参数选择与收敛建议)
- [常见问题与报错](#常见问题与报错)
- [练习](#练习)
- [参考资料](#参考资料)

---

## 1. 本章目标与适用场景

学完本案例，你应当能够：

- 从 PDB 文件中读出金属配位信息（`LINK` 记录），判断哪些残基必须进入 QM 区；
- 用 PROPKA/PDB2PQR 或 H++ 判断质子化态，并正确处理组氨酸互变异构（HID/HIE）；
- 按“Cα–Cβ 切断 + H 封端 + 冻结截断原子”的规则，手工或用脚本截取簇模型；
- 列出簇模型的电荷/多重度账本（Zn(II) d¹⁰，总电荷 +1，单重态）；
- 写出带冻结原子和连续介质（ε = 4）的 Gaussian 16 与 ORCA 6 优化、扫描、TS、IRC、频率和单点输入；
- 正确解读冻结原子带来的小虚频，用 quasi-RRHO 处理低频模式；
- 设计模型尺寸和介电常数的敏感性检验，判断结论是否可靠；
- 知道什么时候改用 ONIOM 或 QM/MM。

| 适用场景 | 是否适合本流程 | 说明 |
|---|---|---|
| 金属酶、活性位点紧凑、反应局域 | 很适合 | hCA II 就是典型例子 |
| 需要比较几条候选机理的能垒 | 适合 | 簇模型便于快速比较多条路径 |
| 长程静电或构象变化主导反应 | 不太适合 | 改用 QM/MM + 采样 |
| 反应涉及远处质子穿梭（如 His64 质子转移） | 只能部分描述 | 需要更大的模型或 QM/MM，见 [Step 10](#step-10-模型收敛性检验) 与 [Step 11](#step-11-进阶oniomb3lypamber-与-cp2k-qmmm-替代方案) |
| 需要绝对速率常数 | 不适合 | 簇模型主要用于机理判别和相对能垒 |

**前置知识**：Gaussian 输入结构（[Gaussian 第2章](../../Gaussian/tutorials/02-input-anatomy.md)）、几何优化与冻结（[第5章](../../Gaussian/tutorials/05-opt.md)）、频率与热化学（[第6章](../../Gaussian/tutorials/06-freq-thermo.md)）、TS 与 IRC（[第7章](../../Gaussian/tutorials/07-ts-irc.md)）、柔性扫描（[第8章](../../Gaussian/tutorials/08-scan.md)）、溶剂模型（[第9章](../../Gaussian/tutorials/09-solvation.md)）。结构建模的通用方法见 [Common 第3章](../../Common/tutorials/03-structure-modeling.md)。

---

## 2. 背景与反应机理

### 2.1 酶与反应

碳酸酐酶催化 CO2 的可逆水合：

$$
\mathrm{CO_2 + H_2O \rightleftharpoons HCO_3^- + H^+}
$$

hCA II 是研究最充分的同工酶之一。它的活性中心位于一个约 15 Å 深的锥形口袋底部。Zn(II) 由三个组氨酸（His94、His96、His119）和一个水分子/氢氧根配位，构成扭曲的四面体。综述见 Lindskog 1997（[DOI: 10.1016/S0163-7258(96)00198-2](https://doi.org/10.1016/S0163-7258(96)00198-2)）和 Christianson & Fierke 1996（[DOI: 10.1021/ar9501232](https://doi.org/10.1021/ar9501232)）。

下面的几何数据是本教程直接从 PDB 文件中读出的（`LINK` 记录或按坐标计算），不是计算结果：

| 量 | 2CBA（天然酶，1.54 Å） | 3D92（CO2 结合态，1.10 Å） | 来源 |
|---|---|---|---|
| His94 NE2–Zn | 2.10 Å | 1.99 Å | PDB `LINK` 记录 |
| His96 NE2–Zn | 2.12 Å | 2.04 Å | PDB `LINK` 记录 |
| His119 ND1–Zn | 2.11 Å | 2.05 Å | PDB `LINK` 记录 |
| Zn–O(水) | 2.05 Å（HOH 263） | 1.92 Å（HOH 401） | PDB `LINK` 记录 |
| O(Zn-水)···OG1(Thr199) | 2.83 Å | 2.65 Å | 由坐标计算 |
| OG1(Thr199)···OE1(Glu106) | 2.56 Å | 2.59 Å | 由坐标计算 |
| O(Zn-水)···C(CO2) | — | 2.79 Å | 由 3D92 坐标计算 |
| Zn···C(CO2) | — | 3.73 Å | 由 3D92 坐标计算 |

注意配位方式不对称：His94 和 His96 用 **Nε2** 配位，His119 用 **Nδ1** 配位。这直接决定了每个咪唑环上哪个氮带 H（见 [Step 1.4](#14-组氨酸互变异构与金属配位)）。

### 2.2 公认的两步机理

教科书上的机理分为两个半反应（Silverman & Lindskog 1988，[DOI: 10.1021/ar00145a005](https://doi.org/10.1021/ar00145a005)；Silverman & McKenna 2007，[DOI: 10.1021/ar7000588](https://doi.org/10.1021/ar7000588)）：

**(1) 水合半反应**：Zn–OH⁻ 的氧亲核进攻 CO2 的碳，生成 Zn 结合的碳酸氢根；随后碳酸氢根被水置换离开。

```text
E-Zn-OH(-) + CO2  ->  E-Zn-HCO3(-)  -(+H2O)->  E-Zn-H2O + HCO3(-)
```

**(2) 再生半反应**：Zn–H2O 去质子化回到 Zn–OH⁻。质子经活性位点的水链传到 His64，再传到溶液。这一步是高 CO2 浓度下 kcat 的限速步骤（见上述两篇 Acc. Chem. Res.；早期理论研究见 Liang & Lipscomb 1988，[DOI: 10.1021/bi00423a025](https://doi.org/10.1021/bi00423a025)）。

```text
E-Zn-H2O + B  ->  E-Zn-OH(-) + BH(+)      (B = His64 via water wire)
```

### 2.3 Lipscomb 与 Lindskog 两条路线

亲核进攻得到的初始碳酸氢根中，Zn 与“原来属于 OH⁻ 的氧”相连，而这个氧上还带着原来的 H。从这一步到产物碳酸氢根的离去构型，文献讨论过两条路线：

| 路线 | 关键步骤 | 代表文献 |
|---|---|---|
| **Lipscomb 机理** | 碳酸氢根内部发生**质子转移**：原 OH⁻ 上的 H 转到另一个氧上，Zn 仍与同一个氧结合 | Liang & Lipscomb 1987，[DOI: 10.1021/bi00391a012](https://doi.org/10.1021/bi00391a012) |
| **Lindskog 机理** | 碳酸氢根**旋转**，由另一个氧与 Zn 配位（不需要质子转移） | Lindskog 1997 综述，[DOI: 10.1016/S0163-7258(96)00198-2](https://doi.org/10.1016/S0163-7258(96)00198-2) |

早期理论研究对两条路线做过比较，例如 Merz 关于 CO2 结合的研究（[DOI: 10.1021/ja00002a004](https://doi.org/10.1021/ja00002a004)）、Zheng & Merz 关于 hCA II 催化 CO2 水合机理的研究（[DOI: 10.1021/ja00052a054](https://doi.org/10.1021/ja00052a054)）、Bottoni 等人的 DFT 簇模型研究（[DOI: 10.1021/ja030336j](https://doi.org/10.1021/ja030336j)）以及 Miscione 等人的后续工作（[DOI: 10.1007/s00214-007-0274-x](https://doi.org/10.1007/s00214-007-0274-x)）。这些工作所用的模型大小、方法和结论各不相同。**本教程不替你下结论**：请阅读原文，再用自己的计算比较这两条路线。这正是 [练习 3](#练习) 的内容。

### 2.4 一个可以用来对照的实验量级

hCA II 的 CO2 水合 kcat 约为 10⁶ s⁻¹ 量级（**文献值**，见 Silverman & Lindskog 1988，[DOI: 10.1021/ar00145a005](https://doi.org/10.1021/ar00145a005)）。用过渡态理论 $k = (k_BT/h)\,e^{-\Delta G^\ddagger/RT}$ 反推，298 K 下对应的表观活化自由能约为 9 kcal/mol 量级。这只是**由文献 kcat 反推的粗略估计**，而且 kcat 的限速步是质子转移（再生半反应），不是本案例计算的亲核进攻。所以计算得到的亲核进攻能垒**不应当**明显高于这个量级；如果你算出来远高于它（例如高出 10 kcal/mol 以上），应当先怀疑模型或构象，而不是去找“新机理”。

---

## 3. 为什么用 QM 簇模型，而不是 QM/MM

QM 簇模型的思路是：从酶中“挖”出活性位点（几十到几百个原子），在截断处加 H 封端并冻结，其余蛋白环境用一个介电常数较小（通常 ε = 4）的连续介质近似。Himo 和 Siegbahn 等人对这一方法做过系统总结（Siegbahn & Himo 2011，[DOI: 10.1002/wcms.13](https://doi.org/10.1002/wcms.13)；Himo 2017，[DOI: 10.1021/jacs.7b02671](https://doi.org/10.1021/jacs.7b02671)；Himo & de Visser 2022，[DOI: 10.1038/s42004-022-00642-2](https://doi.org/10.1038/s42004-022-00642-2)）。

| 比较项 | QM 簇模型 | QM/MM（ONIOM、CP2K 等） |
|---|---|---|
| 准备工作 | 少：只需活性位点坐标、质子化态 | 多：全蛋白加氢、力场参数（金属位点尤其麻烦）、溶剂化、平衡 MD |
| 金属中心参数 | 不需要 MM 参数 | Zn 位点需要 MCPB.py 等工具参数化（[DOI: 10.1021/acs.jcim.5b00674](https://doi.org/10.1021/acs.jcim.5b00674)） |
| 长程静电 | 只能用连续介质近似 | 显式描述（电子嵌入） |
| 构象采样 | 通常只有一个起始结构 | 可以从 MD 抽多个快照 |
| 计算量 | 50–200 原子的 DFT，单机可做 | QM 区相近，再加 MM 部分和微迭代 |
| 结果可重复性 | 较高（模型确定后结果唯一） | 依赖快照和平衡过程 |
| 典型风险 | 模型太小、冻结引起应力、介电常数敏感 | 边界处理、MM 参数质量、采样不足 |

对 hCA II 来说，簇模型是一个好的教学起点，理由如下：

1. 活性位点很小：一个 Zn、三个配位咪唑、一个 OH⁻，加上 Thr199 一个关键氢键；
2. 亲核进攻是局域反应，主要涉及 O–C 成键；
3. Zn(II) 是 d¹⁰ 闭壳层，没有自旋态问题，适合入门；
4. 有 CO2 结合态晶体结构（3D92）可以作为反应物复合物的起点。

但也要清楚它的局限：His64 质子穿梭距离 Zn 约 7–12 Å（2CBA 中 His64 存在 A/B 两个构象，ND1/NE2 到 Zn 的距离为 7.4–12.1 Å，由坐标计算）。这一步需要包含整条水链，模型会大得多，更适合 QM/MM（例如 Riccardi, Yang & Cui 2010，[DOI: 10.1016/j.bbapap.2009.07.026](https://doi.org/10.1016/j.bbapap.2009.07.026)）。

> [!NOTE]
> 模型多大才够？Kulik 等人在 COMT 的 QM/MM 研究中发现，按距离选取 QM 区时，能垒要到数百原子才收敛（[DOI: 10.1021/acs.jpcb.6b07814](https://doi.org/10.1021/acs.jpcb.6b07814)）。Himo 组则用逐步扩大模型的方法检验簇模型的收敛性（如 Liao, Yu & Himo 2011，[DOI: 10.1021/ct200031t](https://doi.org/10.1021/ct200031t)）。本案例的最小模型只是教学起点，[Step 10](#step-10-模型收敛性检验) 必须做。

---

## 4. 计算路线图

```text
 PDB 2CBA (resting state)            PDB 3D92 (CO2 bound)
        |                                   |
        v                                   |
 Step 1  inspect: LINK records, altlocs, waters, His tautomers
         PROPKA / PDB2PQR (pH 7.8)  ->  protonation table
        |                                   |
        v                                   v
 Step 2  cut_cluster.py --add-co2 --co2-pdb 3D92.pdb
         His -> 4-methylimidazole, Thr199 -> ethanol, H caps,
         freeze C-beta, charge +1 / singlet
        |
        v
 Step 3  GFN2-xTB  --opt  ($fix caps, ALPB ether)      [cheap clean-up]
        |
        v
 Step 4  B3LYP-D3(BJ)/def2-SVP, CPCM eps=4: Opt Freq  ->  Reactant R (0 imag)
        |
        v
 Step 5  relaxed scan O(OH-)...C(CO2) 2.8 -> 1.5 A
         -> highest point -> Opt=(TS,CalcFC,NoEigenTest)  ->  TS1 (1 imag)
        |
        v
 Step 6  IRC both directions  ->  ends optimised  ->  R and P1 (Zn-HCO3-)
        |
        v
 Step 7  Freq (frozen-atom caveats, quasi-RRHO)  ->  G_corr
        |
        v
 Step 8  def2-TZVP single points (eps=4, eps=20 check)  ->  E_sp
        |
        v
 Step 9  dG = dE_sp + dG_corr  ->  table + profile (plot_profile.py)
        |
        v
 Step 10 bigger cluster (+Glu106, +waters), eps 4 vs 20, functional check
        |
        v
 Step 11 (optional) ONIOM(B3LYP:Amber) / CP2K QM/MM
 Step 12 NBO charges, IGMH (Multiwfn)
```

后续机理步骤（Lipscomb 质子转移 TS2 或 Lindskog 旋转 TS2'，再到碳酸氢根被水置换）重复 Step 5–9 即可，见 [练习 3](#练习)。

---

## 5. 文件说明

| 路径 | 内容 |
|---|---|
| [`00-structure/cut_cluster.py`](00-structure/cut_cluster.py) | 读取 2CBA（可选 3D92 的 CO2），按 Step 2 的规则截取簇模型，输出 xyz、Gaussian 坐标块、冻结行、原子索引和构建报告。只用 Python 标准库 |
| [`00-structure/fill_template.py`](00-structure/fill_template.py) | 把坐标、电荷多重度和冻结行填入下面各模板的 `@...@` 占位符 |
| [`01-xtb-preopt/xcontrol_fix.inp`](01-xtb-preopt/xcontrol_fix.inp) | xTB `$fix` 冻结截断原子 |
| [`01-xtb-preopt/xcontrol_scan.inp`](01-xtb-preopt/xcontrol_scan.inp) | xTB 弛豫扫描 O···C（预览用） |
| [`01-xtb-preopt/run_xtb.sh`](01-xtb-preopt/run_xtb.sh) | Slurm/本地脚本：静息态与反应物预优化，可选扫描 |
| [`02-opt/g16_opt_reactant.gjf`](02-opt/g16_opt_reactant.gjf)、[`g16_opt_product.gjf`](02-opt/g16_opt_product.gjf) | Gaussian 16：B3LYP-D3(BJ)/def2-SVP，CPCM ε=4，ModRedundant 冻结，Opt Freq |
| [`02-opt/orca_opt_reactant.inp`](02-opt/orca_opt_reactant.inp) | ORCA 6 等价输入（RIJCOSX，`%cpcm epsilon 4.0`，笛卡尔约束） |
| [`03-ts/g16_scan_OC.gjf`](03-ts/g16_scan_OC.gjf)、[`orca_scan_OC.inp`](03-ts/orca_scan_OC.inp) | O···C 弛豫扫描 |
| [`03-ts/g16_ts.gjf`](03-ts/g16_ts.gjf)、[`orca_ts.inp`](03-ts/orca_ts.inp) | TS 优化 + 频率 |
| [`04-irc/g16_irc.gjf`](04-irc/g16_irc.gjf)、[`orca_irc.inp`](04-irc/orca_irc.inp) | IRC |
| [`05-freq/g16_freq_ts.gjf`](05-freq/g16_freq_ts.gjf)、[`orca_freq_ts.inp`](05-freq/orca_freq_ts.inp) | 单独的频率作业 |
| [`06-sp/`](06-sp/g16_sp_tzvp.gjf) | def2-TZVP 单点：CPCM ε=4、ε=20 对照、ORCA SMD 可选 |
| [`07-oniom/`](07-oniom/README.md) | ONIOM(B3LYP:Amber) 骨架模板与说明 |
| [`08-profile/energies_template.csv`](08-profile/energies_template.csv) | **全是占位行（NaN, placeholder=yes）**，由你填入自己的结果 |
| [`08-profile/plot_profile.py`](08-profile/plot_profile.py) | 读 CSV 生成相对自由能表和能量剖面图；遇到占位行会拒绝运行 |
| [`jobs/sub_g16.slurm`](jobs/sub_g16.slurm)、[`jobs/sub_orca.slurm`](jobs/sub_orca.slurm) | Slurm 提交脚本 |

**PDB 文件不随仓库分发**（2CBA.pdb 约 0.2 MB、3D92.pdb 约 0.4 MB，下载命令见 Step 1）。所有输入文件均为 ASCII、LF 换行。

---

## Step 1 获取与处理 PDB 结构

### 1.1 下载

```bash
cd Cases/01-enzyme-carbonic-anhydrase-ii-qm-cluster/00-structure
curl -L -o 2CBA.pdb https://files.rcsb.org/download/2CBA.pdb
curl -L -o 3D92.pdb https://files.rcsb.org/download/3D92.pdb
```

两个条目的基本信息（来自文件头 `TITLE/EXPDTA/REMARK 2/JRNL`，2026-11 核对）：

| PDB | 标题 | 方法 / 分辨率 | 原始文献 |
|---|---|---|---|
| [2CBA](https://www.rcsb.org/structure/2CBA) | Structure of native and apo carbonic anhydrase II and some of its anion-ligand complexes | X 射线，1.54 Å | Håkansson et al., *J. Mol. Biol.* 1992, 227, 1192，[DOI: 10.1016/0022-2836(92)90531-N](https://doi.org/10.1016/0022-2836(92)90531-N) |
| [3D92](https://www.rcsb.org/structure/3D92) | Human carbonic anhydrase II bound with substrate carbon dioxide | X 射线，1.10 Å，pH 7.8，100 K | Domsic et al., *J. Biol. Chem.* 2008, 283, 30766，[DOI: 10.1074/jbc.M805353200](https://doi.org/10.1074/jbc.M805353200) |

为什么用 2CBA 做主体、3D92 只取 CO2？2CBA 是经典的天然酶静息态结构，活性位点没有其他配体，适合演示“从零开始搭模型”。3D92 分辨率更高，并且直接给出了 CO2 的结合位置。两者晶胞相同（P2₁，a = 42.7、b = 41.7、c = 73.0 Å，β = 104.6°），按本教程的检查，活性位点重原子（三个 His 环、Thr199 OG1、Glu106 OE1）在两个文件间的位移约 0.1–0.3 Å（按 Zn 重合后比较坐标）。因此把 3D92 的 CO2 平移到 2CBA 的 Zn 上作为初猜是合理的。你也可以完全用 3D92 搭模型（[练习 1](#练习)）。

### 1.2 检查文件

先看 Zn 的配位和异常记录：

```bash
grep -E "^LINK" 2CBA.pdb              # Zn 与 His94/96/119 和 HOH 263 的配位
grep -E "^HETATM" 2CBA.pdb | grep -v HOH
grep -E "^ATOM.{12}[AB]" 2CBA.pdb | cut -c13-27 | sort -u   # 有替代构象 (altloc) 的原子
grep -E "REMARK 465|REMARK 470" 2CBA.pdb | head               # 缺失残基/原子
```

本教程检查到的情况（你应当能复现）：

- `LINK` 记录给出 His94 NE2、His96 NE2、His119 ND1 与 Zn 配位，以及 Zn–O(HOH 263)；
- 2CBA 中除 Zn 外没有其他 HETATM 配体，共有 220 个结晶水；
- **His64 有 A/B 两个替代构象**（占有率 0.70/0.20 等），这正是著名的 His64 “in/out” 构象。最小簇模型不包含 His64，但做大模型或 QM/MM 时必须选定一个构象；
- 3D92 中除 CO2 A 301（活性位点）外，还有第二个 CO2 A 302（远离 Zn）和两个甘油（GOL 303/304）。GOL 303 的最近原子离 Zn 约 4.6 Å（由坐标计算），搭大模型时要决定删除还是保留。

### 1.3 去除不需要的分子，判断质子化态

簇模型只保留活性位点，所以“删水”在 Step 2 由脚本完成。但质子化态必须在全蛋白中判断，因为 pKa 取决于整个静电环境。常用工具：

| 工具 | 用法 | 说明 |
|---|---|---|
| PROPKA 3 | `propka3 2CBA.pdb` → `2CBA.pka` | 经验 pKa 预测（[DOI: 10.1021/ct100578z](https://doi.org/10.1021/ct100578z)） |
| PDB2PQR | 见下方命令 | 加氢、按 PROPKA 结果指定质子化态、优化氢键网络（[DOI: 10.1002/pro.3280](https://doi.org/10.1002/pro.3280)） |
| H++ 网页服务器 | <http://newbiophysics.cs.vt.edu/H++/> | 基于连续介质静电的 pKa 与加氢（[DOI: 10.1093/nar/gks375](https://doi.org/10.1093/nar/gks375)） |
| PDBFixer | `pdbfixer 2CBA.pdb --add-atoms=all --keep-heterogens=all --ph=7.8 --output=2CBA_fixed.pdb` | 补缺失原子，按 pH 粗略加氢 |

PDB2PQR 的命令（选项名取自 PDB2PQR 源码 `pdb2pqr/main.py`，请用 `pdb2pqr --help` 核对你的版本）：

```bash
pdb2pqr --ff=AMBER --ffout=AMBER --titration-state-method=propka --with-ph=7.8 \
        --keep-chain --pdb-output 2CBA_H.pdb 2CBA.pdb 2CBA.pqr
```

| 选项 | 作用 |
|---|---|
| `--ff=AMBER` | 用 AMBER 参数分配电荷/半径（可选 AMBER、CHARMM、PARSE 等） |
| `--ffout=AMBER` | 输出用 AMBER 的残基命名，可直接看到 HID/HIE/HIP |
| `--titration-state-method=propka` | 调用 PROPKA 计算可滴定残基的 pKa |
| `--with-ph=7.8` | 按该 pH 指定质子化态（这里与 3D92 的结晶 pH 一致；生理条件可用 7.4） |
| `--pdb-output` | 另外输出一个带 H 的 PDB |

> [!WARNING]
> PDB2PQR 和 PROPKA 都**不擅长处理金属位点**：它们会把 Zn 当作普通离子或直接忽略，对 Zn 配位的 His 给出的质子化结果可能不合理，Zn 结合水也会被当成普通水。对金属配体，要以晶体结构的配位几何为准手工指定（见 1.4）。PROPKA 的输出只用来判断**远离 Zn 的**可滴定残基，比如 Glu106、Glu117、His64、Lys/Arg 等。

本案例需要做出的质子化决定：

| 残基/分子 | 决定 | 依据 |
|---|---|---|
| Zn 结合水（HOH 263） | **OH⁻**（去质子化） | 研究的是活性形式 Zn–OH⁻；Zn 结合水的 pKa 约为 7，见 Silverman & Lindskog 1988 综述 |
| His94、His96 | 中性，**Nδ1–H**（HID） | Nε2 与 Zn 配位，H 只能在 Nδ1 |
| His119 | 中性，**Nε2–H**（HIE） | Nδ1 与 Zn 配位，H 只能在 Nε2 |
| Thr199 | 中性 OH | 作为 Zn–OH⁻ 的氢键受体，同时向 Glu106 给出氢键 |
| Glu106 | 去质子化（–1） | 埋藏但与 Thr199 形成氢键网络；用 PROPKA 结果确认 |
| Glu117 | 去质子化（–1） | 与 His119 Nε2–H 形成氢键（第二配位层；2CBA 中 NE2···OE2 = 2.59 Å，由坐标计算） |
| His64 | 视 pH 与构象而定 | 质子穿梭残基，最小模型不含，做大模型时用 PROPKA 判断并与文献对照 |

### 1.4 组氨酸互变异构与金属配位

组氨酸咪唑环有两个氮：Nδ1（ND1，靠近 Cβ）和 Nε2（NE2，远端）。中性 His 只有一个氮带 H：

| Amber 名称 | 质子位置 | 电荷 |
|---|---|---|
| HID | Nδ1–H | 0 |
| HIE | Nε2–H | 0 |
| HIP | 两个 N 都带 H | +1 |

**规则**：与金属配位的那个氮不能带 H，H 必须在另一个氮上。所以 His94/His96 是 HID，His119 是 HIE。这两个 N–H 又分别与第二配位层形成氢键（2CBA 中由坐标计算：His94 ND1···OE1 Gln92 2.69 Å、His96 ND1···O Asn244 2.73 Å、His119 NE2···OE2 Glu117 2.59 Å）。这些“间接配体”会影响 Zn 的 Lewis 酸性，在大模型中可以考虑加入（[Step 10](#step-10-模型收敛性检验)）。

> [!TIP]
> 不要把 X 射线结构中 His 环的“翻转”当作定论：X 射线难以区分 C 和 N。对远离金属的 His，PDB2PQR/Reduce 会根据氢键网络决定是否翻转环；对金属配体，配位键本身就确定了环的朝向。

---

## Step 2 截取活性位点簇模型

### 2.1 选哪些残基，为什么

| 组分 | 截断方式 | 保留原子 | 理由 |
|---|---|---|---|
| Zn(II) | — | Zn | 催化中心 |
| Zn 结合水 → OH⁻ | — | O + 1 H | 亲核试剂 |
| His94、His96、His119 | Cα–Cβ 之间切断 → **4-甲基咪唑** | 环 5 原子 + Cβ（作为甲基 C）+ 3 个甲基 H + 2 个环 C–H + 1 个 N–H | 第一配位层，决定 Zn 的电子结构 |
| Thr199 | Cα–Cβ 之间切断 → 乙醇 | CB、OG1、CG2 + 6 H | 与 Zn–OH⁻ 形成氢键，定位 OH⁻ 的取向（“门卫”残基） |
| Glu106（可选） | Cβ–Cγ 之间切断 → 乙酸根 | CG、CD、OE1、OE2 + 3 H | 接受 Thr199 的氢键，影响氢键网络和总电荷 |
| 第二层水（可选，默认 1 个） | — | O + 2 H | 稳定 OH⁻，后续可参与质子转移 |
| CO2 | — | 3 原子 | 底物，取自 3D92 |

疏水口袋残基（Val121、Val143、Leu198、Trp209 都在 3D92 中 CO2 周围 4 Å 以内，由坐标计算）在最小模型中省略。它们主要提供范德华作用和空间约束，可以在 [Step 10](#step-10-模型收敛性检验) 的大模型中加入。

### 2.2 截断与 H 封端

截断时只切 **C–C 单键**，不切极性键、不切共轭体系，也不切肽键。例如 His 在 Cα–Cβ 处切断（Cβ 变成 4-甲基咪唑的甲基），而不是在 Cβ–Cγ（Cβ 与咪唑环之间）切断，这样咪唑环的电子结构保持不变。具体做法：

1. 保留被切键靠近活性中心一侧的原子（这里是 Cβ），**保留其晶体坐标**；
2. 删掉另一侧的原子（Cα 及以外）；
3. 在 Cβ 上补 H，使它成为甲基。其中一个 H 沿原 Cβ–Cα 键方向放置，C–H 取 1.09 Å，这个 H 就是“封端 H（link H）”；
4. 冻结 Cβ（截断原子），模拟蛋白骨架对侧链位置的约束。

脚本 [`cut_cluster.py`](00-structure/cut_cluster.py) 按这个规则自动完成，并把环上的 C–H、N–H 和 OH 上的 H 都按几何规则放好。**这些 H 都只是初猜**，必须经过 Step 3/4 优化。

```bash
cd 00-structure
python cut_cluster.py -h
# 静息态模型（无 CO2）：Zn-OH- + 3 His + Thr199 + 1 个第二层水
python cut_cluster.py --pdb 2CBA.pdb --outdir build
# 反应物复合物：再加入 3D92 的 CO2
python cut_cluster.py --pdb 2CBA.pdb --outdir build_co2 --add-co2 --co2-pdb 3D92.pdb
# 稍大的模型：加 Glu106（总电荷变为 0）
python cut_cluster.py --pdb 2CBA.pdb --outdir build_glu --with-glu106 --add-co2 --co2-pdb 3D92.pdb
```

主要选项：

| 选项 | 默认 | 说明 |
|---|---|---|
| `--pdb` / `--download` | — | 本地 PDB 文件，或让脚本从 RCSB 下载 2CBA |
| `--add-co2` | 关 | 加入 CO2 |
| `--co2-pdb 3D92.pdb` | 无 | 从 3D92 复制 CO2 A 301（平移使两个 Zn 重合）；不给时按几何规则放置 |
| `--with-glu106` | 关 | 加入 Glu106 乙酸根 |
| `--waters N` | 1 | 最多保留离 OH⁻ 最近的 N 个第二层水 |
| `--water-cutoff R` | 3.4 | 只保留 O···O(H⁻) 距离不超过 R Å 的水；不满足时即使没到 N 个也停止 |
| `--oh-mode thr/anti` | thr | OH⁻ 的 H 朝向 Thr199 OG1，或沿 Zn–O 方向向外 |

> [!NOTE]
> 加入 CO2 时，2CBA 中离 Zn–OH 最近的那个水（HOH 338，O···O 2.60 Å）与 3D92 中 CO2 的位置几乎重合（不到 0.5 Å，由坐标计算）。它就是文献中所说的“deep water”，CO2 结合时会被取代。脚本会自动删除与 CO2 距离小于 2.2 Å 的水，并在 `build_report.txt` 中注明。所以无 CO2 和有 CO2 的模型中，“第二层水”**不是同一个水分子**。比较能量时要注意这一点，最好统一用有 CO2 的模型做反应物。

### 2.3 输出与原子编号

默认选项下（1 个水，无 Glu106），原子顺序是固定的。后面所有输入模板都按这个编号书写：

| 原子编号（1-based） | 内容 |
|---|---|
| 1 | Zn |
| 2、3 | OH⁻ 的 O、H |
| 4–15 | His94（环 4–8，**Cβ 9**，甲基 H 10–12，环 C–H 13–14，Nδ1–H 15） |
| 16–27 | His96（**Cβ 21**） |
| 28–39 | His119（**Cβ 33**，H 在 Nε2 上） |
| 40–48 | Thr199（**CB 40**，OG1 41，CG2 42，H 43–48） |
| 49–51 | 第二层水（**O 49**） |
| 52–54 | CO2（**C 52**，O 53、54），仅 `--add-co2` |

粗体为冻结原子。原子数：静息态 51，反应物复合物 54；加 Glu106 后分别为 58 和 61（Glu106 占 49–55，水和 CO2 依次后移，CO2 的 C 变为 59）。**一定要以 `build_report.txt` 为准**，换了选项后编号会变。示意输出如下（格式示意，坐标与距离以你自己运行的结果为准）：

```text
（示意输出，数值仅供参考）
hCA II QM cluster build report
source PDB      : 2CBA.pdb
atoms           : 54
charge / mult   : 1 / 1
electrons       : 228 (even, closed shell OK)
frozen (1-based): 9, 21, 33, 40, 49

Zn ligands from LINK records: His94 NE2, His96 NE2, His119 ND1, HOH 263 O
...
REACTION COORDINATE: O(OH-) = atom 2, C(CO2) = atom 52 (1-based); ORCA 0-based: 1 and 51; ...
```

生成的文件：

| 文件 | 用途 |
|---|---|
| `cluster.xyz` | 给 xTB、ORCA、可视化软件 |
| `cluster_coords.txt` | 电荷多重度行 + Gaussian 坐标块 |
| `cluster_freeze.txt` | `X 9 F` 等 ModRedundant 冻结行 |
| `cluster_constraints_orca.txt` | `{ C 8 C }` 等 ORCA 约束行（**0-based**） |
| `cluster_fix_xtb.txt` | xTB `$fix` 块 |
| `atoms_index.csv` | 每个原子的残基、原子名、说明 |
| `build_report.txt` | 构建报告：电荷、电子数、冻结列表、反应坐标原子 |

### 2.4 冻结哪些原子

| 冻结原子 | 理由 |
|---|---|
| 三个 His 的 Cβ（9、21、33） | 模拟蛋白骨架对配体位置的约束，防止咪唑环在优化中漂移或聚到一起 |
| Thr199 的 CB（40） | 同上 |
| 第二层水的 O（49） | 防止这个水在小模型中“跑掉”或游走到不真实的位置；H 不冻结，以便调整氢键 |
| （大模型）Glu106 的 CG | 同截断原子 |

冻结的代价是：频率计算中会出现少量与冻结原子相关的小虚频（[Step 7](#step-7-频率与热校正)），熵也带有一定误差。Himo 等人的经验是，只要冻结原子离反应中心足够远，这种误差对相对能垒的影响通常较小，但仍需检验（[DOI: 10.1021/jacs.7b02671](https://doi.org/10.1021/jacs.7b02671)）。

> [!TIP]
> 想更系统地选择残基，可以参考 RINRUS（Residue Interaction Network ResidUe Selector，[GitHub: natedey/RINRUS](https://github.com/natedey/RINRUS)，综述 [DOI: 10.1002/wcms.70078](https://doi.org/10.1002/wcms.70078)）。它根据残基与底物的相互作用网络给残基排序，并生成一系列逐渐增大的模型，适合做模型收敛性检验。

### 2.5 电荷与多重度账本

| 组分 | 形式电荷 | 说明 |
|---|---|---|
| Zn(II) | +2 | d¹⁰ |
| OH⁻ | −1 | 活性形式 |
| 3 × 4-甲基咪唑（中性 His） | 0 | HID/HID/HIE |
| Thr199（乙醇） | 0 | |
| 第二层水 | 0 | |
| CO2 | 0 | |
| **合计（默认模型）** | **+1** | |
| Glu106（乙酸根，可选） | −1 | 加入后合计 **0** |

**多重度**：Zn(II) 为 d¹⁰，所有配体都是闭壳层，体系电子数为偶数（反应物复合物共 228 个电子），所以是**单重态（multiplicity = 1）**，用限制性 RKS 即可。脚本会在报告中检查电子数的奇偶，奇数说明模型搭错了（比如漏了一个 H）。

> [!WARNING]
> 最常见的错误是忘了把 Zn 结合水变成 OH⁻，或者忘了给某个 His 加 N–H。这类错误下，电荷 +1 与电子数的奇偶对不上，Gaussian 会报 `The combination of multiplicity 1 and N electrons is impossible`。见 [常见问题](#常见问题与报错)。

---

## Step 3 xTB 预优化

用几何规则放置的 H 可能有不合理的接触，直接上 DFT 浪费时间。先用 GFN2-xTB（[DOI: 10.1021/acs.jctc.8b01176](https://doi.org/10.1021/acs.jctc.8b01176)）在冻结截断原子的条件下优化一遍。

### 3.1 xcontrol 文件

[`01-xtb-preopt/xcontrol_fix.inp`](01-xtb-preopt/xcontrol_fix.inp)：

```text
$fix
   atoms: 9,21,33,40,49
$end
```

`$fix` 把列出的原子**精确固定**（梯度中对应分量被移除），编号从 1 开始，对应 xyz 文件中的顺序。xTB 还提供 `$constrain`（谐振势约束，可用 `force constant=` 调整强度，原子位置会有少量偏离）。冻结截断原子时用 `$fix` 更直接，与 Gaussian/ORCA 中的笛卡尔冻结一致（[xTB 文档：Detailed Input](https://xtb-docs.readthedocs.io/en/latest/xcontrol.html)）。

### 3.2 运行

```bash
cd 01-xtb-preopt
sbatch run_xtb.sh          # 或者在工作站上：bash run_xtb.sh
```

脚本中的关键命令：

```bash
xtb cluster.xyz --gfn 2 --chrg 1 --uhf 0 --alpb ether --opt tight --input xcontrol.inp > xtb.out
```

| 选项 | 值 | 说明 |
|---|---|---|
| `--gfn 2` | GFN2-xTB | 默认方法，显式写出便于记录 |
| `--chrg 1` | +1 | 加 Glu106 的模型改为 0 |
| `--uhf 0` | 0 | 未成对电子数；Zn(II) d¹⁰ 为 0 |
| `--alpb ether` | 乙醚 | ALPB 隐式溶剂；乙醚 ε 约为 4，粗略模拟蛋白环境。也可用 `water` 做对照 |
| `--opt tight` | tight | 收敛级别（crude、sloppy、loose、lax、normal、tight、vtight、extreme） |
| `--input` | xcontrol.inp | 读取 `$fix` 等详细输入 |

### 3.3 检查

```bash
grep "GEOMETRY OPTIMIZATION CONVERGED" xtb_reactant/xtb.out
grep "TOTAL ENERGY" xtb_reactant/xtb.out | tail -1
ls xtb_reactant/xtbopt.xyz          # 优化后的结构
ls xtb_reactant/NOT_CONVERGED 2>/dev/null && echo "not converged!"
```

成功标准：

- 输出中有 `GEOMETRY OPTIMIZATION CONVERGED`，目录中**没有** `NOT_CONVERGED` 文件；
- 用 VMD/Avogadro/GaussView 打开 `xtbopt.xyz` 目视检查：Zn 仍为四配位（3 N + O），Zn–N 和 Zn–O 距离没有明显拉长（例如超过 2.5 Å）；
- OH⁻ 的 H 仍朝向 Thr199 OG1，Thr199 的 OH 指向 Glu106 方向；
- CO2 没有在 xTB 优化中直接与 OH⁻ 成键（如果成键了，说明起始距离太近或 xTB 势能面不可靠，可在 `$constrain` 中加 `distance: 2, 52, 2.8` 先固定 O···C 距离）。

> [!NOTE]
> xTB 只用于“清理”结构，不用于报告能量。可选的 xTB 弛豫扫描（[`xcontrol_scan.inp`](01-xtb-preopt/xcontrol_scan.inp)，输出 `xtbscan.log`）可以帮你预览 O···C 方向的能量变化趋势，但 xTB 的能垒不可信，最高点只用作参考。需要系统采样构象（例如 OH⁻ 和水的取向）时可以用 CREST（[DOI: 10.1063/5.0197592](https://doi.org/10.1063/5.0197592)），在冻结原子的条件下搜索低能构象。

---

## Step 4 DFT 几何优化

### 4.1 理论水平的选择

| 设置 | 选择 | 理由 |
|---|---|---|
| 泛函 | B3LYP | 簇模型研究中用得最多的泛函，便于与文献比较（[DOI: 10.1063/1.464913](https://doi.org/10.1063/1.464913)、[DOI: 10.1021/j100096a001](https://doi.org/10.1021/j100096a001)） |
| 色散 | D3(BJ) | 簇模型中的氢键、堆积和范德华接触都需要色散校正（[DOI: 10.1002/jcc.21759](https://doi.org/10.1002/jcc.21759)） |
| 基组（优化/频率） | def2-SVP | 速度与精度的折中；Zn 是第四周期元素，def2 系列对 Zn 是全电子基组，**不需要 ECP**（def2 从 Rb 开始才带 ECP，见 [Gaussian 第3章](../../Gaussian/tutorials/03-method-basis.md)） |
| 基组（单点） | def2-TZVP | Step 8 |
| 环境 | CPCM，ε = 4 | 簇模型的标准做法，用低介电常数近似蛋白环境（[DOI: 10.1021/jp9716997](https://doi.org/10.1021/jp9716997)） |
| 冻结 | Cβ 等截断原子 | Step 2.4 |

> [!NOTE]
> 为什么在优化中就加连续介质？对带电的簇（这里总电荷 +1），气相优化有时会得到过分紧凑或电荷过度转移的结构。在 ε = 4 中优化更接近蛋白内环境。也有很多研究在气相中优化、只在单点中加溶剂；两种做法在小模型中差别通常不大，但要在文章中写清楚。

### 4.2 生成输入

```bash
cd 02-opt
python ../00-structure/fill_template.py g16_opt_reactant.gjf opt_reactant.gjf \
    --xyz ../01-xtb-preopt/xtb_reactant/xtbopt.xyz --charge 1 --mult 1 \
    --frozen-file ../00-structure/build_co2/cluster_freeze.txt
```

### 4.3 Gaussian 16 输入

模板 [`02-opt/g16_opt_reactant.gjf`](02-opt/g16_opt_reactant.gjf)：

```text
%nprocshared=16
%mem=24GB
%chk=opt_reactant.chk
#p B3LYP/Def2SVP EmpiricalDispersion=GD3BJ SCRF=(CPCM,Solvent=Water,Read) Opt=(ModRedundant,MaxCycles=200) Freq

hCA II cluster, reactant complex Zn-OH- + CO2, B3LYP-D3(BJ)/def2-SVP, CPCM eps=4, cap atoms frozen

@CHARGE_MULT@
@COORDS@

@FREEZE@

Eps=4.0

```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `B3LYP/Def2SVP` | — | Gaussian 中 def2 基组关键词**不带连字符** |
| `EmpiricalDispersion=GD3BJ` | — | D3(BJ) 色散校正 |
| `SCRF=(CPCM,Solvent=Water,Read)` | — | 以水为基础溶剂，`Read` 读取后面的 PCM 输入段 |
| PCM 输入段 `Eps=4.0` | 4.0 | 把静态介电常数改为 4；其余参数（包括 `EpsInf`）沿用水的值，详见 [Gaussian 第9章 §6.3](../../Gaussian/tutorials/09-solvation.md) |
| `Opt=ModRedundant` | — | 读取坐标后的冻结段 |
| `X 9 F` 等 | — | 冻结 9 号原子的笛卡尔坐标（1 个原子 → 笛卡尔坐标，`F` = freeze，见 [Gaussian 第5章 §6.2](../../Gaussian/tutorials/05-opt.md)） |
| `MaxCycles=200` | 150–300 | 簇模型的软自由度多，默认步数常常不够 |
| `Freq` | — | 优化后立即算频率，确认极小点并得到热校正 |

输入各段的顺序很重要（[Gaussian Input 文档 “Section Ordering”](https://gaussian.com/input/)）：分子坐标 → 空行 → ModRedundant 段 → 空行 → PCM 输入段 → 空行。文件末尾必须有空行。

另一种冻结写法是在坐标行中加冻结码 `-1`（例如 `C  -1  x y z`），效果相同，见 [Gaussian 第5章](../../Gaussian/tutorials/05-opt.md)。两种写法不要混用。

### 4.4 ORCA 6 替代输入

模板 [`02-opt/orca_opt_reactant.inp`](02-opt/orca_opt_reactant.inp)：

```text
! B3LYP D3BJ def2-SVP def2/J RIJCOSX CPCM(water) TightSCF Opt Freq

%pal nprocs 16 end
%maxcore 3000

%cpcm
  epsilon 4.0          # overrides the dielectric constant of water
end

%geom
  MaxIter 200
  Constraints
@ORCA_FIX@
  end
end

%freq
  QuasiRRHO true       # ORCA 6: Grimme quasi-RRHO entropy in the thermochemistry
end

* xyzfile 1 1 reactant.xyz
```

| 关键词 | 说明 |
|---|---|
| `D3BJ` | D3(BJ) 色散 |
| `def2/J RIJCOSX` | RI-J + COSX 近似加速杂化泛函，几乎不损失精度 |
| `CPCM(water)` + `%cpcm epsilon 4.0` | 先指定溶剂，再在 `%cpcm` 块中覆盖介电常数（ORCA 手册 “Solvation” 一节列出了 `epsilon`、`refrac` 等可调参数） |
| `{ C 8 C }` | 笛卡尔约束（固定第 8 号原子）。**ORCA 原子编号从 0 开始**，所以 Gaussian 的 9 号对应 ORCA 的 8 号；`fill_template.py` 会自动减 1 |
| `%freq QuasiRRHO true` | ORCA 6 手册中的 quasi-RRHO 开关（参考频率由 `QRRHORefFreq` 设置） |
| `* xyzfile 1 1 reactant.xyz` | 电荷 1、多重度 1，坐标从外部 xyz 读取 |

### 4.5 运行与检查

```bash
sbatch ../jobs/sub_g16.slurm opt_reactant.gjf
# ORCA:  sbatch ../jobs/sub_orca.slurm orca_opt_reactant.inp
```

```bash
grep -c "Normal termination" opt_reactant.log          # 应为 2（Opt + Freq）
grep "Stationary point found" opt_reactant.log
grep "imaginary frequencies" opt_reactant.log          # 理想情况：无输出
grep -A4 "Converged?" opt_reactant.log | tail -5
# ORCA
grep "THE OPTIMIZATION HAS CONVERGED" orca_opt_reactant.out
grep "ORCA TERMINATED NORMALLY" orca_opt_reactant.out
```

成功标准：

- Gaussian 出现 `Stationary point found`，频率部分没有虚频；或者只有与冻结原子相关、绝对值很小的虚频（见 [Step 7](#step-7-频率与热校正)）；
- Zn 配位保持四面体，CO2 仍是线形（O=C=O 接近 180°），没有与 OH⁻ 成键；
- 冻结原子的坐标与输入完全一致（可用 GaussView 测量，或者比较首尾两个结构）。

产物（Zn–HCO3⁻）的优化用 [`g16_opt_product.gjf`](02-opt/g16_opt_product.gjf)，起始结构取 Step 6 IRC 的产物端点。

---

## Step 5 过渡态搜索

### 5.1 反应坐标

亲核进攻的主要反应坐标是 O(OH⁻)···C(CO2) 距离：在反应物中约 2.8 Å（3D92 晶体结构中为 2.79 Å），在碳酸氢根中约为 C–O 单键长度。默认模型中 O = 2 号原子、C = 52 号原子（ORCA 中为 1 和 51）。

### 5.2 弛豫扫描

[`03-ts/g16_scan_OC.gjf`](03-ts/g16_scan_OC.gjf) 从优化好的反应物出发：

```text
%nprocshared=16
%mem=24GB
%oldchk=opt_reactant.chk
%chk=scan_OC.chk
#p B3LYP/Def2SVP EmpiricalDispersion=GD3BJ SCRF=(CPCM,Solvent=Water,Read) Opt=(ModRedundant,MaxCycles=100) Geom=Check Guess=Read

hCA II cluster, relaxed scan of O(OH-)...C(CO2) from the optimised reactant, 13 steps of -0.10 A

1 1

@FREEZE@
B 2 52 S 13 -0.10

Eps=4.0

```

`B 2 52 S 13 -0.10` 表示从当前值开始，把 2–52 距离每步缩短 0.10 Å，共 13 步，每一点都做受限优化（[Gaussian 第8章](../../Gaussian/tutorials/08-scan.md)）。冻结行仍然要写，否则截断原子会在扫描中被放开。

ORCA 写法（[`03-ts/orca_scan_OC.inp`](03-ts/orca_scan_OC.inp)，0-based 编号）：

```text
%geom
  Scan
    B 1 51 = 2.80, 1.50, 14
  end
  Constraints
@ORCA_FIX@
  end
end
```

ORCA 会把每个扫描点的结构写入 `orca_scan_OC.allxyz`，可以从中取出能量最高的结构作为 TS 初猜。ORCA 模板之间按默认输出文件名衔接：`orca_opt_reactant.inp` 的最终结构写入 `orca_opt_reactant.xyz`（扫描和单点读取它），`orca_ts.inp` 写出 `orca_ts.xyz` 与 `orca_ts.hess`（IRC 和单独频率作业读取它们）。

> [!NOTE]
> 从 chk 读取结构的 Gaussian 模板（`g16_scan_OC.gjf`、`g16_sp_tzvp*.gjf`）中直接写了电荷多重度行 `1 1`，对应默认模型。用 `--with-glu106` 构建的模型要改成 `0 1`；ORCA 模板的 `* xyzfile 1 1 ...` 同理。

提取扫描曲线：

```bash
python ../../../Gaussian/scripts/scan_irc_profile.py scan_OC.log --mode scan -o scan_OC.png
```

（脚本说明见 [Gaussian/scripts](../../Gaussian/scripts/README.md)，用 `-h` 查看所有选项。）

（示意输出，数值仅供参考）你应该看到一条先上升、在某个 O···C 距离出现极大、然后下降的曲线。极大点所在的距离和相对能量要以你自己的计算为准。如果曲线一路单调下降，说明在这个模型中亲核进攻几乎无垒，或者 OH⁻ 的朝向不利于反应（例如 OH⁻ 的 H 正对着 CO2），需要回头检查反应物结构。

### 5.3 TS 优化

把扫描极大点的结构填入 [`03-ts/g16_ts.gjf`](03-ts/g16_ts.gjf)：

```text
#p B3LYP/Def2SVP EmpiricalDispersion=GD3BJ SCRF=(CPCM,Solvent=Water,Read) Opt=(TS,CalcFC,NoEigenTest,ModRedundant,MaxCycles=200) Freq
```

| 选项 | 说明 |
|---|---|
| `TS` | 优化到一阶鞍点 |
| `CalcFC` | 第一步计算精确 Hessian；对 54 原子的体系可以承受 |
| `NoEigenTest` | 不在 Hessian 负本征值个数不为 1 时中止。簇模型中冻结原子常带来额外的小负本征值，这个选项几乎是必需的 |
| `ModRedundant` | 继续冻结截断原子 |
| `Freq` | 优化后立刻验证 |

推荐把扫描极大点的结构导出成 xyz（例如用 GaussView 打开 `scan_OC.log` 后另存，或用仓库脚本 `scan_irc_profile.py --write-max`），再用 `fill_template.py --xyz` 填入模板，这样能清楚地知道自己用的是哪一个点。注意 [Geom 关键词文档](https://gaussian.com/geom/) 中的 `Step=N` 只适用于失败或未完成的优化，不要用它从已完成的扫描中取点。

ORCA（[`03-ts/orca_ts.inp`](03-ts/orca_ts.inp)）用 `! OptTS Freq`，并在 `%geom` 中设置 `Calc_Hess true`。想更便宜，可以用 `Hybrid_Hess { ... } end` 只对反应中心原子精确计算 Hessian（ORCA 手册 “Transition States” 一节）。ORCA 也提供 `! ScanTS`，一步完成扫描和 TS 优化。

### 5.4 检查：唯一虚频

```bash
grep "imaginary frequencies" ts_attack.log
grep -A12 "Harmonic frequencies" ts_attack.log | head -14
```

成功标准：

- 恰好 **1 个明显的虚频**，振动模式是 O(OH⁻) 向 C(CO2) 靠近、同时 CO2 弯曲（在 GaussView/Chemcraft 中播放动画确认）；
- 如果另外还有几个绝对值很小（例如几十 cm⁻¹ 以内）的虚频，并且动画显示它们集中在冻结的 Cβ/甲基附近，通常是冻结造成的（[Step 7](#step-7-频率与热校正)）；
- 如果第二个虚频涉及反应中心，说明还不是正确的 TS，需要沿该模式微扰结构后重新优化。

（示意输出，数值仅供参考）

```text
 ******    1 imaginary frequencies (negative Signs) ******
 Frequencies --  -XXX.XXXX   YY.YYYY   ZZ.ZZZZ
```

---

## Step 6 IRC 验证

TS 必须连接正确的反应物和产物。[`04-irc/g16_irc.gjf`](04-irc/g16_irc.gjf)：

```text
%nprocshared=16
%mem=24GB
%oldchk=ts_attack.chk
%chk=irc_attack.chk
#p B3LYP/Def2SVP EmpiricalDispersion=GD3BJ SCRF=(CPCM,Solvent=Water,Read) IRC=(RCFC,MaxPoints=40,StepSize=8,MaxCycle=50) Geom=AllCheck Guess=Read

Eps=4.0

```

| 选项 | 说明 |
|---|---|
| `RCFC` | 从 chk 读取 TS 频率作业的笛卡尔力常数（IRC 必须给出 `RCFC` 或 `CalcFC`，见 [IRC 关键词文档](https://gaussian.com/irc/)） |
| `MaxPoints=40` | 每个方向最多 40 个点（默认 10，对于有软模式的簇模型通常不够） |
| `StepSize=8` | 步长，单位 0.01 Bohr（默认 10） |
| `Geom=AllCheck` | 结构、电荷多重度、标题都从 chk 读，所以没有标题和坐标段 |
| PCM 输入段 | 仍需要 `Eps=4.0`，否则会回到水的介电常数 |

ORCA（[`04-irc/orca_irc.inp`](04-irc/orca_irc.inp)）使用 `! IRC`，在 `%irc` 中设置 `Direction both`、`InitHess read`、`Hess_Filename "orca_ts.hess"`，并用 `Monitor_Internals { B 1 51 } end` 在输出中跟踪 O···C 距离。

> [!WARNING]
> 冻结原子与 IRC：IRC 在质量加权坐标中积分，冻结原子的处理方式在不同程序/版本中可能不同。**一定要检查** IRC 各点中 Cβ 是否保持不动。如果 Cβ 移动了，得到的路径只能定性参考：取 IRC 两端结构，按 Step 4 的方法（加冻结）重新优化，再用重新优化的端点计算能量。这也是标准做法：**IRC 端点永远要再优化**。

检查：

```bash
grep -A60 "Summary of reaction path following" irc_attack.log | head -70
python ../../../Gaussian/scripts/scan_irc_profile.py irc_attack.log --mode irc -o irc_attack.png
```

成功标准：一端回到 Zn–OH⁻ + CO2（O···C 拉长到约 2.5 Å 以上，CO2 接近线形）；另一端得到 Zn 结合的碳酸氢根，O–C 成键，CO2 部分弯成三角平面。把产物端点拿去做 Step 4 的 `g16_opt_product.gjf`，得到 P1。

---

## Step 7 频率与热校正

### 7.1 冻结原子带来的问题

ModRedundant 冻结只作用于**优化**。频率计算会对全部 3N 个自由度构建 Hessian。冻结原子所在的位置并不是这些自由度上的驻点，所以：

- 会出现少量**小虚频**（通常与冻结的甲基/Cβ 平动或扭转有关）；
- 投影掉平动和转动的假设不完全成立；
- 低频模式的熵对结构很敏感。

处理原则（Siegbahn & Himo 综述中的常规做法，[DOI: 10.1002/wcms.13](https://doi.org/10.1002/wcms.13)）：

| 情况 | 处理 |
|---|---|
| 极小点出现 1–3 个绝对值很小的虚频，且只涉及冻结原子 | 可接受，在文章中说明；计算熵时把它们当作低频处理，或者忽略 |
| 虚频涉及反应中心 | 不是真正的驻点，必须重新优化 |
| TS 有 1 个反应虚频 + 若干冻结相关小虚频 | 可接受 |
| 想严格一些 | 用部分 Hessian（PHVA）只对非冻结原子做振动分析；ORCA 的 `%freq` 中有 `Partial_Hess` 选项，具体用法请查阅 ORCA 手册 “Frequency calculations” 一节 |

### 7.2 quasi-RRHO

簇模型中有大量低于 100 cm⁻¹ 的软模式。谐振子近似会严重高估它们的熵贡献，常用 Grimme 的 quasi-RRHO 处理（[DOI: 10.1002/chem.201200497](https://doi.org/10.1002/chem.201200497)）：

- **Gaussian**：自身不提供 quasi-RRHO。可以用仓库中的 [`thermo_qrrho.py`](../../Gaussian/scripts/thermo_qrrho.py)，或者 GoodVibes（[patonlab/GoodVibes](https://github.com/patonlab/GoodVibes)，[DOI: 10.12688/f1000research.22758.1](https://doi.org/10.12688/f1000research.22758.1)）对 log 做后处理：

  ```bash
  python ../../../Gaussian/scripts/thermo_qrrho.py opt_reactant.log ts_attack.log opt_product.log \
         --qs grimme --cutoff 100 --csv thermo.csv
  ```

- **ORCA 6**：`%freq QuasiRRHO true end`（模板中已经打开）。

> [!NOTE]
> 一些簇模型研究只报告 ΔE + ΔZPE，而不报告 ΔG，理由是冻结模型的熵不可靠。对单分子步骤（例如本例中“CO2 已在口袋中”的反应物 → TS），平动/转动熵的变化很小，ΔE+ZPE 与 ΔG 往往相近；但对底物结合这类双分子步骤，熵就非常重要，必须单独讨论。你的文章中要写清楚报告的是哪一个量。

### 7.3 从输出中取热校正

```bash
grep "Thermal correction to Gibbs Free Energy" opt_reactant.log ts_attack.log opt_product.log
grep "Zero-point correction" opt_reactant.log ts_attack.log opt_product.log
# ORCA
grep "G-E(el)" orca_*.out
```

把这些数（Hartree）填入 [`08-profile/energies_template.csv`](08-profile/energies_template.csv) 的 `G_corr_hartree` 列。

---

## Step 8 高精度单点

def2-SVP 的能量误差较大，用 def2-TZVP 在 SVP 结构上做单点：

$$
G = E_{\rm sp}^{\rm TZVP,\,\varepsilon=4} + G_{\rm corr}^{\rm SVP}
$$

[`06-sp/g16_sp_tzvp.gjf`](06-sp/g16_sp_tzvp.gjf)：

```text
%nprocshared=16
%mem=32GB
%oldchk=opt_reactant.chk
%chk=sp_tzvp_reactant.chk
#p B3LYP/Def2TZVP EmpiricalDispersion=GD3BJ SCRF=(CPCM,Solvent=Water,Read) Geom=Check Guess=Read

hCA II cluster, reactant, single point B3LYP-D3(BJ)/def2-TZVP, CPCM eps=4, at the def2-SVP geometry

1 1

Eps=4.0

```

`Guess=Read` 从 SVP 波函数投影到 TZVP 作为初猜，可以加快 SCF 收敛。对 TS 和 P1，复制该文件并修改 `%oldchk`、`%chk` 和标题。

需要准备的单点（每个驻点都要做）：

| 文件 | 目的 |
|---|---|
| `g16_sp_tzvp.gjf` / `orca_sp_tzvp.inp` | 主结果：def2-TZVP，CPCM ε=4 |
| `g16_sp_tzvp_eps20.gjf` | 敏感性检验：ε=20（Step 10） |
| `orca_sp_tzvp_smd.inp` | 可选：SMD 乙醚（ε 约 4.2，[SMD 原文 DOI: 10.1021/jp810292n](https://doi.org/10.1021/jp810292n)），与 CPCM ε=4 对照 |

> [!TIP]
> 也可以在单点中换用其他泛函（例如 ωB97X-D、M06-2X 或双杂化泛函）检验泛函依赖性，见 [Gaussian 第3章](../../Gaussian/tutorials/03-method-basis.md)。如果不同泛函给出的能垒相差很大（例如超过 5 kcal/mol），结论的可信度就要打折扣。

```bash
grep "SCF Done" sp_tzvp_*.log | awk '{print $1, $5}'
# ORCA
grep "FINAL SINGLE POINT ENERGY" orca_sp_tzvp*.out
```

---

## Step 9 能垒与能量剖面

### 9.1 填表

把 [`08-profile/energies_template.csv`](08-profile/energies_template.csv) 复制为 `energies.csv`，用你自己的结果替换每一个 `NaN`，再把 `placeholder` 改为 `no`：

```text
label,state_type,E_sp_hartree,G_corr_hartree,n_imag,placeholder,source_files
R,minimum,<your E_sp>,<your G_corr>,0,no,sp_tzvp_reactant.log; opt_reactant.log
TS1,ts,<your E_sp>,<your G_corr>,1,no,sp_tzvp_ts.log; ts_attack.log
P1,minimum,<your E_sp>,<your G_corr>,0,no,sp_tzvp_product.log; opt_product.log
```

> [!IMPORTANT]
> 仓库自带的 CSV **全部是占位行**，不是结果。`plot_profile.py` 会检查 `placeholder` 列和缺失值，只要还有一行没填好就拒绝作图。同时它会核对虚频个数（极小点为 0、TS 为 1）。如果你决定接受冻结相关的小虚频，请在 `n_imag` 中填“主要虚频”的个数，并在 `source_files` 列或你的笔记中说明。

### 9.2 作图

```bash
cd 08-profile
python plot_profile.py energies.csv                  # 只打印表格
python plot_profile.py energies.csv -o profile.png   # 表格 + 图
python plot_profile.py energies.csv --ref R --unit kj
```

输出表格的格式（示意输出，数值仅供参考；`x.x` 处将显示你的数值）：

```text
label    type          G (Hartree)  dG (kcal/mol)       n_imag
R        minimum      -XXXX.XXXXXX            0.0            0
TS1      ts           -XXXX.XXXXXX            x.x            1
P1       minimum      -XXXX.XXXXXX            x.x            0
```

### 9.3 报告时写清楚

一张合格的能量表至少要注明：

| 项目 | 例子 |
|---|---|
| 几何与频率的理论水平 | B3LYP-D3(BJ)/def2-SVP，CPCM(ε=4) |
| 单点的理论水平 | B3LYP-D3(BJ)/def2-TZVP，CPCM(ε=4) |
| 报告的量 | ΔG（quasi-RRHO，100 cm⁻¹）或 ΔE+ZPE |
| 模型 | 原子数、总电荷、包含哪些残基、冻结哪些原子 |
| 参考态 | 以哪个结构为零点（这里是 R：CO2 已在口袋中的复合物） |
| 敏感性 | ε=4 与 ε=20 的差值；大模型与小模型的差值 |

把你的结果与 [2.4 节](#24-一个可以用来对照的实验量级) 的实验量级和第 2.3 节列出的计算文献对照。差别大时要讨论原因（模型大小、泛函、OH⁻ 和水的构象等），而不是直接下结论。

---

## Step 10 模型收敛性检验

簇模型最受质疑的地方就是“模型够不够大”。至少要做以下检验：

### 10.1 模型尺寸

| 模型 | 命令 | 原子数 | 总电荷 |
|---|---|---|---|
| M1（最小，本教程主线） | `cut_cluster.py --add-co2 --co2-pdb 3D92.pdb` | 54 | +1 |
| M2：+ Glu106 | `... --with-glu106` | 61 | 0 |
| M3：+ 第 2 个水 | `... --with-glu106 --waters 2 --water-cutoff 3.7` | 64（以 `build_report.txt` 为准） | 0 |
| M4：+ 第二配位层（Gln92、Asn244 侧链、Glu117 等） | 需要扩展脚本或用 RINRUS | 约 100 以上 | 视加入的残基而定 |
| M5：+ 疏水口袋（Val121、Val143、Leu198、Trp209） | 同上 | 约 150–200 | 同上 |

对每个模型重复 Step 4–8（至少重复 TS 和 R），比较 ΔG‡ 的变化。当继续扩大模型后能垒变化小于约 1–2 kcal/mol 时，可以认为基本收敛（这是簇模型文献中常用的经验判据，并非严格标准）。

> [!NOTE]
> 加入 Glu106 后总电荷从 +1 变为 0，静电环境变化明显，能垒的变化可能比加入一个中性残基大。这本身就是很有教学意义的观察：带电残基对 Zn–OH⁻ 亲核性的影响。

### 10.2 介电常数

用 `g16_sp_tzvp_eps20.gjf` 在相同结构上把 ε 从 4 改为 20，比较能垒变化：

| 观察 | 解释 |
|---|---|
| ε=4 与 ε=20 的能垒差别很小（例如小于 1–2 kcal/mol） | 模型已经足够大，环境的静电屏蔽影响已经基本被显式残基包含（Himo 等人强调过这一判据，[DOI: 10.1021/jacs.7b02671](https://doi.org/10.1021/jacs.7b02671)） |
| 差别大 | 模型对长程静电敏感，需要扩大模型，或者改用 QM/MM |

### 10.3 冻结策略

把冻结的原子从 Cβ 改为 Cα（即保留到 Cα，在 Cα 上封端，模型会稍大），或者只冻结 Cβ 的一个方向，比较能垒变化。如果冻结方式改变后结果变化很大，说明模型中存在不真实的应力。

---

## Step 11 进阶：ONIOM(B3LYP:Amber) 与 CP2K QM/MM 替代方案

当簇模型的收敛性检验不理想，或者你要研究的是 His64 质子穿梭这类长程过程时，就需要 QM/MM。

### 11.1 Gaussian ONIOM

ONIOM(B3LYP/def2-SVP:Amber) 把 QM 区（与簇模型相同）放在完整的、经过 MD 平衡的蛋白中。骨架模板见 [`07-oniom/g16_oniom_template.gjf`](07-oniom/g16_oniom_template.gjf) 和 [`07-oniom/README.md`](07-oniom/README.md)，完整语法见 [Gaussian 第15章 ONIOM](../../Gaussian/tutorials/15-oniom.md)。准备流程：

1. 用 PDB2PQR/H++ 指定质子化态（Step 1）；
2. Zn 位点参数化：用 AmberTools 的 MCPB.py（[教程](https://ambermd.org/tutorials/advanced/tutorial20/mcpbpy.php)，[DOI: 10.1021/acs.jcim.5b00674](https://doi.org/10.1021/acs.jcim.5b00674)）生成键合模型参数；
3. tleap 溶剂化、加反离子，经过 MM 最小化和 MD 平衡；
4. 从轨迹中抽取若干快照，把 QM 区之外 10–15 Å 以上的原子冻结；
5. 写 ONIOM 输入：`ONIOM(B3LYP/Def2SVP:Amber)=EmbedCharge`，先用机械嵌入（ME）优化，再用电子嵌入（EE）优化（[第15章](../../Gaussian/tutorials/15-oniom.md) 推荐的策略）；
6. 对多个快照重复 TS 搜索，报告平均值和离散程度。

ONIOM 的综述见 Chung et al. 2015（[DOI: 10.1021/cr5004419](https://doi.org/10.1021/cr5004419)），QM/MM 方法的一般综述见 Senn & Thiel 2009（[DOI: 10.1002/anie.200802019](https://doi.org/10.1002/anie.200802019)）。

### 11.2 CP2K QM/MM

CP2K 的 `&QMMM` 适合在周期性体系中做 QM/MM 分子动力学和自由能计算，也可以用 GFN-xTB 或 DFTB 作为 QM 方法做长时间采样。设置方法见 [CP2K 第16章 QM/MM 与半经验方法](../../CP2K/tutorials/16-qmmm-xtb.md)。CP2K 自带的 SCC-DFTB 参数包含 Zn（见该章），可以作为 Zn 酶 QM/MM 采样的便宜起点，但定量能垒仍要用 DFT/MM 精修。

| 方案 | 适合 | 主要成本 |
|---|---|---|
| 簇模型（本案例） | 局域反应、机理筛选 | 模型收敛性检验 |
| ONIOM(QM:Amber) | 单个或少量快照上的 TS 搜索 | MM 准备、Zn 参数、多快照 |
| CP2K QM/MM MD + 增强采样 | 自由能面、质子穿梭、构象涨落 | 采样时间 |

其他值得参考的工具和教程：ASH（Python 多尺度建模框架，[GitHub: RagnarB83/ash](https://github.com/RagnarB83/ash)，[文档](https://ash.readthedocs.io/)）、CCPBioSim QM/MM workshop（[GitHub: CCPBioSim/qmmm-workshop](https://github.com/CCPBioSim/qmmm-workshop)），以及 Cui 等人对 CA 质子转移的 QM/MM 系列研究（[DOI: 10.1021/jp021931v](https://doi.org/10.1021/jp021931v)、[DOI: 10.1021/bi701950j](https://doi.org/10.1021/bi701950j)）。

---

## Step 12 结果分析

### 12.1 NBO 电荷与二阶微扰

在优化结构上做 NBO 分析，比较 R、TS1、P1 中 Zn、O(OH⁻)、C(CO2) 的自然电荷，以及 O 孤对 → CO2 π\* 的二阶微扰能 E(2)：

```text
%oldchk=ts_attack.chk
%chk=nbo_ts.chk
#p B3LYP/Def2SVP EmpiricalDispersion=GD3BJ SCRF=(CPCM,Solvent=Water,Read) Pop=NBO Geom=AllCheck Guess=Read

Eps=4.0

```

Gaussian 内置的是 NBO 3.1；需要 NBO 7 功能时用 `Pop=NBO7`（需要单独安装 NBO 7，[DOI: 10.1002/jcc.25873](https://doi.org/10.1002/jcc.25873)）。电荷分析的更多方法（CM5、Hirshfeld、ESP）见 [Gaussian 第12章](../../Gaussian/tutorials/12-population-charges.md)。

可以关注的问题：

- 从 R 到 TS1，CO2 上的负电荷增加了多少？Zn 的电荷是否基本不变（Zn 作为 Lewis 酸稳定 OH⁻，而不是直接参与电子转移）？
- Thr199 与 OH⁻ 之间的氢键在 TS1 中变强还是变弱？

### 12.2 IGMH/NCI 弱相互作用

用 Multiwfn（[DOI: 10.1002/jcc.22885](https://doi.org/10.1002/jcc.22885)、[DOI: 10.1063/5.0216272](https://doi.org/10.1063/5.0216272)）做 IGMH（[DOI: 10.1002/jcc.26812](https://doi.org/10.1002/jcc.26812)），把 CO2 定义为一个片段、其余簇为另一个片段，就能直观看到 CO2 在疏水口袋中的结合方式，以及 TS 中 O···C 正在形成的键：

```bash
formchk ts_attack.chk ts_attack.fchk
Multiwfn ts_attack.fchk        # 主功能 20 -> 11 (IGMH)，按提示定义两个片段
```

完整的操作步骤、VMD 渲染和批处理命令流见 [Gaussian 第16章 后处理与 Multiwfn](../../Gaussian/tutorials/16-postprocess-multiwfn.md)。

### 12.3 几何参数表

最后整理一张关键几何参数表（R、TS1、P1 的 Zn–O、O···C、O–C–O 角、Zn–N 平均值、O(H)···OG1(Thr199)），与 2CBA/3D92 的晶体数据（[2.1 节](#21-酶与反应)）对照。冻结原子附近结构的偏差能反映模型应力。

---

## 参数选择与收敛建议

| 项目 | 建议 | 备注 |
|---|---|---|
| 模型大小 | 先做 50–60 原子的最小模型，再至少做一个更大的模型 | Step 10 |
| 泛函 | B3LYP-D3(BJ) 为主，至少用一个其他泛函做单点对照 | |
| 优化基组 | def2-SVP | Zn 不需要 ECP |
| 单点基组 | def2-TZVP | 有余力时可试 def2-TZVPP |
| 介电常数 | ε = 4，并检验 ε = 20 | Step 10.2 |
| 积分格点 | Gaussian 16 默认 UltraFine；ORCA 6 默认格点对本体系一般够用 | 频率出现奇怪的低频时可加密格点 |
| SCF | 带电簇用默认设置通常可以收敛；困难时用 `SCF=XQC`（Gaussian）或 `SlowConv`（ORCA） | |
| 优化步数 | MaxCycles 150–300 | 软模式多 |
| TS | 扫描 → `Opt=(TS,CalcFC,NoEigenTest)`；必要时 `RecalcFC=N` | |
| 热化学 | quasi-RRHO（截断 100 cm⁻¹） | Step 7 |
| 构象 | OH⁻ 的 H 朝向、第二层水取向都要试；CREST 可以辅助 | |

---

## 常见问题与报错

**1. `The combination of multiplicity 1 and 229 electrons is impossible.`**

- 原因：电子数为奇数。常见原因是 Zn 结合水没有去掉一个 H（变成了 H2O 而不是 OH⁻，但电荷仍写 +1），或者某个 His 少了 N–H。
- 处理：数一数 H 原子数，与 `build_report.txt` 中的电子数对照；核对 Step 2.5 的电荷账本。

**2. 优化中 OH⁻ 的 H 转移到了 Thr199 或水上，Zn 上变成了水。**

- 原因：在过小的模型和低介电常数中，Zn–OH⁻ 与 Thr199–OH 或水之间的质子分布可能不同于酶中的情况；也可能是初始 H 朝向不合理。
- 处理：检查质子化态是否合理；尝试 `--oh-mode anti`；加入 Glu106（它接受 Thr199 的氢键，使 Thr199 更难被去质子化）；增大模型。在文章中要报告这种行为，它本身就说明模型过小。

**3. CO2 在优化中飘出口袋，或者被拉到 Zn 上。**

- 原因：最小模型中没有疏水口袋残基（Val121、Val143、Leu198、Trp209），CO2 缺乏范德华约束。
- 处理：加入这些残基的侧链（并冻结其截断原子），或者先在 O···C 约束下优化（ModRedundant `B 2 52 F`），再放开约束。

**4. `Error termination via Lnk1e in .../l103.exe` 并伴随 `Linear angle in Tors` 或 `FormBX had a problem`。**

- 原因：CO2 的 O=C=O 接近 180°，冗余内坐标在线性角附近出现奇异。
- 处理：先从最新结构重启（`Geom=Check Guess=Read`），让程序重建内坐标；或者在 ModRedundant 段中用 `L` 类型为 O=C=O 定义线性弯曲（[Gaussian 第5章](../../Gaussian/tutorials/05-opt.md)）。不要直接改用 `Opt=Cartesian`：纯笛卡尔优化不支持冻结变量，截断原子会被放开。更多排查方法见 [Gaussian 第18章](../../Gaussian/tutorials/18-errors-capstone.md)。

**5. TS 优化后有 2 个或更多明显的虚频。**

- 原因：初猜离 TS 太远，或者第二个虚频对应 OH⁻/甲基转动等软模式。
- 处理：先看动画。若第二个模式只涉及冻结原子附近，可以接受（Step 7）；若涉及反应中心，用 GaussView 沿该模式微扰结构后重新做 `Opt=(TS,CalcFC,NoEigenTest)`，或者从扫描中更接近极大值的点重新开始。

**6. `Number of steps exceeded, NStep= ...`（L9999）。**

- 原因：软自由度多，优化慢。
- 处理：`Opt=(ModRedundant,MaxCycles=300)`，从最后一步重启（`Geom=Check Guess=Read`）；ORCA 中加大 `MaxIter`。如果能量在反复振荡，检查是否有某个水或甲基在两个取向间来回翻转。

**7. 频率中出现与冻结原子相关的小虚频。**

- 原因：见 Step 7.1，这是冻结簇模型的固有现象。
- 处理：确认它们只涉及冻结原子附近、绝对值小；在文章中说明；用 quasi-RRHO 处理低频模式。不要为了“消掉”这些虚频而解冻所有原子，否则模型会塌缩。

**8. 改了 `Eps` 却发现结果和水中完全一样。**

- 原因：`SCRF=(...)` 中忘了写 `Read`，或者 PCM 输入段前后少了空行，导致 `Eps=4.0` 没有被读入。
- 处理：在输出中查找 `Eps=`（[Gaussian 第9章](../../Gaussian/tutorials/09-solvation.md) 给出了输出格式），确认数值是 4.0。

**9. ORCA 报错或约束没有生效。**

- 原因：ORCA 原子编号从 **0** 开始，直接照抄 Gaussian 的编号会冻结错误的原子。
- 处理：用 `fill_template.py` 生成约束（它会自动减 1），或手工减 1；在输出开头查看约束列表进行核对。

**10. xTB 输出 `FAILED TO CONVERGE GEOMETRY OPTIMIZATION` 并生成 `NOT_CONVERGED` 文件。**

- 原因：初始 H 有严重接触，或者 SCC 不收敛。
- 处理：先用 `--opt crude` 或 `--opt loose` 粗优化再收紧；SCC 不收敛时可适当提高电子温度（`--etemp`，默认 300 K），或先用不需要 SCC 迭代的 GFN0-xTB（`--gfn 0`）或 GFN-FF（`--gfnff`）粗清理结构，再回到 GFN2-xTB；同时目视检查是否有原子重叠。

---

## 练习

1. **换一个起点**：完全用 3D92 搭模型（Zn、His、Thr199、CO2 都取自 3D92），与“2CBA + 3D92 的 CO2”模型比较 xTB 预优化后的关键距离。提示：`cut_cluster.py --pdb 3D92.pdb` 时 Zn 结合水为 HOH 401，脚本会从 `LINK` 记录自动识别；`--co2-pdb 3D92.pdb` 仍然可用。
2. **质子化态的影响**：把 Zn–OH⁻ 改成 Zn–H2O（总电荷 +2），优化反应物复合物，扫描 O···C 距离，体会为什么 Zn 结合水必须先去质子化才有亲核性。
3. **Lipscomb vs Lindskog**：从 P1 出发分别搜索 (a) 碳酸氢根内部质子转移的 TS（Lipscomb），(b) 碳酸氢根旋转、换另一个氧与 Zn 配位的 TS（Lindskog）。比较两个能垒，再与 2.3 节的文献对照。
4. **介电常数敏感性**：对 M1 和 M2 模型都做 ε=4 与 ε=20 的单点，看哪个模型对 ε 更不敏感。
5. **泛函依赖性**：在 B3LYP-D3(BJ)/def2-SVP 结构上，用另外两个泛函做 def2-TZVP 单点，比较 ΔE‡ 的离散程度。
6. **IGMH 可视化**：对 R 和 TS1 做 IGMH，把 CO2 定义为片段 1，比较 CO2 与 Thr199、疏水残基之间的相互作用等值面。

---

## 参考资料

**结构与实验**

- PDB [2CBA](https://www.rcsb.org/structure/2CBA)：K. Håkansson, M. Carlsson, L. A. Svensson, A. Liljas, *J. Mol. Biol.* 1992, 227, 1192. [DOI: 10.1016/0022-2836(92)90531-N](https://doi.org/10.1016/0022-2836(92)90531-N)
- PDB [3D92](https://www.rcsb.org/structure/3D92)：J. F. Domsic et al., *J. Biol. Chem.* 2008, 283, 30766. [DOI: 10.1074/jbc.M805353200](https://doi.org/10.1074/jbc.M805353200)（CO2 结合态）
- S. Lindskog, *Pharmacol. Ther.* 1997, 74, 1. [DOI: 10.1016/S0163-7258(96)00198-2](https://doi.org/10.1016/S0163-7258(96)00198-2)（结构与机理综述）
- D. N. Silverman, S. Lindskog, *Acc. Chem. Res.* 1988, 21, 30. [DOI: 10.1021/ar00145a005](https://doi.org/10.1021/ar00145a005)（限速质子转移）
- D. N. Silverman, R. McKenna, *Acc. Chem. Res.* 2007, 40, 669. [DOI: 10.1021/ar7000588](https://doi.org/10.1021/ar7000588)（溶剂介导的质子转移）
- D. W. Christianson, C. A. Fierke, *Acc. Chem. Res.* 1996, 29, 331. [DOI: 10.1021/ar9501232](https://doi.org/10.1021/ar9501232)（Zn 结合位点）

**CA 机理的计算研究**

- J.-Y. Liang, W. N. Lipscomb, *Biochemistry* 1987, 26, 5293. [DOI: 10.1021/bi00391a012](https://doi.org/10.1021/bi00391a012)（Zn 结合碳酸氢根的内部质子转移，即 Lipscomb 路线）
- J.-Y. Liang, W. N. Lipscomb, *Biochemistry* 1988, 27, 8676. [DOI: 10.1021/bi00423a025](https://doi.org/10.1021/bi00423a025)（Zn 结合水与 His64 之间的质子转移，对应再生半反应）
- K. M. Merz Jr., *J. Am. Chem. Soc.* 1991, 113, 406. [DOI: 10.1021/ja00002a004](https://doi.org/10.1021/ja00002a004)
- Y.-J. Zheng, K. M. Merz Jr., *J. Am. Chem. Soc.* 1992, 114, 10498. [DOI: 10.1021/ja00052a054](https://doi.org/10.1021/ja00052a054)
- A. Bottoni et al., *J. Am. Chem. Soc.* 2004, 126, 1542. [DOI: 10.1021/ja030336j](https://doi.org/10.1021/ja030336j)
- G. P. Miscione et al., *Theor. Chem. Acc.* 2007, 118, 193. [DOI: 10.1007/s00214-007-0274-x](https://doi.org/10.1007/s00214-007-0274-x)
- Q. Cui, M. Karplus, *J. Phys. Chem. B* 2003, 107, 1071. [DOI: 10.1021/jp021931v](https://doi.org/10.1021/jp021931v)；D. Riccardi et al., *Biochemistry* 2008, 47, 2369. [DOI: 10.1021/bi701950j](https://doi.org/10.1021/bi701950j)；D. Riccardi, S. Yang, Q. Cui, *BBA Proteins Proteom.* 2010, 1804, 342. [DOI: 10.1016/j.bbapap.2009.07.026](https://doi.org/10.1016/j.bbapap.2009.07.026)（QM/MM 质子转移）

**簇模型方法学**

- P. E. M. Siegbahn, F. Himo, *WIREs Comput. Mol. Sci.* 2011, 1, 323. [DOI: 10.1002/wcms.13](https://doi.org/10.1002/wcms.13)
- P. E. M. Siegbahn, F. Himo, *J. Biol. Inorg. Chem.* 2009, 14, 643. [DOI: 10.1007/s00775-009-0511-y](https://doi.org/10.1007/s00775-009-0511-y)
- F. Himo, *J. Am. Chem. Soc.* 2017, 139, 6780. [DOI: 10.1021/jacs.7b02671](https://doi.org/10.1021/jacs.7b02671)
- F. Himo, S. P. de Visser, *Commun. Chem.* 2022, 5, 29. [DOI: 10.1038/s42004-022-00642-2](https://doi.org/10.1038/s42004-022-00642-2)
- X. Sheng, F. Himo, *Acc. Chem. Res.* 2023, 56, 938. [DOI: 10.1021/acs.accounts.2c00795](https://doi.org/10.1021/acs.accounts.2c00795)
- M. R. A. Blomberg, T. Borowski, F. Himo, R.-Z. Liao, P. E. M. Siegbahn, *Chem. Rev.* 2014, 114, 3601. [DOI: 10.1021/cr400388t](https://doi.org/10.1021/cr400388t)（金属酶）
- R.-Z. Liao, J.-G. Yu, F. Himo, *J. Chem. Theory Comput.* 2011, 7, 1494. [DOI: 10.1021/ct200031t](https://doi.org/10.1021/ct200031t)（模型尺寸与介电常数）
- H. J. Kulik et al., *J. Phys. Chem. B* 2016, 120, 11381. [DOI: 10.1021/acs.jpcb.6b07814](https://doi.org/10.1021/acs.jpcb.6b07814)（QM 区大小）
- RINRUS：D. A. Wappett et al., *WIREs Comput. Mol. Sci.* 2026, 16, e70078. [DOI: 10.1002/wcms.70078](https://doi.org/10.1002/wcms.70078)；代码 [natedey/RINRUS](https://github.com/natedey/RINRUS)

**方法与软件**

- GFN2-xTB：C. Bannwarth, S. Ehlert, S. Grimme, *J. Chem. Theory Comput.* 2019, 15, 1652. [DOI: 10.1021/acs.jctc.8b01176](https://doi.org/10.1021/acs.jctc.8b01176)；ALPB：[DOI: 10.1021/acs.jctc.1c00471](https://doi.org/10.1021/acs.jctc.1c00471)；CREST：[DOI: 10.1063/5.0197592](https://doi.org/10.1063/5.0197592)
- B3LYP：[DOI: 10.1063/1.464913](https://doi.org/10.1063/1.464913)、[DOI: 10.1021/j100096a001](https://doi.org/10.1021/j100096a001)；D3(BJ)：[DOI: 10.1002/jcc.21759](https://doi.org/10.1002/jcc.21759)；def2 基组：[DOI: 10.1039/b508541a](https://doi.org/10.1039/b508541a)
- CPCM：V. Barone, M. Cossi, *J. Phys. Chem. A* 1998, 102, 1995. [DOI: 10.1021/jp9716997](https://doi.org/10.1021/jp9716997)；SMD：[DOI: 10.1021/jp810292n](https://doi.org/10.1021/jp810292n)
- quasi-RRHO：S. Grimme, *Chem. Eur. J.* 2012, 18, 9955. [DOI: 10.1002/chem.201200497](https://doi.org/10.1002/chem.201200497)；GoodVibes：[DOI: 10.12688/f1000research.22758.1](https://doi.org/10.12688/f1000research.22758.1)
- PROPKA3：[DOI: 10.1021/ct100578z](https://doi.org/10.1021/ct100578z)，[jensengroup/propka](https://github.com/jensengroup/propka)；PDB2PQR/APBS：[DOI: 10.1002/pro.3280](https://doi.org/10.1002/pro.3280)，[Electrostatics/pdb2pqr](https://github.com/Electrostatics/pdb2pqr)，[文档](https://pdb2pqr.readthedocs.io/)；H++：[DOI: 10.1093/nar/gks375](https://doi.org/10.1093/nar/gks375)；PDBFixer：[openmm/pdbfixer](https://github.com/openmm/pdbfixer)
- MCPB.py：[DOI: 10.1021/acs.jcim.5b00674](https://doi.org/10.1021/acs.jcim.5b00674)；ONIOM 综述：[DOI: 10.1021/cr5004419](https://doi.org/10.1021/cr5004419)；QM/MM 综述：[DOI: 10.1002/anie.200802019](https://doi.org/10.1002/anie.200802019)
- NBO 7：[DOI: 10.1002/jcc.25873](https://doi.org/10.1002/jcc.25873)；Multiwfn：[DOI: 10.1002/jcc.22885](https://doi.org/10.1002/jcc.22885)、[DOI: 10.1063/5.0216272](https://doi.org/10.1063/5.0216272)；IGMH：[DOI: 10.1002/jcc.26812](https://doi.org/10.1002/jcc.26812)
- ORCA 6：F. Neese, *WIREs Comput. Mol. Sci.* 2025, 15, e70019. [DOI: 10.1002/wcms.70019](https://doi.org/10.1002/wcms.70019)

**官方文档**

- Gaussian 16：[Opt](https://gaussian.com/opt/)、[Geom](https://gaussian.com/geom/)、[IRC](https://gaussian.com/irc/)、[SCRF](https://gaussian.com/scrf/)、[Input Section Ordering](https://gaussian.com/input/)、[ONIOM](https://gaussian.com/oniom/)
- ORCA 6.0 手册：[Geometry optimizations, scans, TS, IRC](https://www.faccts.de/docs/orca/6.0/manual/contents/typical/optimizations.html)、[Frequencies（含 QuasiRRHO、IRC 选项）](https://www.faccts.de/docs/orca/6.0/manual/contents/detailed/frequencies.html)、[Solvation（CPCM/SMD）](https://www.faccts.de/docs/orca/6.0/manual/contents/detailed/solvationmodels.html)、[ORCA tutorials](https://www.faccts.de/docs/orca/6.0/tutorials/)
- xTB：[文档首页](https://xtb-docs.readthedocs.io/en/latest/)、[xcontrol 详细输入](https://xtb-docs.readthedocs.io/en/latest/xcontrol.html)、[PES 扫描](https://xtb-docs.readthedocs.io/en/latest/scan.html)、[隐式溶剂](https://xtb-docs.readthedocs.io/en/latest/gbsa.html)

**相关教程与资源（只链接，不复制）**

- Sobereva《要善用簇模型，不要盲目用ONIOM算蛋白质-小分子相互作用问题》：[sobereva.com/597](http://sobereva.com/597)（中文博文；簇模型与 ONIOM 的取舍）
- CCPBioSim QM/MM workshop：[CCPBioSim/qmmm-workshop](https://github.com/CCPBioSim/qmmm-workshop)（酶 QM/MM 的动手教程）
- ASH 多尺度框架：[RagnarB83/ash](https://github.com/RagnarB83/ash)（用 Python 统一调度 ORCA/xTB/OpenMM 做 QM/MM）
- AmberTools MCPB.py 教程：[ambermd.org tutorial 20](https://ambermd.org/tutorials/advanced/tutorial20/mcpbpy.php)（金属位点力场参数化）
- 本仓库：[Gaussian 第7章 TS 与 IRC](../../Gaussian/tutorials/07-ts-irc.md)、[第15章 ONIOM](../../Gaussian/tutorials/15-oniom.md)、[第16章 Multiwfn](../../Gaussian/tutorials/16-postprocess-multiwfn.md)、[CP2K 第16章 QM/MM](../../CP2K/tutorials/16-qmmm-xtb.md)

---

[返回目录](../README.md)（本板块目前只有一个案例，暂无上一章/下一章）
