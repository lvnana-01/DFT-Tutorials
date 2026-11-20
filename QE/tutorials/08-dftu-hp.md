# 第8章 DFT+U/U+V 与 hp.x 线性响应（DFT+U, DFT+U+V & Linear-Response Hubbard Parameters with hp.x）

本章讲 Quantum ESPRESSO 中的 Hubbard 修正：Dudarev DFT+U、DFT+U+J0、扩展的 DFT+U+V，以及 QE 7.1 起取代旧关键词的 `HUBBARD` 卡（`atomic` / `ortho-atomic` 投影）。重点是用 hp.x（基于 DFPT 的线性响应，Timrov–Marzari–Cococcioni）从第一性原理计算 U 和 V，并迭代到自洽 U。算例为反铁磁 NiO（U）和 LiCoO2（U+V）（Quantum ESPRESSO DFT+U / Hubbard U / hp.x linear response tutorial）。

> [!NOTE]
> 全部输入在示例目录 [08-dftu-hp](../examples/08-dftu-hp/README.md)。`HUBBARD` 卡语法对照官方 [Hubbard_input.pdf](https://www.quantum-espresso.org/Doc/user_guide_PDF/Hubbard_input.pdf) 与 `INPUT_PW`，hp.x 关键词对照 `INPUT_HP`（QE 7.5）。本仓库**不运行任何计算**；示意输出只演示格式，数值仅供参考。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 写出 QE ≥ 7.1 的 `HUBBARD` 卡：`U`、`J0`、`V`、多个流形、`Hubbard_occ`；
- 把 7.0 以前的旧输入（`lda_plus_u`、`Hubbard_U(i)`、`U_projection_type`）翻译成新语法；
- 理解 `atomic` 与 `ortho-atomic` 投影的差别，以及为什么 U 值离开具体投影/赝势就没有意义；
- 用 hp.x 计算磁性绝缘体 NiO 的 U（含"假金属 → 绝缘体"两步 SCF），并迭代到自洽；
- 用 hp.x 计算 DFT+U+V 并把 `HUBBARD.dat` 直接粘贴回 pw.x 输入；
- 读懂 `HUBBARD OCCUPATIONS` 块、`Hubbard energy` 与 hp.x 的 U 表。

| 场景 | 推荐做法 | 本章小节 |
|---|---|---|
| 过渡金属氧化物、带隙/磁矩被 PBE 严重低估 | Dudarev DFT+U，`ortho-atomic` | §3、§4 |
| 需要第一性原理 U | hp.x 线性响应（DFPT），必要时自洽迭代 | §5 |
| 共价性强的氧化物（Co–O、Mn–O、Fe–O） | DFT+U+V | §6 |
| 用他人论文中的 U | 只在赝势、泛函、投影都一致时才可比较 | §2.3 |
| 非共线 / SOC + U | 7.3 起支持 Dudarev 非共线 DFT+U(+V) 及其力、应力与 hp.x | §3.5 |

### 与 VASP 的对应关系

| VASP | QE（≥ 7.1） | 说明 |
|---|---|---|
| `LDAU = .TRUE.` | 出现 `HUBBARD` 卡即开启 | 不再需要 `lda_plus_u` |
| `LDAUTYPE = 2`（Dudarev） | 只写 `U`（不写 `J`） | 默认 Dudarev |
| `LDAUTYPE = 1`（Liechtenstein） | 写 `U` + `J`（d 壳层可加 `B`） | 程序自动识别 |
| `LDAUL` / `LDAUU` / `LDAUJ` | `U Fe-3d 4.0`（按 species 标签 + 流形） | QE 必须显式写流形 |
| PAW 投影子 | `atomic` / `ortho-atomic` / `norm-atomic` / `wf` / `pseudo` | **投影不同，U 不可互换** |
| `LDAUTYPE = 3`（超胞线性响应） | hp.x（原胞 + q 网格的 DFPT） | hp.x 不需要超胞 |
| `LMAXMIX = 4/6` | 无对应（QE 中占据矩阵随密度一起保存） | — |

VASP 侧的 DFT+U 与线性响应讲解见 [VASP 第7章](../../VASP/tutorials/07-magnetism-u-soc.md)。

---

## 2. 原理速览

### 2.1 DFT+U（Dudarev）

半局域泛函对局域 d/f 电子存在自相互作用误差，倾向于分数占据、过度离域。Dudarev 形式在选定的 Hubbard 流形上加一个只依赖 $U_{\rm eff}$ 的修正：

$$
E_{\rm DFT+U}=E_{\rm DFT}+\frac{U}{2}\sum_{I,\sigma}\mathrm{Tr}\!\left[\mathbf n^{I\sigma}\left(1-\mathbf n^{I\sigma}\right)\right]
$$

其中 $\mathbf n^{I\sigma}$ 是原子 I 上 Hubbard 流形的占据矩阵，$n^{I\sigma}_{mm'}=\sum_{k,v} f_{kv}\langle\psi_{kv}|P^I_{m'm}|\psi_{kv}\rangle$，由**投影函数**定义。修正项在整数占据（0 或 1）时为零，把占据态往下推、空态往上推，从而打开带隙。

QE 中 `U` 后再写 `J0`（Himmetoglu 等 2011 的 DFT+U+J0）可额外加入 Hund 修正；写 `J`（及 d 壳层的 `B`、f 壳层的 `E2`/`E3`）则切换到 Liechtenstein 形式。Hubbard_input.pdf 说明：**HUBBARD 卡中不出现 `J` 时自动使用 Dudarev 形式**。

### 2.2 DFT+U+V

Campo 与 Cococcioni（2010）把 Hubbard 修正推广到不同原子之间：

$$
E_{UV}=\sum_{I}\frac{U^I}{2}\mathrm{Tr}\!\left[\mathbf n^{II}(1-\mathbf n^{II})\right]-\sum_{I\neq J}\frac{V^{IJ}}{2}\mathrm{Tr}\!\left[\mathbf n^{IJ}\mathbf n^{JI}\right]
$$

$V^{IJ}$ 作用在原子间的广义占据矩阵上，有利于杂化（例如 Co-3d 与 O-2p），对共价性强的氧化物通常比单纯 +U 更能同时描述带隙与结构。

### 2.3 投影函数决定 U 的含义

| `HUBBARD {...}` 选项 | 含义（`INPUT_PW`） | hp.x 支持 | 力/应力 |
|---|---|---|---|
| `atomic` | 赝势中的原子轨道（非正交） | 是 | 是 |
| `ortho-atomic` | Löwdin 正交化的原子轨道，**官方推荐** | 是 | 是 |
| `norm-atomic` | 只做 Löwdin 归一化（未正交） | 否 | 否 |
| `wf` | 由 pmw.x 生成的 Wannier 函数（7.5 起另有 Wannier90 接口） | 否 | 否 |
| `pseudo` | 赝势投影子（需含全电子原子轨道的赝势） | 否 | 是 |

Hubbard_input.pdf 强调：**没有默认投影，必须显式写出**；`ortho-atomic` 避免了在轨道重叠区重复施加 U，通常给出更准确的占据。同一材料在 `atomic` 与 `ortho-atomic` 下的线性响应 U 可以差 1 eV 以上，赝势生成时的电子组态不同也会改变 U。结论：**U 不是元素的普适常数**，文献 U（尤其是 VASP 的 PAW 投影下的 U，或 Materials Project 拟合的 U）不能直接搬到 QE。

### 2.4 线性响应与 hp.x

Cococcioni 与 de Gironcoli（2005）的线性响应方法：对原子 J 的 Hubbard 流形施加势移 α，测量各原子占据的变化，得到裸响应 χ₀（非自洽）与屏蔽响应 χ（自洽），

$$
U^{I}=\left(\chi_0^{-1}-\chi^{-1}\right)_{II},\qquad V^{IJ}=\left(\chi_0^{-1}-\chi^{-1}\right)_{IJ}
$$

旧做法要在超胞中做一系列有限差分计算；hp.x 用 DFPT 把局域微扰分解为原胞中的一组单色 q 微扰（Timrov 等 2018、2021、2022），**q 网格 `nq1×nq2×nq3` 扮演超胞尺寸的角色**，无需超胞、可自动化。该方法针对的是"过度离域/过度杂化"的开壳层流形；Hubbard_input.pdf 指出它**不适用于 Hubbard 流形全满的闭壳层体系**（例如 Zn²⁺ 的 3d¹⁰）。

自洽 U：hp.x 得到的 U 依赖于基态（PBE 还是 PBE+U）。常见做法是 U_in → 计算 U_out → 令 U_in = U_out，直到二者之差小于阈值（例如 0.01–0.1 eV）。

---

## 3. HUBBARD 卡语法

### 3.1 基本写法（Dudarev DFT+U）

```fortran
HUBBARD {ortho-atomic}
U Ni1-3d 6.0
U Ni2-3d 6.0
```

- 卡名后的大括号（或圆括号）里写投影类型；
- 每行 `U  标签-流形  数值(eV)`；**标签必须与 `ATOMIC_SPECIES` 完全一致**（拆成 `Ni1`/`Ni2` 的 AFM 两个 species 都要写）；
- 流形（`3d`、`4f`、`2p`…）必须是赝势中存在的原子波函数；7.3.1 起要求赝势的 `PP_CHI` 带 `label`（多数库已具备，旧版 ATOMPAW 生成的 PAW 可能缺失，需手工补 label，见 Hubbard_input.pdf 第 5 节）；
- 官方示例（`HP/examples`、Hubbard_input.pdf）都把 HUBBARD 卡放在输入末尾（`K_POINTS` 之后），照此书写最稳妥；
- 初始占据默认从赝势读取，可用 `&SYSTEM` 中的 `Hubbard_occ(ityp,i)` 覆盖（i = 1–3 对应第 1–3 个流形）。

### 3.2 更多写法（摘自 Hubbard_input.pdf 的语法，数值为演示值）

```fortran
! DFT+U+J0 (Dudarev with Hund J0)
HUBBARD {ortho-atomic}
U  Mn-3d 5.0
J0 Mn-3d 1.0

! two Hubbard manifolds for one species (2nd channel = "background")
HUBBARD {ortho-atomic}
U Mn-3d 5.0
U Mn-3p 3.0

! Liechtenstein DFT+U+J (J present -> Liechtenstein; B for d shells)
HUBBARD {ortho-atomic}
U Mn-3d 5.0
J Mn-3d 1.0
B Mn-3d 1.1

! DFT+U+V: on-site U may also be written as "V Co-3d Co-3d 1 1 <U>"
HUBBARD {ortho-atomic}
U Co-3d 7.70
V Co-3d O-2p 1 19 0.75
```

`V` 行的两个整数是原子对的索引 I、J：I 是 `ATOMIC_POSITIONS` 中的原子序号，J 是程序内部构造的 3×3×3 虚拟超胞中的近邻序号，所以会出现比 `nat` 大得多的数字（如 19、46）。**这些索引不要手算**，用 hp.x 生成（`determine_num_pert_only = .true.` 可只输出索引而不做响应计算）。

### 3.3 旧语法 → 新语法对照

| 7.0 及以前（`&SYSTEM`） | 7.1+ | 备注 |
|---|---|---|
| `lda_plus_u = .true.` | 删除；写 `HUBBARD` 卡即开启 | 旧关键词会触发报错 |
| `lda_plus_u_kind = 0` | 只写 `U`（Dudarev） | |
| `lda_plus_u_kind = 1` | 写 `U` + `J`（Liechtenstein） | 要 Liechtenstein 而 J=0 时，写极小 J（如 1.d-12） |
| `lda_plus_u_kind = 2` | 出现 `V` 行即为 DFT+U+V | |
| `U_projection_type = 'ortho-atomic'` | `HUBBARD {ortho-atomic}` | 现在必须显式写 |
| `Hubbard_U(1) = 5.0` | `U Mn-3d 5.0` | 流形不再硬编码，需显式写 |
| `Hubbard_J0(1) = 1.0` | `J0 Mn-3d 1.0` | |
| `Hubbard_V(1,19,1) = 0.75` | `V Co-3d O-2p 1 19 0.75` | k=1（standard–standard） |
| `Hubbard_alpha(i)` | `ALPHA Ni-3d 0.05` | 旧式超胞线性响应的势移 |
| `starting_ns_eigenvalue(m,is,ityp)` | 仍在 `&SYSTEM` 中 | 引导轨道序（§4.4） |

只要 `&SYSTEM` 中出现任何旧关键词，pw.x 7.1+ 打印类似

```text
     WARNING!!! The input parameter lda_plus_u is obsolete.
     WARNING!!! The input syntax for DFT+Hubbard codes has changed since v7.1
```

并以 `DFT+Hubbard input syntax has changed since v7.1` 退出。网上大量 QE 6.x 教程（包括一些热门 YouTube 系列的配套输入）仍是旧语法，照抄前务必翻译。

### 3.4 关键词表

| 关键词 / 卡 | 推荐值 | 说明 |
|---|---|---|
| `HUBBARD {ortho-atomic}` | 首选 `ortho-atomic` | 无默认值，必须写 |
| `U label-manifold value` | 来自 hp.x | 单位 eV |
| `J0 label-manifold value` | 一般不用 | hp.x 不支持 `J0` ≠ 0 |
| `V labelI-mI labelJ-mJ I J value` | 来自 hp.x 的 `HUBBARD.dat` | 索引由程序生成 |
| `Hubbard_occ(ityp,i)` | 默认从赝势读 | 赝势按离子态生成时可能需要覆盖 |
| `starting_ns_eigenvalue(m,ispin,ityp)` | 默认 −1（不设） | 第一步覆盖占据矩阵本征值，引导轨道序 |
| `Hubbard_beta(i)` | 默认 0 | 旧式线性响应求 J0 的势移 |

### 3.5 非共线与 SOC

Hubbard_input.pdf 第 6 节：7.3 起 Dudarev DFT+U 与 DFT+U+V 扩展到非共线框架（含 Hubbard 力与应力，hp.x 亦支持；Binci & Marzari 2023）。不写 `J` 时使用 Dudarev，写非零 `J` 时切换为 Liechtenstein；**非共线 DFT+U+V 只有 Dudarev 形式**。非共线/SOC 的设置方法见 [第7章](07-magnetism-soc.md)。

---

## 4. 实战一：NiO 的 DFT+U

### 4.1 输入（`nio.dftu.scf.in`）

```fortran
&CONTROL
  calculation = 'scf'
  prefix      = 'nio_u'
  outdir      = './tmp'
  pseudo_dir  = './pseudo'
  tprnfor     = .true.
  tstress     = .true.
/
&SYSTEM
  ibrav       = 0
  nat         = 4
  ntyp        = 3
  ecutwfc     = 50.0
  ecutrho     = 400.0
  occupations = 'smearing'
  smearing    = 'gaussian'
  degauss     = 0.005
  nspin       = 2
  starting_magnetization(1) =  0.5
  starting_magnetization(2) = -0.5
  nbnd        = 30
/
&ELECTRONS
  conv_thr    = 1.0d-9
  mixing_beta = 0.3
/
ATOMIC_SPECIES
  Ni1  58.6934  ni_pbe_v1.4.uspp.F.UPF
  Ni2  58.6934  ni_pbe_v1.4.uspp.F.UPF
  O    15.999   O.pbe-n-kjpaw_psl.0.1.UPF
CELL_PARAMETERS angstrom
  4.170  2.085  2.085
  2.085  4.170  2.085
  2.085  2.085  4.170
ATOMIC_POSITIONS angstrom
  Ni1  0.000  0.000  0.000
  Ni2  2.085  2.085  0.000
  O    2.085  0.000  0.000
  O    4.170  2.085  0.000
K_POINTS automatic
  6 6 6 0 0 0
HUBBARD {ortho-atomic}
U Ni1-3d 6.0
U Ni2-3d 6.0
```

结构与 AFM 设置沿用 [第7章 §4](07-magnetism-soc.md)。`U = 6.0` 只是演示值——Hubbard_input.pdf 本身也声明其示例 U "must not be used for production calculations"。

### 4.2 运行与检查

```bash
cd QE/examples/08-dftu-hp
mpirun -np 4 pw.x -nk 2 -in nio.dftu.scf.in > nio.dftu.scf.out

grep -e '^!' -e 'Hubbard energy' nio.dftu.scf.out | tail -2
grep -e 'highest occupied' -e 'the Fermi energy' nio.dftu.scf.out | tail -1
grep 'Atomic magnetic moment for atom' nio.dftu.scf.out | tail -2
grep -A3 'Tr\[ns' nio.dftu.scf.out | tail -8
```

输出中的 `HUBBARD OCCUPATIONS` 块（每个 Hubbard 原子一段）给出：

（示意输出，数值仅供参考）

```text
     =================== HUBBARD OCCUPATIONS ===================
     ------------------------ ATOM    1 ------------------------
     Tr[ns(  1)] (up, down, total) =   x.xxxxx  x.xxxxx  x.xxxxx
     Atomic magnetic moment for atom   1 =   x.xxxxx
     SPIN  1
     eigenvalues:
     ...
     occupation matrix ns (before diag.):
     ...
```

成功判据：

- `Hubbard energy` 出现在总能分解中（`The total energy is ...` 之后）；
- 两个 Ni 的 `Atomic magnetic moment` 大小相同、符号相反，比纯 PBE 更大；
- 占据矩阵本征值接近 0 或 1（Ni²⁺ d⁸：多数自旋 5 个 ≈1，少数自旋 3 个 ≈1、2 个 ≈0 的 e_g）；
- 绝缘体的 `highest occupied, lowest unoccupied level` 之差（近似带隙）明显大于纯 PBE。若使用 smearing，输出的是 Fermi 能，带隙需从 DOS 读（[第6章](06-band-dos-pdos.md)）。

### 4.3 关于"不同 U 的总能不可比较"

U 改变了能量泛函本身，**不同 U 的总能不能直接相减**（例如画 E–U 曲线找"最佳 U"是错误做法）；形成能、反应能必须在**同一 U**下计算。混合 GGA/GGA+U 的能量修正方案（如 Materials Project 的 MP2020）是针对 VASP 设置拟合的，不适用于 QE。

### 4.4 多个自洽解与轨道序

DFT+U 能量面常有多个局部极小（不同轨道序、不同自旋态）。Hubbard_input.pdf 也提醒，磁性体系的 DFT+U 计算会依据初始占据收敛到不同基态。可用 `starting_ns_eigenvalue(m, ispin, ityp)` 在第一步覆盖占据矩阵本征值以引导：

```fortran
&SYSTEM
  ...
  ! push the 5th eigenvalue (in ascending order) of the spin-down
  ! occupation matrix of species 1 to zero in the first iteration
  starting_ns_eigenvalue(5,2,1) = 0.0
/
```

做法是多试几个起点，取能量最低者，并检查各解的占据矩阵。

---

## 5. 实战二：hp.x 计算 NiO 的 U

### 5.1 为什么磁性绝缘体要两步 SCF

官方 `HP/examples/example02/README` 解释：未知磁态时必须用 smearing（把体系当作"假金属"）才能让磁化自由演化；但 DFPT 对金属在 q=0 有一项正比于费米能级 DOS 倒数的修正，假金属的 DOS 极小，会导致数值发散。官方给出两种方案：

| 方案 | 步骤 | 说明 |
|---|---|---|
| 1（本例） | ① 小 smearing 假金属 SCF → ② 用整数 `tot_magnetization`、`occupations='fixed'`、相同 `nbnd`，从 ① 的势与波函数重启 SCF → ③ hp.x | 更干净 |
| 2 | ① 之后直接 hp.x，并设 `no_metq0 = .true.` | 跳过 q=0 金属响应项；仅适用于磁性绝缘体 |

第 ② 步必须：

- `tot_magnetization` 取 ① 中总磁化四舍五入后的整数（AFM NiO 为 0）；
- `nbnd` 等于 ① 的 `number of Kohn-Sham states`（smearing 计算默认多 20% 的带，不一致会出问题）；
- `startingpot = 'file'`、`startingwfc = 'file'`，否则仅凭约束磁化很难收敛到正确的磁态。

### 5.2 hp.x 输入（`nio.hp.in`）

```fortran
&INPUTHP
  prefix       = 'nio_lr'
  outdir       = './tmp_lr'
  nq1 = 2, nq2 = 2, nq3 = 2
  conv_thr_chi = 1.0d-8
  iverbosity   = 2
/
```

| 关键词 | 推荐值 | 说明（`INPUT_HP`） |
|---|---|---|
| `prefix`, `outdir` | 与 pw.x 一致 | 默认 `'pwscf'` 与 `ESPRESSO_TMPDIR` |
| `nq1, nq2, nq3` | 2×2×2 起，测试收敛 | 无偏移的 MP q 网格，默认 1,1,1 |
| `conv_thr_chi` | 1.0d-8 | χ 的收敛阈值（默认 1.0d-5，官方示例用 1.0d-8） |
| `iverbosity` | 2 | 2 会额外打印 χ₀、χ 及其逆矩阵 |
| `find_atpert` | 默认 1 | 1 按未扰动占据判断不等价 Hubbard 原子；3 按对称性；4 扰动全部（最贵） |
| `skip_type(i)` / `equiv_type(i)` | AFM 可用 | 跳过某类型并声明与另一类型等价；不能与 `find_atpert=1` 同用 |
| `perturb_only_atom(i)`, `start_q`, `last_q`, `sum_pertq`, `compute_hp` | 拆分大计算 | 按原子/按 q 点并行后再汇总 |
| `no_metq0` | `.false.` | 方案 2 时设 `.true.` |
| `num_neigh` | 6 | U+V：写入 `HUBBARD.dat` 的近邻数 |
| `lmin` | 2 | U+V：只为 l ≥ lmin 的 Hubbard 原子写 V |
| `rmax` | 100 Bohr | U+V：搜索近邻的最大距离 |
| `max_seconds` | 按队列时长 | 平滑停止 |

### 5.3 运行（脚本 `run_hp.sh`）

```bash
cd QE/examples/08-dftu-hp
bash run_hp.sh                         # one cycle from PBE (U_in = 0.0001 eV)
NCYCLE=5 TOL=0.01 bash run_hp.sh       # self-consistent U
PW_CMD="srun pw.x -nk 4" HP_CMD="srun hp.x -nk 4" NCYCLE=5 bash run_hp.sh
```

脚本的核心逻辑（与官方 example02 一致）：

```bash
nbnd=$(grep 'number of Kohn-Sham states' step1.out | awk '{print $5}')
totmag=$(grep 'total magnetization' step1.out | tail -1 | awk '{print $4}')
# round to integer, fill @NBND@ / @TOTMAG@ / @U_IN@ in nio.lr_step2.in
hp.x -in nio.hp.in > hp.out          # -> nio_lr.Hubbard_parameters.dat
```

> [!NOTE]
> 第 1 步的 `HUBBARD` 卡用 `U ... 0.0001`：hp.x 要求 pw.x 计算是 DFT+Hubbard 计算（否则报 `The HP code can be used only on top of DFT+Hubbard (i.e. when the HUBBARD card is used in pw.x)`），极小 U 等价于从 PBE 出发。官方示例同样这样做。

### 5.4 读结果

hp.x 把结果写在当前目录的 `<prefix>.Hubbard_parameters.dat`：

```bash
grep -A6 'Hubbard U parameters' nio_lr.Hubbard_parameters.dat
cat u_history.dat
```

（示意输出，数值仅供参考；QE 7.5 的表头含 `manifold` 列）

```text
                                 Hubbard U parameters:

       site n.  type  label  spin  new_type  new_label  manifold  Hubbard U (eV)
         1        1    Ni1     1      1         Ni1        3d        x.xxxx
         2        2    Ni2    -1      1         Ni1        3d        x.xxxx
```

- `spin` 列 1 / −1 表示两个子晶格自旋相反；`new_type`/`new_label` 显示 hp.x 认定 Ni2 与 Ni1 等价，二者 U 必须相同；
- `iverbosity = 2` 时还打印 `chi0 matrix`、`chi matrix`、各自逆矩阵和 `Hubbard matrix`（对角元即 U，非对角元为原子间相互作用）；
- 官方 example02 参考输出（PBEsol、pslibrary 超软赝势、k 与 q 均为 2×2×2）的 U 约 7.9 eV；本例换成 PBE + GBRV/pslibrary 0.1 组合，**数值会不同**，不要拿它当"正确答案"。

自洽判据：`u_history.dat` 中相邻两轮 |U_out − U_in| < `TOL`。实践中常同时在每一轮做结构弛豫（"自洽 U + 结构"），这一类完整流程可用 aiida-hubbard 自动化（见参考资料）。

### 5.5 收敛测试

| 参数 | 测试方法 |
|---|---|
| q 网格 | 2×2×2 → 3×3×3 → 4×4×4，U 变化 < 0.1 eV |
| k 网格（SCF） | k 网格必须与 q 网格相容（k 为 q 的整数倍最稳妥），一起加密 |
| `ecutwfc` | 比普通 SCF 稍高更保险，U 对截断通常不敏感 |
| `conv_thr`（SCF） | 第 ② 步 1.0d-14 量级，太松会污染响应 |
| `conv_thr_chi` | 1.0d-8 起 |

---

## 6. 实战三：LiCoO2 的 DFT+U+V

### 6.1 输入

`licoo2.uv.scf.in`：菱方原胞（`ibrav=5`），HUBBARD 卡同时给 U 与一条 V（索引 1 2 只是占位，表示需要 Co–O 的 V）：

```fortran
ATOMIC_POSITIONS crystal
  Co  0.0000000  0.0000000  0.0000000
  O   0.2604885  0.2604885  0.2604885
  O   0.7395115  0.7395115  0.7395115
  Li  0.5000000  0.5000000  0.5000000
K_POINTS automatic
  4 4 4 0 0 0
HUBBARD {ortho-atomic}
U Co-3d 0.0001
V Co-3d O-2p 1 2 0.0001
```

LiCoO2 是非磁绝缘体（低自旋 Co³⁺），不需要两步 SCF。`licoo2.hp.in` 额外给 `num_neigh = 6`、`lmin = 2`。

### 6.2 运行与 `HUBBARD.dat`

```bash
bash run_uv.sh
cat HUBBARD.dat
grep -A14 'Hubbard V parameters' licoo2.Hubbard_parameters.dat
```

hp.x 在 DFT+U+V 模式下额外写出 `HUBBARD.dat`（文件首行注释 `# Copy this data in the pw.x input file for DFT+Hubbard calculations`），可直接替换下一轮 pw.x 输入中的 HUBBARD 卡。格式（示意输出，数值仅供参考）：

```text
# Copy this data in the pw.x input file for DFT+Hubbard calculations
HUBBARD {ortho-atomic}
V     Co-3d     Co-3d    1     1   x.xxxx
V     Co-3d      O-2p    1    19   y.yyyy
V     Co-3d      O-2p    1    46   y.yyyy
...   (six Co-O nearest neighbours in total)
```

第一行 `V Co-3d Co-3d 1 1` 就是 Co 的在位 U（等价写法）；其余 6 行对应八面体配位的 6 个 O 近邻。`Hubbard_parameters.dat` 中的 V 表还列出更远近邻及其距离（Bohr），可用于判断 V 随距离的衰减。

> [!IMPORTANT]
> Hubbard_input.pdf：hp.x 目前每个原子类型只支持**一个** Hubbard 流形；跨流形相互作用（standard–background）等高级写法 pw.x 能读，但 hp.x 算不了。另外 hp.x 不支持带 background 的 DFT+U+V。

### 6.3 自洽 U+V

把 `HUBBARD.dat` 写回 pw.x 输入 → 重新 SCF（并可做 relax）→ 再跑 hp.x，直到 U 与 V 都收敛。注意 V 的索引在结构改变后可能变化，每一轮都应使用 hp.x 新写出的 `HUBBARD.dat`，不要只改数值。

---

## 7. hp.x 的限制速查

以下限制直接来自 QE 7.5 `HP/src/hp_readin.f90` 的检查：

| 条件 | 报错信息 |
|---|---|
| pw.x 没有 HUBBARD 卡 | `The HP code can be used only on top of DFT+Hubbard (i.e. when the HUBBARD card is used in pw.x)` |
| Liechtenstein（写了 `J`） | `The HP code does not support the Liechtenstein formulation of DFT+U` |
| 投影不是 atomic/ortho-atomic | `The HP code for this Hubbard_projectors type is not implemented` |
| `J0` ≠ 0 | `Hubbard_J0 /= 0 is not allowed.` |
| 前一步是 vc-relax | `The HP code is not working after vc-relax` |
| 约束磁化 | `The HP code with constrained magnetization is not yet available` |
| `tot_magnetization` + smearing | `The HP code with two Fermi energies is not available for metals` |
| tetrahedra | `HP with tetrahedra is not supported` |
| `K_POINTS gamma` | `Cannot start from pw.x data file using Gamma-point tricks` |
| meta-GGA / 杂化泛函 | `The HP code with meta-GGA functionals is not yet available` / `... hybrid functionals ...` |
| `assume_isolated='2D'` | `The HP code does not support the 2D cutoff` |
| 带组并行 `-nb > 1` | `band parallelization is not implemented in HP` |
| 一个 species 写了两个及以上 Hubbard 流形 | `Two (or more) Hubbard channels per atomic type is not implemented` |
| `V` 行用了 background 通道 | `The HP code does not support DFT+U+V with the background` |
| `occupations='from_input'` | `The HP code with arbitrary occupations not tested` |
| `tqr = .true.` | `The HP code with Q in real space is not supported` |

---

## 参数选择与收敛建议

| 项目 | 建议 |
|---|---|
| 投影 | 默认 `ortho-atomic`；需要与文献比较时，确认对方的投影与赝势 |
| U 的来源 | 优先 hp.x 线性响应（必要时自洽）；经验 U 只用于定性趋势，并在论文中注明来源 |
| 自洽 U | `|U_out − U_in|` < 0.05–0.1 eV；结构敏感的体系同时迭代结构 |
| q / k 网格 | q 从 2×2×2 起；k 与 q 相容（整数倍）；U 变化 < 0.1 eV |
| 磁性绝缘体 | 两步 SCF（官方方案 1），或 `no_metq0 = .true.`（方案 2） |
| 原子顺序 | Hubbard 原子在 `ATOMIC_POSITIONS` 中排在最前 |
| 能量比较 | 只比较同一 U（同一 U+V 集合）下的能量 |
| 金属 | 弱关联金属一般不加 U；确需时 hp.x 支持金属（`lgauss`）但收敛更慢 |
| 闭壳层流形 | 3d¹⁰、4f¹⁴ 等不适合线性响应 U |

---

## 常见问题与报错

| 报错 / 现象 | 原因 | 解决 |
|---|---|---|
| `DFT+Hubbard input syntax has changed since v7.1`（前面有 `WARNING!!! The input parameter lda_plus_u is obsolete.` 等） | 使用了旧关键词 | 按 §3.3 改写为 HUBBARD 卡 |
| HUBBARD 卡读取出错 / 流形找不到 | 标签与 `ATOMIC_SPECIES` 不一致（如写 `Ni-3d` 但 species 是 `Ni1`）；流形不在赝势中；赝势 `PP_CHI` 缺 label | 逐一核对标签；`grep label *.UPF` 检查轨道；按 Hubbard_input.pdf 第 5 节补 label |
| U 加上了但带隙没变化 | U 加在了错误的流形或 species（AFM 只写了 `Ni1`） | 每个 species 都要写 |
| hp.x 在初始化阶段停止（检查原子顺序） | Hubbard 原子没排在 `ATOMIC_POSITIONS` 最前（官方 HP examples README 的要求） | 调整 `ATOMIC_POSITIONS` 顺序后重跑 pw.x |
| hp.x 响应发散 / 不收敛（磁性绝缘体） | 假金属 DOS 极小导致 q=0 金属项发散 | 两步 SCF（方案 1）或 `no_metq0 = .true.` |
| `The HP code with two Fermi energies is not available for metals` | 第 ② 步还保留 smearing 同时设了 `tot_magnetization` | 第 ② 步改 `occupations='fixed'` |
| 第 ② 步 SCF 不收敛 | 未从第 ① 步重启，或 `nbnd` 不一致 | `startingpot/startingwfc='file'`，`nbnd` 取第 ① 步的值 |
| `fixed occupations requires integer tot_magnetization` | 第 ② 步直接抄了第 ① 步的非整数总磁化 | 四舍五入为整数（`run_hp.sh` 已自动处理） |
| `The HP code is not working after vc-relax` | 直接对 vc-relax 的 save 目录运行 hp.x | 用优化后的结构重新做 scf |
| 不同次 DFT+U 计算能量差几十 meV | 收敛到不同轨道序 | 多起点 + `starting_ns_eigenvalue`，取最低能量解 |
| 线性响应 U 与文献差 1–2 eV | 投影、赝势、泛函、q 网格不同 | 逐项对照；只在相同设置下比较 |

---

## 练习

1. 用 `nio.dftu.scf.in` 把 U 从 0 扫到 8 eV（步长 1 eV），记录 Ni 磁矩和带隙随 U 的变化；思考为什么不能用这些总能画 E–U 曲线。
2. 把 `HUBBARD {ortho-atomic}` 改成 `HUBBARD {atomic}`，用 hp.x 重新计算 U，比较两种投影下的 U、磁矩和带隙。
3. 运行 `NCYCLE=5 bash run_hp.sh`，画出 U_in–U_out 收敛曲线；再把 q 网格从 2×2×2 改为 3×3×3（k 网格同步加密），评估 U 的收敛。
4. 对 NiO 试用官方"方案 2"：只做第 ① 步，然后在 `nio.hp.in` 中加 `no_metq0 = .true.`，比较与方案 1 的 U。
5. 对 LiCoO2 完成 2–3 轮 U+V 自洽，每轮用新 `HUBBARD.dat`；比较 DFT+U 与 DFT+U+V 的带隙。
6. （阅读）对比 [VASP 第7章](../../VASP/tutorials/07-magnetism-u-soc.md) 中 `LDAUTYPE = 3` 的超胞线性响应流程，列出与 hp.x 的三点主要差异。

---

## 参考资料

**官方文档**

- [New DFT+Hubbard input（Hubbard_input.pdf）](https://www.quantum-espresso.org/Doc/user_guide_PDF/Hubbard_input.pdf)：HUBBARD 卡全部语法、旧→新对照、赝势 label、非共线 DFT+U
- [INPUT_PW](https://www.quantum-espresso.org/Doc/INPUT_PW.html)：`HUBBARD` 卡、`Hubbard_occ`、`Hubbard_beta`、`starting_ns_eigenvalue`
- [INPUT_HP](https://www.quantum-espresso.org/Doc/INPUT_HP.html)：hp.x 全部关键词
- QE 源码 `HP/examples`（example01–example10：LiCoO2、NiO、CrI3、Ni、Ni2MnGa 拆分并行、2D 体系、U+V 等）：<https://gitlab.com/QEF/q-e>
- [QE release notes](https://gitlab.com/QEF/q-e/-/blob/develop/Doc/release-notes)：7.1 新语法、7.3 非共线 DFT+U、7.5 轨道分辨 DFT+U 与 Wannier 投影

**课程与教程（只链接）**

- [materialscloud-org/QuantumESPRESSO-school-2023](https://github.com/materialscloud-org/QuantumESPRESSO-school-2023)：EPFL/Materials Cloud 高级 QE 学校，含 DFT+U+V、hp.x 与 Koopmans 上机，是本章最直接的进阶材料
- [materialscloud-org/hubbard-koopmans-2022](https://github.com/materialscloud-org/hubbard-koopmans-2022)：2022 版 Hubbard/Koopmans 教程
- [pranabdas/espresso](https://pranabdas.github.io/espresso/)：入门教程网站，DFT+U 部分已使用新 HUBBARD 语法
- [Imperial College MSE404-MM](https://imperialcollegelondon.github.io/MSE404-MM/)：QE 实验课程，适合补自旋极化基础
- [Cornell PARADIM Computational Tutorials](https://www.paradim.org/toolbox/computation/tutorials)：QE 基础操作（版本较旧，未涉及新 HUBBARD 语法）

**相关仓库（只链接）**

- [aiidateam/aiida-hubbard](https://github.com/aiidateam/aiida-hubbard)：基于 AiiDA 的自洽 U/V + 结构弛豫自动工作流
- [aiidateam/aiida-quantumespresso](https://github.com/aiidateam/aiida-quantumespresso)：pw.x 工作流基础
- [quantumNerd/Quantum-Espresso-Tutorial-2019-Projects](https://github.com/quantumNerd/Quantum-Espresso-Tutorial-2019-Projects)：其 DFT+U 输入是 7.1 以前的旧语法，阅读时请对照 §3.3 翻译

**关键文献**

- V. I. Anisimov, J. Zaanen, O. K. Andersen, *Phys. Rev. B* 44, 943 (1991). DOI: [10.1103/PhysRevB.44.943](https://doi.org/10.1103/PhysRevB.44.943)
- S. L. Dudarev et al., *Phys. Rev. B* 57, 1505 (1998). DOI: [10.1103/PhysRevB.57.1505](https://doi.org/10.1103/PhysRevB.57.1505)
- M. Cococcioni, S. de Gironcoli, *Phys. Rev. B* 71, 035105 (2005). DOI: [10.1103/PhysRevB.71.035105](https://doi.org/10.1103/PhysRevB.71.035105)
- V. L. Campo Jr, M. Cococcioni, *J. Phys.: Condens. Matter* 22, 055602 (2010)（DFT+U+V）. DOI: [10.1088/0953-8984/22/5/055602](https://doi.org/10.1088/0953-8984/22/5/055602)
- B. Himmetoglu, R. M. Wentzcovitch, M. Cococcioni, *Phys. Rev. B* 84, 115108 (2011)（DFT+U+J0）. DOI: [10.1103/PhysRevB.84.115108](https://doi.org/10.1103/PhysRevB.84.115108)
- H. J. Kulik et al., *Phys. Rev. Lett.* 97, 103001 (2006)（自洽 U）. DOI: [10.1103/PhysRevLett.97.103001](https://doi.org/10.1103/PhysRevLett.97.103001)
- I. Timrov, N. Marzari, M. Cococcioni, *Phys. Rev. B* 98, 085127 (2018)（DFPT 求 U）. DOI: [10.1103/PhysRevB.98.085127](https://doi.org/10.1103/PhysRevB.98.085127)
- I. Timrov, N. Marzari, M. Cococcioni, *Phys. Rev. B* 103, 045141 (2021)（USPP/PAW 下的自洽 U）. DOI: [10.1103/PhysRevB.103.045141](https://doi.org/10.1103/PhysRevB.103.045141)
- I. Timrov, N. Marzari, M. Cococcioni, *Comput. Phys. Commun.* 279, 108455 (2022)（HP 程序）. DOI: [10.1016/j.cpc.2022.108455](https://doi.org/10.1016/j.cpc.2022.108455)
- I. Timrov et al., *Phys. Rev. B* 102, 235159 (2020)（ortho-atomic 投影下的 Pulay 力）. DOI: [10.1103/PhysRevB.102.235159](https://doi.org/10.1103/PhysRevB.102.235159)
- L. Binci, N. Marzari, *Phys. Rev. B* 108, 115157 (2023)（非共线 DFT+U）. DOI: [10.1103/PhysRevB.108.115157](https://doi.org/10.1103/PhysRevB.108.115157)
- E. Macke et al., *J. Chem. Theory Comput.* 20, 4824 (2024)（轨道分辨 DFT+U）. DOI: [10.1021/acs.jctc.3c01403](https://doi.org/10.1021/acs.jctc.3c01403)
- L. Bastonero et al., *npj Comput. Mater.* 11, 183 (2025)（aiida-hubbard 自动化工作流）. DOI: [10.1038/s41524-025-01685-4](https://doi.org/10.1038/s41524-025-01685-4)

---

[← 上一章](07-magnetism-soc.md) | [返回目录](../README.md) | [下一章 →](09-vdw-surface-workfunction.md)
