# 第7章 磁性、非共线与自旋轨道耦合（Spin Polarization, Noncollinear Magnetism & SOC）

本章讲 Quantum ESPRESSO（pw.x）中的磁性计算：共线自旋极化（`nspin=2`、`starting_magnetization`、`tot_magnetization`）、铁磁/反铁磁构型（把同一元素拆成多个 species）、非共线磁性（`noncolin`、`angle1`/`angle2`）、约束磁化（`constrained_magnetization`、`lambda`）、自旋轨道耦合（`lspinorb` + 全相对论赝势）以及用力定理（`lforcet`）估算磁晶各向异性能（MAE）。算例覆盖 bcc Fe、fcc Ni、AFM NiO 和 GaAs 的 SOC 价带劈裂（Quantum ESPRESSO magnetism / noncollinear / spin-orbit coupling tutorial）。

> [!NOTE]
> 全部输入在示例目录 [07-magnetism-soc](../examples/07-magnetism-soc/README.md)。关键词已对照 QE 7.5 的 `INPUT_PW` 与源码核实；本仓库**不运行任何计算**，文中"示意输出"只演示格式，数值仅供参考。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 正确写出 `nspin = 2` 与 `starting_magnetization(i)`，理解"不给初始磁化就得到非磁解"；
- 从输出中读出 `total magnetization`、`absolute magnetization` 和 `Magnetic moment per site`；
- 用"拆 species"的方式建立反铁磁（AFM）构型，比较 FM / AFM / NM 能量；
- 知道什么时候用 `tot_magnetization` 固定总磁矩（例如磁性绝缘体、O₂ 三重态）；
- 做非共线计算（`noncolin = .true.`）并用 `angle1`/`angle2` 指定磁矩方向；
- 用全相对论（FR）赝势 + `lspinorb = .true.` 计算 SOC 能带，读懂 GaAs 价带顶的 Γ8/Γ7 劈裂；
- 理解力定理 MAE 的流程（LSDA 自洽 → `lforcet` 非自洽 SOC）及其精度要求。

| 场景 | 推荐做法 | 本章小节 |
|---|---|---|
| Fe/Co/Ni 等金属、合金 | `nspin=2` + `starting_magnetization` + smearing | §3 |
| 反铁磁氧化物（NiO、MnO） | 拆 species（`Ni1`/`Ni2`）+ 相反初始磁化；通常还需 DFT+U（第8章） | §4 |
| 已知总自旋的分子/绝缘体 | `tot_magnetization`（可用 `occupations='fixed'`） | §3.4 |
| 磁矩方向、自旋螺旋、阻挫 | `noncolin = .true.` + `angle1/angle2` | §5 |
| 固定磁矩方向/大小求能量 | `constrained_magnetization` + `lambda` | §5.3 |
| 重元素能带、Rashba、拓扑 | `noncolin` + `lspinorb` + FR 赝势 | §6 |
| 易磁化轴、MAE | 力定理：LSDA → `lforcet` 非自洽 | §7 |

> [!TIP]
> 本章默认你已掌握 [第2章](02-input-anatomy.md) 的输入结构、[第4章](04-convergence.md) 的收敛测试和 [第6章](06-band-dos-pdos.md) 的能带流程。磁性计算对 k 点和 smearing 更敏感，收敛测试要在**自旋极化**下重做。

### 与 VASP 的对应关系

从 VASP 转过来的读者可以对照下表（概念对应，并非逐字等价）；VASP 侧的完整讲解见 [VASP 第7章](../../VASP/tutorials/07-magnetism-u-soc.md)。

| VASP | QE（pw.x） | 说明 |
|---|---|---|
| `ISPIN = 2` | `nspin = 2` | 共线自旋极化 |
| `MAGMOM`（每原子一个数） | `starting_magnetization(i)`（**每个 species 一个数**） | QE 按 species 给值，AFM 需拆 species |
| `NUPDOWN` | `tot_magnetization` | 固定 N↑ − N↓ |
| `LNONCOLLINEAR = .TRUE.` | `noncolin = .true.` | 两分量旋量 |
| `MAGMOM`（每原子 3 个数） | `starting_magnetization` + `angle1` + `angle2` | 极角/方位角，单位为度 |
| `LSORBIT = .TRUE.` | `lspinorb = .true.` + FR 赝势 | QE 的 SOC 来自赝势中的 j 相关投影子 |
| `SAXIS` | `angle1`/`angle2` 直接给方向 | QE 中没有单独的量子化轴参数 |
| `I_CONSTRAINED_M` + `LAMBDA` | `constrained_magnetization` + `lambda` | 罚函数形式不同，`lambda` 数值不可互换 |
| `ICHARG = 11` + SOC（力定理） | `lforcet = .true.`（`calculation='nscf'`） | 均为非自洽 SOC |

---

## 2. 原理速览

### 2.1 共线自旋极化（LSDA / σ-GGA）

把电子密度分成两种自旋分量：

$$
n(\mathbf r)=n_\uparrow(\mathbf r)+n_\downarrow(\mathbf r),\qquad
m(\mathbf r)=n_\uparrow(\mathbf r)-n_\downarrow(\mathbf r)
$$

两套 Kohn–Sham 方程分别求解，磁化只有 z 分量（`nspin = 2`）。输出中的两个量：

- **total magnetization** $=\int m(\mathbf r)\,d\mathbf r$，单位 μB/cell；
- **absolute magnetization** $=\int |m(\mathbf r)|\,d\mathbf r$。AFM 态 total ≈ 0 而 absolute 明显不为零，这是判断 AFM 是否成立的最快指标。

自旋极化 DFT 的能量面常有多个自洽解（FM、各种 AFM、高/低自旋、NM），**最终落入哪个解取决于初始磁化和对称性**。因此磁性研究的标准做法是从多个初始构型出发、比较能量。

### 2.2 `starting_magnetization` 的两种含义

`INPUT_PW` 对 `starting_magnetization(i)` 的定义（QE ≥ 7.3）：

- 绝对值 **< 1**：解释为"每个价电子的自旋极化"，范围 −1 到 1；
- 绝对值 **≥ 1**：解释为该位点的磁矩（μB）——源码中会除以该 species 的价电子数 `zv`。QE 7.2 及更早版本只接受第一种写法。

> [!WARNING]
> 不同赝势的价电子数不同（例如 SSSP 中的 Fe PAW 为 16 个价电子），所以 `starting_magnetization = 0.5` 对应的初始磁矩随赝势而变。它只是**初始猜测**，最终磁矩由自洽决定；但初始值太小可能收敛到 NM 解，太大可能落入高自旋亚稳态。

### 2.3 非共线磁性

`noncolin = .true.` 时每个 Kohn–Sham 态是两分量旋量，密度矩阵为 2×2，磁化 $\mathbf m(\mathbf r)$ 是矢量（内部 `nspin = 4`，但**不要**在输入中写 `nspin = 4`）。初始磁矩方向由每个 species 的两个角度给出：

- `angle1(i)`：磁矩与 z 轴的夹角（度）；
- `angle2(i)`：磁矩在 xy 面投影与 x 轴的夹角（度）。

没有 SOC 时，自旋空间与实空间不耦合，**整体旋转所有磁矩不改变能量**——这是一个很好的自检。

### 2.4 自旋轨道耦合

SOC 是相对论修正，近似写作 $H_{\rm SOC}\propto \xi(r)\,\mathbf L\cdot\boldsymbol\sigma$。在 QE 中 SOC 通过**全相对论（fully relativistic, FR）赝势**中按总角动量 j 区分的投影子引入（Dal Corso 等的 US/PAW 实现，见参考资料）。因此：

- `lspinorb = .true.` 必须配合 FR 赝势（UPF 头部 `relativistic="full"`、`has_so="true"`）；
- 若全部赝势都不含 SO 信息，pw.x 报错 `Spin-orbit calculations require at least one spin-orbit pseudo`；
- 反过来，在**不开** `lspinorb` 的计算中使用 FR 赝势时，pw.x 会把它平均成标量相对论赝势（源码 `average_pp`），因此同一个 FR 赝势可以同时用于 LSDA 与 SOC 两步——力定理流程正是利用这一点。

SOC 有两种常见情形：

| 情形 | 输入 | 输出提示 |
|---|---|---|
| 非磁 + SOC（GaAs、Pt、Bi₂Se₃） | `noncolin`、`lspinorb`，**不给** `starting_magnetization` | `Non magnetic calculation with spin-orbit` |
| 磁性 + SOC（Fe、CrI₃） | 另外给 `starting_magnetization`、`angle1/angle2` | `Noncollinear calculation with spin-orbit` |

`INPUT_PW` 明确说明：非共线/SOC 情形下所有原子的 `starting_magnetization` 都为 0 时，程序施加时间反演对称性，**磁化根本不计算、恒为零**。

### 2.5 磁晶各向异性与力定理

SOC 把自旋方向与晶格耦合，磁化沿不同晶向时能量不同，差值即 MAE。完全自洽地比较两个方向需要极高精度（MAE 常为 μeV–meV/atom）。力定理近似：固定共线自洽得到的密度，只改变磁化方向做一次非自洽 SOC 计算，两个方向的能量差 ≈ 占据态本征值之和（band energy）之差。

---

## 3. 共线自旋极化：bcc Fe 与 fcc Ni

### 3.1 输入

`fe.scf.in`（完整文件见示例目录）：

```fortran
&CONTROL
  calculation = 'scf'
  prefix      = 'fe'
  outdir      = './tmp'
  pseudo_dir  = './pseudo'
  tprnfor     = .true.
  tstress     = .true.
/
&SYSTEM
  ibrav       = 3                ! bcc, primitive cell with 1 atom
  celldm(1)   = 5.42             ! a = 2.87 Angstrom (exp.); use your relaxed value
  nat         = 1
  ntyp        = 1
  ecutwfc     = 90.0
  ecutrho     = 1080.0
  occupations = 'smearing'
  smearing    = 'mv'             ! Marzari-Vanderbilt cold smearing
  degauss     = 0.02             ! Ry; test together with the k mesh
  nspin       = 2                ! collinear spin polarization
  starting_magnetization(1) = 0.5  ! per valence electron, range [-1,1]
/
&ELECTRONS
  conv_thr    = 1.0d-9
  mixing_beta = 0.3
/
ATOMIC_SPECIES
  Fe  55.845  Fe.pbe-spn-kjpaw_psl.0.2.1.UPF
ATOMIC_POSITIONS crystal
  Fe  0.0  0.0  0.0
K_POINTS automatic
  14 14 14 0 0 0
```

截断能取自 SSSP 1.3 PBE Efficiency 对 Fe 的推荐值（90 / 1080 Ry）。`fe.nm.scf.in` 与之相同，只是去掉了 `nspin` 和 `starting_magnetization`，用于比较磁性稳定化能 $E_{\rm FM}-E_{\rm NM}$。

### 3.2 关键词表

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `nspin` | 2 | 共线自旋极化；默认 1（非极化） |
| `starting_magnetization(i)` | 0.2–0.6（磁性 species），其余 0 | 默认 0。**至少一个 species 非零**，否则 SCF 得到 NM 态 |
| `occupations` | `'smearing'` | 默认 `'fixed'`；`nspin=2` 时若既不用 smearing/tetrahedra 又不给 `tot_magnetization`，pw.x 报 `fixed occupations and lsda need tot_magnetization` |
| `smearing` / `degauss` | `'mv'` / 0.01–0.02 Ry | 金属磁体；与 k 网格一起收敛 |
| `mixing_beta` | 0.2–0.3 | 磁性体系比默认 0.7 更稳 |
| `report` | 默认 −1 | 控制局域磁矩的打印频率：0 从不，−1 开始与收敛时，N 每 N 步 |
| `tot_magnetization` | 不设 | 固定总磁化（见 §3.4） |

### 3.3 运行与检查

```bash
cd QE/examples/07-magnetism-soc
mpirun -np 4 pw.x -nk 4 -in fe.scf.in > fe.scf.out
mpirun -np 4 pw.x -nk 4 -in fe.nm.scf.in > fe.nm.scf.out

grep '^!' fe.scf.out fe.nm.scf.out
grep -e 'total magnetization' -e 'absolute magnetization' fe.scf.out | tail -2
grep -A2 'Magnetic moment per site' fe.scf.out | tail -2
grep 'the Fermi energy' fe.scf.out
```

（示意输出，数值仅供参考）

```text
     total magnetization       =     2.25 Bohr mag/cell
     absolute magnetization    =     2.40 Bohr mag/cell

     Magnetic moment per site  (integrated on atomic sphere of radius R)
     atom    1 (R=0.xxx)  charge=  x.xxxx  magn=  x.xxxx
```

成功判据：

1. 每个 SCF 步都打印 `total magnetization`，最后几步稳定；`convergence has been achieved` 与 `JOB DONE` 均出现；
2. `total magnetization` 与文献中 PBE 的 bcc Fe 磁矩（约 2.2 μB）同一量级；
3. $E_{\rm FM}<E_{\rm NM}$，两者用同样的截断与 k 点；
4. `Magnetic moment per site` 是在半径 R 的原子球内积分得到的，**依赖 R**，只能作定性比较；总磁化才是严格量。

> [!NOTE]
> 自旋极化金属的输出会给 `the Fermi energy is ... ev`（一个费米能）；若使用 `tot_magnetization`，自旋向上/向下各有一个费米能，输出为 `the spin up/dw Fermi energies are ...`。

fcc Ni（`ni.scf.in`）完全同理，只把 `ibrav=2`、`celldm(1)=6.65`，赝势换成 SSSP 的 `ni_pbe_v1.4.uspp.F.UPF`（45 / 360 Ry），`starting_magnetization(1) = 0.3`。Ni 的磁矩约 0.6 μB，对 k 点和 `degauss` 比 Fe 更敏感，建议至少 16×16×16。

### 3.4 固定总磁化：`tot_magnetization`

`tot_magnetization` 规定 N↑ − N↓（单位 μB/cell），用于：

- 分子（如 O₂ 三重态 `tot_magnetization = 2`）、已知自旋的缺陷；
- 磁性绝缘体：先用 smearing 得到磁态，再取整数总磁矩改用 `occupations = 'fixed'`（第8章 hp.x 的两步 SCF 正是这样做的）；
- 固定自旋矩扫描（fixed spin moment），画 E(M) 曲线。

规则（来自 `INPUT_PW` 与源码检查）：

- 只能用于 `nspin = 2`，否则报 `tot_magnetization requires nspin=2`；
- 文档注明：**固定 `tot_magnetization` 时不要再给 `starting_magnetization`**；
- 非共线情形改用 `constrained_magnetization = 'total'`（§5.3）。

```fortran
&SYSTEM
  ...
  nspin             = 2
  occupations       = 'fixed'
  tot_magnetization = 2        ! e.g. O2 triplet in a box
/
```

---

## 4. 反铁磁：拆分 species（NiO AFM-II）

### 4.1 为什么必须拆 species

QE 的 `starting_magnetization` 按 **species** 而不是按原子给值。要让两个 Ni 磁矩相反，就必须在 `ATOMIC_SPECIES` 中定义 `Ni1`、`Ni2` 两个 species（可指向同一个 UPF 文件），这也是 `INPUT_PW` 推荐的做法。

NiO 是岩盐结构，实验基态为 AFM-II：(111) 面内铁磁、相邻 (111) 面反平行。最小磁性胞是含 2 个 Ni、2 个 O 的菱方胞。

### 4.2 输入（`nio_afm.scf.in`）

```fortran
&CONTROL
  calculation = 'scf'
  prefix      = 'nio'
  outdir      = './tmp'
  pseudo_dir  = './pseudo'
/
&SYSTEM
  ibrav       = 0
  nat         = 4
  ntyp        = 3
  ecutwfc     = 50.0
  ecutrho     = 400.0
  occupations = 'smearing'       ! safe start for a magnetic system;
  smearing    = 'gaussian'       ! switch to 'fixed' + tot_magnetization
  degauss     = 0.005            ! once the magnetic state is known
  nspin       = 2
  starting_magnetization(1) =  0.5   ! Ni1 spin up
  starting_magnetization(2) = -0.5   ! Ni2 spin down
  starting_magnetization(3) =  0.0   ! O
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
```

截断取 Ni（45/360 Ry）与 O（50/400 Ry）SSSP 推荐值中的较大者。

### 4.3 检查 AFM 是否成立

```bash
mpirun -np 4 pw.x -nk 2 -in nio_afm.scf.in > nio_afm.scf.out
grep -e 'total magnetization' -e 'absolute magnetization' nio_afm.scf.out | tail -2
grep -A5 'Magnetic moment per site' nio_afm.scf.out | tail -5
```

（示意输出，数值仅供参考）

```text
     total magnetization       =     0.00 Bohr mag/cell
     absolute magnetization    =     3.xx Bohr mag/cell
```

- `total` ≈ 0、`absolute` 明显大于 0，且 atom 1 与 atom 2 的 `magn` 大小相等、符号相反 → AFM 成立；
- 若两个 Ni 磁矩同号，说明 SCF 跑到了 FM 解：检查 `starting_magnetization` 的符号与 species 对应关系；
- 纯 PBE 下 NiO 的带隙严重偏小（远低于实验约 4 eV），磁矩也偏小——这正是 [第8章](08-dftu-hp.md) DFT+U 的出发点。

### 4.4 FM / AFM / NM 能量比较

用**同一个超胞、同一套截断与 k 点**分别计算三种初始构型：

| 构型 | `starting_magnetization(1)` | `(2)` | 备注 |
|---|---|---|---|
| AFM-II | 0.5 | −0.5 | 实验基态 |
| FM | 0.5 | 0.5 | 对比 |
| NM | 去掉 `nspin` | — | 参考 |

能量差 $E_{\rm FM}-E_{\rm AFM}$ 可映射到 Heisenberg 模型求交换常数 J（映射公式与约定见 [VASP 第7章 §7](../../VASP/tutorials/07-magnetism-u-soc.md)，物理与代码无关）。注意：**AFM 与 FM 的能量必须来自同一超胞**，不能用 1 原子 FM 原胞直接比较 4 原子 AFM 胞（k 点采样不等价）。

> [!TIP]
> QE 7.4 起新增按 species 标签识别对称性（让 `Ni1`/`Ni2` 这类带编号的 species 参与对称性分析）以及可选的自旋翻转对称性（7.4 release notes）。注意：release notes 中写的是 `symmetries_with_labels`，但 7.5 源码（`Modules/input_parameters.f90`）中实际的变量名是 **`symmetry_with_labels`**（与 `use_spinflip` 一样位于 `&CONTROL`，默认 `.false.`），且这两个变量尚未写入 `INPUT_PW.html`。初学阶段保持默认即可；确需使用时以你所用版本的源码为准。

---

## 5. 非共线磁性

### 5.1 输入要点

`fe_noncolin.scf.in` 把 bcc Fe 的磁矩放在 x 方向（`angle1 = 90`、`angle2 = 0`），不开 SOC：

```fortran
&SYSTEM
  ibrav       = 3
  celldm(1)   = 5.42
  nat         = 1
  ntyp        = 1
  ecutwfc     = 90.0
  ecutrho     = 1080.0
  occupations = 'smearing'
  smearing    = 'mv'
  degauss     = 0.02
  noncolin    = .true.
  starting_magnetization(1) = 0.5
  angle1(1)   = 90.0
  angle2(1)   = 0.0
  report      = 1                ! print local moments every scf iteration
/
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `noncolin` | `.true.` | 默认 `.false.`；打开后**不要**写 `nspin`（`noncolin` 与 `nspin=2` 同时出现会报 `noncolin .and. nspin==2 are conflicting flags`） |
| `starting_magnetization(i)` | 同共线 | 全为 0 时不计算磁化（时间反演对称） |
| `angle1(i)` | 0–180 | 与 z 轴夹角（度） |
| `angle2(i)` | 0–360 | xy 投影与 x 轴夹角（度） |
| `report` | 1 | 每步打印各原子磁矩矢量，便于观察方向演化 |
| `mixing_beta` | 0.1–0.2 | 非共线比共线更难收敛 |

### 5.2 输出

非共线计算输出每个原子的磁化矢量与极坐标（`report` 控制频率）：

```bash
grep -e 'total magnetization' fe_noncolin.scf.out | tail -1
grep -A1 'magnetization :' fe_noncolin.scf.out | tail -4
grep 'polar coord' fe_noncolin.scf.out | tail -1
```

（示意输出，数值仅供参考）

```text
     total magnetization       =     2.25     0.00     0.00 Bohr mag/cell
     magnetization :          x.xxxxxx    0.000000    0.000000
     polar coord.: r, theta, phi [deg] :     x.xxxxxx   90.000000    0.000000
```

自检：同样截断与 k 点下，`fe_noncolin.scf.out` 的总能应与共线 `fe.scf.out` 一致（差异在 SCF 噪声量级）。若相差显著，先检查 k 点对称性是否因方向改变而被降低。

### 5.3 约束磁化：`constrained_magnetization`

| 取值 | 约束对象 | 罚函数（`INPUT_PW` 原式） |
|---|---|---|
| `'none'` | 不约束（默认） | — |
| `'total'` | 总磁化矢量 = `fixed_magnetization(1:3)` | λ Σᵢ (Mᵢ − Mᵢ⁰)²，**仅非共线**；共线请用 `tot_magnetization` |
| `'atomic'` | 每个 species 的原子磁矩 = 由 `starting_magnetization`、`angle1/angle2` 定义的 `mcons` | λ Σ (m − mcons)² |
| `'total direction'` | 总磁化与 z 轴夹角 θ = `fixed_magnetization(3)` | λ (arccos(M_z/|M|) − θ)² |
| `'atomic direction'` | 只约束 cos(angle1) | λ Σ (m_z/|m| − cos angle1)² |

`lambda` 默认 1.0。`INPUT_PW` 的建议：若 SCF 不收敛，先用小 `lambda` 收敛，再以更大的 `lambda` 重启。官方测试集 `test-suite/pw_noncolin/noncolin-constrain_atomic.in` 中 bcc Fe 使用 `lambda = 0.005`，可作为起点参考。

> [!WARNING]
> **QE 7.4 不兼容变更**：`'atomic'` 与 `'atomic direction'` 约束中 `lambda` 的含义被修改，与其他约束统一（7.4 release notes）。7.4 之前的论文/输入中的 `lambda` 数值不能直接照搬到 7.4+，需要重新测试。罚能会进入总能，比较能量时应确认罚项足够小（磁矩已贴近目标值）。

`INPUT_PW` 还提醒：对称化可能阻止磁化转到目标方向，尽量不要从高对称构型出发；`nosym` 只作最后手段。

---

## 6. 自旋轨道耦合：GaAs 价带劈裂

### 6.1 全相对论赝势

SSSP 中多数赝势是标量相对论的，SOC 需要另行下载 FR 赝势。pslibrary 1.0.0 的 FR 文件名带 `rel-` 前缀，可从 QE 官方赝势页直接下载：

```bash
mkdir -p pseudo && cd pseudo
wget https://pseudopotentials.quantum-espresso.org/upf_files/Ga.rel-pbe-dn-kjpaw_psl.1.0.0.UPF
wget https://pseudopotentials.quantum-espresso.org/upf_files/As.rel-pbe-n-kjpaw_psl.1.0.0.UPF
grep -o 'relativistic="[a-z]*"\|has_so="[a-z]*"' *.UPF
grep 'Suggested minimum cutoff' *.UPF
cd ..
```

两个文件头部均为 `relativistic="full" has_so="true"`；Ga 的建议最小截断为 60 / 255 Ry，As 为 20 / 103 Ry（来自 UPF 头部 `Suggested minimum cutoff`）。PseudoDojo 也提供 FR 模守恒赝势（`fr` 系列），做 DFPT 或 GW 时更常用，见 [第3章](03-pseudopotentials.md)。

### 6.2 输入（`gaas_soc.scf.in`）

```fortran
&CONTROL
  calculation = 'scf'
  prefix      = 'gaas_soc'
  outdir      = './tmp'
  pseudo_dir  = './pseudo'
/
&SYSTEM
  ibrav       = 2
  celldm(1)   = 10.68            ! a = 5.65 Angstrom (exp.)
  nat         = 2
  ntyp        = 2
  ecutwfc     = 60.0
  ecutrho     = 480.0
  occupations = 'fixed'
  noncolin    = .true.           ! two-component spinors (do NOT set nspin)
  lspinorb    = .true.           ! use the j-dependent (FR) projectors
/
&ELECTRONS
  conv_thr    = 1.0d-10
  mixing_beta = 0.5
/
ATOMIC_SPECIES
  Ga  69.723  Ga.rel-pbe-dn-kjpaw_psl.1.0.0.UPF
  As  74.922  As.rel-pbe-n-kjpaw_psl.1.0.0.UPF
ATOMIC_POSITIONS crystal
  Ga  0.00  0.00  0.00
  As  0.25  0.25  0.25
K_POINTS automatic
  8 8 8 0 0 0
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `noncolin` | `.true.` | SOC 只能在非共线框架下做；共线 + `lspinorb` 报 `spin orbit requires a non collinear calculation` |
| `lspinorb` | `.true.` | 使用赝势中的 SO 信息 |
| `starting_magnetization` | 不给 | 非磁 SOC；给了则变成磁性 SOC |
| `nbnd` | 默认或 2× 标量计算 | 旋量下每条带容纳 1 个电子，占据带数 = 价电子数 |
| `starting_spin_angle` | 默认 `.false.` | 只影响磁性 SOC 的初始波函数构造，一般不用改 |

### 6.3 能带与读数

```bash
mpirun -np 4 pw.x -nk 4 -in gaas_soc.scf.in   > gaas_soc.scf.out
mpirun -np 4 pw.x -nk 4 -in gaas_soc.bands.in > gaas_soc.bands.out
mpirun -np 4 bands.x -in gaas_soc.bands_pp.in > gaas_soc.bands_pp.out

grep 'Non magnetic calculation with spin-orbit' gaas_soc.scf.out
grep 'highest occupied' gaas_soc.scf.out
```

`gaas_soc.bands.in` 使用 `K_POINTS tpiba_b` 的 L–Γ–X 路径；`bands.x` 中 `lsigma(3) = .true.` 额外写出每个旋量态的 ⟨σz⟩ 文件 `gaas_soc.bands.dat.3`（非磁 SOC 中因 Kramers 简并而接近 0，磁性 SOC 中才有意义）。绘图可用 `python ../../scripts/qe_plot_bands.py -h` 查看用法（在示例目录中执行）；k 路径与 `bands.x` 的细节见 [第6章](06-band-dos-pdos.md)。

读取劈裂：在 Γ 点（路径第 41 个点）的价带顶附近，SOC 前三重简并（不计自旋）的 Γ15 态劈裂为四重 Γ8 与二重 Γ7，二者之差就是自旋轨道劈裂 Δ_SO。GaAs 实验值约 0.34 eV，PBE 的 Δ_SO 通常与之接近，但 PBE 带隙严重偏小，这一点与 SOC 无关（需要杂化泛函，见 [第14章](14-hybrid-metagga.md)）。

（示意输出，数值仅供参考；路径共 81 个 k 点，少于 100，`gaas_soc.bands.out` 默认就会打印各 k 点本征值，以下为 Γ 点片段）

```text
          k = 0.0000 0.0000 0.0000 (  xxxx PWs)   bands (ev):
   ...   x.xxxx  x.xxxx  y.yyyy  y.yyyy  y.yyyy  y.yyyy  z.zzzz  z.zzzz ...
```

其中两重的 `x.xxxx` 是 Γ7，四重的 `y.yyyy` 是 Γ8（价带顶）。

### 6.4 SOC 的计算量

- 旋量使波函数数组翻倍、对称操作减少（磁性 SOC 中时间反演被破坏），k 点数往往成倍增加；
- 大体系 SOC 常见策略：先用标量相对论做结构优化，再在固定结构上做 SOC 单点/能带；
- 7.x 中 SOC 与 PAW、USPP、DFT+U（第8章）均可组合；非共线 DFT+U+V 在 7.3 起支持（Dudarev 形式）。

---

## 7. 磁晶各向异性能（MAE）：力定理流程

### 7.1 流程

```text
fe_mae.lsda.in   (scf, nspin=2, FR USPP averaged to scalar)  -> tmp/fe_mae.save
       |
       +--> copy to tmp_001/, tmp_100/, tmp_110/
       |
fe_mae.nscf_<dir>.in  (nscf, noncolin, lspinorb, lforcet, angle1/angle2)
       |
grep "eband, Ef (eV)"  ->  MAE = eband[100] - eband[001]
```

`lforcet` 的 `INPUT_PW` 描述：从已有的共线 LSDA 密度文件开始非共线计算时，假设旧密度的磁化沿 z，并按 type 1 的 `angle1`/`angle2` 旋转。源码（`potinit.f90`）中这一分支调用 `nc_magnetization_from_lsda`，`print_ks_energies.f90` 在 `lforcet` 时打印：

```text
------
 eband, Ef (eV) =   <band energy>   <Fermi energy>
------
```

> [!IMPORTANT]
> QE 7.5 源码对 PAW 直接报错 `force theorem with PAW not implemented`，因此示例使用超软 FR 赝势 `Fe.rel-pbe-spn-rrkjus_psl.1.0.0.UPF`（建议最小截断 71 / 496 Ry，示例取 75 / 600 Ry）。另外，在 k 点 ≥ 100 时 `eband` 所在的打印块需要 `verbosity = 'high'`，示例已设置。

### 7.2 关键输入（`fe_mae.nscf_100.in` 的 `&SYSTEM`）

```fortran
&SYSTEM
  ibrav       = 3
  celldm(1)   = 5.42
  nat         = 1
  ntyp        = 1
  ecutwfc     = 75.0
  ecutrho     = 600.0
  occupations = 'smearing'
  smearing    = 'mv'
  degauss     = 0.01
  noncolin    = .true.
  lspinorb    = .true.
  lforcet     = .true.
  nosym       = .true.           ! same full k-point set for every direction
  starting_magnetization(1) = 0.5
  angle1(1)   = 90.0
  angle2(1)   = 0.0
/
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `calculation` | `'nscf'` | 非自洽 |
| `lforcet` | `.true.` | 读 LSDA 密度并旋转磁化 |
| `angle1(1)`, `angle2(1)` | 方向 | [001]: (0, 0)；[100]: (90, 0)；[110]: (90, 45) |
| `degauss` | 0.005–0.01 Ry | MAE 对 smearing 敏感，需测试 |
| `nosym` | `.true.` | 官方 `PP/examples/ForceTheorem_example` 的非自洽 SOC 步骤即如此设置，保证各方向使用同一套（全布里渊区）k 点 |
| K 网格 | 24³ 起，逐步加密 | bcc Fe 的 MAE 极小，需要非常密的 k 网格；`nosym` 下 24³ 即 13824 个 k 点，计算量很大 |

### 7.3 运行

```bash
bash run_all.sh mae
# or manually:
mpirun -np 8 pw.x -nk 8 -in fe_mae.lsda.in > fe_mae.lsda.out
for d in 001 100 110; do
  rm -rf tmp_$d; mkdir -p tmp_$d; cp -r tmp/fe_mae.save tmp_$d/
  mpirun -np 8 pw.x -nk 8 -in fe_mae.nscf_$d.in > fe_mae.nscf_$d.out
done
grep 'eband, Ef (eV)' fe_mae.nscf_*.out
```

MAE（以 [001] 为参考）：$E_{\rm MAE}=E_{\rm band}[100]-E_{\rm band}[001]$。bcc Fe 是立方体系，[100] 与 [001] 在对称上等价，**理论上差值应为 0**——这正好是一个估计数值噪声的基准；真正有意义的是 [001]/[110]/[111] 之间的差，以及在 CoPt、FePt 这类单轴体系中的 [001]/[100] 差。

### 7.4 精度要求

- bcc Fe 的 MAE 文献值在 μeV/atom 量级，远小于常规 k 点误差；做 MAE 时用"等价方向之差"估计噪声，噪声必须远小于你关心的 MAE；
- 所有方向必须用**同一套 k 网格和 smearing**；磁化方向不同时对称性约简得到的不可约 k 点集也不同，因此示例与官方 `ForceTheorem_example` 一样在 nscf 步设 `nosym = .true.`（代价是 k 点数大增）；
- 对单轴材料（L1₀ FePt、CrI₃ 单层），MAE 可达 meV 量级，力定理与自洽 SOC 差值的比较是很好的练习。

---

## 参数选择与收敛建议

| 项目 | 建议 |
|---|---|
| 截断能 | 按 SSSP 推荐值起步；FR 赝势参照 UPF 头部 `Suggested minimum cutoff`，并在**自旋极化**下复测 |
| k 点 | 磁性金属比非磁金属更密；Ni 至少 16³（原胞），MAE 24³ 以上并逐步加密 |
| smearing | 磁性金属用 `mv`，`degauss` 0.01–0.02 Ry；磁矩随 `degauss` 变化 < 0.01 μB 再认为收敛 |
| `starting_magnetization` | 磁性 species 0.2–0.6；对可能有高/低自旋的离子，分别从小值与大值出发比较能量 |
| AFM | 拆 species；能量对比在同一超胞内进行 |
| 混合 | `mixing_beta` 0.1–0.3；非均匀体系（表面、slab）可试 `mixing_mode = 'local-TF'` |
| 非共线/SOC | `mixing_beta` 0.1–0.2；先用共线结果判断磁态，再开 SOC |
| MAE | 力定理 + 极密 k 点；用等价方向之差估计噪声 |
| 约束 | `lambda` 从小值开始逐步增大，最终确认罚能可忽略 |

---

## 常见问题与报错

| 报错 / 现象 | 原因 | 解决 |
|---|---|---|
| `some starting_magnetization MUST be set` | `nspin=2`、smearing、未设 `tot_magnetization`，且所有 `starting_magnetization` 都没写 | 至少给一个 species 非零初始磁化 |
| `fixed occupations and lsda need tot_magnetization` | `nspin=2`，`occupations` 保持默认 `'fixed'` 且没给 `tot_magnetization` | 加 `occupations='smearing'`，或（绝缘体）给整数 `tot_magnetization` |
| `fixed occupations requires integer tot_magnetization` | `occupations='fixed'` 时 `tot_magnetization` 不是整数 | 四舍五入为整数，或改用 smearing |
| `tot_magnetization requires nspin=2` | 在 `nspin=1` 或 `noncolin` 下用了 `tot_magnetization` | 共线用 `nspin=2`；非共线改用 `constrained_magnetization='total'` |
| `noncolin .and. nspin==2 are conflicting flags` | 同时写了 `noncolin=.true.` 与 `nspin=2` | 删除 `nspin` |
| `spin orbit requires a non collinear calculation` | 设置了 `lspinorb` 却没开 `noncolin` | 加 `noncolin = .true.` |
| `Spin-orbit calculations require at least one spin-orbit pseudo` | 所有赝势都是标量相对论的 | 下载 FR 赝势（`*.rel-*` 或 PseudoDojo FR） |
| 信息 `Not all pseudopotentials have spin-orbit data` | 部分 species 用了标量赝势 | 轻元素（O、H）可以接受，重元素必须用 FR |
| `force theorem with PAW not implemented` | `lforcet` 与 PAW 赝势组合 | 改用 FR 超软或 NC 赝势 |
| 结果总磁化为 0（期望铁磁） | 初始磁化太小或为 0；对称性约束 | 增大 `starting_magnetization`，检查 species 设定 |
| AFM 计算两个原子磁矩同号 | 未拆 species 或初始符号写反 | 拆成 `X1`/`X2` 并给相反符号 |
| `convergence NOT achieved after 100 iterations: stopping` | 磁性 SCF 振荡 | 降 `mixing_beta` 到 0.1–0.2，增大 `electron_maxstep`，增大 `degauss` 先收敛再减小 |
| 非共线 + 无 SOC 时不同方向能量不同 | k 点约简随方向变化、SCF 精度不足 | 同一 k 点集（必要时 `nosym`），收紧 `conv_thr` |
| 7.4 前后约束计算结果不同 | `lambda` 含义在 7.4 改变 | 重新测试 `lambda` 并注明 QE 版本 |

---

## 练习

1. 对 bcc Fe 扫描 `degauss`（0.005、0.01、0.02、0.04 Ry）与 k 网格（10³–20³），画出总磁化随两者的变化，确定磁矩收敛到 0.01 μB 所需的参数。
2. 用 `nio_afm.scf.in` 构造 FM 版本（两个 Ni 都给 +0.5），在同一胞内比较 $E_{\rm FM}-E_{\rm AFM}$；思考为什么纯 PBE 的交换常数比实验偏大。
3. 对 O₂ 分子（立方盒子，Γ 点）分别做 `tot_magnetization = 0` 与 `= 2` 的计算，比较能量，验证三重态基态。
4. 计算 GaAs 不开 SOC（同一组 FR 赝势，但去掉 `noncolin/lspinorb`）的能带，与 SOC 能带叠加在一张图中，读取 Δ_SO。
5. 在 MAE 示例中把 k 网格从 16³ 增加到 32³，观察 [100] 与 [001] 的 `eband` 之差（理论为 0）如何减小，据此估计 [110] 相对 [001] 结果的可信度。
6. （进阶）用 `constrained_magnetization = 'atomic direction'` 让 bcc Fe 磁矩偏离 z 轴 30°，比较不同 `lambda` 下的罚能与磁矩方向。

---

## 参考资料

**官方文档**

- [INPUT_PW（pw.x 全部关键词）](https://www.quantum-espresso.org/Doc/INPUT_PW.html)：`nspin`、`starting_magnetization`、`tot_magnetization`、`noncolin`、`angle1/angle2`、`lforcet`、`constrained_magnetization`、`lambda`、`lspinorb`、`report`
- [INPUT_BANDS（bands.x）](https://www.quantum-espresso.org/Doc/INPUT_BANDS.html)：`lsigma`、`spin_component`
- [PWscf User's Guide](https://www.quantum-espresso.org/Doc/pw_user_guide/)
- QE 源码中的官方示例：`PW/examples/example06`（非共线磁性的总能与能带）、`example07`（Pt 的 FR 超软赝势 SOC）、`example11`（PAW 非共线/SOC 的 bcc Fe），`PP/examples/ForceTheorem_example`（Co slab 的力定理 MAE，`lforcet` + `nosym`），以及 `test-suite/pw_noncolin`、`test-suite/pw_spinorbit`：<https://gitlab.com/QEF/q-e>
- [QE 7.4 / 7.5 release notes（`lambda` 含义变化等）](https://gitlab.com/QEF/q-e/-/blob/develop/Doc/release-notes)
- 赝势：[SSSP 1.3 Efficiency](https://www.materialscloud.org/discover/sssp/table/efficiency) / [QE 官方赝势页（pslibrary FR 文件）](https://pseudopotentials.quantum-espresso.org/) / [PseudoDojo](https://www.pseudo-dojo.org/)

**课程与教程（只链接）**

- [Imperial College MSE404-MM](https://imperialcollegelondon.github.io/MSE404-MM/)：其中的自旋极化实验（Fe、O₂ 等）适合作为本章的配套练习，持续维护
- [pranabdas/espresso](https://pranabdas.github.io/espresso/)：从安装到磁性、SOC 的一站式教程网站
- [nguyen-group/QE-SSP](https://github.com/nguyen-group/QE-SSP)：《QE Course for Solid-State Physics》配套算例
- [Cornell PARADIM Summer School 2018](https://www.paradim.org/cu_ss_2018)：讲 DFT 原理与 QE 基本操作（QE 6.x、LDA），本章磁性内容可作延伸；下载/编译部分已过时
- [MaX e-School 2021（ICTP）](https://indico.ictp.it/event/9616/)：含视频讲座

**相关仓库（只链接）**

- [QEF/q-e](https://github.com/QEF/q-e)：官方镜像，`PW/examples` 与 `test-suite` 是最权威的输入示例
- [mailhexu/TB2J](https://github.com/mailhexu/TB2J)：基于 Wannier 函数/格林函数计算 Heisenberg 交换参数（可与 QE + Wannier90 结合，见 [第15章](15-wannier90.md)）
- [materialsproject/pymatgen](https://github.com/materialsproject/pymatgen)：磁构型枚举（`MagneticStructureEnumerator`）

**关键文献**

- A. Dal Corso, A. Mosca Conte, *Phys. Rev. B* 71, 115106 (2005)（超软赝势的 SOC）. DOI: [10.1103/PhysRevB.71.115106](https://doi.org/10.1103/PhysRevB.71.115106)
- A. Dal Corso, *Phys. Rev. B* 82, 075116 (2010)（PAW 的相对论自旋 DFT）. DOI: [10.1103/PhysRevB.82.075116](https://doi.org/10.1103/PhysRevB.82.075116)
- A. Dal Corso, *Comput. Mater. Sci.* 95, 337 (2014)（pslibrary）. DOI: [10.1016/j.commatsci.2014.07.043](https://doi.org/10.1016/j.commatsci.2014.07.043)
- D. Hobbs, G. Kresse, J. Hafner, *Phys. Rev. B* 62, 11556 (2000)（非共线磁性）. DOI: [10.1103/PhysRevB.62.11556](https://doi.org/10.1103/PhysRevB.62.11556)
- G. Prandini et al., *npj Comput. Mater.* 4, 72 (2018)（SSSP）. DOI: [10.1038/s41524-018-0127-2](https://doi.org/10.1038/s41524-018-0127-2)
- P. Giannozzi et al., *J. Phys.: Condens. Matter* 29, 465901 (2017). DOI: [10.1088/1361-648X/aa8f79](https://doi.org/10.1088/1361-648X/aa8f79)

---

[← 上一章](06-band-dos-pdos.md) | [返回目录](../README.md) | [下一章 →](08-dftu-hp.md)
