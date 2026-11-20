# 第5章 结构优化与状态方程（Structure Relaxation & Equation of State: relax, vc-relax, ev.x）

本章讲 Quantum ESPRESSO pw.x 的结构优化（structure relaxation）：固定晶胞的 `calculation='relax'` 与变胞的 `calculation='vc-relax'`，涉及 `&IONS`（`ion_dynamics='bfgs'`）、`&CELL`（`cell_dynamics`、`press`、`press_conv_thr`、`cell_dofree`）、收敛阈值 `etot_conv_thr`/`forc_conv_thr` 与 `if_pos` 冻结原子。第二部分用手动扫描晶格常数 + `ev.x`（Birch/Murnaghan 状态方程，Equation of State, EOS）或 Python 拟合得到平衡晶格常数与体模量，并解释 vc-relax 末尾 final scf 与 Pulay 应力的关系。

> [!NOTE]
> 本章全部输入在示例目录 [05-relax-eos](../examples/05-relax-eos/README.md)。截断能、k 点和 smearing 应先按 [第4章](04-convergence.md) 测好。文中"示意输出"只演示格式，**数值仅供参考，不是真实计算结果**。

---

## 1. 本章目标与适用场景

| 目标 | 你将学会 |
|---|---|
| relax | 固定晶胞优化原子位置：分子、表面吸附、缺陷超胞 |
| vc-relax | 同时优化晶胞与原子位置；在目标压强 `press` 下求平衡结构 |
| 收敛阈值 | `etot_conv_thr`、`forc_conv_thr`、`press_conv_thr` 的含义与单位 |
| 约束 | `if_pos` 冻结原子/分量；`cell_dofree` 限制晶胞自由度（`ibrav`、`2Dxy`、`shape` 等） |
| 结果读取 | `Begin final coordinates` 块、`bfgs converged`、`P=`、`Final enthalpy` |
| final scf | 理解 vc-relax 末尾以新晶胞重做 scf 的意义（Pulay 应力） |
| EOS | E(a) 或 E(V) 扫描 + `ev.x`/Python 拟合 Birch-Murnaghan，得到 a0、B0、B0' |

与 VASP 用户的对应关系：

| VASP | QE（pw.x） | 说明 |
|---|---|---|
| `IBRION=2`/`1`，`ISIF=2` | `calculation='relax'`，`ion_dynamics='bfgs'` | 固定晶胞 |
| `ISIF=3` | `calculation='vc-relax'`，`cell_dofree='all'` | 晶胞形状+体积+原子 |
| `ISIF=4` | `cell_dofree='shape'` | 固定体积 |
| `ISIF=7` | `cell_dofree='volume'` | 只变体积；源码限制只能用于 `ibrav=1`（简单立方） |
| `EDIFFG<0`（eV/Å） | `forc_conv_thr`（Ry/Bohr） | 1e-3 Ry/Bohr ≈ 0.026 eV/Å |
| `EDIFFG>0` | `etot_conv_thr`（Ry） | QE 中**两者都要满足** |
| `PSTRESS`（kbar） | `press`（kbar） | 目标压强 |
| `NSW` | `nstep` | 最大离子步数 |
| Selective dynamics `F F F` | `ATOMIC_POSITIONS` 行尾 `0 0 0` | `if_pos` |
| CONTCAR | 输出中的 `Begin final coordinates` 块 | 需自己拷回输入 |

VASP 的同类内容见 [VASP 第4章](../../VASP/tutorials/04-relax-eos.md)。

---

## 2. 原理速览

### 2.1 BFGS 准牛顿法

pw.x 的 relax/vc-relax 默认 `ion_dynamics='bfgs'`：用力（能量梯度）和逐步更新的近似逆 Hessian 预测下一步位置，并用**信赖半径**（trust radius）限制步长。相关关键词（全部只对 bfgs 有效）：

| 关键词 | 默认值 | 含义 |
|---|---|---|
| `trust_radius_ini` | `0.5` Bohr | 初始离子步长 |
| `trust_radius_max` | `0.8` Bohr | 最大步长 |
| `trust_radius_min` | `1.D-3` Bohr | 步长小于它时重置 BFGS 历史 |
| `bfgs_ndim` | `1` | =1 为标准 BFGS；>1 使用 GDIIS 混合 |
| `upscale` | `100.D0` | 接近收敛时 `conv_thr` 最多被自动缩小到 `conv_thr/upscale` |

其他可选：relax 下还有 `'damp'`（阻尼动力学，可配合 `CONSTRAINTS` 卡）和 `'fire'`（FIRE 算法）；vc-relax 下 `ion_dynamics` 可选 `'bfgs'` 或 `'damp'`，且 `'bfgs'` 要求 `cell_dynamics` 也是 `'bfgs'`。

### 2.2 收敛判据：能量和力**同时**满足

| 关键词 | 默认值 | 单位 | 说明 |
|---|---|---|---|
| `etot_conv_thr` | `1.0D-4` | Ry | 相邻两次 scf 总能变化；**广延量** |
| `forc_conv_thr` | `1.0D-3` | Ry/Bohr | 所有原子所有力分量都要小于它 |
| `press_conv_thr` | `0.5D0` | kbar | 仅 vc-relax；其他判据同样适用 |
| `nstep` | relax 类为 50 | — | 最大离子步数；`nstep=0` 只做初始化"干跑" |

官方文档明确：`etot_conv_thr` 与 `forc_conv_thr` **两个判据都必须满足**。pw.x 不会替你检查 `conv_thr` 与这两个判据是否匹配（BFGS 接近收敛时只会把 `conv_thr` 自动收紧，最多到 `conv_thr/upscale`）。实践上让 `conv_thr` 比 `etot_conv_thr` 小 2–3 个数量级，否则力里的 SCF 噪声会让 BFGS 来回振荡（见常见问题）。

单位换算速记：1 Ry/Bohr ≈ 25.71 eV/Å；默认 `forc_conv_thr = 1e-3` Ry/Bohr ≈ 0.026 eV/Å。精细优化（声子前、表面能）常用 `1e-4` Ry/Bohr。

### 2.3 变胞优化与 Pulay 应力

vc-relax 在晶胞变化时**固定平面波的数目**（由初始晶胞和 `ecutwfc` 决定）。晶胞缩小/膨胀后，固定数目的平面波对应的有效截断能变了，基组不完备产生虚假应力，即 **Pulay 应力**。

pw.x 的处理是：BFGS 收敛后，用最终晶胞**重新生成 G 矢量**再做一次 scf，输出中会看到：

```text
     Final scf calculation at the relaxed structure.
     The G-vectors are recalculated for the final unit cell
     Results may differ from those at the preceding step.
```

如果 final scf 给出的压强与最后一个 BFGS 步差很多（例如几 kbar），说明截断能对应力不够收敛。对策：

1. 提高 `ecutwfc`（和 `ecutrho`），以**应力收敛**为准（见 [第4章](04-convergence.md) 4.4 节）；
2. 用最终结构作为新起点再跑一次 vc-relax，直到 final scf 与最后一步基本一致；
3. `cell_factor`（变胞计算默认 2.0）用于构造赝势插值表，应大于晶胞在优化中的最大线性收缩；一般不用改，晶胞剧烈收缩时报错再调大。

> [!TIP]
> 官方文档中的 `ecfixed`、`qcutz`、`q2sigma` 是为变胞 MD 或应力计算设计的"修正动能泛函"参数（保持有效截断不变），一般结构优化不需要用；先用足够高的 `ecutwfc` 是更简单可靠的方案。

### 2.4 状态方程

对立方晶体，把总能 E 对体积 V 作图并拟合 EOS，就能得到平衡体积 V0（或晶格常数 a0）、体模量 B0 和其压力导数 B0'。常用的三阶 Birch-Murnaghan：

$$
E(V) = E_0 + \frac{9V_0B_0}{16}\left\{\left[\left(\frac{V_0}{V}\right)^{2/3}-1\right]^3 B_0' + \left[\left(\frac{V_0}{V}\right)^{2/3}-1\right]^2\left[6-4\left(\frac{V_0}{V}\right)^{2/3}\right]\right\}
$$

EOS 扫描相比 vc-relax 的优点：

- 每个点是固定晶胞的单点，**不受 Pulay 应力影响**（每个体积的基组都按自己的晶胞重新生成）；
- 同时得到体模量 B0，可与实验对比；
- 对多个相（fcc/bcc/hcp）做 E(V) 曲线可以比较相对稳定性与相变压强。

---

## 3. 实战一：Si 的 vc-relax

### 3.1 输入

[`si.vc-relax.in`](../examples/05-relax-eos/si.vc-relax.in)：

```fortran
&CONTROL
! Variable-cell relaxation of diamond Si (lattice parameter only, by symmetry).
  calculation   = 'vc-relax'
  prefix        = 'si'
  outdir        = './tmp'
  pseudo_dir    = './pseudo'
  etot_conv_thr = 1.0d-5      ! Ry, total (extensive)
  forc_conv_thr = 1.0d-4      ! Ry/Bohr, every force component
  nstep         = 100
/
&SYSTEM
  ibrav     = 2
  celldm(1) = 10.40
  nat       = 2
  ntyp      = 1
  ecutwfc   = 45.0
  ecutrho   = 360.0
  occupations = 'fixed'
/
&ELECTRONS
  conv_thr    = 1.0d-10
  mixing_beta = 0.7
/
&IONS
  ion_dynamics = 'bfgs'
/
&CELL
  cell_dynamics  = 'bfgs'
  press          = 0.0        ! target pressure in kbar
  press_conv_thr = 0.2        ! kbar (default 0.5)
  cell_dofree    = 'ibrav'    ! keep the fcc lattice of ibrav=2
/
ATOMIC_SPECIES
  Si  28.0855  Si.pbe-n-rrkjus_psl.1.0.0.UPF
ATOMIC_POSITIONS crystal
  Si  0.00  0.00  0.00
  Si  0.25  0.25  0.25
K_POINTS automatic
  8 8 8 0 0 0
```

| 关键词 | 取值 | 说明 |
|---|---|---|
| `calculation` | `'vc-relax'` | 自动打开 `tstress`；relax/vc-relax 也会自动计算力 |
| `etot_conv_thr` | `1.0d-5` Ry | 比默认严 10 倍；注意它是总量，2 原子体系没问题 |
| `forc_conv_thr` | `1.0d-4` Ry/Bohr | Si 原子在高对称位置，力本来为零 |
| `&IONS` | 必须出现 | relax/vc-relax 必须有 `&IONS` 段（即使用默认值） |
| `&CELL` | 必须出现 | vc-relax/vc-md 必须有 `&CELL` 段 |
| `cell_dynamics` | `'bfgs'` | 与 `ion_dynamics='bfgs'` 配套 |
| `press` | `0.0` kbar | 目标压强；高压研究可设为如 `100.0`（=10 GPa） |
| `press_conv_thr` | `0.2` kbar | 默认 0.5 kbar |
| `cell_dofree` | `'ibrav'` | 晶格保持 `ibrav=2`（fcc）的形式，只允许对称性允许的变化 |
| `ecutwfc`/`ecutrho` | 45/360 Ry | 高于能量收敛值，照顾应力 |

> [!NOTE]
> `cell_dofree` 的常用取值：`'all'`（默认，全部轴和角）、`'ibrav'`（保持初始 ibrav 晶格类型，可与其它选项组合为 `'ibrav+...'`）、`'shape'`（固定体积）、`'volume'`（只变体积、角度不变；pw.x 源码中限制为 `ibrav=1`，否则报 `Isotropic expansion is only allowed for ibrav=1`）、`'2Dxy'`（只变 x、y 分量，二维材料用，见 [第9章](09-vdw-surface-workfunction.md)）、`'z'`、`'xyz'`、`'epitaxial_ab'` 等。轴不正交时部分选项会破坏对称性，官方文档有提醒。

### 3.2 运行

```bash
cd QE/examples/05-relax-eos
mpirun -np 4 pw.x -nk 4 -in si.vc-relax.in > si.vc-relax.out
```

`-nk 4` 表示 4 个 k 点池（pool）并行；Si 8×8×8 有足够的不可约 k 点。并行参数详见 [第18章](18-automation-parallel-troubleshooting.md)。集群上使用 [`job_relax_slurm.sh`](../examples/05-relax-eos/job_relax_slurm.sh)（Slurm 基础见 [Common 第2章](../../Common/tutorials/02-job-scheduler.md)）。

### 3.3 读输出

```bash
# 每个离子步的压强与总能
grep -e '^!' -e 'P=' si.vc-relax.out
# 每步后的晶胞与体积
grep -e 'new unit-cell volume' si.vc-relax.out
# BFGS 是否收敛
grep -A1 -e 'bfgs converged' si.vc-relax.out
grep -e 'Final enthalpy' si.vc-relax.out
# 最终结构
grep -A12 'Begin final coordinates' si.vc-relax.out
grep 'JOB DONE' si.vc-relax.out
```

典型片段（示意输出，数值仅供参考）：

```text
     bfgs converged in   x scf cycles and   x bfgs steps
     (criteria: energy <  1.0E-05 Ry, force <  1.0E-04 Ry/Bohr, cell <  2.0E-01 kbar)

     End of BFGS Geometry Optimization

     Final enthalpy =     -15.8xxxxxxx Ry
Begin final coordinates
     new unit-cell volume =    xxx.xxxxx a.u.^3 (    xx.xxxxx Ang^3 )
     density =      x.xxxxx g/cm^3

CELL_PARAMETERS (alat= 10.40000000)
  -0.49xxxxxxx   0.000000000   0.49xxxxxxx
   0.000000000   0.49xxxxxxx   0.49xxxxxxx
  -0.49xxxxxxx   0.49xxxxxxx   0.000000000

ATOMIC_POSITIONS (crystal)
Si            0.0000000000        0.0000000000        0.0000000000
Si            0.2500000000        0.2500000000        0.2500000000
End final coordinates
```

要点：

- vc-relax 优化的目标是**焓** $H = E + PV$，所以收敛后打印 `Final enthalpy`；relax 打印 `Final energy`；
- `CELL_PARAMETERS (alat= 10.40000000)` 中的数值以**初始** `celldm(1)` 为单位。最终晶格常数 = 10.40 × (2 × 0.49xxx) Bohr（fcc 原胞矢量的分量为 a/2）。更省事的办法是从体积换算：fcc 原胞 $V = a^3/4$；
- 收敛块之后还有 final scf；最后一次 `P=` 是 final scf 的压强。

```bash
grep 'P=' si.vc-relax.out | tail -2
```

如果最后两个值相差显著，回到 2.3 节处理 Pulay 应力。

> [!TIP]
> 在示例目录下用 `python ../../scripts/qe_parse_out.py -h` 查看解析脚本用法，它的 `--final-structure` 选项可以把最终结构写成新的输入片段；ASE 也能读 pw.x 输出：`ase.io.read('si.vc-relax.out', index=-1)`（见 [Common 第4章](../../Common/tutorials/04-ase-pymatgen.md)）。

### 3.4 把结果用于后续计算

后续 scf/能带计算应使用优化后的晶格常数。两种做法：

1. 把新的 a（Bohr）写回 `celldm(1)`，保持 `ibrav=2`；
2. 用 `ibrav=0` + 输出中的 `CELL_PARAMETERS`：若保留 `(alat=...)` 单位，则 `celldm(1)` 必须等于括号中的 alat；若换算成 `bohr` 或 `angstrom`，则**不能**再给 `celldm(1)`/`A`，否则报 `lattice parameter specified twice`。

> [!WARNING]
> 不写单位的 `CELL_PARAMETERS`/`ATOMIC_POSITIONS` 已被标记为 DEPRECATED。拷贝结构时**务必显式写单位**，并确认 `celldm(1)` 与 `A` 不同时出现。

---

## 4. 实战二：CO 分子的 relax（冻结原子）

### 4.1 输入

[`co.relax.in`](../examples/05-relax-eos/co.relax.in)：

```fortran
&CONTROL
! Fixed-cell relaxation of a CO molecule in a 12 Angstrom cubic box.
! The C atom is frozen with if_pos = 0 0 0; only O moves.
  calculation   = 'relax'
  prefix        = 'co'
  outdir        = './tmp'
  pseudo_dir    = './pseudo'
  etot_conv_thr = 1.0d-5
  forc_conv_thr = 1.0d-4
  nstep         = 50
/
&SYSTEM
  ibrav     = 1
  A         = 12.0            ! cubic box edge in Angstrom
  nat       = 2
  ntyp      = 2
  ecutwfc   = 50.0            ! max of the C and O hints
  ecutrho   = 400.0
  occupations = 'fixed'       ! closed-shell molecule
/
&ELECTRONS
  conv_thr    = 1.0d-9
  mixing_beta = 0.5
/
&IONS
  ion_dynamics = 'bfgs'
/
ATOMIC_SPECIES
  C  12.011  C.pbe-n-kjpaw_psl.1.0.0.UPF
  O  15.999  O.pbe-n-kjpaw_psl.0.1.UPF
ATOMIC_POSITIONS angstrom
  C  6.000  6.000  6.000   0 0 0
  O  6.000  6.000  7.200
K_POINTS gamma
```

| 关键词 | 取值 | 说明 |
|---|---|---|
| `calculation` | `'relax'` | 晶胞固定，只动原子 |
| `ibrav=1` + `A` | 12.0 Å | 简单立方盒子；`A` 单位为 Å（`celldm(1)` 为 Bohr，二者只能给一个） |
| `if_pos` | C：`0 0 0` | 该原子三个力分量乘 0，即冻结；O 不写默认 `1 1 1` |
| `K_POINTS gamma` | — | 孤立分子只需 Γ 点，并启用 Γ 专用算法（省内存与时间） |
| `ecutwfc`/`ecutrho` | 50/400 Ry | 取 C（45/360）与 O（50/400）SSSP 提示值的最大值 |

为什么冻结 C？孤立分子在盒子里可以整体平移，总力为零时平移不影响能量；冻结一个原子去掉平移自由度，也让输出更易读。对于表面 slab，常把底部几层原子用 `0 0 0` 冻结来模拟体相（[第9章](09-vdw-surface-workfunction.md)）。也可以只冻结某个方向，例如 `1 1 0` 只允许在 xy 平面内移动。

> [!NOTE]
> 分子在周期盒子中会与镜像相互作用。12 Å 对 CO 这样的非极性小分子一般够用，但偶极分子（如 H2O）可以考虑 `assume_isolated='mt'`（Martyna-Tuckerman 修正），盒子大小也需要测试（见 [第2章](02-input-anatomy.md)）。

### 4.2 运行与检查

```bash
mpirun -np 4 pw.x -in co.relax.in > co.relax.out
grep -e '^!' co.relax.out                     # 每个离子步的总能
grep -A3 'Forces acting on atoms' co.relax.out | tail -3
grep -e 'bfgs converged' -e 'Final energy' co.relax.out
grep -A4 'Begin final coordinates' co.relax.out
```

C–O 键长 = 两原子 z 坐标之差。与实验值（约 1.128 Å）比较时请记住 PBE 对键长有系统偏差，并且结果依赖赝势与截断。

---

## 5. 实战三：fcc Cu 的状态方程与 ev.x

### 5.1 为什么用 Cu

Cu 是金属（需要 smearing 与较密 k 点），SSSP 1.3 Efficiency 中 Cu 的提示截断较高（90/720 Ry），适合演示"金属 + 高截断"的 EOS 工作流。

### 5.2 模板与扫描

[`cu_eos/cu_eos_template.in`](../examples/05-relax-eos/cu_eos/cu_eos_template.in) 用 `A = @ALAT@`（Å）给晶格常数：

```fortran
&CONTROL
  calculation = 'scf'
  prefix      = 'cu'
  outdir      = './tmp'
  pseudo_dir  = './pseudo'
  tstress     = .true.
/
&SYSTEM
  ibrav     = 2
  A         = @ALAT@
  nat       = 1
  ntyp      = 1
  ecutwfc   = 90.0
  ecutrho   = 720.0
  occupations = 'smearing'
  smearing    = 'mv'
  degauss     = 0.02
/
&ELECTRONS
  conv_thr    = 1.0d-10
  mixing_beta = 0.7
/
ATOMIC_SPECIES
  Cu  63.546  Cu.paw.z_11.ld1.psl.v1.0.0-low.upf
ATOMIC_POSITIONS crystal
  Cu  0.00  0.00  0.00
K_POINTS automatic
  16 16 16 0 0 0
```

扫描 7 个晶格常数（约 ±3% 范围，覆盖预期极小值两侧）：

```bash
cd QE/examples/05-relax-eos/cu_eos
ln -s ../pseudo pseudo
RUN=1 PW_CMD="mpirun -np 8 pw.x -nk 4" bash run_eos.sh
cat cu_eos.dat          # 两列: a (Angstrom)  E (Ry)
```

扫描设计要点：

- 点数 ≥ 7，范围约 ±3%–5%（体积约 ±10%–15%），极小值应落在**区间内部**；
- 所有点使用**完全相同**的截断、k 网格、smearing；
- 范围太大会进入非谐区、使 B0' 不可靠；范围太小则拟合对噪声敏感；
- `tstress=.true.` 让每个点也输出压强，可作为一致性检查：E 极小值附近压强应接近 0。

### 5.3 用 ev.x 拟合

`ev.x` 是 QE 自带的交互式 EOS 拟合程序（源码 `PW/tools/ev.f90`），依次读入 5 个回答：

| 提示 | 回答示例 | 说明 |
|---|---|---|
| `Lattice parameter or Volume are in (au, Ang) >` | `Ang` | 输入数据的长度单位 |
| `Enter type of bravais lattice (fcc, bcc, sc, noncubic) >` | `fcc` | 立方晶格时第一列是晶格常数；`noncubic` 时第一列是体积 |
| `Enter type of equation of state : 1=birch1, 2=birch2, 3=keane, 4=murnaghan >` | `2` | 2 = 三阶 Birch（4 参数）；1、4 为 3 参数 |
| `Input file >` | `cu_eos.dat` | 两列：a（或 V）与 E（Ry） |
| `Output file >` | `cu_eos.ev.out` | 拟合结果 |

非交互运行：

```bash
ev.x << EOF
Ang
fcc
2
cu_eos.dat
cu_eos.ev.out
EOF
head -3 cu_eos.ev.out
```

输出文件开头的格式（示意输出，数值仅供参考）：

```text
# equation of state: birch 3rd order.  chisq =   0.xxxxD-xx
# a0 = x.xxxx a.u., k0 =  xxxx kbar, dk0 =  x.xx d2k0 = -x.xxx emin =  -xxx.xxxxx
# a0 =  x.xxxxx Ang, k0 = xxx.x GPa,  V0 =    xx.xx (a.u.)^3,  V0 =    xx.xx A^3
```

随后是每个点的 `Lat.Par / E_calc / E_fit / E_diff / Pressure / Enthalpy` 表。检查 `E_diff` 是否都很小、没有某一点明显偏离（若有，检查该点是否 SCF 未收敛或 k 网格不同）。

### 5.4 用 Python 拟合（可选）

仓库脚本 `qe_eos_fit.py` 读取 (V,E) 或 (a,E) 两列数据拟合 Birch-Murnaghan，在 `cu_eos/` 目录下用 `python ../../../scripts/qe_eos_fit.py -h` 查看用法（fcc 的 E(a) 数据配合 `--x alat --vfactor 0.25`）。也可以直接用 ASE：

```python
import numpy as np
from ase.eos import EquationOfState
from ase.units import Ry, GPa

a, e = np.loadtxt("cu_eos.dat", unpack=True)   # Angstrom, Ry
v = a**3 / 4.0                                 # fcc primitive-cell volume, A^3
eos = EquationOfState(v, e * Ry, eos="birchmurnaghan")
v0, e0, b = eos.fit()
print("a0 = %.4f A, B0 = %.1f GPa" % ((4 * v0) ** (1 / 3), b / GPa))
eos.plot("cu_eos.png")
```

ASE 的 `EquationOfState` 用法见 [Common 第4章](../../Common/tutorials/04-ase-pymatgen.md)。

### 5.5 EOS 与 vc-relax 的交叉验证

| 方法 | 优点 | 缺点 |
|---|---|---|
| vc-relax | 一次运行；可处理非立方、内坐标同时优化 | 受 Pulay 应力影响，需足够高截断 |
| E(V) 扫描 + EOS | 不受 Pulay 应力影响；同时得到 B0、B0' | 非立方晶体每个体积还要优化 c/a 与内坐标 |

两种方法得到的 a0 应在约 0.1% 量级内一致；差得多时，优先怀疑 vc-relax 的截断。

---

## 6. relax/vc-relax 的实用技巧

- **先 relax 后 vc-relax**：初始结构很差时（如手工搭建的表面吸附构型），先固定晶胞 relax，避免晶胞和原子同时大幅移动导致 SCF 崩溃；
- **步子太大**：若中途出现 `S matrix not positive definite` 或 `problems computing cholesky`，减小 `trust_radius_ini`（如 0.1–0.2 Bohr），检查原子是否靠得太近；
- **重启**：最简单的方法是把 `Begin final coordinates`（未收敛时用最后一个 `ATOMIC_POSITIONS` 块）拷回输入重新运行；`restart_mode='restart'` 与 `max_seconds` 的用法见 [第18章](18-automation-parallel-troubleshooting.md)；
- **对称性**：pw.x 会保持输入结构的对称性。若想让体系自发破缺对称（如 Jahn-Teller 畸变），需要在初始结构中加入小扰动；
- **磁性体系**：优化时 `nspin=2` 与 `starting_magnetization` 必须保留，见 [第7章](07-magnetism-soc.md)；
- **范德华体系**：层状材料的 vc-relax 必须加色散修正（`vdw_corr`），见 [第9章](09-vdw-surface-workfunction.md)；
- **快速检查输入**：`nstep = 0` 会在初始化后立即停止，用于检查对称性、k 点数、内存估计。

---

## 参数选择与收敛建议

| 参数 | 粗优化 | 精细优化 / 声子前 | 说明 |
|---|---|---|---|
| `forc_conv_thr` | `1.0d-3`（默认） | `1.0d-4` 或更小 | Ry/Bohr |
| `etot_conv_thr` | `1.0d-4`（默认） | `1.0d-5`～`1.0d-6` | Ry，广延量；大体系按原子数放宽 |
| `press_conv_thr` | `0.5`（默认） | `0.1`～`0.2` | kbar |
| `conv_thr` | `1.0d-6`～`1.0d-8` | `1.0d-10` 量级 | 太大时力有噪声，BFGS 来回振荡 |
| `ecutwfc` | 能量收敛值 | 应力收敛值（通常更高） | vc-relax/EOS 以应力为准 |
| `trust_radius_ini` | `0.5`（默认） | `0.1`～`0.2`（结构差/易崩时） | Bohr |
| `nstep` | 50（默认） | 100–200 | 软模式体系需要更多步 |
| EOS 范围 | ±3%–5% a | 7–11 个点 | 极小值在区间内部 |

---

## 常见问题与报错

| 报错 / 现象 | 原因 | 解决办法 |
|---|---|---|
| 反复出现 `uphill step: resetting bfgs history`，能量来回振荡 | `conv_thr` 相对 `etot_conv_thr`/`forc_conv_thr` 太松，力有 SCF 噪声 | 收紧 `conv_thr`（比 `etot_conv_thr` 小 2–3 个数量级，如 1e-8～1e-10）；检查 `ecutrho` 是否足够 |
| `The maximum number of steps has been reached.` | `nstep` 用完仍未收敛 | 用最后的结构继续；检查是否有软模式或振荡；增大 `nstep` |
| `trust_radius < trust_radius_min` / `resetting bfgs history` | BFGS 步长不断缩小，常因力有噪声 | 收紧 `conv_thr`、提高截断；必要时换 `ion_dynamics='fire'`（relax） |
| `bfgs failed after ... scf cycles and ... bfgs steps, convergence not achieved` | BFGS 无法进一步降低能量 | 检查力的噪声（`conv_thr`、`ecutrho`）；从当前结构重启 |
| 中途出现 `S matrix not positive definite` / `problems computing cholesky` | 离子步过大导致原子过近 | 减小 `trust_radius_ini`；先 relax 再 vc-relax；检查单位 |
| final scf 的 `P=` 与最后一个 BFGS 步差很多 | Pulay 应力，截断不足 | 提高 `ecutwfc`/`ecutrho`；从新结构重跑 vc-relax |
| vc-relax 后晶格对称性降低 | `cell_dofree='all'` 下数值噪声累积 | 用 `cell_dofree='ibrav'`；或用对称化工具（spglib）处理后重跑 |
| `ibrav=0: must read cell parameters` | `ibrav=0` 但未给 `CELL_PARAMETERS` | 拷贝结构时补全 `CELL_PARAMETERS {angstrom}` |
| `lattice parameter specified twice` | `CELL_PARAMETERS` 用 `bohr`/`angstrom` 单位，同时又给了 `celldm(1)` 或 `A` | 删掉 `celldm(1)`/`A`，或改用 `alat` 单位 |
| EOS 拟合的 a0 落在扫描端点 | 扫描范围没有包住极小值 | 平移/扩大扫描区间 |
| EOS 中某一点明显偏离曲线 | 该点 SCF 未收敛、用了不同参数或金属 k 点不足 | 检查该点 `pw.out`，统一参数后重算 |

---

## 练习

1. 把 `si.vc-relax.in` 的 `ecutwfc`/`ecutrho` 降到 30/240 Ry 再跑一次，比较 final scf 前后的压强差，体会 Pulay 应力。
2. 对 Si 做 E(a) 扫描（7 点，±3%），用 `ev.x` 拟合，与 vc-relax 的 a0 对比。
3. 在 `co.relax.in` 中去掉 C 原子的 `0 0 0`，比较收敛步数和最终键长；再把盒子从 12 Å 改为 10 Å 与 15 Å，观察键长和总能的变化。
4. 用 `ev.x` 的 4 种 EOS（1–4）分别拟合同一组 Cu 数据，比较 a0、k0 的差异。
5. 设 `press = 100.0`（10 GPa）对 Si 做 vc-relax，用得到的体积检查是否与 EOS 拟合给出的 P(V) 一致。
6. 用 fcc 与 bcc 两种结构分别做 Cu 的 E(V) 扫描（bcc 用 `ibrav=3`，`ev.x` 输入 `bcc`），比较每原子能量，确认 fcc 更稳定。

---

## 参考资料

**官方文档**

- [INPUT_PW（pw.x）](https://www.quantum-espresso.org/Doc/INPUT_PW.html) — `&IONS`、`&CELL`、`etot_conv_thr`、`forc_conv_thr`、`press_conv_thr`、`cell_dofree`、`if_pos`
- [PWscf User's Guide](https://www.quantum-espresso.org/Doc/pw_user_guide/) — 结构优化与变胞计算说明
- [ev.x 源码 `PW/tools/ev.f90`（GitLab）](https://gitlab.com/QEF/q-e/-/blob/develop/PW/tools/ev.f90) — 交互输入与输出格式以源码为准
- [QE GitLab 主仓库](https://gitlab.com/QEF/q-e) — `PW/examples` 中有 relax/vc-relax 官方示例

**课程与教程（推荐，链接而非复制）**

- [Imperial College MSE404-MM Lab 7: Structural Optimisation](https://imperialcollegelondon.github.io/MSE404-MM/labs/lab07/readme/) — relax/vc-relax 的系统讲解
- [Cornell PARADIM Summer School 2016](https://www.paradim.org/cu_ss_2016) — "分子与体相平衡结构""自动优化与弹性常数"上机（QE 5.x/6.x 时代，适合理解原理）
- [Cornell PARADIM Summer School 2018](https://www.paradim.org/cu_ss_2018) — Cl2 结合能、金刚石弹性常数，有视频
- [pranabdas/espresso: Structure optimization](https://pranabdas.github.io/espresso/hands-on/structure-optimization) — 带 Python 后处理的 relax/vc-relax 示例
- [MIT OCW 3.320（2005）](https://ocw.mit.edu/courses/3-320-atomistic-computer-modeling-of-materials-sma-5107-spring-2005/) — 早期 PWscf 晶格常数实验，历史参考

**相关仓库（仅链接）**

- [jkitchin/dft-book-espresso](https://github.com/jkitchin/dft-book-espresso) — 用 ASE 驱动 QE 做 EOS、表面与吸附
- [QEF/postqe](https://github.com/QEF/postqe) — 官方 Python 后处理，含 EOS
- [dalcorso/thermo_pw](https://github.com/dalcorso/thermo_pw) — 自动化 EOS、弹性常数、准谐近似
- [AlexanderKvashnin/QE_tutorial](https://github.com/AlexanderKvashnin/QE_tutorial) — Skoltech 课程 Lab：优化、弹性等
- [spglib/spglib](https://github.com/spglib/spglib) — 对称性判断与结构对称化

**关键文献**

- F. Birch, *Phys. Rev.* 71, 809 (1947)（Birch 状态方程）. DOI: [10.1103/PhysRev.71.809](https://doi.org/10.1103/PhysRev.71.809)
- F. D. Murnaghan, *Proc. Natl. Acad. Sci. USA* 30, 244 (1944). DOI: [10.1073/pnas.30.9.244](https://doi.org/10.1073/pnas.30.9.244)
- P. Giannozzi et al., *J. Phys.: Condens. Matter* 21, 395502 (2009). DOI: [10.1088/0953-8984/21/39/395502](https://doi.org/10.1088/0953-8984/21/39/395502)
- P. Giannozzi et al., *J. Phys.: Condens. Matter* 29, 465901 (2017). DOI: [10.1088/1361-648X/aa8f79](https://doi.org/10.1088/1361-648X/aa8f79)
- G. Prandini et al., *npj Comput. Mater.* 4, 72 (2018)（SSSP，含应力/EOS 验证）. DOI: [10.1038/s41524-018-0127-2](https://doi.org/10.1038/s41524-018-0127-2)
- E. Bosoni et al., *Nat. Rev. Phys.* 6, 45 (2024)（全电子与赝势代码的 EOS 对比基准）. DOI: [10.1038/s42254-023-00655-3](https://doi.org/10.1038/s42254-023-00655-3)

---

[← 上一章](04-convergence.md) | [返回目录](../README.md) | [下一章 →](06-band-dos-pdos.md)
