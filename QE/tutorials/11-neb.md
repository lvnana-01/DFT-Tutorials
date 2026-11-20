# 第11章 过渡态搜索：neb.x（Nudged Elastic Band with neb.x）

本章讲解用 Quantum ESPRESSO 的 neb.x 搜索最小能量路径（MEP）和过渡态（transition state, TS），从而得到反应或扩散的能垒。内容包括：neb.x 的超级卡片输入结构（`BEGIN`/`BEGIN_PATH_INPUT`/`&PATH`/`BEGIN_ENGINE_INPUT`/`BEGIN_POSITIONS`/`FIRST_IMAGE`/`INTERMEDIATE_IMAGE`/`LAST_IMAGE`），`&PATH` 中 `num_of_images`、`opt_scheme='broyden'`、`CI_scheme='auto'`、`path_thr`、`k_max`/`k_min` 等关键词，`-inp` 运行方式与 `-ni` image 并行，`prefix.dat`/`.int`/`.xyz`/`.axsf`/`.crd`/`.path` 输出文件的含义，以及用 ph.x Γ 点计算验证鞍点只有一个虚频。全章以 H 原子在 Al(100) 表面两个相邻四重空位之间的扩散为例。

> [!NOTE]
> 本章全部输入在示例目录 [11-neb](../examples/11-neb/README.md)。`&PATH` 关键词及默认值对照 QE 7.5/7.6 的 [INPUT_NEB](https://www.quantum-espresso.org/Doc/INPUT_NEB.html) 和 NEB 源码核实。仓库不提交 UPF 文件；文中输出片段均为**示意输出，数值仅供参考**。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 说清 NEB 与 climbing-image NEB（CI-NEB）的原理，以及能垒为什么会被普通 NEB 低估；
- 在同一套设置下严格弛豫初态和末态，保证原子顺序与固定约束一致；
- 写出完整的 neb.x 输入，理解每个 `&PATH` 关键词的作用和默认值（尤其是默认 `nstep_path=1` 这个坑）；
- 用 `neb.x -ni N -inp file` 运行，并按 image 数分配 MPI 进程；
- 读懂 neb.x 主输出中的 `activation energy`、每个 image 的 `energy`/`error`/`frozen` 表，以及 `.dat`/`.int`/`.xyz`/`.axsf` 文件；
- 中断后用 `restart_mode='restart'` 续算；
- 用 ph.x（`nat_todo`）+ dynmat.x 快速确认鞍点只有一个虚频，并知道更严格的验证方法。

| 想做的事 | 方法 | 本章小节 |
|---|---|---|
| 已知初末态，求 MEP 和能垒 | NEB（`CI_scheme='no-CI'`） | §3–§6 |
| 需要更准确的鞍点能量和结构 | CI-NEB（`CI_scheme='auto'`） | §2.2、§4 |
| 手动指定爬坡 image（多峰路径） | `CI_scheme='manual'` + `CLIMBING_IMAGES` | §4.3 |
| 路径初猜很差、需要中间结构 | `INTERMEDIATE_IMAGE` | §4.4 |
| 鞍点验证 | ph.x Γ 点 + dynmat.x | §7 |
| 有限温度自由能垒、溶剂效应 | 元动力学/约束 MD（本节不覆盖），见 [第17章](17-md-cp.md) 的 MD 基础 | — |

**前置知识**：结构优化（[第5章](05-relax-eos.md)）、slab 建模与固定原子（[第9章](09-vdw-surface-workfunction.md)）、收敛测试（[第4章](04-convergence.md)）。

### 与 VASP 的对应关系

| VASP（VTST） | QE neb.x | 说明 |
|---|---|---|
| `IMAGES`（不含端点） | `num_of_images`（**含**两个端点） | VASP 中 IMAGES=5 相当于 QE 中 `num_of_images=7` |
| `00/`、`01/` … 目录 + POSCAR | 单个输入文件中的 `FIRST_IMAGE`/`INTERMEDIATE_IMAGE`/`LAST_IMAGE` | QE 自动线性插值中间 image |
| `LCLIMB = .TRUE.` | `CI_scheme = 'auto'` | |
| `SPRING = -5` | `k_max`、`k_min`（Hartree 原子单位） | QE 可以用可变弹簧常数 |
| `IBRION=3`、`POTIM=0` + `IOPT=1`（LBFGS）/ `IOPT=3`（QuickMin） | `opt_scheme = 'broyden'` / `'quick-min'` | 官方推荐 broyden（拟牛顿）；QE neb.x 没有 FIRE（VTST 的 `IOPT=7`） |
| `EDIFFG = -0.05` | `path_thr = 0.05`（eV/Å） | QE 只看垂直于路径的力 |
| `nebresults.pl`、`neb.dat` | `prefix.dat`、`prefix.int`、主输出中的 `activation energy` | |
| 按 image 数整除的 MPI 进程 | `-ni`（`-nimage`）image 并行 | |

VASP 的 NEB/Dimer 流程见 [VASP 第11章](../../VASP/tutorials/11-neb-dimer.md)，CP2K 见 [CP2K 第13章](../../CP2K/tutorials/13-neb-dimer.md)。QE 的 neb.x **没有** Dimer 方法；需要 Dimer 时可以借助 ASE 等外部工具驱动 pw.x（见 [Common 第4章](../../Common/tutorials/04-ase-pymatgen.md)）。

---

## 2. 原理速览

### 2.1 最小能量路径与鞍点

势能面 E(R) 上，初态（IS）和末态（FS）是两个局部极小点。连接它们、每一点上垂直于路径的力都为零的路径称为**最小能量路径**（MEP）。MEP 上能量最高的点是一阶鞍点，即过渡态，它满足：

1. 力为零：∇E(R‡) = 0；
2. Hessian 恰有**一个**负本征值，对应的本征矢就是反应坐标方向。

正向能垒 E_a = E(R‡) − E(IS)，逆向能垒 E_a' = E(R‡) − E(FS)。谐振子过渡态理论（hTST）给出速率常数

$$
k = \nu^{*}\exp\!\left(-\frac{E_a}{k_\mathrm B T}\right),\qquad
\nu^{*} = \frac{\prod_{i=1}^{3N}\nu_i^\mathrm{IS}}{\prod_{j=1}^{3N-1}\nu_j^{\ddagger}}
$$

指前因子 ν* 需要初态和 TS 的全部实频（TS 去掉虚频），这就是 §7 做频率计算的原因之一。

### 2.2 NEB 与 CI-NEB

NEB 把路径离散成 M 个 image R₁, …, R_M（QE 中 M = `num_of_images`，**包括**两个端点）。端点固定，中间 image 受到的力为

$$
\mathbf F_i = -\nabla E(\mathbf R_i)\big|_\perp + k\left(|\mathbf R_{i+1}-\mathbf R_i| - |\mathbf R_i-\mathbf R_{i-1}|\right)\hat\tau_i
$$

真实力只保留**垂直**于路径切向 τ̂ 的分量，弹簧力只保留**沿**路径的分量。"nudged"指的就是这种投影：弹簧不会把路径拉直抄近道，真实力也不会让 image 滑向两端。

普通 NEB 中能量最高的 image 只是离鞍点最近的离散点，能垒会被低估。CI-NEB 对能量最高的 image 去掉弹簧力，并把真实力沿切向的分量反号：

$$
\mathbf F_{\max} = -\nabla E(\mathbf R_{\max}) + 2\left(\nabla E\cdot\hat\tau\right)\hat\tau
$$

这个 image 在垂直方向往下走、沿路径往上爬，收敛后正好停在鞍点上。INPUT_NEB 中 `CI_scheme='auto'` 的说明是："能量最高的 image 不受弹簧作用，可以沿路径爬升"。

### 2.3 收敛判据

neb.x 的收敛判据是 **`path_thr`**：INPUT_NEB 定义为"垂直于路径的力的范数"小于 `path_thr`（默认 0.05 eV/Å）时停止。源码 `NEB/src/path_base.f90` 中每个 image 的 error 取该 image 梯度分量的最大绝对值（换算成 eV/Å），整条路径的 error 取所有 image 中的最大值。

### 2.4 可变弹簧常数

`k_max`、`k_min`（默认都是 0.1 Hartree 原子单位）设为不同值时，neb.x 使用可变弹簧常数：能量越高的 image 弹簧越硬，鞍点附近的 image 分布更密。

---

## 3. 准备端点：同一套设置下的严格弛豫

NEB 对端点有三条硬性要求：

1. **相同的计算设置**：`&SYSTEM`、`&ELECTRONS`、赝势、k 点、晶胞都必须与 NEB 中的 `BEGIN_ENGINE_INPUT` 完全一致，否则端点能量与中间 image 不可比；
2. **相同的固定约束**：弛豫时固定了哪些原子，NEB 中也固定同样的原子；
3. **相同的原子顺序**：第 k 个原子在初态和末态中必须是"同一个原子"，neb.x 按原子序号逐个线性插值。

### 3.1 体系与模型

Al(100) 2×2 表面晶胞（a = 4.04 Å，表面晶格常数 a/√2 ≈ 2.857 Å），3 层共 12 个 Al，底层 4 个固定；顶层 Al 到上方周期像底层 Al 的真空约 16 Å。H 从 (1/4, 1/4) 处的四重空位 A 跳到相邻的 (3/4, 1/4) 空位 B，中间经过桥位。

> [!WARNING]
> 3 层、2×2、6×6×1 k 点只是为了教学压缩的设置。正式研究需要用自己优化的晶格常数（第5章），检查 slab 层数、表面晶胞大小（H 与其周期像之间的相互作用）、k 点和真空厚度的收敛。

### 3.2 初态弛豫输入

`h_al100.initial.relax.in`：

```fortran
&CONTROL
  calculation   = 'relax'
  prefix        = 'h_al100_initial'
  outdir        = './tmp'
  pseudo_dir    = './pseudo'
  etot_conv_thr = 1.0d-5
  forc_conv_thr = 5.0d-4
  nstep         = 100
/
&SYSTEM
  ibrav       = 0
  nat         = 13
  ntyp        = 2
  ecutwfc     = 60.0
  ecutrho     = 480.0
  occupations = 'smearing'
  smearing    = 'mv'
  degauss     = 0.02
/
&ELECTRONS
  conv_thr    = 1.0d-8
  mixing_beta = 0.3
/
&IONS
  ion_dynamics = 'bfgs'
/
ATOMIC_SPECIES
  Al  26.9815  Al.pbe-n-kjpaw_psl.1.0.0.UPF
  H    1.008   H.pbe-rrkjus_psl.1.0.0.UPF
CELL_PARAMETERS angstrom
  5.71342  0.00000   0.00000
  0.00000  5.71342   0.00000
  0.00000  0.00000  20.00000
ATOMIC_POSITIONS angstrom
  Al  0.00000  0.00000  2.00000   0 0 0
  Al  0.00000  2.85671  2.00000   0 0 0
  Al  2.85671  0.00000  2.00000   0 0 0
  Al  2.85671  2.85671  2.00000   0 0 0
  Al  1.42836  1.42836  4.02000
  Al  1.42836  4.28507  4.02000
  Al  4.28507  1.42836  4.02000
  Al  4.28507  4.28507  4.02000
  Al  0.00000  0.00000  6.04000
  Al  0.00000  2.85671  6.04000
  Al  2.85671  0.00000  6.04000
  Al  2.85671  2.85671  6.04000
  H   1.42836  1.42836  6.80000
K_POINTS automatic
  6 6 1 0 0 0
```

`h_al100.final.relax.in` 只把 `prefix` 改为 `'h_al100_final'`、H 的 x 坐标改为 4.28507，其余逐字相同。

| 关键词 | 本例取值 | 说明 |
|---|---|---|
| `forc_conv_thr` | 5.0d-4 Ry/Bohr（约 0.013 eV/Å） | 比 NEB 的 `path_thr` 更严，保证端点是真正的极小点（默认 1.0d-3） |
| `etot_conv_thr` | 1.0d-5 Ry | 默认 1.0d-4；注意它是广延量 |
| `if_pos` | 底层 `0 0 0` | 冻结 slab 底部，模拟体相 |
| `ecutwfc`/`ecutrho` | 60 / 480 Ry | Al 与 H 的 SSSP 提示值取大者（H 决定） |
| `conv_thr` | 1.0d-8 | 力的精度依赖 SCF 精度，NEB 中尤其如此 |

### 3.3 运行与取出弛豫结构

```bash
cd QE/examples/11-neb
mpirun -np 16 pw.x -nk 4 -in h_al100.initial.relax.in > h_al100.initial.relax.out
mpirun -np 16 pw.x -nk 4 -in h_al100.final.relax.in   > h_al100.final.relax.out
grep -e 'bfgs converged' -e '^!' h_al100.initial.relax.out | tail -2
grep -A15 'Begin final coordinates' h_al100.initial.relax.out
```

把 `Begin final coordinates` 与 `End final coordinates` 之间的 `ATOMIC_POSITIONS` 块分别粘贴到 NEB 输入的 `FIRST_IMAGE` 和 `LAST_IMAGE` 下。

> [!TIP]
> 两个端点等价时（本例的两个四重空位由平移对称性联系），两端弛豫后的能量应几乎相同，这是检查设置是否一致的好办法。若相差超过几个 meV，先检查 k 点网格是否对两个位置同样对称。

---

## 4. neb.x 输入结构

### 4.1 超级卡片

neb.x 的输入由若干"超级卡片"（supercard）嵌套组成：

```text
BEGIN
  BEGIN_PATH_INPUT
    &PATH ... /                  <- NEB-specific namelist
    [CLIMBING_IMAGES card]       <- only for CI_scheme='manual'
  END_PATH_INPUT
  BEGIN_ENGINE_INPUT
    &CONTROL ... /               <- normal pw.x namelists
    &SYSTEM ... /
    &ELECTRONS ... /
    ATOMIC_SPECIES ...
    BEGIN_POSITIONS
      FIRST_IMAGE
        ATOMIC_POSITIONS ...
      [INTERMEDIATE_IMAGE
        ATOMIC_POSITIONS ...]    <- zero or more
      LAST_IMAGE
        ATOMIC_POSITIONS ...
    END_POSITIONS
    K_POINTS ...                 <- other pw.x cards
    CELL_PARAMETERS ...
  END_ENGINE_INPUT
END
```

INPUT_NEB 的要点：

- `BEGIN_ENGINE_INPUT` 中是普通的 pw.x namelist 和卡片，**只有** `ATOMIC_POSITIONS` 改为写在 `BEGIN_POSITIONS ... END_POSITIONS` 中，每个 `ATOMIC_POSITIONS` 前面加 `FIRST_IMAGE`、`INTERMEDIATE_IMAGE` 或 `LAST_IMAGE`；
- `INTERMEDIATE_IMAGE` 可以没有，也可以有任意多个；
- `K_POINTS`、`CELL_PARAMETERS` 等其他卡片写在 `END_POSITIONS` 之后；
- 文档特别警告：**制表符（TAB）和 DOS 的 CR/LF 换行**都可能导致问题，在 Windows 上编辑后要转换成 LF；
- neb.x 解析输入后会生成 `neb.dat`（`&PATH` 部分）和每个位置一个 `pw_1.in`、`pw_2.in` …… 文件。

### 4.2 完整输入

`h_al100.neb.in`（端点坐标需替换为弛豫后的结果）：

```fortran
BEGIN
BEGIN_PATH_INPUT
&PATH
  string_method = 'neb'
  restart_mode  = 'from_scratch'
  nstep_path    = 200
  num_of_images = 7
  opt_scheme    = 'broyden'
  CI_scheme     = 'auto'
  path_thr      = 0.05
  ds            = 1.0
  k_max         = 0.3
  k_min         = 0.2
/
END_PATH_INPUT
BEGIN_ENGINE_INPUT
&CONTROL
  prefix     = 'h_al100_neb'
  outdir     = './tmp'
  pseudo_dir = './pseudo'
/
&SYSTEM
  ibrav       = 0
  nat         = 13
  ntyp        = 2
  ecutwfc     = 60.0
  ecutrho     = 480.0
  occupations = 'smearing'
  smearing    = 'mv'
  degauss     = 0.02
/
&ELECTRONS
  conv_thr    = 1.0d-8
  mixing_beta = 0.3
/
ATOMIC_SPECIES
  Al  26.9815  Al.pbe-n-kjpaw_psl.1.0.0.UPF
  H    1.008   H.pbe-rrkjus_psl.1.0.0.UPF
BEGIN_POSITIONS
FIRST_IMAGE
ATOMIC_POSITIONS angstrom
  Al  0.00000  0.00000  2.00000   0 0 0
  Al  0.00000  2.85671  2.00000   0 0 0
  Al  2.85671  0.00000  2.00000   0 0 0
  Al  2.85671  2.85671  2.00000   0 0 0
  Al  1.42836  1.42836  4.02000
  Al  1.42836  4.28507  4.02000
  Al  4.28507  1.42836  4.02000
  Al  4.28507  4.28507  4.02000
  Al  0.00000  0.00000  6.04000
  Al  0.00000  2.85671  6.04000
  Al  2.85671  0.00000  6.04000
  Al  2.85671  2.85671  6.04000
  H   1.42836  1.42836  6.80000
LAST_IMAGE
ATOMIC_POSITIONS angstrom
  Al  0.00000  0.00000  2.00000
  Al  0.00000  2.85671  2.00000
  Al  2.85671  0.00000  2.00000
  Al  2.85671  2.85671  2.00000
  Al  1.42836  1.42836  4.02000
  Al  1.42836  4.28507  4.02000
  Al  4.28507  1.42836  4.02000
  Al  4.28507  4.28507  4.02000
  Al  0.00000  0.00000  6.04000
  Al  0.00000  2.85671  6.04000
  Al  2.85671  0.00000  6.04000
  Al  2.85671  2.85671  6.04000
  H   4.28507  1.42836  6.80000
END_POSITIONS
K_POINTS automatic
  6 6 1 0 0 0
CELL_PARAMETERS angstrom
  5.71342  0.00000   0.00000
  0.00000  5.71342   0.00000
  0.00000  0.00000  20.00000
END_ENGINE_INPUT
END
```

注意 `&CONTROL` 里**不写** `calculation`：路径计算由 neb.x 控制，pw.x 只为每个 image 计算能量和力。固定约束 `if_pos` 只需写在 `FIRST_IMAGE` 中，neb.x 从第一个 image 读取固定原子信息并应用于所有 image（官方 example01 也是这样写的）。

### 4.3 `&PATH` 关键词

| 关键词 | 默认值 | 推荐值 | 说明 |
|---|---|---|---|
| `string_method` | `'neb'` | `'neb'` | 另一选项 `'smd'`（string method dynamics） |
| `restart_mode` | `'from_scratch'` | 续算时 `'restart'` | 续算读取 `prefix.path` |
| `nstep_path` | **1** | 100–300 | 路径优化的最大步数。**默认只有 1 步，一定要设** |
| `num_of_images` | 0 | 5–9 | 路径上的点数，**包含两个端点**，必须大于 3 |
| `opt_scheme` | `'quick-min'` | `'broyden'` | 官方文档标注 broyden 为"suggested"；其他：`'broyden2'`、`'sd'`、`'langevin'`（仅 smd） |
| `CI_scheme` | `'no-CI'` | `'auto'` | `'manual'` 时用 `CLIMBING_IMAGES` 卡指定 image 序号（逗号分隔） |
| `first_last_opt` | `.false.` | `.false.` | 设为 `.true.` 时端点也在过程中优化（不受弹簧作用）；端点已严格弛豫时不需要 |
| `minimum_image` | `.false.` | 跨越周期边界时 `.true.` | 若某原子在相邻输入 image 之间移动超过半个晶轴，选择合适的周期像，避免路径"跳跃" |
| `ds` | 1.0 | 1.0 | 优化步长（Hartree 原子单位）；broyden 时作为 Jacobian 对角部分的初猜 |
| `k_max`、`k_min` | 0.1、0.1 | 如 0.3、0.2 | 弹簧常数范围（Hartree 原子单位）；两者不同时启用可变弹簧 |
| `path_thr` | 0.05 eV/Å | 0.05；精修时 0.02–0.03 | 垂直于路径的力的阈值 |
| `use_masses` | `.false.` | `.false.` | 质量加权坐标，配合 quick-min 使用 |
| `use_freezing` | `.false.` | 视情况 | 只优化误差大于最大误差一半的 image，节省计算量 |
| `lfcp`、`fcp_mu` … | `.false.` | — | 恒电势（FCP）NEB，需要 ESM，超出本章范围 |

手动指定爬坡 image 的写法（放在 `&PATH` 之后、`END_PATH_INPUT` 之前）：

```text
CLIMBING_IMAGES
  4
```

### 4.4 何时需要中间 image

neb.x 默认在 `FIRST_IMAGE` 与 `LAST_IMAGE` 之间线性插值。下列情况线性插值会给出很差甚至不物理的初始路径：

- 分子在表面上转动（线性插值会把键长压短）；
- 交换机制（两个原子互换位置，线性插值会让它们在中点重合）；
- 路径明显弯曲（例如 H 经过顶位而不是桥位）。

这时可以在中间插入一个或几个 `INTERMEDIATE_IMAGE`（例如手工摆出的桥位构型），neb.x 会在相邻给定结构之间分段插值。所有给定结构的原子数和顺序必须相同。

---

## 5. 运行与并行

### 5.1 命令行

INPUT_NEB 明确写着 **"neb.x DOES NOT READ FROM STANDARD INPUT"**，必须用 `-inp`（或 `-input`）给出输入文件：

```bash
mpirun -np 40 neb.x -ni 5 -inp h_al100.neb.in > h_al100.neb.out
```

另一种方式是事先准备好 `neb.dat` 和 `pw_1.in`、`pw_2.in` ……，然后用 `-input_images N` 运行；初学者用 `-inp` 即可。

Slurm 示例（示例目录中的 `job_neb_slurm.sh`）：

```bash
STEP=relax sbatch job_neb_slurm.sh          # relax both end points
STEP=neb NI=5 sbatch job_neb_slurm.sh       # neb.x -ni 5 -inp h_al100.neb.in
```

通用的 Slurm/PBS 模板见 `job_qe_slurm.sh`、`job_qe_pbs.sh`（用 `bash ../scripts/job_qe_slurm.sh -h` 查看用法），作业调度的基础见 [Common 第2章](../../Common/tutorials/02-job-scheduler.md)。

### 5.2 image 并行（`-ni`）

`-ni N`（`-nimage N`）把 MPI 进程分成 N 组，每组同时计算一个 image 的能量和力。经验做法：

- 端点在 `first_last_opt=.false.` 时不需要重复计算，真正每步要算的是 `num_of_images − 2` 个中间 image；本例 7 个 image 中有 5 个中间 image，`-ni 5` 让每组算一个；
- 总进程数应能被 N 整除；每组内部还可以再用 `-nk` 做 k 点并行，例如 `mpirun -np 40 neb.x -ni 5 -nk 2 -inp ...`；
- image 之间通信很少，适合跨节点；但所有组要等最慢的 image 算完才能进入下一步。

### 5.3 目录与文件

运行时，neb.x 为每个 image 在 `outdir` 下建立 `prefix_1/`、`prefix_2/` …… 子目录（源码 `NEB/src/compute_scf.f90`），每个 image 的 pw.x 输出追加写入其中的 `PW.out`：

```bash
ls tmp/
#  h_al100_neb_1/  h_al100_neb_2/  ...  h_al100_neb_7/
tail -30 tmp/h_al100_neb_4/PW.out           # SCF details of image 4
```

某个 image 的 SCF 不收敛时，先看这里。

---

## 6. 读懂输出

### 6.1 主输出

每个路径优化步结束后，主输出中会打印能垒和各 image 的状态（格式取自 `NEB/src/path_io_routines.f90`）：

```text
     activation energy (->) =   0.1...  eV
     activation energy (<-) =   0.1...  eV

     image        energy (eV)        error (eV/A)        frozen

         1      -....            0.0...            T
         2      -....            0.0...            F
         3      -....            0.0...            F
         4      -....            0.0...            F
         5      -....            0.0...            F
         6      -....            0.0...            F
         7      -....            0.0...            T

     climbing image =  4
```

（示意输出，数值仅供参考。）

- `activation energy (->)`：正向能垒 E_max − E(image 1)；`(<-)`：逆向能垒；
- `error (eV/A)`：每个 image 垂直于路径的力；
- `frozen`：该 image 本步是否被冻结（端点和 `use_freezing` 时被冻结的 image 显示 T）；
- `climbing image`：使用 CI 时，当前作为爬坡 image 的序号。

收敛时打印：

```text
     neb: convergence achieved in  NN iterations
```

达到 `nstep_path` 仍未收敛时打印 `neb: reached the maximum number of steps`。常用 grep：

```bash
grep 'activation energy' h_al100.neb.out | tail -2
grep -e 'neb: convergence achieved' -e 'reached the maximum number' h_al100.neb.out
grep 'climbing image' h_al100.neb.out | tail -1
```

### 6.2 输出文件

文件名的前缀取自 `&CONTROL` 中的 `prefix`（源码 `path_base.f90`），本例为 `h_al100_neb`：

| 文件 | 内容 |
|---|---|
| `prefix.dat` | 每个 image 一行：归一化反应坐标、E − E(image 1)（eV）、error（eV/Å） |
| `prefix.int` | 用三次插值（利用沿路径的力）得到的平滑能量曲线，两列 |
| `prefix.xyz` | 所有 image 的结构，xyz 格式（Å） |
| `prefix.axsf` | 所有 image 的结构，XCrySDen 动画格式 |
| `prefix.crd` | 所有 image 的坐标，pw.x `ATOMIC_POSITIONS` 格式（与输入单位相同） |
| `prefix.path` | 续算所需的完整信息 |
| `neb.dat`、`pw_*.in` | 输入解析后生成的文件 |

画能量剖面：

```bash
gnuplot -e "set terminal png; set output 'mep.png'; \
            set xlabel 'reaction coordinate'; set ylabel 'E - E(IS) (eV)'; \
            plot 'h_al100_neb.int' w l t 'interpolated', 'h_al100_neb.dat' u 1:2 w p pt 7 t 'images'"
```

在 XCrySDen 中打开 `h_al100_neb.axsf` 可以逐帧检查路径是否合理（H 是否经过桥位、Al 原子是否有异常移动）。爬坡 image 的坐标可以从 `prefix.crd` 中直接拷出，作为 §7 的鞍点结构。

### 6.3 成功标准

1. 出现 `neb: convergence achieved`，最后一步所有 image 的 error < `path_thr`；
2. 能量剖面平滑、只有一个峰（本例峰在 image 4）；出现多个峰时，说明存在中间态，应拆成两段分别计算；
3. 正向与逆向能垒之差等于 E(FS) − E(IS)（本例两端等价，二者应接近）；
4. 端点能量与单独弛豫得到的能量一致（不一致说明设置不同）；
5. 鞍点有且只有一个虚频（§7）。

### 6.4 续算

作业因时间限制中断时，保留 `outdir` 和 `prefix.path`，把 `restart_mode` 改为 `'restart'` 重新提交：

```fortran
&PATH
  restart_mode  = 'restart'
  ...
/
```

neb.x 会从 `prefix.path` 读入最后的路径继续优化。也可以先用较松的 `path_thr` 和 `CI_scheme='no-CI'` 得到大致路径，再续算时打开 `CI_scheme='auto'` 并收紧 `path_thr`。

---

## 7. 验证鞍点：ph.x Γ 点频率

### 7.1 思路

在爬坡 image 的结构上计算振动频率：鞍点应恰好有一个虚频（ph.x/dynmat.x 中以**负频率**表示 ω² < 0），其方向沿 IS→FS。对表面上的小吸附物，可以只位移吸附原子作为快速近似：ph.x 的 `nat_todo` 指定参与线性响应的原子。INPUT_PH 对此的原话是"近似，使用时要小心，风险自负"，所以它只用来确认虚频个数和方向。

### 7.2 输入

`h_ts.scf.in`：与 NEB 设置相同，`calculation='scf'`，坐标替换为爬坡 image 的坐标，`conv_thr` 收紧到 1.0d-12（DFPT 需要高精度的基态）。

`h_ts.ph_gamma.in`：

```text
Gamma-point modes of H only at the NEB saddle point (nat_todo approximation)
&INPUTPH
  prefix   = 'h_ts'
  outdir   = './tmp'
  tr2_ph   = 1.0d-14
  fildyn   = 'h_ts.dyn'
  nat_todo = 1
/
0.0 0.0 0.0
13
```

- 第一行是标题（ph.x 必需）；
- `0.0 0.0 0.0` 是 q 点（Γ）；
- 最后一行是要位移的原子序号列表（第 13 个原子是 H）。

`h_ts.dynmat.in`：

```fortran
&INPUT
  fildyn = 'h_ts.dyn'
  asr    = 'no'
  remove_interaction_blocks = .true.
  filout = 'h_ts.modes'
/
```

INPUT_DYNMAT 说明，`remove_interaction_blocks=.true.` 会把未位移原子对应的行和列置零，用于配合 `nat_todo`，去掉被位移部分与固定部分之间的相互作用。表面吸附物没有平移不变性，所以 `asr='no'`。

| 关键词 | 本例取值 | 说明 |
|---|---|---|
| `fildyn` | `'h_ts.dyn'` | ph.x 写出的 Γ 点动力学矩阵 |
| `asr` | `'no'`（默认） | 吸附物 + 固定 slab，不施加平移 ASR |
| `remove_interaction_blocks` | `.true.` | 默认 `.false.`；与 `nat_todo` 配合 |
| `filout` | `'h_ts.modes'` | 频率与归一化位移写到这个文件（默认 `'dynmat.out'`），**不是**标准输出 |

### 7.3 运行与检查

```bash
mpirun -np 16 pw.x -nk 4 -in h_ts.scf.in > h_ts.scf.out
mpirun -np 16 ph.x -nk 4 -in h_ts.ph_gamma.in > h_ts.ph_gamma.out
dynmat.x < h_ts.dynmat.in > h_ts.dynmat.out
grep 'freq' h_ts.modes
```

```text
     freq (    1) =      -....  [THz] =     -....  [cm-1]
     freq (    2) =       0.0...  [THz] =       0.0...  [cm-1]
     ...
     freq (   37) =       0.0...  [THz] =       0.0...  [cm-1]
     freq (   38) =       ....  [THz] =      ....  [cm-1]
     freq (   39) =       ....  [THz] =      ....  [cm-1]
```

（示意输出，数值仅供参考。）dynmat.x 总是输出 3·nat = 39 个模式：36 个来自被置零的 Al 行列，频率为 0（没有物理意义）；剩下 3 个是 H 的模式，其中应恰好有一个负值（虚频），其位移矢量应主要沿 x 方向（A→B），另外两个为正。ph.x 的输出 `h_ts.ph_gamma.out` 末尾也会打印对角化结果，但那是未去掉相互作用块的矩阵，以 dynmat.x 的结果为准。

> [!NOTE]
> 更严格的验证：位移所有可动原子（不用 `nat_todo`，计算量大很多），或者用 phonopy 等有限位移工具（见 [第12章](12-phonon-dfpt.md)）。还可以做"伪 IRC"：把鞍点结构沿虚频方向正负各推 0.05–0.1 Å，分别 `relax`，确认分别回到初态和末态。

### 7.4 零点能修正（可选）

如果对初态和 TS 都做了同样的频率计算，可以对能垒做零点能（ZPE）修正：

$$
E_a^{\mathrm{ZPE}} = E_a + \sum_{j}\tfrac{1}{2}h\nu_j^{\ddagger} - \sum_{i}\tfrac{1}{2}h\nu_i^{\mathrm{IS}}
$$

TS 的求和不包括虚频。H 这样的轻原子 ZPE 修正可能相当可观，而且同位素（D）的修正不同，这是同位素效应的来源。

---

## 参数选择与收敛建议

| 项目 | 建议 |
|---|---|
| `num_of_images` | 简单扩散 5–7 个；路径长或弯曲时 9–11 个；先粗后细，image 越多每步越贵 |
| `nstep_path` | 100–300；**不要依赖默认值 1** |
| `opt_scheme` | 首选 `'broyden'`；broyden 不稳定时换 `'quick-min'`（可配合 `use_masses`） |
| `CI_scheme` | 先用 `'no-CI'` 或直接 `'auto'`；路径很差时先跑若干步普通 NEB 再打开 CI |
| `path_thr` | 默认 0.05 eV/Å 足够得到能垒；需要精确鞍点结构和频率时收紧到 0.02–0.03 |
| `k_max`、`k_min` | 默认 0.1；想让鞍点附近 image 更密时用可变弹簧，例如 0.3/0.2 |
| SCF 精度 | `conv_thr` ≤ 1e-8；金属体系注意 smearing 和 k 点对力的影响 |
| 端点 | `forc_conv_thr` 比 `path_thr` 更严；端点能量差应与 NEB 输出一致 |
| 体系尺寸 | 吸附物与其周期像之间的距离、slab 厚度、底层固定的层数都要测试 |
| 能垒精度 | 能垒对泛函非常敏感；比较不同路径时保持完全相同的设置，论文中注明泛函与 vdW 修正 |

---

## 常见问题与报错

**1. neb.x 启动后什么也不做，或报 `Neither a file to parse nor input files for each image found`。**
原因：用 `neb.x < neb.in` 从标准输入读入（neb.x 不读标准输入），或输入文件路径不对。解决：`neb.x -inp neb.in > neb.out`。

**2. 只走了一步就结束，输出 `neb: reached the maximum number of steps`。**
原因：`nstep_path` 默认值是 1。解决：在 `&PATH` 中显式设置 `nstep_path = 200`。

**3. 读取输入时报 namelist/卡片错误，但内容看起来没有问题。**
原因：输入中有 TAB 或 Windows CR/LF 换行（INPUT_NEB 特别警告）。解决：`dos2unix neb.in`，把 TAB 换成空格；检查超级卡片名拼写（`BEGIN_PATH_INPUT`、`END_PATH_INPUT`、`BEGIN_ENGINE_INPUT`、`END_ENGINE_INPUT`、`BEGIN_POSITIONS`、`END_POSITIONS`、`END`）。

**4. `num_of_images` 报错或 image 数不对。**
原因：`num_of_images` 必须大于 3，且**包括**两个端点；从 VASP 转来的用户常把它当作中间 image 数。解决：按"中间 image 数 + 2"设置。

**5. 初始路径中原子"穿越"或彼此重叠，第一步 SCF 就不收敛。**
原因：初态和末态的原子顺序不一致，或原子跨越了周期边界。解决：检查两组坐标的原子顺序；跨边界时设 `minimum_image = .true.` 或手动把末态坐标平移到同一个周期像；必要时加 `INTERMEDIATE_IMAGE`。

**6. 某个 image 的 SCF 不收敛，整个计算停止。**
原因：插值得到的中间结构原子过近，或金属/磁性体系的 SCF 本身就难收敛。解决：查看 `outdir/prefix_N/PW.out`；降低 `mixing_beta`、增大 `electron_maxstep`；加入合理的 `INTERMEDIATE_IMAGE`；检查初末态是否合理。

**7. 能量剖面有两个或更多峰。**
原因：路径上存在中间态（局部极小）。解决：从 `prefix.crd` 中取出中间极小附近的 image，单独弛豫，把路径拆成两段分别计算 NEB。

**8. 爬坡 image 收敛很慢，error 在某个值附近振荡。**
原因：路径还很粗糙时就打开了 CI，或者 broyden 的历史累积了不好的信息。解决：先用 `CI_scheme='no-CI'` 收敛到一个较松的 `path_thr`，再 `restart_mode='restart'` 并打开 CI；或者换 `opt_scheme='quick-min'`；减小 `ds`。

**9. 续算报错读取不到路径文件。**
原因：`prefix.path` 被删除，或 `prefix` 与之前不同。解决：保持 `prefix`、`outdir` 不变，并保留 `prefix.path`；如果文件丢失，只能把最近的 `prefix.crd` 坐标作为 `INTERMEDIATE_IMAGE` 重新开始。

**10. 端点能量与单独弛豫得到的能量不一致。**
原因：NEB 的 `&SYSTEM`/`&ELECTRONS`/k 点/晶胞与弛豫时不同，或端点坐标不是弛豫后的坐标。解决：从弛豫输入复制设置，用 `Begin final coordinates` 中的坐标。

**11. ph.x 验证时出现 2 个或更多虚频，或者没有虚频。**
原因：NEB 没有真正收敛到鞍点（`path_thr` 太松），或者 `nat_todo` 近似遗漏了重要的耦合，也可能是数值噪声。解决：收紧 `path_thr` 继续优化；提高 `conv_thr`、`tr2_ph`；对小的虚频（几十 cm⁻¹ 以下）检查截断能和 SCF 精度；必要时位移更多原子。

**12. ph.x 报 `cannot start from pw.x data file using Gamma-point tricks`。**
原因：SCF 使用了 `K_POINTS gamma`。解决：改用 `K_POINTS automatic 1 1 1 0 0 0`（或正常的 k 网格）重新 SCF。

---

## 练习

1. 把 `num_of_images` 改为 5 和 9，比较能垒和收敛步数，体会 image 密度的影响。
2. 关掉 CI（`CI_scheme='no-CI'`），比较最高 image 的能量与 CI-NEB 的能垒，验证"普通 NEB 低估能垒"。
3. 分别用 `opt_scheme='broyden'` 和 `'quick-min'` 计算，记录达到相同 `path_thr` 所需的步数。
4. 把 slab 加厚到 4 层（底部 2 层固定），评估能垒随层数的变化。
5. 设计 H 经过**顶位**的路径：在 `FIRST_IMAGE` 和 `LAST_IMAGE` 之间加一个 H 位于 Al 顶位的 `INTERMEDIATE_IMAGE`，看 NEB 最终是否仍然走桥位。
6. 对初态也做 `nat_todo=1` 的 Γ 点频率计算，计算 ZPE 修正后的能垒；把 H 的质量改成 D（`ATOMIC_SPECIES` 中质量改为 2.014，ph.x 中也可以用 `amass`），比较同位素效应。
7. 把 `-ni` 分别设为 1、5，记录每个路径步的墙钟时间，理解 image 并行的效率。

---

## 参考资料

**官方文档**

- neb.x 输入说明 INPUT_NEB：<https://www.quantum-espresso.org/Doc/INPUT_NEB.html>
- pw.x 输入说明 INPUT_PW：<https://www.quantum-espresso.org/Doc/INPUT_PW.html>
- ph.x 输入说明 INPUT_PH：<https://www.quantum-espresso.org/Doc/INPUT_PH.html>；dynmat.x：<https://www.quantum-espresso.org/Doc/INPUT_DYNMAT.html>
- PWscf 用户手册（含并行与 image 说明）：<https://www.quantum-espresso.org/Doc/pw_user_guide/>
- 官方 NEB 示例：<https://gitlab.com/QEF/q-e/-/tree/develop/NEB/examples>（`example01` 为 H2 + H 交换反应，依次演示 CI-NEB、不带 CI 的 NEB，以及利用 `.path` 文件续算；另有 `ESM_example`）
- NEB 源码（输出格式与文件名）：`NEB/src/path_base.f90`、`NEB/src/path_io_routines.f90`，见 <https://gitlab.com/QEF/q-e>

**课程与教程（只链接，不复制）**

- pranabdas/espresso：<https://pranabdas.github.io/espresso/>，结构优化（`structure-optimization`）一节可作为 §3 端点弛豫的补充；站点目前没有 neb.x 专题。
- jkitchin/dft-book-espresso：<https://github.com/jkitchin/dft-book-espresso>，CMU dft-book 的 QE 版，目前只有 ASE 入门和分子两章（不含 NEB），适合想用 ASE 驱动 pw.x、再用 ASE 的 NEB 工具时入门。
- QE-SSP（《Quantum ESPRESSO Course for Solid-State Physics》配套）：<https://github.com/nguyen-group/QE-SSP>，不含 NEB，其石墨烯、MoS2 等表面/二维算例可作为 slab 建模的参考。
- Cornell PARADIM 的公开 QE 教程（2016/2018 暑校）不包含 NEB，其平衡结构与振动部分可作为 §3、§7 的背景阅读：<https://www.paradim.org/cu_ss_2018>
- MaX e-School 2021（ICTP/MaX/QEF）：<https://indico.ictp.it/event/9616/>，含 QE 高级功能与并行讲座视频。

**本仓库**

- 示例：[11-neb](../examples/11-neb/README.md)
- 脚本：`qe_parse_out.py`（解析弛豫输出与最终结构）、`job_qe_slurm.sh`；用 `python ../scripts/<name> -h` 查看用法。
- 对照阅读：[VASP 第11章 NEB 与 Dimer](../../VASP/tutorials/11-neb-dimer.md)、[CP2K 第13章 NEB 与 Dimer](../../CP2K/tutorials/13-neb-dimer.md)

**关键文献**

- G. Henkelman, B. P. Uberuaga, H. Jónsson, *J. Chem. Phys.* **113**, 9901 (2000)，CI-NEB：<https://doi.org/10.1063/1.1329672>
- G. Henkelman, H. Jónsson, *J. Chem. Phys.* **113**, 9978 (2000)，改进切向估计：<https://doi.org/10.1063/1.1323224>
- D. Sheppard, R. Terrell, G. Henkelman, *J. Chem. Phys.* **128**, 134106 (2008)，MEP 优化算法比较：<https://doi.org/10.1063/1.2841941>
- S. Smidstrup, A. Pedersen, K. Stokbro, H. Jónsson, *J. Chem. Phys.* **140**, 214106 (2014)，IDPP 初始路径：<https://doi.org/10.1063/1.4878664>
- G. H. Vineyard, *J. Phys. Chem. Solids* **3**, 121 (1957)，谐振子过渡态理论：<https://doi.org/10.1016/0022-3697(57)90059-8>
- P. Giannozzi et al., *J. Phys.: Condens. Matter* **21**, 395502 (2009)，QE：<https://doi.org/10.1088/0953-8984/21/39/395502>
- P. Giannozzi et al., *J. Phys.: Condens. Matter* **29**, 465901 (2017)：<https://doi.org/10.1088/1361-648X/aa8f79>

---

[← 上一章](10-charge-elf-bader.md) | [返回目录](../README.md) | [下一章 →](12-phonon-dfpt.md)
