# 第17章 分子动力学：pw.x MD 与 Car-Parrinello cp.x（Ab Initio MD: Born-Oppenheimer with pw.x and Car-Parrinello with cp.x）

本章介绍 Quantum ESPRESSO 中两种从头算分子动力学（ab initio molecular dynamics, AIMD）。一种是 `pw.x` 的 Born-Oppenheimer MD（`calculation='md'`/`'vc-md'`），涉及时间步 `dt`（Rydberg 原子单位）、积分器 `ion_dynamics='verlet'`/`'velocity-verlet'`、恒温器 `ion_temperature`（`'svr'`、`'berendsen'`、`'andersen'`、`'rescaling'`、`'nose'` 等）及其参数 `tempw`、`nraise`、`fnosep`，以及外推选项 `pot_extrapolation`/`wfc_extrapolation`。另一种是 `cp.x` 的 Car-Parrinello 动力学，涉及 `electron_dynamics='damp'` → `'verlet'`、虚拟电子质量 `emass`/`emass_cutoff`、Hartree 原子单位的 `dt`、`ndr`/`ndw` 重启目录，以及 `.evp` 文件里的 `ekinc`/`econs`/`econt`。最后介绍如何用 ASE、dpdata 读取轨迹，为机器学习势（DeePMD、NEP 等）准备训练数据。

> [!NOTE]
> 本章输入在 [examples/17-md-cp](../examples/17-md-cp/README.md)。仓库**没有运行**这些计算，文中输出片段均标注"（示意输出，数值仅供参考）"。关键词对照的是官方 [INPUT_PW](https://www.quantum-espresso.org/Doc/INPUT_PW.html)、[INPUT_CP](https://www.quantum-espresso.org/Doc/INPUT_CP.html)、CP 用户指南（`CPV/Doc/user_guide.md`）和 QE 7.6 源码。

---

## 1. 本章目标与适用场景

| 目标 | 你将学会 |
|---|---|
| BOMD 输入 | 写出 `pw.x` 的 NVE/NVT MD 输入，正确换算 `dt` 的单位 |
| 恒温器 | 区分速度重标定类（`rescaling`、`rescale-v`、`berendsen`）、随机类（`svr`、`andersen`）和 Nose-Hoover（`nose`）恒温器，知道哪些给出正则系综 |
| 效率 | 用 `pot_extrapolation`、`wfc_extrapolation` 减少每步 SCF 迭代 |
| 变胞 MD | 了解 `vc-md` 的 `cell_dynamics='pr'`/`'w'`、`press`、`wmass` |
| 检查 | 读 `temperature`、`Ekin + Etot (const)`，判断能量漂移是否可接受 |
| CP 动力学 | 走通 `cp.x` 的"阻尼电子基态 → Verlet CP 动力学"流程，选择 `emass` 与 `dt`，读 `.evp` |
| 重启 | `max_seconds` + `restart_mode='restart'`，`prefix.EXIT` 文件，`ndr`/`ndw` |
| 后处理 | 用 ASE 读 pw.x 轨迹，用 dpdata 读 cp.x 轨迹，为 MLFF 采集数据 |

适用场景：液体与熔体（液态 Si、水）、有限温度下的结构与扩散、非谐性显著的材料、反应体系的短时间采样，以及为机器学习势生成构型和 DFT 标签。

与 VASP、CP2K 的对应关系（概念对照，数值不能直接照搬）：

| 概念 | QE `pw.x` | QE `cp.x` | VASP | CP2K |
|---|---|---|---|---|
| MD 开关 | `calculation='md'` | `calculation='cp'` | `IBRION=0` | `RUN_TYPE MD` |
| 时间步 | `dt`（Ry a.u.，4.8378×10⁻¹⁷ s） | `dt`（Ha a.u.，2.4189×10⁻¹⁷ s） | `POTIM`（fs） | `TIMESTEP`（fs） |
| 步数 | `nstep` | `nstep` | `NSW` | `STEPS` |
| 温度 | `tempw` | `tempw` | `TEBEG`/`TEEND` | `TEMPERATURE` |
| 恒温器 | `ion_temperature` | `ion_temperature='nose'`/`'rescaling'` | `MDALGO` + `SMASS` 等 | `&THERMOSTAT` |
| 变胞 | `calculation='vc-md'` | `calculation='vc-cp'` | `ISIF=3` + `MDALGO=3` | `ENSEMBLE NPT_I/NPT_F` |
| 波函数外推 | `wfc_extrapolation` | 不需要（电子随离子一起演化） | 默认的外推 | `EXTRAPOLATION ASPC` |

VASP 的 MD 见 [VASP 第15章](../../VASP/tutorials/15-md-mlff.md)，CP2K 的 AIMD 见 [CP2K 第12章](../../CP2K/tutorials/12-aimd.md)。如果主要做大体系、长时间 AIMD，CP2K 的 Gaussian + 平面波混合基组通常更快；QE 的优势在于和本系列其它章节使用同一套平面波赝势和设置。

---

## 2. 原理速览

### 2.1 Born-Oppenheimer MD

BOMD 在每个离子构型 R(t) 上把电子完全收敛到基态，再用 Hellmann-Feynman 力推进离子：

$$
M_I\ddot{\mathbf R}_I=-\nabla_{\mathbf R_I}\min_{\{\psi\}}E_{KS}[\{\psi\},\{\mathbf R\}].
$$

`pw.x` 默认用 Verlet 算法（`ion_dynamics='verlet'`），也可以选 `'velocity-verlet'`。积分器本身是时间可逆、辛的；实际计算中 SCF 不可能完全收敛，力里的噪声会导致总能量漂移。因此在 NVE 中，**`conv_thr` 决定能量守恒的质量**，这一点比时间步的影响更容易被忽视。

### 2.2 Car-Parrinello MD

CP 方法不在每一步做 SCF，而是给 Kohn-Sham 轨道一个虚拟质量 μ（`emass`），把轨道当作经典动力学变量与离子一起演化：

$$
\mathcal L=\sum_i\mu\langle\dot\psi_i|\dot\psi_i\rangle+\sum_I\frac12M_I\dot{\mathbf R}_I^2-E_{KS}[\{\psi\},\{\mathbf R\}]+\sum_{ij}\Lambda_{ij}(\langle\psi_i|\psi_j\rangle-\delta_{ij}).
$$

只要电子的特征频率远高于离子的振动频率（**绝热分离**），电子就会紧贴 Born-Oppenheimer 面运动，虚拟动能 K_e（`.evp` 中的 `ekinc`）保持很小。要点：

- `emass` 越小，绝热性越好、力的系统误差越小，但电子的最高频率 ∝ 1/√μ 越高，`dt` 必须越小；
- CP 用户指南指出，`emass` 取约 50 a.u. 一般不会明显影响动力学性质；对静态性质（RDF、晶格常数）影响更小；
- CP 中电子带着虚拟质量"拖"着离子运动，离子的有效质量偏大；算扩散系数这类动力学量时需要注意（参见 P. Tangney 2006 的分析）；
- 有能隙体系最适合 CP。金属的能隙为零，绝热分离会被破坏，这时需要电子恒温器或者改用 BOMD。

### 2.3 时间步与单位

| 程序 | `dt` 单位 | 1 a.u. | 示例 |
|---|---|---|---|
| `pw.x` | Rydberg 原子单位 | 4.8378×10⁻¹⁷ s | `dt=20` ≈ 0.97 fs，`dt=40` ≈ 1.94 fs |
| `cp.x` | Hartree 原子单位 | 2.4189×10⁻¹⁷ s | `dt=5` ≈ 0.12 fs，`dt=10` ≈ 0.24 fs |

官方文档特地提醒两者相差一倍（INPUT_PW 中 `dt` 的说明：*beware, the CP code uses Hartree atomic units, half that much!!!*）。经验上 BOMD 的步长取体系最快振动周期的 1/10–1/20：含 H 体系约 0.5 fs（`pw.x` 中 `dt≈10`），只含重元素的固体 1–2 fs（`dt≈20–40`）。CP 的步长受电子虚拟动力学限制，通常比 BOMD 小 5–10 倍。

### 2.4 恒温器概览

| `ion_temperature` | 参数 | 系综 | 说明 |
|---|---|---|---|
| `'not_controlled'`（默认） | — | NVE | 不控温 |
| `'initial'` | `tempw` | NVE | 只在开始时按 `tempw` 给随机速度 |
| `'rescaling'` | `tempw`、`tolp` | 非正则 | 平均温度偏离 `tempw` 超过 `tolp` 时重标定速度；适合升温/平衡 |
| `'rescale-v'` | `tempw`、`nraise` | 非正则 | 每 `nraise` 步把最近 `nraise` 步的平均温度重标定到 `tempw` |
| `'rescale-T'`、`'reduce-T'` | `tempw`、`delta_t`、`nraise` | 非正则 | 程序化升温/降温（模拟退火、淬火） |
| `'berendsen'` | `tempw`、`nraise` | 非正则 | τ = `nraise`×`dt`；平均温度正确，但温度涨落被压低，不产生正则分布 |
| `'andersen'` | `tempw`、`nraise` | 正则 | 随机碰撞，ν·dt = 1/`nraise`；扰乱动力学，不适合算扩散 |
| `'svr'` | `tempw`、`nraise` | 正则 | Bussi-Donadio-Parrinello 随机速度重标定，τ = `nraise`×`dt`；**推荐的默认选择** |
| `'nose'` | `tempw`、`fnosep`、`nhpcl`、`nhptyp` 等 | 正则 | Nose-Hoover（链）恒温器；QE 7.5 起 pw.x 正式支持（与 Verlet / velocity-Verlet 配合） |

> [!NOTE]
> QE 7.5 的 release notes 写明：pw.x 新增了基于 Nose-Hoover 恒温器的 NVT（Verlet 和 velocity-Verlet 引擎），以及 Nose-Hoover 恒温器 + Parrinello-Rahman/Wentzcovitch 压浴的 NPT（Beeman 引擎）。若你用的是 7.4 或更早的版本，`pw.x` 的 `'nose'` 选项可能不可用或行为不同，请以对应版本的 INPUT_PW 为准。

---

## 3. 实战一：64 原子 Si 的 BOMD（pw.x）

### 3.1 体系与输入

2×2×2 惯用胞的金刚石结构 Si 超胞（64 原子，简立方 `ibrav=1`），1000 K NVT。完整输入见 [si64.md.in](../examples/17-md-cp/si64.md.in)，坐标部分省略：

```fortran
&CONTROL
  calculation = 'md'
  prefix      = 'si64'
  outdir      = './tmp'
  pseudo_dir  = './pseudo'
  dt          = 40.0          ! Ry a.u. (~1.94 fs)
  nstep       = 1000          ! ~1.9 ps; first run a short test (e.g. 20 steps)
  iprint      = 1
  tprnfor     = .true.
  max_seconds = 82800         ! stop cleanly before a 24 h wall-time limit
  disk_io     = 'low'
/
&SYSTEM
  ibrav       = 1
  celldm(1)   = 20.52         ! Bohr = 2 x 10.26 (2x2x2 conventional cells); use your relaxed a
  nat         = 64
  ntyp        = 1
  ecutwfc     = 30.0
  ecutrho     = 240.0
  occupations = 'smearing'    ! safe if the system becomes metallic (hot/liquid Si)
  smearing    = 'gaussian'
  degauss     = 0.01
/
&ELECTRONS
  conv_thr    = 1.0d-8        ! about 1e-10 Ry/atom; loose conv_thr gives energy drift
  mixing_beta = 0.5
/
&IONS
  ion_dynamics      = 'verlet'
  ion_temperature   = 'svr'          ! Bussi-Donadio-Parrinello thermostat
  tempw             = 1000.0         ! K
  nraise            = 25             ! tau = nraise*dt (~48 fs)
  pot_extrapolation = 'second_order'
  wfc_extrapolation = 'second_order'
/
ATOMIC_SPECIES
  Si  28.0855  Si.pbe-n-rrkjus_psl.1.0.0.UPF
ATOMIC_POSITIONS crystal
  Si  0.000000  0.000000  0.000000
  ...                                  ! 64 lines in total, see the example file
K_POINTS gamma
```

### 3.2 关键词表

| 关键词 | 默认值 | 推荐值 | 说明 |
|---|---|---|---|
| `calculation` | `'scf'` | `'md'` | `'vc-md'` 为变胞 MD |
| `dt` | 20.0 | 重元素 20–40；含 H 约 10 | **Ry 原子单位** |
| `nstep` | 50（md） | 视需要 | 设为 0 时只做初始化，可用来检查输入（"dry run"） |
| `ion_dynamics` | `'verlet'`（md） | `'verlet'` 或 `'velocity-verlet'` | 另有过阻尼 `'langevin'`、`'langevin-smc'`；vc-md 只能用 `'beeman'` |
| `ion_temperature` | `'not_controlled'` | NVT 用 `'svr'` | 见 2.4 节 |
| `tempw` | 300.0 | 目标温度（K） | 也是初始速度分布的温度 |
| `nraise` | 1 | SVR/Berendsen 取 τ = `nraise`×`dt` ≈ 50–200 fs | 含义随恒温器而变 |
| `ion_velocities` | `'default'` | 一般不改 | `'from_input'` 时从 `ATOMIC_VELOCITIES` 卡读入速度（单位 a.u.） |
| `pot_extrapolation` | `'atomic'` | `'second_order'` | `'first_order'`/`'second_order'` 只对 MD 有意义 |
| `wfc_extrapolation` | `'none'` | `'second_order'` | 同上；能明显减少每步 SCF 次数 |
| `conv_thr` | 1e-6 | 每原子 1e-9–1e-11 Ry（按体系大小换算） | 广延量；决定力的噪声和 NVE 漂移 |
| `refold_pos` | `.false.` | 视后处理需要 | `.true.` 时每步把原子折回晶胞；算 MSD 时保持 `.false.` |
| `remove_rigid_rot` | `.false.` | 孤立分子可设 `.true.` | 去掉周期像的力矩；注意此时总能量**不再守恒** |
| `max_seconds` | 1e7 | 墙钟时间减去 10–30 分钟 | 到时间后程序正常退出，可以续算 |
| `K_POINTS gamma` | — | 大超胞 | 只用 Γ 点时调用专门优化的例程，内存和时间约减半 |

> [!TIP]
> 64 原子 Si 用 Γ 点采样只是演示。液态 Si 是金属，需要 `occupations='smearing'`；晶体 Si 有能隙。示例从晶体出发并直接加热到 1000 K（仍低于 Si 的熔点约 1687 K），所以加了一个较小的 `degauss` 以防万一。熔化液态 Si 需要更高的温度和更长的模拟时间。

### 3.3 运行

```bash
# dry run: stop after initialisation, check memory and the parallel setup
sed 's/nstep       = 1000/nstep       = 0/' si64.md.in > dry.in
mpirun -np 32 pw.x -in dry.in > dry.out
grep -e "Estimated max dynamical RAM" -e "Parallelization info" -A3 dry.out

# production (Gamma point only: no k-point pools)
mpirun -np 32 pw.x -in si64.md.in > si64.md.out
```

Γ 点计算只有一个 k 点，不能用 `-nk` pool 并行（`Gamma-only calculations not allowed with pools`），只能用平面波/FFT 并行（必要时加 `-nt`、`-nd`），见 [第18章](18-automation-parallel-troubleshooting.md)。

### 3.4 读输出

每个 MD 步，`pw.x` 打印（格式取自 `PW/src/dynamics_module.f90`；示意输出，数值仅供参考）：

```text
     Entering Dynamics:    iteration =    25
                           time      =   0.0484 pico-seconds

     kinetic energy (Ekin) =           0.xxxxxxxx Ry
     temperature           =         xxx.xxxxxxxx K
     Ekin + Etot (const)   =        -xxx.xxxxxxxx Ry
```

使用 `ion_temperature='nose'` 时还会多两行：`Ions Nose Energy` 和 `Ekin + Etot + Ions Nose`，后者才是 Nose-Hoover 动力学的守恒量。

```bash
grep "temperature           =" si64.md.out | awk '{print $3}' > T.dat
grep "Ekin + Etot (const)"     si64.md.out | awk '{print $6}' > Econst.dat
grep -c "Entering Dynamics" si64.md.out
grep "convergence has been achieved" si64.md.out | awk '{print $6}' | sort -n | uniq -c   # SCF iterations per step
```

也可以用示例目录中的只读分析脚本：

```bash
python ../examples/17-md-cp/md_analyze.py pw si64.md.out --skip 200 --natoms 64
```

（示意输出，数值仅供参考）

```text
MD steps parsed      : 1000 (skipping first 200)
T mean +- std        : xxxx.x +- xx.x K
conserved-E drift    : x.xxxx meV/atom/ps
conserved-E range    : xx.xxx meV (whole cell)
```

注意：NVT 中 `Ekin + Etot (const)` 不守恒（恒温器在交换能量），漂移只有在 NVE 中才有诊断意义。

### 3.5 NVE 检查

[si64.nve.in](../examples/17-md-cp/si64.nve.in) 把 `ion_temperature` 改为 `'initial'`，并把 `conv_thr` 收紧到 1e-10。判据：

- `Ekin + Etot (const)` 在一个常数附近振荡，没有线性漂移；经验上漂移小于约 1 meV/atom/ps 可以接受（取决于用途）；
- 漂移明显时，先收紧 `conv_thr`，再减小 `dt`，最后检查截断能（Pulay 力在 MD 中不是主要问题，但 `ecutrho` 太低会给出噪声）；
- 从平衡好的 NVT 构型开始 NVE 更合理。可以用 ASE 读 NVT 输出的最后一帧，写成新的 `ATOMIC_POSITIONS`；速度用 `ion_velocities='from_input'` + `ATOMIC_VELOCITIES` 卡传入，或者直接让 `'initial'` 重新生成。

### 3.6 变胞 MD（vc-md）简介

| 关键词 | 选项/默认值 | 说明 |
|---|---|---|
| `calculation` | `'vc-md'` | 需要 `&CELL` namelist |
| `ion_dynamics` | `'beeman'`（vc-md 的默认值，也是唯一选项） | |
| `cell_dynamics` | `'none'`、`'pr'`（Parrinello-Rahman）、`'w'`（Wentzcovitch） | Beeman 积分 |
| `press` | 0.0（kbar） | 目标压强 |
| `wmass` | PR：0.75×总质量/π²；W：0.75×总质量/π²/Ω^(2/3) | 虚拟晶胞质量（amu） |
| `cell_factor` | 变胞计算默认 2.0 | 晶胞线性收缩的上限，用于构建赝势插值表 |
| `cell_dofree` | `'all'` | 可以限制可动的晶胞自由度，如 `'2Dxy'` |

vc-md 必须用较高的截断能（Pulay 应力），并且要先单独测试压强的收敛。QE 7.5 起可以把 Nose-Hoover 恒温器与 PR/Wentzcovitch 压浴组合成 NPT。

### 3.7 重启

`pw.x` 只有在**正常停止**之后才能续算（PW 用户指南 "Restarting" 一节）。正常停止的方式有：

1. 到达 `max_seconds`；
2. 用户在工作目录或 `outdir` 中创建 `prefix.EXIT` 文件（例如 `touch tmp/si64.EXIT`）；
3. （实验性）编译时加了 `-D__TERMINATE_GRACEFULLY`，并收到 SIGTERM/SIGUSR1 等信号。

续算时把 `restart_mode` 设为 `'restart'`，**使用相同的进程数和并行参数**，其它输入不变。示例中的 [job_md_slurm.sh](../examples/17-md-cp/job_md_slurm.sh) 在时间到之前 10 分钟收到 `USR1` 信号时会创建 `EXIT` 文件。Slurm 作业脚本的写法见 [Common 第2章](../../Common/tutorials/02-job-scheduler.md)。

> [!WARNING]
> 不要在 `outdir` 里残留旧的 `prefix.EXIT` 文件：程序启动时如果发现它，会在初始化后立刻停止。示例作业脚本在开头会删除这些文件。

---

## 4. 实战二：H₂O 分子的 Car-Parrinello 动力学（cp.x）

### 4.1 CP 模拟是一串运行

CP 用户指南强调，一次 CP 模拟由多次运行组成：

```text
h2o.cp_damp.in    electrons -> ground state (damped dynamics, ions fixed)   writes tmp/h2o_50.save
       |
h2o.cp_verlet.in  CP dynamics: electrons and ions with Verlet (NVE)        reads 50, writes 51
       |
(next runs)       restart from 51, optionally switch on a Nose thermostat, collect statistics
```

`cp.x` 的重启目录为 `outdir/prefix_ndw.save/`（写）和 `outdir/prefix_ndr.save/`（读），`ndr`、`ndw` 的默认值都是 50。

### 4.2 第一步：阻尼动力学求电子基态

```fortran
&CONTROL
  calculation   = 'cp'
  title         = 'H2O molecule, electronic minimisation'
  restart_mode  = 'from_scratch'
  prefix        = 'h2o'
  outdir        = './tmp'
  pseudo_dir    = './pseudo'
  ndw           = 50
  nstep         = 300           ! number of damped electronic steps (upper bound)
  iprint        = 10
  isave         = 100
  dt            = 5.0           ! Hartree a.u. (~0.12 fs)
  tprnfor       = .true.
  etot_conv_thr = 1.0d-6        ! Ha
  ekin_conv_thr = 1.0d-6        ! fictitious electronic kinetic energy, Ha
/
&SYSTEM
  ibrav     = 1
  celldm(1) = 20.0              ! Bohr (about 10.6 Angstrom) cubic box
  nat       = 3
  ntyp      = 2
  ecutwfc   = 85.0              ! Ry; PseudoDojo "normal" hint for O is 42 Ha = 84 Ry
/
&ELECTRONS
  electron_dynamics = 'damp'
  electron_damping  = 0.2
  emass             = 400.0     ! fictitious electron mass (a.u.)
  emass_cutoff      = 2.5       ! Ry, Fourier acceleration
/
&IONS
  ion_dynamics = 'none'
/
ATOMIC_SPECIES
  O  15.999  O.upf
  H   2.014  H.upf              ! deuterium mass: allows a larger dt
ATOMIC_POSITIONS angstrom
  O   5.290  5.290  5.410
  H   5.290  6.090  4.790      ! O-H = 1.01 Angstrom (slightly stretched, so the molecule vibrates)
  H   5.290  4.530  4.820      ! O-H = 0.96 Angstrom
```

| 关键词（cp.x） | 默认值 | 说明 |
|---|---|---|
| `calculation` | `'cp'` | 其它：`'scf'`、`'relax'`、`'vc-relax'`、`'vc-cp'`、`'cp-wf'` 等 |
| `restart_mode` | **`'restart'`** | 注意 cp.x 的默认值与 pw.x 不同；第一次运行必须写 `'from_scratch'` |
| `ndr` / `ndw` | 50 / 50 | 读/写重启目录 `prefix_ndr.save`、`prefix_ndw.save` |
| `nstep` | 50 | 阻尼/最速下降时是电子步数；CP 动力学时是 MD 步数 |
| `iprint` | 10 | 每隔多少步写一次轨迹文件；标准输出每 10×`iprint` 步写一次 |
| `isave` | 100 | 每隔多少步保存一次重启信息 |
| `dt` | 1.0 | **Hartree 原子单位** |
| `ekin_conv_thr` | 1e-6 | 电子虚拟动能小于该值认为到达基态（与 `etot_conv_thr` 同时满足） |
| `electron_dynamics` | `'none'` | `'sd'`、`'damp'`、`'verlet'`、`'cg'` |
| `electron_damping` | 0.1 | 阻尼系数，常用 0.1–0.5 |
| `emass` | 400.0 | 虚拟电子质量（a.u.） |
| `emass_cutoff` | 2.5 | Fourier 加速的截断（Ry），高于此动能的 G 分量使用重标定后的质量 |
| `nr1b, nr2b, nr3b` | — | 超软赝势的"盒子"网格，**有 US 赝势时必须设置**；本例用模守恒赝势，不需要 |

> [!IMPORTANT]
> `cp.x` 不支持 PAW 赝势（读入时会报 `PAW not implemented`）。CP 计算请用模守恒或超软赝势。

CP 用户指南推荐用 `electron_dynamics='cg'` 初始化，它还能给出与离子速度一致的波函数"速度"，开始 CP 动力学时更平稳（目前 CG 只支持平面波并行）。本例用阻尼动力学，是为了展示经典的"damp → verlet"流程。

检查（示意输出，数值仅供参考）：

```bash
grep -A3 "nfi" h2o.cp_damp.out | head      # columns include ekinc and etot
tail -3 tmp/h2o.evp
ls tmp/h2o_50.save/
```

`ekinc` 应单调下降，最后低于 `ekin_conv_thr`。

### 4.3 第二步：Car-Parrinello 动力学

```fortran
&CONTROL
  calculation  = 'cp'
  title        = 'H2O molecule, CP dynamics'
  restart_mode = 'restart'
  prefix       = 'h2o'
  outdir       = './tmp'
  pseudo_dir   = './pseudo'
  ndr          = 50
  ndw          = 51
  nstep        = 2000           ! 2000 x 5 a.u. = 0.24 ps
  iprint       = 10             ! trajectory files every 10 steps
  isave        = 500
  dt           = 5.0
  tprnfor      = .true.
/
&SYSTEM
  ibrav     = 1
  celldm(1) = 20.0
  nat       = 3
  ntyp      = 2
  ecutwfc   = 85.0
/
&ELECTRONS
  electron_dynamics   = 'verlet'
  electron_velocities = 'zero'  ! first CP run only
  emass               = 400.0
  emass_cutoff        = 2.5
/
&IONS
  ion_dynamics   = 'verlet'
  ion_velocities = 'zero'       ! first CP run only; the stretched O-H bond starts the motion
/
ATOMIC_SPECIES
  O  15.999  O.upf
  H   2.014  H.upf
ATOMIC_POSITIONS angstrom
  O   5.290  5.290  5.410
  H   5.290  6.090  4.790
  H   5.290  4.530  4.820
```

| 关键词 | 选项 | 说明 |
|---|---|---|
| `electron_velocities` | `'zero'`、`'default'`、`'change_step'` | 阻尼/最速下降之后的第一次 CP 运行用 `'zero'`；之后的续算**删掉这一行**（默认读入上次的速度），否则每次重启都会"淬火" |
| `ion_velocities` | `'default'`、`'zero'`、`'random'`、`'from_input'`、`'change_step'` | `'random'` 配合 `tempw` 从 Maxwell-Boltzmann 分布取初速度 |
| `ion_temperature` | `'not_controlled'`、`'nose'`、`'rescaling'` | Nose 恒温器配合 `fnosep`（THz）、`tempw`、`nhpcl` |
| `electron_temperature` | `'not_controlled'`、`'nose'`、`'rescaling'` | 电子恒温器（`ekincw`、`fnosee`），用法微妙，见 CP 用户指南 |

`emass=400` 配合氘（2.014 amu）和 `dt=5` 是常见的保守组合。CP 用户指南提到 `emass≈50` 时对动力学性质影响很小，但此时 `dt` 也要相应减小。改变 `dt` 续算时，从 QE 7.3 起只需设置 `electron_velocities='change_step'`，代码会自动处理（不必再用 `tolp` 指定旧步长）。

### 4.4 读 `.evp` 和轨迹文件

`cp.x` 每 `iprint` 步把数据写入 `outdir` 下的若干文件（CP 用户指南 "Output files"）：

| 文件 | 内容与单位 |
|---|---|
| `prefix.evp` | 热力学量，一行一个输出步 |
| `prefix.pos` | 每个输出步：第一行为步数和时间（ps），其后为坐标（Bohr，未折回晶胞）；v6.6 起按输入顺序排列 |
| `prefix.vel` | 速度（Hartree a.u.；pw.x 的速度用 Rydberg a.u.，相差 2 倍） |
| `prefix.for` | 力（Hartree a.u.），MD 时或 `tprnfor=.true.` 时写出 |
| `prefix.cel` | 晶胞矩阵（**转置**：晶胞矢量在列上） |
| `prefix.str` | 应力张量（GPa，含离子动能部分） |

`.evp` 的列：

```text
#   nfi  time(ps)  ekinc  Tcell(K)  Tion(K)  etot  enthal  econs  econt  Volume  Pressure(GPa)
```

| 列 | 含义 | 期望行为 |
|---|---|---|
| `ekinc` | 电子虚拟动能 K_e | CP 动力学中保持很小，没有持续增长 |
| `Tion` | 离子温度 | 在平衡值附近涨落 |
| `etot` | DFT 势能 E_DFT | 与离子动能来回交换 |
| `econs` | E_DFT + K_ions（零虚拟质量极限下的物理守恒量） | 小幅振荡；平均值无明显漂移 |
| `econt` | E_DFT + K_ions + K_e（CP 拉格朗日量的守恒量） | 在 `dt` 足够小时几乎严格不变 |

```bash
awk '!/^#/ {print $2, $3, $8, $9}' tmp/h2o.evp > evp_cols.dat     # time, ekinc, econs, econt
python ../examples/17-md-cp/md_analyze.py evp tmp/h2o.evp
```

（示意输出，数值仅供参考）

```text
rows parsed          : 200
Tion mean +- std     : xxx.x +- xxx.x K
ekinc mean (max)     : x.xxxe-0x (x.xxxe-0x) Ha
econs drift          : x.xxxx meV/ps (whole cell)
econt drift          : x.xxxx meV/ps (whole cell)
```

判断：`econt` 漂移大，说明 `dt` 太大；`ekinc` 持续增长、`econs` 与 `econt` 越来越偏离，说明绝热性被破坏（`emass` 太大、体系能隙太小或金属）。处理办法是减小 `emass`/`dt`，或者用几步 `electron_dynamics='cg'` 把电子拉回 BO 面后再继续 `'verlet'`（INPUT_CP 对 `'cg'` 选项的说明）。

### 4.5 从单分子到液态水（概念）

单个 H₂O 分子只用来演示流程。液态水的 CP 模拟通常用 32–64 个水分子的立方盒子（密度约 1 g/cm³），流程不变：先用 `electron_dynamics='damp'`（或 `'cg'`）得到电子基态，再用 `'verlet'` 做 CP 动力学，然后打开 Nose 恒温器平衡，最后在 NVE 或 NVT 下采样。与单分子相比要注意：体系变大后 `nr1b/nr2b/nr3b`（若用超软赝势）、`emass`、`dt` 的选择方法相同，但计算量随电子数快速增长，应先用短轨迹检查 `ekinc` 与 `econt`；GGA（如 PBE）描述的液态水结构偏"过度结构化"，这是泛函本身的问题，不是 CP 参数的问题。本仓库没有给出液态水的输入文件。

### 4.6 控温（Nose 恒温器）

在 `&IONS` 中加入（CP 用户指南的写法）：

```fortran
&IONS
  ion_dynamics    = 'verlet'
  ion_temperature = 'nose'
  fnosep          = 60.0      ! THz; comparable to the centre of the vibrational spectrum
  tempw           = 300.0
/
```

CP 用户指南的建议：不要从完全弛豫的构型直接打开 Nose 恒温器（先加一点随机位移，例如 `tranp`/`amprp`）；先用几步 NVE 看平均温度；分子体系有多个模式需要热化，最好用恒温器链（`nhpcl`）。

---

## 5. 轨迹后处理与 MLFF 数据采集

### 5.1 用 ASE 读 pw.x 轨迹

ASE 的 `espresso-out` 读取器可以从 `pw.x` 输出中读出每一步的结构、能量和力：

```python
from ase.io import read, write

frames = read("si64.md.out", index=":", format="espresso-out")
print(len(frames), frames[-1].get_potential_energy())   # eV
write("si64_md.extxyz", frames, format="extxyz")         # positions, energy, forces
```

示例脚本 `md_analyze.py pw ... --xyz out.extxyz` 做的就是这件事。之后可以用 ASE、[MDAnalysis](https://github.com/MDAnalysis/mdanalysis) 或 OVITO 计算 RDF、MSD 等（ASE/pymatgen 的基础用法见 [Common 第4章](../../Common/tutorials/04-ase-pymatgen.md)）。

> [!TIP]
> 算 MSD 时要用"未折回"的坐标。`pw.x` 默认 `refold_pos=.false.`，不会折回；`cp.x` 的 `.pos` 也是未折回的坐标。

### 5.2 为机器学习势准备数据

AIMD 是生成 MLFF 训练数据（构型 + DFT 能量/力/维里）的常用方式。注意事项：

- **dpdata**（DeePMD-kit 的数据工具）有 QE 的格式插件：`qe/pw/scf`（单点 pw.x 输出）和 `qe/cp/traj`（cp.x 轨迹）。`pw.x` 的 MD 输出可以先用 ASE 转成 extxyz，再用 dpdata 的 ASE/extxyz 接口读入。格式名以 [dpdata 文档](https://docs.deepmodeling.com/projects/dpdata/) 为准；
- **标签一致性**：训练集所有构型必须用同一套泛函、赝势、截断能、smearing，最好用更严格的 `conv_thr` 单独重算单点（`calculation='scf'`，`tprnfor`、`tstress`），而不是直接用 MD 过程中的能量和力；
- **CP 轨迹的力**含有虚拟电子质量带来的系统误差，用作标签前应在 BO 面上重算；
- **去相关**：相邻 MD 步高度相关，按自相关时间抽样（例如每 10–50 步取一帧），或用主动学习（DP-GEN 等）挑选构型；
- 需要应力标签时，要在 SCF 单点中设置 `tstress=.true.`，并确认截断能足以收敛应力（[第4章](04-convergence.md)）。

VASP 的 MLFF 流程对照见 [VASP 第15章](../../VASP/tutorials/15-md-mlff.md)。

---

## 参数选择与收敛建议

| 参数 | 建议 | 检查方法 |
|---|---|---|
| `dt`（pw.x） | 含 H 约 10（≈0.5 fs）；重元素 20–40 | NVE 漂移随 `dt` 的变化 |
| `dt`（cp.x） | 3–10（Ha a.u.），与 `emass` 一起选 | `econt` 守恒 |
| `emass` | 300–400 起步，追求精度时降到约 50–100（同时减小 `dt`） | `ekinc` 小且稳定；力与 BO 力对比 |
| `conv_thr`（BOMD） | 每原子约 1e-10 Ry | NVE 漂移 |
| 外推 | `pot_extrapolation='second_order'`、`wfc_extrapolation='second_order'` | 每步 SCF 迭代数 |
| 恒温器 | 平衡与生产都可用 `'svr'`，τ ≈ 50–200 fs | 温度分布宽度符合 √(2/(3N))·T |
| 初始结构 | 先弛豫，或从经典 MD/MLFF 预平衡 | 前几步温度不出现暴涨 |
| 超胞 | 关联长度之上；液体至少 50–100 原子 | 改变超胞大小，比较 RDF |
| 模拟时长 | 平衡 ≥ 1 ps，生产按待测性质的关联时间 | 分块平均的误差 |
| k 点 | 大超胞用 Γ（`K_POINTS gamma`）；金属/小胞要测 | 与 2×2×2 k 网格对比力 |

---

## 常见问题与报错

| 报错/现象 | 原因 | 解决 |
|---|---|---|
| 温度在开头几步暴涨 | 初始结构力很大（未弛豫），或 `dt` 太大 | 先 `relax`；减小 `dt`；先用 `'rescaling'` 平衡 |
| NVE 下 `Ekin + Etot (const)` 线性漂移 | `conv_thr` 太松或 `dt` 太大 | 收紧 `conv_thr`（每原子 1e-10 Ry 量级），再减小 `dt` |
| 每步 SCF 迭代很多 | 没有外推；金属体系 smearing 太小 | 设置 `wfc_extrapolation`/`pot_extrapolation='second_order'`；调整 `mixing_beta`、`degauss` |
| `Error in routine setup: Gamma-only calculations not allowed with pools` | `K_POINTS gamma` 时用了 `-nk` | 去掉 `-nk` |
| `convergence NOT achieved after 100 iterations: stopping` 出现在 MD 中途 | 构型变得很差（原子过近）或体系变成金属 | 检查轨迹；`occupations='smearing'`；减小 `dt`；增大 `electron_maxstep` |
| `Error in routine iosys: unknown ion_temperature ...` | 拼写错误，或版本不支持该恒温器 | 对照当前版本的 INPUT_PW |
| `restart_mode='restart'` 没有接上上次的轨迹（输出中出现 `Message from routine iosys: restart disabled: needed files not found`，或读重启数据时报错） | 上次没有正常停止（被 kill）、`outdir`/`prefix` 变了，或换了进程数/并行方式 | 用 `max_seconds` 或 `prefix.EXIT` 正常停止；保持 `prefix`、`outdir` 和并行参数不变；否则从最后一帧重新开始 |
| 作业刚开始就停止，输出中没有 MD 步 | `outdir` 或工作目录中残留了 `prefix.EXIT` | 删除 EXIT 文件 |
| cp.x：`Error in routine readpp: PAW not implemented` | 用了 PAW 赝势 | 换成模守恒或超软赝势 |
| cp.x：超软赝势运行出错或结果异常 | 没有设置盒子网格 `nr1b, nr2b, nr3b` | 按 INPUT_CP 设置盒子网格 |
| cp.x：第一次就要从 `prefix_50.save` 读数据而失败 | `restart_mode` 默认是 `'restart'` | 第一次运行写 `restart_mode='from_scratch'` |
| cp.x：`ekinc` 持续上升 | 绝热性被破坏（`emass` 太大、`dt` 太大、金属/小能隙） | 减小 `emass` 和 `dt`；插入几步 `electron_dynamics='cg'`；金属改用 BOMD |
| cp.x：orthogonalization 不收敛 | `dt` 太大，或初始波函数太差 | 减小 `dt`；调整 `ortho_eps`/`ortho_max`；开头几步用 `orthogonalization='Gram-Schmidt'`（INPUT_CP 说明它只能用于最初几步，能量和本征值不正确） |
| cp.x：重启后系统被"冻住" | 续算时仍保留了 `electron_velocities='zero'` / `ion_velocities='zero'` | 第一次 CP 运行之后删掉这两行 |
| pw.x 与 cp.x 的步长或速度对不上 | Ry 与 Ha 原子单位相差 2 倍 | 按 2.3 节换算 |

---

## 练习

1. 在 `si64.nve.in` 中分别取 `dt=20`、`40`、`60`，各跑 200 步，比较 `Ekin + Etot (const)` 的振荡幅度和漂移。
2. 固定 `dt=40`，把 `conv_thr` 从 1e-6 依次收紧到 1e-10，比较 NVE 漂移和每步平均 SCF 迭代数，找出性价比最好的设置。
3. 比较 `wfc_extrapolation='none'` 与 `'second_order'` 的每步 SCF 迭代数和总耗时。
4. 用 `'berendsen'` 和 `'svr'`（相同 τ）各跑一段 NVT，比较温度分布的宽度，并与正则系综的理论值 √(2/(3N))·T 对照。
5. 在 `h2o.cp_verlet.in` 中把 `emass` 分别设为 100、400、800（相应调整 `dt`），比较 `ekinc` 的大小和 O-H 伸缩振动的频率（从 `.pos` 中的键长时间序列做 Fourier 变换）。
6. 用 ASE 把 Si 的 NVT 轨迹转成 extxyz，计算 Si-Si 径向分布函数，并与理想晶体的近邻距离对比。

---

## 参考资料

**官方文档**

- [INPUT_PW](https://www.quantum-espresso.org/Doc/INPUT_PW.html)：`&CONTROL`（`calculation`、`dt`、`nstep`、`max_seconds`、`restart_mode`）、`&IONS`（`ion_temperature`、`tempw`、`nraise`、`fnosep`、`pot_extrapolation`、`wfc_extrapolation`）、`&CELL`
- [INPUT_CP](https://www.quantum-espresso.org/Doc/INPUT_CP.html)：`electron_dynamics`、`emass`、`emass_cutoff`、`ndr`/`ndw`、`ion_velocities` 等
- CP 用户指南：[CPV/Doc/user_guide.md](https://github.com/QEF/q-e/blob/develop/CPV/Doc/user_guide.md)（输出文件格式、阻尼 → Verlet 流程、`emass` 的选择、Nose 恒温器）；autopilot 动态改参数见同目录 `autopilot_guide.md`
- [PWscf User's Guide](https://www.quantum-espresso.org/Doc/pw_user_guide/)（重启、EXIT 文件、信号处理）
- [QE release notes](https://gitlab.com/QEF/q-e/-/blob/develop/Doc/release-notes)：7.5 新增 pw.x 的 Nose-Hoover NVT/NPT
- 官方算例：`PW/examples/`（MD 示例）与 `CPV/examples/`（[GitHub 镜像](https://github.com/QEF/q-e)）

**课程与教程（只链接）**

- [MaX e-School 2021](https://indico.ictp.it/event/9616/)：QE 官方暑校，含 MD 相关讲座视频
- [Imperial MSE404-MM](https://imperialcollegelondon.github.io/MSE404-MM/)：Lab 9（Finite Temperature Properties）讲的是基于声子的有限温度性质（准谐近似），不是 AIMD，可作为与 MD 方法对照的背景阅读
- [ENCCS: Efficient materials modelling on HPC](https://enccs.github.io/efficient-materials-modelling-on-hpc/)：大体系 QE 并行设置，对长时间 AIMD 很有用
- Cornell PARADIM 暑校（[2016](https://www.paradim.org/cu_ss_2016)、[2018](https://www.paradim.org/cu_ss_2018)）：以 pw.x/ph.x 为主，没有 MD 上机，可作为 DFT 基础补充

**相关仓库（只链接）**

- [QEF/q-e](https://github.com/QEF/q-e)：`PW/examples`、`CPV/examples`
- [aiidateam/qe-tools](https://github.com/aiidateam/qe-tools)：纯 Python 解析 pw.x/cp.x 输入输出
- [elcorto/pwtools](https://github.com/elcorto/pwtools)：解析 QE MD 轨迹、计算声子 DOS 等
- [deepmodeling/dpdata](https://github.com/deepmodeling/dpdata)：`qe/pw/scf`、`qe/cp/traj` 格式，转换为 DeePMD 训练集
- [MDAnalysis/mdanalysis](https://github.com/MDAnalysis/mdanalysis)：通用轨迹分析

**本仓库脚本**：`qe_parse_out.py`（解析 pw.x 输出）、`job_qe_slurm.sh`、`job_qe_pbs.sh`（作业模板），用 `python ../scripts/<name> -h` 查看用法；示例目录中的 `md_analyze.py` 用 `-h` 查看用法。

**关键文献**

- R. Car, M. Parrinello, *Phys. Rev. Lett.* **55**, 2471 (1985)，Car-Parrinello 方法. DOI: [10.1103/PhysRevLett.55.2471](https://doi.org/10.1103/PhysRevLett.55.2471)
- G. Bussi, D. Donadio, M. Parrinello, *J. Chem. Phys.* **126**, 014101 (2007)，SVR 恒温器. DOI: [10.1063/1.2408420](https://doi.org/10.1063/1.2408420)
- S. Nosé, *J. Chem. Phys.* **81**, 511 (1984). DOI: [10.1063/1.447334](https://doi.org/10.1063/1.447334)
- W. G. Hoover, *Phys. Rev. A* **31**, 1695 (1985). DOI: [10.1103/PhysRevA.31.1695](https://doi.org/10.1103/PhysRevA.31.1695)
- M. Parrinello, A. Rahman, *Phys. Rev. Lett.* **45**, 1196 (1980)，变胞 MD. DOI: [10.1103/PhysRevLett.45.1196](https://doi.org/10.1103/PhysRevLett.45.1196)
- P. Tangney, *J. Chem. Phys.* **124**, 044111 (2006)，CP 方法的理论基础与虚拟质量的作用. DOI: [10.1063/1.2162893](https://doi.org/10.1063/1.2162893)
- P. Giannozzi et al., *J. Chem. Phys.* **152**, 154105 (2020)，QE 面向百亿亿次计算的进展. DOI: [10.1063/5.0005082](https://doi.org/10.1063/5.0005082)
- D. Marx, J. Hutter, *Ab Initio Molecular Dynamics: Basic Theory and Advanced Methods*, Cambridge University Press (2009)

---

[← 上一章](16-elph-epw.md) | [返回目录](../README.md) | [下一章 →](18-automation-parallel-troubleshooting.md)
