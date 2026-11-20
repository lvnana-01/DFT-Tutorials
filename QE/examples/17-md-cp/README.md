# 示例 17：pw.x 分子动力学与 cp.x Car-Parrinello 动力学

配套教程：[第17章 分子动力学：pw.x MD 与 Car-Parrinello cp.x](../../tutorials/17-md-cp.md)

> [!IMPORTANT]
> 本目录只提供输入、作业脚本和一个文本分析脚本，**仓库作者没有运行过这些计算**。README 与教程中的输出片段都是"（示意输出，数值仅供参考）"。

## 目的

1. `si64.md.in`：64 原子 Si 超胞（2×2×2 惯用胞，金刚石结构），`pw.x` Born-Oppenheimer MD，SVR 恒温器，1000 K，NVT。
2. `si64.nve.in`：同一体系的 NVE 运行（`ion_temperature='initial'`），用来检查 `Ekin + Etot (const)` 守恒。
3. `h2o.cp_damp.in` → `h2o.cp_verlet.in`：单个 H₂O 分子的 `cp.x` 流程：先用阻尼动力学把电子弛豫到基态，再做 Car-Parrinello（Verlet）动力学，检查 `.evp` 文件中的 `ekinc`、`econs`、`econt`。

## 文件说明

| 文件 | 说明 |
|---|---|
| `si64.md.in` | pw.x，`calculation='md'`，`dt=40`（Ry 原子单位，约 1.94 fs），`ion_temperature='svr'`、`tempw=1000`、`nraise=25`，二阶势/波函数外推，`K_POINTS gamma` |
| `si64.nve.in` | 同上，但 `ion_temperature='initial'`（只在开始时按 `tempw` 给随机速度），`conv_thr=1e-10` |
| `h2o.cp_damp.in` | cp.x，`electron_dynamics='damp'`，`ion_dynamics='none'`，写出 `tmp/h2o_50.save` |
| `h2o.cp_verlet.in` | cp.x，`restart_mode='restart'`，`ndr=50`/`ndw=51`，电子与离子都用 `'verlet'`，NVE |
| `md_analyze.py` | 只读文本的分析脚本：pw.x 输出的温度统计与守恒量漂移（可选用 ASE 导出 extxyz），cp.x `.evp` 文件的漂移与 `ekinc` |
| `job_md_slurm.sh` | Slurm 作业脚本，参数 `pw` / `pwnve` / `cp`；时间到之前 10 分钟写 `prefix.EXIT`，让程序正常退出 |

## 赝势准备（不提交 UPF 文件）

- Si（pw.x）：[SSSP 1.3 PBE Efficiency](https://www.materialscloud.org/discover/sssp/table/efficiency) 中的 `Si.pbe-n-rrkjus_psl.1.0.0.UPF`（超软），放到 `pseudo/`。
- H、O（cp.x）：PseudoDojo 模守恒赝势（cp.x **不支持 PAW**，会报 `PAW not implemented`；模守恒赝势也不需要设置 `nr1b/nr2b/nr3b`）：

```bash
mkdir -p pseudo && cd pseudo
wget https://www.pseudo-dojo.org/pseudos/nc-sr-05_pbe_standard/H.upf.gz && gunzip H.upf.gz
wget https://www.pseudo-dojo.org/pseudos/nc-sr-05_pbe_standard/O.upf.gz && gunzip O.upf.gz
# SSSP file for Si: extract Si.pbe-n-rrkjus_psl.1.0.0.UPF from the SSSP 1.3 PBE efficiency tarball
cd ..
```

赝势的选择与核对方法见 [第3章](../../tutorials/03-pseudopotentials.md)。

## 运行步骤

```bash
cd QE/examples/17-md-cp

# 0) quick checks: nstep = 0 makes pw.x stop right after initialisation
sed 's/nstep       = 1000/nstep       = 0/' si64.md.in > dry.in
mpirun -np 4 pw.x -in dry.in > dry.out
grep "Estimated max dynamical RAM" dry.out

# 1) pw.x NVT (and later NVE)
sbatch job_md_slurm.sh pw
sbatch job_md_slurm.sh pwnve

# 2) cp.x damped minimisation followed by CP dynamics
sbatch job_md_slurm.sh cp
```

如果作业因时间限制在中途停下（`max_seconds` 或 `prefix.EXIT`），把对应输入中的 `restart_mode` 改为 `'restart'`（pw.x 的默认值是 `'from_scratch'`），**不要改进程数和并行参数**，然后重新提交。

分析（只读文本，不运行任何计算）：

```bash
python md_analyze.py -h
python md_analyze.py pw si64.md.out --skip 200 --natoms 64
python md_analyze.py pw si64.nve.out --natoms 64 --xyz si64_nve.extxyz   # needs ASE
python md_analyze.py evp tmp/h2o.evp
```

## 预期结果 / 如何检查

pw.x 每个 MD 步都会打印（格式见 `PW/src/dynamics_module.f90`；示意输出，数值仅供参考）：

```text
     Entering Dynamics:    iteration =    10
                           time      =   0.0194 pico-seconds

     kinetic energy (Ekin) =           0.xxxxxxxx Ry
     temperature           =         xxx.xxxxxxxx K
     Ekin + Etot (const)   =        -xxx.xxxxxxxx Ry
```

- NVT（SVR）：平衡后平均温度应接近 `tempw`；64 原子体系的温度涨落相对值约为 √(2/(3N)) ≈ 10%，即 1000 K 时 ±100 K 左右，这是正常的。
- NVE：`Ekin + Etot (const)` 应在一个常数附近小幅振荡，没有系统漂移。若漂移明显，先收紧 `conv_thr`，再减小 `dt`。
- 每个 MD 步的 SCF 迭代数（`convergence has been achieved in N iterations`）在外推生效后应明显少于第一步。

cp.x 的 `tmp/h2o.evp`（列名见 CP 用户指南）：

```text
#   nfi  time(ps)  ekinc  Tcell(K)  Tion(K)  etot  enthal  econs  econt  Volume  Pressure(GPa)
```

- 阻尼步：`ekinc` 应单调下降到 `ekin_conv_thr`（1e-6 Ha）以下；
- Verlet 步：`econt` 应几乎不变；`econs` 允许有小的振荡；`ekinc` 应保持很小且没有持续增长（绝热条件）。`ekinc` 持续增长说明电子从离子"吸热"，要减小 `emass` 或 `dt`。
- 轨迹文件 `tmp/h2o.pos`、`h2o.vel`、`h2o.for`、`h2o.cel` 的单位是 Bohr 和 Hartree 原子单位。

## 注意事项

- **时间单位不同**：pw.x 的 `dt` 用 Rydberg 原子单位（4.8378×10⁻¹⁷ s），cp.x 用 Hartree 原子单位（2.4189×10⁻¹⁷ s）。同样写 `dt=20`，pw.x 的步长是 cp.x 的两倍。
- `si64.md.in` 的晶格常数 `celldm(1)=20.52` Bohr 是占位值（2 × 10.26），请换成你用同一赝势优化得到的值（[第5章](../../tutorials/05-relax-eos.md)）。
- `h2o.cp_*.in` 中 H 用了氘的质量（2.014），这样可以用更大的 `dt`；若关心动力学性质，请改回 1.008 并相应减小 `dt`。
- 64 原子 Si 的 Γ 点采样只适合演示；做定量的液态/非晶 Si 需要更大的超胞或多个 k 点，并测试 `ecutwfc`。
- 第一次 CP 运行之后，下一次续算前一定要删掉 `electron_velocities='zero'` 和 `ion_velocities='zero'`，否则每次重启都会把体系"冻住"。
