# 示例 16：电声耦合与超导（ph.x + lambda.x 路线 / EPW 路线）

配套教程：[第16章 电声耦合与超导：ph.x 与 EPW 入门](../../tutorials/16-elph-epw.md)

> [!IMPORTANT]
> 本目录只提供输入文件和作业脚本，**仓库作者没有运行过这些计算**。文中出现的所有输出片段都是"（示意输出，数值仅供参考）"，请以你自己的计算结果为准。

## 目的

用两条路线计算金属的电声耦合常数 λ、对数平均频率 ω_log 和超导转变温度 T_c：

1. `phx-al/`：fcc Al，纯 QE 路线：`pw.x`（`la2F=.true.`）→ `ph.x`（`electron_phonon='interpolated'`）→ `q2r.x`/`matdyn.x`（`la2F=.true.`）→ `lambda.x`（Allen-Dynes 公式）。结构与官方 `PHonon/examples/example03` 一致，但输入是重新写的：改用 PBE 模守恒赝势，k/q 网格也加密了。
2. `epw-pb/`：fcc Pb，EPW 路线：`pw.x` → `ph.x`（保存 dV_SCF）→ `pp.py` → `pw.x` nscf（完整 k 网格）→ `epw.x`（Wannier 化）→ `epw.x`（细网格插值 + 各向同性 Migdal-Eliashberg 方程）。流程参照 [EPW 官方教程 "Superconducting properties"](https://docs.epw-code.org/tutorials/tutorial_04/index.html)，参数做了缩减，便于入门。

## 文件说明

| 文件 | 说明 |
|---|---|
| `phx-al/al.scf.fit.in` | 稠密 k 网格（24×24×24）SCF，`la2F=.true.`，在 `outdir` 写出 `a2Fsave` |
| `phx-al/al.scf.in` | 声子用 k 网格（12×12×12）SCF，`prefix`/`outdir` 与上一步相同 |
| `phx-al/al.elph.in` | `ph.x`：`ldisp` 4×4×4 q 网格，`fildvscf`，`electron_phonon='interpolated'`，`el_ph_sigma=0.005`，`el_ph_nsigma=10` |
| `phx-al/al.q2r.in` | `q2r.x`，`la2F=.true.` 同时变换线宽 |
| `phx-al/al.matdyn.disp.in` | `matdyn.x` 沿 Γ–X–W–L–Γ 输出频率与声子线宽 |
| `phx-al/al.matdyn.dos.in` | `matdyn.x` 在 24³ q 网格上算声子 DOS 与 α²F（`a2F.dos1`…`a2F.dos10`、`lambda` 文件） |
| `phx-al/al.lambda.in` | `lambda.x` 输入：q 点与权重、`elph_dir/elph.inp_lambda.*` 文件名、μ* |
| `phx-al/job_al_elph_slurm.sh` | Slurm 串行执行上述各步 |
| `epw-pb/pb.scf.in` | Pb SCF（8×8×8 k） |
| `epw-pb/pb.ph.in` | Pb DFPT，6×6×6 q 网格，`fildvscf='dvscf'` |
| `epw-pb/pb.q2r.in` | 力常数（供检查声子色散用） |
| `epw-pb/pb.nscf.in.head` | nscf 输入**头部**，没有 K_POINTS，需用 `kmesh.pl 6 6 6` 追加 216 个 k 点 |
| `epw-pb/pb.epw1.in` | EPW 第 1 步：Wannier 化（4 个 sp3 WF，跳过 5 条 5d 带），写出 Wannier 表象数据 |
| `epw-pb/pb.epw2.in` | EPW 第 2 步：`epwread`，24³/12³ 细网格，`a2f_iso` + `liso` + `limag` + `lpade` |
| `epw-pb/job_pb_epw_slurm.sh` | Slurm 串行执行 EPW 全流程 |

## 赝势准备（不提交 UPF 文件）

两条路线都使用 PseudoDojo 的模守恒赝势（`nc-sr-05_pbe_standard`，标量相对论，PBE）。下载后放在各自的 `pseudo/` 目录：

```bash
# Al route
mkdir -p phx-al/pseudo && cd phx-al/pseudo
wget https://www.pseudo-dojo.org/pseudos/nc-sr-05_pbe_standard/Al.upf.gz && gunzip Al.upf.gz
cd ../..

# Pb route (the inputs use pseudo_dir = '../pseudo', i.e. a shared directory next to epw-pb/)
mkdir -p pseudo && cd pseudo
wget https://www.pseudo-dojo.org/pseudos/nc-sr-05_pbe_standard/Pb.upf.gz && gunzip Pb.upf.gz
cd ..
```

下载后用 `grep -i "z_valence\|functional\|pseudo_type" phx-al/pseudo/Al.upf pseudo/Pb.upf` 核对价电子数（Al 3，Pb 14）和泛函。截断能请参考 PseudoDojo 网页给出的 hint，并按 [第4章](../../tutorials/04-convergence.md) 做收敛测试。赝势的选择与下载见 [第3章](../../tutorials/03-pseudopotentials.md)。

## 运行步骤

### 路线 A：ph.x + lambda.x（Al）

```bash
cd QE/examples/16-elph-epw/phx-al
sbatch job_al_elph_slurm.sh         # or run the commands inside the script by hand
```

运行到 `ph.x` 之后，先核对 `al.lambda.in` 里的 q 点表：

```bash
ls al.dyn*                                    # al.dyn0 + al.dyn1 ... al.dynN (N = number of irreducible q)
head -2 al.dyn0                               # q grid and N
for i in elph_dir/elph.inp_lambda.*; do head -1 $i; done   # q (2pi/a) of each file
grep -c "Dynamical  Matrix in cartesian axes" al.dyn3      # number of q in the star = weight
```

`lambda.x` 需要的是"q 点 + 星中 q 点数作为权重"，顺序必须和 `elph.inp_lambda.1…N` 一致。本例给出的 8 个 q 点和权重（1, 8, 4, 6, 24, 12, 3, 6，总和 64 = 4³）与官方 `PHonon/examples/example03` 中 fcc Al 4×4×4 网格的列表相同，但仍请用你的 `elph.inp_lambda.*` 第一行和 `al.dynN` 核对顺序；`lambda.x` 只用权重做求和，不检查坐标。

### 路线 B：EPW（Pb）

`epw-pb/` 中所有输入都用 `outdir='./'`：EPW 自带的 `pp.py` 只在当前目录里寻找 `_ph0/`、`pb.save/` 和 `pb.dyn*`，所以 `pw.x`、`ph.x`、`pp.py`、`epw.x` 都必须在同一个目录里运行（作业脚本就是这样做的）。

```bash
cd QE/examples/16-elph-epw/epw-pb
export QE_SRC=/path/to/q-e-qe-7.6               # contains EPW/bin/pp.py
export KMESH=/path/to/wannier90/utility/kmesh.pl
sbatch job_pb_epw_slurm.sh
```

`epw.x` 需要单独编译：在 QE 源码目录执行 `make epw`（会同时编译 pw.x、ph.x 和 wannier90）。安装细节见 [第1章](../../tutorials/01-install.md) 和 [EPW 下载安装页](https://docs.epw-code.org/DownloadAndInstall.html)。

## 预期结果 / 如何检查

路线 A：

```bash
grep "lambda(" al.elph.out | head              # per-mode lambda and gamma for each q
grep "Broadening" lambda                        # matdyn: lambda, dos(Ef), omega_ln for each sigma
tail -12 al.lambda.out                          # lambda.x: lambda, <log w>, Tc for each sigma
```

`al.lambda.out` 末尾的格式（示意输出，数值仅供参考）：

```text
     lambda = 0.xxxxxx (   0.xxxxxx )  <log w>=  xxx.xxx K  N(Ef)=  x.xxxxxx at degauss= 0.005
     ...
lambda        omega_log          T_c
   0.xxxxx       xxx.xxx            x.xxx
   ...
```

检查要点：

- 10 个 `degauss`（0.005…0.050 Ry）给出的 λ 应在中间一段形成平台；太小的 σ 噪声大，太大的 σ 会抹平费米面结构。取平台值，不要取单个 σ 的结果。
- `lambda.x` 括号内的数值（对 α²F 积分得到的 λ）应与前面的 λ 接近，否则检查 `emax` 是否覆盖全部声子频率、`degaussq` 是否合理。
- 文献中 Al 的 λ 约为 0.4，T_c 实验值约 1.2 K。你的结果对 k/q 网格和 μ* 都很敏感，量级对即可，**不要把本例当作收敛结果**。

路线 B：

```bash
grep -A6 "Wannier Function centers" pb.epw1.out      # 4 sp3-like WFs, spreads
grep "Estimated Tc using McMillan" pb.epw2.out         # quick Allen-Dynes estimate
ls pb.a2f pb.imag_iso_*                                # alpha^2F and gap on the imaginary axis
```

- `pb.a2f`：第 1 列频率（meV），第 2 列 α²F(ω)，第 3 列累积耦合常数 λ(ω)（EPW Tutorial 4 的说明）。
- `pb.imag_iso_XXX.XX`：每个温度一个文件，各列依次为 Matsubara 频率、重整化函数 Z、能隙 Δ（eV）。最低 Matsubara 频率处的 Δ 随温度降到 0 的温度就是 T_c。
- Wannier 插值能带（EPW 第 1 步可用 `wdata` 让 Wannier90 输出 `pb_band.dat`）要与 DFT 能带在冻结窗口内重合，这是后续所有结果可信的前提。

## 注意事项

- `celldm(1)`、`ecutwfc` 都是占位值：先用同一赝势做 `vc-relax`（第5章）和截断能收敛（第4章）。
- `la2F=.true.` 的 SCF 和声子 SCF 必须用相同的 `prefix`、`outdir`，并且先跑稠密网格那一步；顺序反了会覆盖 `a2Fsave` 以外的数据，或者 ph.x 读不到 `a2Fsave`。
- `ph.x` 的电声计算**不支持** `-ni`（image 并行），只能用 `-nk`/`-nd`；EPW 则可以用 `-nk` 与 `-ni`（nprocs = npool × nimage）。
- 半导体（如 Si）不能直接用 `electron_phonon='interpolated'`，ph.x 会报 `Electron-phonon only for metals`。
- EPW 的输入关键字随版本变化：本例按 EPW 6.1（随 QE 7.6 发布）文档写，旧版 EPW 4/5 教程中的部分写法已不适用。
