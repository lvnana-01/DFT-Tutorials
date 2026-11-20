# QE 辅助脚本（Quantum ESPRESSO helper scripts）

本目录收录 Quantum ESPRESSO（QE）计算的前处理、收敛测试、作业提交和后处理脚本：从结构文件生成 pw.x 输入（[`qe_inputgen.py`](qe_inputgen.py)）、输入静态检查（[`qe_check_input.py`](qe_check_input.py)）、ecutwfc/k 点/degauss 收敛扫描（[`qe_conv_test.sh`](qe_conv_test.sh)）、pw.x 输出解析（[`qe_parse_out.py`](qe_parse_out.py)）、能带/DOS/PDOS 作图（[`qe_plot_bands.py`](qe_plot_bands.py)、[`qe_plot_dos.py`](qe_plot_dos.py)）、声子色散（[`qe_phonon_plot.py`](qe_phonon_plot.py)）、状态方程拟合（[`qe_eos_fit.py`](qe_eos_fit.py)）、功函数（[`qe_workfunction.py`](qe_workfunction.py)）、cube 差分电荷（[`qe_cube_diff.py`](qe_cube_diff.py)），以及 Slurm/PBS 作业模板。解析的输出格式按 QE GitLab 仓库 develop 分支（7.x）源码中的 `FORMAT` 语句核对，涉及的源文件在各节列出。

> [!NOTE]
> 所有 Python 脚本都支持 `-h`（[`qe_cube_diff.py`](qe_cube_diff.py) 用 `python qe_cube_diff.py <子命令> -h`），bash 脚本用 `bash <脚本> -h` 显示文件头部的用法说明。本 README 中的参数均取自各脚本的 argparse 定义或 bash 选项解析；命令行示例只演示用法，没有附带运行结果，脚本未在真实计算上回归测试，第一次使用时请先在小算例上核对。

## 总览

| 脚本 | 功能 | 依赖 | 对应章节 |
|---|---|---|---|
| [qe_inputgen.py](qe_inputgen.py) | 从 CIF/POSCAR 等生成 pw.x 输入：`ibrav=0`、显式单位、按 SSSP JSON 设置赝势与截断、按 k 间距设网格、可选自旋与 crystal_b 路径 | ase（numpy）；可选 seekpath | [第2章](../tutorials/02-input-anatomy.md)、[第3章](../tutorials/03-pseudopotentials.md)、[第6章](../tutorials/06-band-dos-pdos.md)、[第18章](../tutorials/18-automation-parallel-troubleshooting.md) |
| [qe_check_input.py](qe_check_input.py) | 静态检查 pw.x 输入：namelist 的 `/`、CRLF、单位缺失、7.1 起删除的 DFT+U 旧关键字、废弃关键字、ecutrho/ecutwfc、nat/ntyp 行数、赝势文件 | 仅标准库 | [第2章](../tutorials/02-input-anatomy.md)、[第4章](../tutorials/04-convergence.md)、[第8章](../tutorials/08-dftu-hp.md)、[第18章](../tutorials/18-automation-parallel-troubleshooting.md) |
| [qe_conv_test.sh](qe_conv_test.sh) | 模板占位符替换，批量生成 ecut / kpt / smear 扫描输入，可选运行 pw.x，汇总 `!` 总能与压强 | bash、sed、awk；`--run` 需 pw.x | [第4章](../tutorials/04-convergence.md) |
| [qe_parse_out.py](qe_parse_out.py) | 解析 pw.x 输出：总能、费米能/HOMO-LUMO、SCF 历史、力、应力、磁矩、优化最终结构、JOB DONE，导出 JSON/CSV | 仅标准库 | [第2章](../tutorials/02-input-anatomy.md)、[第4章](../tutorials/04-convergence.md)、[第5章](../tutorials/05-relax-eos.md)、[第7章](../tutorials/07-magnetism-soc.md)、[第18章](../tutorials/18-automation-parallel-troubleshooting.md) |
| [qe_plot_bands.py](qe_plot_bands.py) | 读取 bands.x 的 `.gnu` 或 `data-file-schema.xml`，费米能对齐、标注高对称点、给出直接/间接带隙 | numpy、matplotlib | [第6章](../tutorials/06-band-dos-pdos.md)、[第7章](../tutorials/07-magnetism-soc.md)、[第14章](../tutorials/14-hybrid-metagga.md)、[第15章](../tutorials/15-wannier90.md) |
| [qe_plot_dos.py](qe_plot_dos.py) | 读取 dos.x 与 projwfc.x 的 `pdos_atm#` 文件，按元素/轨道/原子求和，自旋分辨作图 | numpy、matplotlib | [第6章](../tutorials/06-band-dos-pdos.md)、[第7章](../tutorials/07-magnetism-soc.md)、[第8章](../tutorials/08-dftu-hp.md) |
| [qe_phonon_plot.py](qe_phonon_plot.py) | 读取 matdyn.x 的 `.gp` 色散与声子 DOS，换算 cm⁻¹/THz/meV，标记虚频 | numpy、matplotlib | [第12章](../tutorials/12-phonon-dfpt.md)、[第13章](../tutorials/13-dielectric-ir-raman.md)、[第16章](../tutorials/16-elph-epw.md) |
| [qe_eos_fit.py](qe_eos_fit.py) | 由 (V,E) 或 (a,E) 或一组 pw.x 输出拟合三阶 Birch-Murnaghan / Murnaghan，给出 V0/a0、B0、B0' | numpy、scipy、matplotlib | [第5章](../tutorials/05-relax-eos.md) |
| [qe_workfunction.py](qe_workfunction.py) | 读取 average.x 的 `avg.dat` 或 pp.x cube，求真空能级与功函数 Φ = V_vac − E_F | numpy、matplotlib | [第9章](../tutorials/09-vdw-surface-workfunction.md) |
| [qe_cube_diff.py](qe_cube_diff.py) | cube 文件信息与积分、差分电荷（任意线性组合）、平面平均与 ΔQ(z) | numpy；作图需 matplotlib | [第10章](../tutorials/10-charge-elf-bader.md) |
| [job_qe_slurm.sh](job_qe_slurm.sh) | Slurm 模板：环境加载、OpenMP、`srun pw.x -nk/-nd/-nt/-ni`、`ESPRESSO_TMPDIR`，结束后检查 JOB DONE | bash、Slurm | [第1章](../tutorials/01-install.md)、[第18章](../tutorials/18-automation-parallel-troubleshooting.md) |
| [job_qe_pbs.sh](job_qe_pbs.sh) | PBS Pro / OpenPBS / Torque 模板，功能同上（`mpirun -machinefile $PBS_NODEFILE`） | bash、PBS/Torque | [第18章](../tutorials/18-automation-parallel-troubleshooting.md) |

## 安装依赖与通用约定

```bash
# Python >= 3.9
pip install numpy matplotlib          # plot_bands / plot_dos / phonon_plot / workfunction / cube_diff
pip install scipy                     # qe_eos_fit.py
pip install ase                       # qe_inputgen.py
pip install seekpath                  # qe_inputgen.py --calculation bands --kpath seekpath (pulls spglib)
# qe_parse_out.py and qe_check_input.py use only the standard library
```

或者用 conda：`conda install -c conda-forge numpy scipy matplotlib ase seekpath`。

通用约定：

- **作图**：matplotlib 只在需要作图时才导入，默认用 Agg 后端直接保存图片（集群上没有图形界面也能用），加 `--show` 才弹出窗口。`-o/--output none` 表示不出图，只打印结果或写 CSV。
- **单位**：pw.x 能量为 Ry、长度为 bohr；bands.x、dos.x、projwfc.x 输出的能量为 eV；matdyn.x 频率为 cm⁻¹。脚本内部用 CODATA 2018 常数：1 Ry = 13.605693122994 eV，1 bohr = 0.529177210903 Å。各节写明了输入输出单位。
- **费米能**：画能带/DOS 时，参考能量应取自 scf（或 nscf）的输出，**不要**用 `calculation='bands'` 运行里的费米能。各绘图脚本都提供 `--efermi` 和 `--scf-out` 两种给法。
- **换行符**：在 Windows 下编辑的输入需要先 `dos2unix`。[`qe_check_input.py`](qe_check_input.py) 和作业模板都会检查 CRLF。

> [!IMPORTANT]
> 脚本生成的输入（[`qe_inputgen.py`](qe_inputgen.py)、[`qe_conv_test.sh`](qe_conv_test.sh)）在提交前都要人工读一遍，并用 [`qe_check_input.py`](qe_check_input.py) 检查。赝势文件不随仓库分发，下载方法见 [第3章](../tutorials/03-pseudopotentials.md)。

---

## 1. qe_inputgen.py — 从结构文件生成 pw.x 输入

**用途**：读入 ASE 能识别的结构（cif、POSCAR、带晶胞的 extxyz 等），写出 `ibrav = 0` + `CELL_PARAMETERS angstrom` + `ATOMIC_POSITIONS crystal` 的完整 pw.x 输入。给定 SSSP JSON（如 `SSSP_1.3.0_PBE_efficiency.json`）时，按其中每个元素的 `filename`、`cutoff_wfc`、`cutoff_rho` 填写 `ATOMIC_SPECIES`，并取结构中所有元素截断的**最大值**作为 `ecutwfc`/`ecutrho`（这是 SSSP 对化合物的推荐做法）。

**依赖**：`pip install ase`；`--kpath seekpath` 另需 `pip install seekpath`。两者都在运行时才导入，缺失时会提示安装命令。

**用法示例**：

```bash
python qe_inputgen.py Si.cif --sssp SSSP_1.3.0_PBE_efficiency.json --kspacing 0.2 -o si.scf.in
python qe_inputgen.py POSCAR --format vasp --sssp sssp.json --calculation vc-relax --smearing mv --degauss 0.02
python qe_inputgen.py Fe.cif --sssp sssp.json --magmom Fe=0.5 --smearing mv --degauss 0.02
python qe_inputgen.py Si.cif --sssp sssp.json --calculation scf --kpath-cell -o si.scf.in      # SeeK-path cell
python qe_inputgen.py Si.cif --sssp sssp.json --calculation bands --kpath seekpath --npoints 40 --nbnd 16
python qe_inputgen.py MoS2.vasp --format vasp --ecutwfc 50 --dual 8 --kgrid 12 12 1 --pseudo-dir ./pseudo
python qe_check_input.py si.scf.in          # always check the result
```

**参数**（括号内为默认值）：

| 类别 | 参数 | 说明 |
|---|---|---|
| 通用 | `structure` | 结构文件（位置参数），必须带晶胞 |
| | `--format FMT` | ASE 格式名，无法从文件名猜出时指定（如 `vasp`） |
| | `-o/--output` | 输出文件（`<prefix>.<calculation>.in`） |
| | `--force` | 覆盖已存在的文件 |
| | `--calculation` | `scf`（默认）/ `nscf` / `bands` / `relax` / `vc-relax` |
| | `--prefix`、`--outdir`、`--pseudo-dir` | 写入 `&CONTROL`（`pwscf`、`./tmp`、`./pseudo`） |
| 赝势与截断 | `--sssp JSON` | SSSP JSON：赝势文件名 + 推荐截断 |
| | `--pseudo El=FILE ...` | 覆盖个别元素的赝势文件名；不给 `--sssp` 时默认写 `<El>.UPF`，必须手动改 |
| | `--ecutwfc`、`--ecutrho` | 覆盖截断（Ry） |
| | `--dual` | `ecutrho = dual × ecutwfc`；未给 `--sssp`/`--ecutrho` 时默认 4（模守恒），超软/PAW 通常 8–12 |
| k 点 | `--kspacing` | k 间距（Å⁻¹，含 2π，与 VASP KSPACING 约定相同；默认 0.25），`n_i = max(1, ceil(2π|b_i|/kspacing))` |
| | `--kgrid N1 N2 N3` | 直接给网格（与 `--kspacing` 互斥） |
| | `--kshift S1 S2 S3` | 网格偏移（`0 0 0`） |
| | `--kpath` | bands 路径来源：`seekpath`（默认）或 `ase` |
| | `--kpath-cell` | scf/nscf 时也写 SeeK-path 标准原胞，保证与后续 bands 运行同一晶胞 |
| | `--npoints` | crystal_b 每段点数（40） |
| | `--nbnd` | 能带数（bands/nscf 时需要包含空带） |
| 占据与自旋 | `--smearing` | `gaussian` / `mv` / `m-p` / `fd`，写 `occupations='smearing'` |
| | `--degauss` | 展宽（Ry，0.02） |
| | `--tetra` | nscf 时用 `occupations='tetrahedra_opt'`（算 DOS） |
| | `--magmom El=VAL ...` | 打开 `nspin = 2` 并写 `starting_magnetization(i)`（取值 −1 到 1） |
| | `--vdw` | `vdw_corr`：`grimme-d2` / `grimme-d3` / `dft-d3` / `ts` / `xdm` / `mbd` |
| 收敛 | `--conv-thr-per-atom` | `conv_thr = 值 × nat`（1e-9 Ry） |
| | `--mixing-beta` | 0.4 |
| | `--forc-conv-thr`、`--etot-conv-thr` | relax 判据（1e-4 Ry/bohr、1e-5 Ry） |
| | `--press-conv-thr`、`--cell-dofree` | vc-relax 的 `&CELL`（0.5 kbar；如 `ibrav`、`shape`、`2Dxy`） |

> [!WARNING]
> `--calculation bands --kpath seekpath` 写出的是 SeeK-path **标准化原胞**，它可能与原始结构的取向或原胞选择不同。crystal_b 坐标只对这个晶胞有效，因此 scf 也要用同一晶胞：先用 `--kpath-cell` 生成 scf 输入，再用相同参数生成 bands 输入。

## 2. qe_check_input.py — pw.x 输入静态检查

**用途**：在提交作业前发现常见低级错误，不需要 pw.x。每条信息分为 `ERROR`（pw.x 会停止）、`WARN`（多半有问题或已过时）、`INFO`（提示）三级。

主要检查项：

- namelist（`&CONTROL`、`&SYSTEM`、`&ELECTRONS` ……）没有以 `/` 结束，或卡片写进了 namelist；
- QE 7.1 起删除的旧 DFT+U 关键字：`lda_plus_u`、`lda_plus_u_kind`、`Hubbard_U(i)`、`Hubbard_J0`、`Hubbard_V`、`Hubbard_J`、`U_projection_type`（应改用 `HUBBARD` 卡，详见 [第8章](../tutorials/08-dftu-hp.md)）；
- 已废弃、会被忽略的 `wf_collect`、`lkpoint_dir`；
- `CELL_PARAMETERS` / `ATOMIC_POSITIONS` 没写单位（已不推荐）；`ibrav = 0` 却没有 `CELL_PARAMETERS`；同时给了 `celldm(1)` 和 `A`；
- `ecutrho/ecutwfc < 4`；赝势为超软/PAW（读 UPF 头部的 `pseudo_type`）但 `ecutrho = 4 × ecutwfc`；
- `ntyp`、`nat` 与 `ATOMIC_SPECIES`、`ATOMIC_POSITIONS` 行数不一致，原子标签不在 `ATOMIC_SPECIES` 中；
- `occupations='smearing'` 没给 `degauss`（默认 0）；`nspin=2` 没给 `starting_magnetization`/`tot_magnetization`；
- `calculation='bands'` 却用 `K_POINTS automatic`；`K_POINTS automatic` 不是 6 个整数；
- CRLF 换行、非 ASCII 字符（全角引号/空格）；赝势文件是否存在（目录依次取 `--pseudo-dir`、输入里的 `pseudo_dir`、`$ESPRESSO_PSEUDO`、`~/espresso/pseudo`）。

**依赖**：仅标准库。

**用法示例**：

```bash
python qe_check_input.py si.scf.in
python qe_check_input.py *.in --pseudo-dir ./pseudo
python qe_check_input.py fe.scf.in --strict       # WARN also gives exit status 1
```

| 参数 | 说明 |
|---|---|
| `inputs` | 一个或多个输入文件 |
| `--pseudo-dir` | 覆盖输入中的 `pseudo_dir` |
| `--strict` | 出现 WARN 也以状态 1 退出（适合放进批量脚本） |

退出状态：无 ERROR（`--strict` 时无 WARN）为 0，否则为 1。

## 3. qe_conv_test.sh — 收敛测试批量生成与汇总

**用途**：对模板输入中的占位符做 `sed` 替换，为每个取值生成 `DIR/<值>/pw.in`。默认**只生成不运行**；加 `--run` 才在各子目录运行 `$PW_CMD -in pw.in > pw.out`。最后汇总到 `DIR/summary.dat`。原理和判据见 [第4章](../tutorials/04-convergence.md)，该章的示例模板与本脚本的占位符相同。

占位符：

| 占位符 | 模式 | 替换为 |
|---|---|---|
| `@ECUT@` | `ecut` | ecutwfc（Ry） |
| `@ECUTRHO@` | `ecut` | `dual × ecutwfc`，dual 由 `-d` 给出 |
| `@K@` | `kpt` | `N N N` |
| `@K1@ @K2@ @K3@` | `kpt` | `N N N3`，N3 由 `-z` 固定（板/二维体系） |
| `@DEGAUSS@` | `smear` | 展宽（Ry） |
| `@KEY@`（任意） | 任意 | 由 `-s KEY=VALUE` 固定，例如 `-s "K=8 8 8"`、`-s SMEAR=mv` |

生成后如果还有未替换的 `@...@`，脚本会报错退出。

**用法示例**：

```bash
bash qe_conv_test.sh -t si_ecut_template.in -m ecut -v "25 30 35 40 50 60" -d 8 -s "K=8 8 8"
bash qe_conv_test.sh -t al_kpt_template.in -m kpt -v "8 12 16 20" -s DEGAUSS=0.02 -s SMEAR=mv
bash qe_conv_test.sh -t slab_template.in -m kpt -v "4 6 8" -z 1
PW_CMD="mpirun -np 4 pw.x -nk 2" bash qe_conv_test.sh -t si_ecut_template.in -m ecut -v "30 40 50" -d 8 -s "K=8 8 8" --run
bash qe_conv_test.sh --collect conv_ecut
```

| 选项 | 说明 |
|---|---|
| `-t FILE` | 模板输入（`--collect` 时不需要） |
| `-m MODE` | `ecut` / `kpt` / `smear` |
| `-v "LIST"` | 空格分隔的取值（加引号） |
| `-d DUAL` | `@ECUTRHO@` 的倍数（默认 4；超软/PAW 用 8–12） |
| `-z N3` | 第三个 k 方向固定为 N3 |
| `-s KEY=VAL` | 其它占位符的固定值，可重复 |
| `-o DIR` | 工作目录（默认 `conv_<mode>`） |
| `--run` | 运行 pw.x；命令取自环境变量 `PW_CMD`（默认 `pw.x`） |
| `--collect DIR` | 只汇总已有目录 |

`summary.dat` 各列：取值、最后一个 `!    total energy` 行的能量（Ry）、相对**最后一个（最收敛）**取值的能量差（meV/atom，原子数取自 `number of atoms/cell`）、压强（`P=`，kbar，需 `tstress=.true.`）、SCF 迭代总数、状态（`done` / `unfinished` / `scf_not_converged` / `not_run`）。

## 4. qe_parse_out.py — pw.x 输出解析

**用途**：一次性读出 pw.x 输出中最常用的信息，适用于 scf、nscf、relax、vc-relax、md。可批量处理多个文件并导出 CSV（每个文件一行）或 JSON（全部解析结果）。

解析的行（QE 7.x `PW/src` 中的格式）：

| 内容 | 输出中的行 |
|---|---|
| 收敛总能 | `!    total energy              =     ... Ry`（每个离子步一行） |
| SCF 历史 | `iteration #`、`total energy =`、`estimated scf accuracy    <` |
| SCF 是否收敛 | `convergence has been achieved in` / `convergence NOT achieved after` |
| 费米能 / 带边 | `the Fermi energy is`、`the spin up/dw Fermi energies are`、`highest occupied, lowest unoccupied level (ev):`、`highest occupied level (ev):` |
| 力 | `Forces acting on atoms (cartesian axes, Ry/au):` 下的 `atom ... force =`，以及 `Total force =` |
| 应力 | `total   stress  (Ry/bohr**3) ... P=`（取 kbar 部分） |
| 磁矩 | `total magnetization`、`absolute magnetization`、`Magnetic moment per site` 下的 `charge= ... magn=` |
| 优化 | `bfgs converged in`、`bfgs failed after`、`The maximum number of steps has been reached.`、`Final energy =` / `Final enthalpy =` |
| 最终结构 | `Begin final coordinates` … `End final coordinates` |
| 结束 / 报错 | `JOB DONE.`；`%%%%%%` 包围的错误信息 |

**依赖**：仅标准库。

**用法示例**：

```bash
python qe_parse_out.py si.scf.out
python qe_parse_out.py relax.out --scf                          # SCF history of every ionic step
python qe_parse_out.py relax.out --final-structure final.in     # CELL_PARAMETERS + ATOMIC_POSITIONS
python qe_parse_out.py conv_ecut/*/pw.out --csv summary.csv -q
python qe_parse_out.py fe.scf.out --json fe.json
```

| 参数 | 说明 |
|---|---|
| `outputs` | 一个或多个 pw.x 输出 |
| `--scf` | 打印 SCF 历史 |
| `--json FILE` | 写出全部解析数据 |
| `--csv FILE` | 每个文件一行：能量（Ry、eV、eV/atom）、费米能、HOMO/LUMO、压强、合力、磁矩、离子步数、未收敛 SCF 数 |
| `--final-structure FILE` | 把（最后一个文件的）最终结构块写入 FILE，可粘贴进新输入 |
| `-q/--quiet` | 不打印文本报告 |

退出状态：所有文件都有 `JOB DONE.` 且没有未收敛的 SCF 时为 0，否则为 1。HOMO/LUMO 只针对**当前 k 网格**，不等于能带图上的带隙；绝缘体要输出 LUMO，需要设置 `nbnd` 包含空带。

## 5. qe_plot_bands.py — 能带作图与带隙

**用途**：读取 bands.x 写出的 `<filband>.gnu`（两列：路径坐标（2π/alat）与能量（eV），每条带一块，块之间空行分隔），或 `calculation='bands'` 运行后的 `<outdir>/<prefix>.save/data-file-schema.xml`（`<ks_energies>` 中的本征值为 Hartree，脚本换算成 eV）。高对称点位置取自 bands.x 标准输出中的 `high-symmetry point: ... x coordinate` 行（`--bandsx-out`），或手动给 `--ticks`；读 XML 时由 k 路径方向的变化自动判断。流程见 [第6章](../tutorials/06-band-dos-pdos.md)。

**依赖**：numpy、matplotlib。

**用法示例**：

```bash
python qe_plot_bands.py si.bands.dat.gnu --bandsx-out si.bands_pp.out --scf-out si.scf.out \
       --labels G X W K G L --ylim -12 6 -o si_bands.png
python qe_plot_bands.py si.bands.dat.gnu --efermi 6.2 --labels "G X|U K G L" --csv bands.csv
python qe_plot_bands.py tmp/si.save/data-file-schema.xml --labels G X W K G L --scf-out si.scf.out
python qe_plot_bands.py fe_up.gnu --spin2 fe_dw.gnu --scf-out fe.scf.out        # bands.x spin_component=1/2
```

| 参数 | 说明 |
|---|---|
| `data` | `.gnu` 文件或 `data-file-schema.xml` |
| `--spin2 FILE` | 自旋向下的 `.gnu`（bands.x `spin_component=2`）；XML 中 `lsda=true` 时自动分开 |
| `--bandsx-out FILE` | bands.x 的标准输出，读取高对称点 x 坐标 |
| `--ticks X ...` | 手动给刻度位置 |
| `--labels L ...` | 每个刻度的标签；`G` 显示为 Γ；路径跳跃处写成 `X\|U`；整串加引号也可以 |
| `--efermi EV` / `--scf-out FILE` | 参考能量；否则用 XML 中的值（会给出提示），都没有则不平移 |
| `--no-shift` | 画绝对能量 |
| `--tol` | 判定占据/空态的容差（eV，默认 1e-3） |
| `--ylim E1 E2` | 能量范围（−6 6） |
| `--title`、`--csv FILE` | 标题；导出平移后的能量 |
| `-o/--output`、`--dpi`、`--show` | 图片（`bands.png`，`none` 不出图）、分辨率（200）、交互窗口 |

脚本会打印 VBM、CBM 所在的 k 点序号与路径坐标，判断直接/间接带隙，并给出最小直接带隙；如果有能带穿过参考能量，则判定为金属。带隙只针对所画路径，且 LDA/GGA 会低估带隙（见 [第14章](../tutorials/14-hybrid-metagga.md)）。

## 6. qe_plot_dos.py — DOS 与 PDOS

**用途**：读取 dos.x 输出（`<fildos>`）和 projwfc.x 输出（`<filpdos>.pdos_tot`、`<filpdos>.pdos_atm#N(El)_wfc#M(l)`），按元素、元素+角动量或原子求和后作图。自旋极化时自旋向下画成负值。格式依据 `PP/src/dos.f90` 与 `PP/src/partialdos.f90`：

- dos.x：表头 `#  E (eV)   dos(E)     Int dos(E) EFermi = ... eV`，nspin=2 时列为 `E dosup dosdw Int`。脚本可以直接从表头读取费米能；
- `pdos_atm#…`：第 2 列（nspin=2 时为第 2、3 列）是对 m 求和的 `ldos`，脚本对这一列求和；自旋轨道耦合时文件名为 `..._wfc#2(p_j1.5)`，同一 l 的不同 j 会被加在一起；
- `kresolveddos=.true.` 的 k 分辨文件（多一列 `ik`）不支持。

**依赖**：numpy、matplotlib。

**用法示例**：

```bash
python qe_plot_dos.py --dos si.dos --xlim -12 6 -o si_dos.png
python qe_plot_dos.py --pdos-prefix si.pdos --by element --scf-out si.nscf.out
python qe_plot_dos.py --pdos-prefix fe.pdos --by orbital --efermi 17.8 --xlim -8 4
python qe_plot_dos.py --dos tio2.dos --pdos-prefix tio2.pdos --by orbital --select Ti-d O-p --csv pdos.csv
python qe_plot_dos.py --pdos-prefix nio.pdos --by atom --select 1 2 --pdos-tot
```

| 参数 | 说明 |
|---|---|
| `--dos FILE` | dos.x 输出 |
| `--pdos-prefix P` | projwfc.x 的 `filpdos`（读取 `P.pdos_atm#*`） |
| `--by` | `element`（默认，标签如 `Fe`）/ `orbital`（`Fe-d`）/ `atom`（`1-Fe`） |
| `--select ...` | 只画指定分组；`--by atom` 时给原子序号 |
| `--pdos-tot` | 同时画 `P.pdos_tot` 中的投影总和（与总 DOS 的差反映投影的完整程度） |
| `--efermi` / `--scf-out` | 费米能；否则用 dos.x 表头的 `EFermi` |
| `--no-shift`、`--xlim`（−10 5）、`--ymax` | 能量轴与 DOS 上限 |
| `--csv FILE` | 导出所有曲线（能量网格不同的曲线会被跳过） |
| `-o/--output`、`--dpi`、`--show` | `dos.png`；`none` 不出图 |

脚本还会打印每条曲线在 `--xlim` 窗口内的积分（电子数的粗略估计），方便检查 PDOS 是否合理。

## 7. qe_phonon_plot.py — 声子色散与 DOS

**用途**：读取 matdyn.x 的 `<flfrq>.gp`（每行：累计路径长度（2π/alat）+ 3·nat 个频率（cm⁻¹））和 `<fldos>`（`# Frequency[cm^-1] DOS PDOS`），画色散和 DOS。matdyn.x 把 ω² < 0 的模式写成**负频率**，脚本会报告低于 `−imag_tol` 的虚频及其所在的 q 点序号。格式依据 `PHonon/PH/matdyn.f90`。完整流程见 [第12章](../tutorials/12-phonon-dfpt.md)。

刻度位置：

- `--matdyn-in`（推荐）：读取 `q_in_band_form=.true.` 时 namelist 后面的 q 点列表（`qx qy qz n`），第 i 个高对称点位于前面各段点数之和处。行尾的 `! X` 注释会作为标签；
- `--freq-file`：读取 `<flfrq>` 中的 q 坐标，由方向变化判断；
- 都不给时只标出路径两端。

`.gp` 文件在路径跳跃处不会断开，跳跃表现为两个刻度之间的一小段。

**依赖**：numpy、matplotlib。

**用法示例**：

```bash
python qe_phonon_plot.py --gp si.freq.gp --matdyn-in matdyn.in --dos si.phdos -o si_phonon.png
python qe_phonon_plot.py --gp si.freq.gp --freq-file si.freq --labels G X U K G L
python qe_phonon_plot.py --gp mos2.freq.gp --matdyn-in matdyn.in --unit THz --imag-tol 5
python qe_phonon_plot.py --dos si.phdos --unit meV --pdos
```

| 参数 | 说明 |
|---|---|
| `--gp FILE` | `<flfrq>.gp` |
| `--dos FILE` | matdyn.x `dos=.true.` 写出的 `fldos`（默认文件名 `matdyn.dos`） |
| `--matdyn-in FILE` / `--freq-file FILE` | 确定高对称点位置（见上） |
| `--labels ...` | 刻度标签，`G` 显示为 Γ |
| `--unit` | `cm-1`（默认）/ `THz` / `meV` |
| `--imag-tol` | 虚频阈值（cm⁻¹，默认 1.0）；Γ 附近的小负值多为声学求和规则或收敛问题 |
| `--ylim`、`--pdos`、`--csv FILE` | 纵轴范围；画各原子 PDOS 列；导出换算后的色散 |
| `-o/--output`、`--dpi`、`--show` | `phonon.png` |

## 8. qe_eos_fit.py — 状态方程拟合

**用途**：对一组 E(V) 或 E(a) 拟合三阶 Birch-Murnaghan（默认）或 Murnaghan 方程，给出 E0、V0（或 a0）、B0（GPa）、B0' 及参数误差，并检查最小值是否落在数据范围内部。QE 自带的 `ev.x` 用交互方式做同样的拟合，用法见 [第5章](../tutorials/05-relax-eos.md)。

输入可以是两列数据文件，也可以用 `--from-out` 直接读一组 pw.x 输出：脚本取每个文件的 `unit-cell volume = ... (a.u.)^3`（vc-relax 时取 `new unit-cell volume`）和最后一个 `!` 总能。

**依赖**：numpy、scipy、matplotlib。

**用法示例**：

```bash
python qe_eos_fit.py ev.dat                                         # V (A^3)  E (Ry)
python qe_eos_fit.py ea.dat --x alat --vfactor 0.25 -o si_eos.png   # fcc / diamond primitive cell
python qe_eos_fit.py ea.dat --x alat --xunit bohr --vfactor 0.25    # a in bohr (celldm(1))
python qe_eos_fit.py --from-out scan/*/pw.out --eos murnaghan --natoms 2
```

| 参数 | 说明 |
|---|---|
| `data` | 两列数据文件（`#` 为注释） |
| `--from-out OUT ...` | 改为读取 pw.x 输出 |
| `--x` | 第一列是 `volume`（默认）还是 `alat` |
| `--xunit` | `angstrom`（默认）或 `bohr` |
| `--vfactor` | `V = vfactor × a³`：简单立方 1，bcc 原胞 1/2，fcc/金刚石原胞 1/4 |
| `--eunit` | `Ry`（默认）/ `eV` / `Ha` |
| `--natoms` | 按原子数归一，输出每原子的 E、V（B0 不变） |
| `--eos` | `bm3`（默认）或 `murnaghan` |
| `-o/--output`、`--dpi`、`--show` | `eos.png` |

至少需要 5 个点，建议在平衡晶格常数附近 ±5% 范围内取 7–9 个点。

## 9. qe_workfunction.py — 平面平均势与功函数

**用途**：从 average.x 的 `avg.dat`，或 pp.x 以 `output_format=6` 写出的 cube 文件，得到垂直于表面的平面平均势，自动寻找真空平台（或用 `--window` 指定），计算 Φ = V_vac − E_F。流程（见 [第9章](../tutorials/09-vdw-surface-workfunction.md)）：pw.x scf → pp.x `plot_num=11`（V_bare + V_H）→ average.x → 本脚本。

- `avg.dat`（`PP/src/average.f90`）：三列 `z p(z) m(z)`，z 单位为 bohr，势的单位与 pp.x 量一致（plot_num=11 为 Ry）。脚本换算为 Å 与 eV；
- cube：长度单位为 bohr，数值单位为 Ry，沿 `--axis` 做平面平均；
- 自动判断真空：在 |dV/dz| < `--flat-tol`（eV/Å）的最宽区域中取中间 20% 的平均值。非对称板加偶极校正时，两侧真空能级不同，脚本会分别报告。

**依赖**：numpy、matplotlib。

**用法示例**：

```bash
python qe_workfunction.py avg.dat --scf-out pt111.scf.out -o pt111_wf.png
python qe_workfunction.py avg.dat --efermi 5.123 --window 22 26
python qe_workfunction.py slab_pot.cube --axis 3 --scf-out slab.scf.out --csv vz.csv
python qe_workfunction.py avg.dat --scf-out slab.scf.out --column macro
```

| 参数 | 说明 |
|---|---|
| `data` | `avg.dat` 或 `.cube` |
| `--axis` | cube 的表面法向（1/2/3，默认 3） |
| `--column` | avg.dat 用平面平均 `planar`（第 2 列，默认）或宏观平均 `macro`（第 3 列） |
| `--efermi` / `--scf-out` | 费米能（取自板模型 scf 输出） |
| `--window Z1 Z2` | 真空窗口（Å） |
| `--flat-tol` | 自动判断平台的斜率阈值（eV/Å，默认 0.05） |
| `--csv FILE`、`-o/--output`、`--dpi`、`--show` | 导出 z(Å)、V(eV)；图片 `workfunction.png` |

## 10. qe_cube_diff.py — 差分电荷与平面平均

**用途**：处理 pp.x 写出的 Gaussian cube（`plot_num=0`、`output_format=6`）。三个子命令：

- `info`：网格、晶胞、数值范围与积分（Σ值 × 体素体积）。对电荷密度 cube，积分应等于价电子数，建议先用它核对单位；
- `diff`：`OUT = Σ cᵢ · cubeᵢ`，默认系数为 `1 −1 −1 …`，即 Δρ = ρ(AB) − ρ(A) − ρ(B)。所有输入必须使用同一晶胞和同一 FFT 网格（相同晶胞与 `ecutrho`，必要时显式设定 `nr1/nr2/nr3`），否则脚本报错；
- `planar`：沿某一晶轴求平面平均、单位长度电荷 λ(z) 和累积电荷 ΔQ(z) = ∫₀ᶻ λ dz'，可作双轴图并导出 CSV。

差分电荷的完整流程、可视化（VESTA）和 Bader 分析见 [第10章](../tutorials/10-charge-elf-bader.md)。

**依赖**：numpy；`planar` 作图需要 matplotlib。

**用法示例**：

```bash
python qe_cube_diff.py info ab_rho.cube
python qe_cube_diff.py diff drho.cube ab_rho.cube a_rho.cube b_rho.cube
python qe_cube_diff.py diff spin.cube up.cube dn.cube --coeffs 1 -1
python qe_cube_diff.py planar drho.cube --axis 3 --csv drho_z.csv -o drho_z.png
```

| 子命令 | 参数 | 说明 |
|---|---|---|
| `info` | `cubes ...` | 一个或多个 cube |
| `diff` | `output inputs ...` | 输出文件在前，输入在后 |
| | `--coeffs c1 c2 ...` | 线性组合系数（默认 `1 -1 -1 ...`） |
| `planar` | `cube` | 输入 cube |
| | `--axis` | 1/2/3（默认 3） |
| | `--csv FILE` | 列：`z_A, rho_avg_e_per_A3, lambda_e_per_A, DeltaQ_e` |
| | `-o/--output`、`--dpi`、`--show` | 图片 `planar.png`（`none` 不出图） |

## 11. job_qe_slurm.sh / job_qe_pbs.sh — 作业模板

**用途**：通用的 QE 作业脚本，适用于 pw.x、ph.x、pp.x、neb.x、hp.x 等。主要步骤：

1. 加载环境：在脚本的 environment 块中**只保留一种**方式（module、conda 或自编译路径），按需设置 `ESPRESSO_PSEUDO`；
2. 设置 `ESPRESSO_TMPDIR`：输入里没写 `outdir` 时，QE 把临时文件放在这里。Slurm 版默认 `$SLURM_TMPDIR`，PBS 版默认 `$TMPDIR`，都没有时用提交目录下的 `tmp`；
3. 检查输入是否存在、是否为 CRLF 换行，可执行文件是否在 PATH 中；
4. 设置 `OMP_NUM_THREADS`；已有的输出备份为 `<out>.<n>.bak`；
5. 运行 `srun EXE <flags> -in INP > OUT`（PBS 版为 `mpirun -np N -machinefile $PBS_NODEFILE ...`）；
6. 检查 `JOB DONE.`、`convergence NOT achieved` 和 `%%%%%%` 错误块，打印最终能量、费米能、BFGS 状态，并提示用 `max_seconds` + `restart_mode='restart'` 续算。

**依赖**：bash；Slurm 或 PBS Pro / OpenPBS / Torque。

**用法示例**：

```bash
# Slurm
sbatch job_qe_slurm.sh si.scf.in
sbatch --export=ALL,INP=si.scf.in,NK=4 job_qe_slurm.sh
sbatch --export=ALL,EXE=ph.x,INP=si.ph.in,NK=4 job_qe_slurm.sh
sbatch --export=ALL,EXE=neb.x,INP=neb.in,NI=4 job_qe_slurm.sh
# PBS (script arguments are not passed, use -v)
qsub -v INP=si.scf.in,NK=4 job_qe_pbs.sh
qsub -v EXE=ph.x,INP=si.ph.in,NK=4 job_qe_pbs.sh
bash job_qe_slurm.sh -h
```

| 变量 | 说明 |
|---|---|
| `INP` / `OUT` | 输入（Slurm 版也可作为第一个参数；默认 `pw.in`）/ 输出（默认 `<INP 去扩展名>.out`） |
| `EXE` | 可执行程序（默认 `pw.x`） |
| `NK` | `-nk`，k 点池数：应整除 MPI 进程数，且不超过不可约 k 点数 |
| `ND` | `-nd`，并行对角化（ScaLAPACK）进程数，须为平方数 |
| `NT` | `-nt`，task groups |
| `NI` | `-ni`，images（neb.x、ph.x） |
| `EXTRA` | 其它命令行参数 |
| `LAUNCHER` | Slurm 版默认 `srun`，可改为 `mpirun`；PBS 版默认 `mpirun` |
| `NTHREADS`、`NRANKS` | 仅 PBS 版：OpenMP 线程数（默认取 `OMP_NUM_THREADS`，否则 1）；MPI 进程数（默认为 `$PBS_NODEFILE` 行数 ÷ NTHREADS） |

不指定 `NK`/`ND` 时，QE 7.1 及以后会自动选择并行参数。并行参数的选取、GPU（每张 GPU 一个 MPI 进程）和性能分析见 [第18章](../tutorials/18-automation-parallel-troubleshooting.md)，安装与环境配置见 [第1章](../tutorials/01-install.md)。

---

## 参考资料

- QE 输入文档：[INPUT_PW](https://www.quantum-espresso.org/Doc/INPUT_PW.html)、[INPUT_PP](https://www.quantum-espresso.org/Doc/INPUT_PP.html)、[INPUT_PROJWFC](https://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html)、[INPUT_DOS](https://www.quantum-espresso.org/Doc/INPUT_DOS.html)、[INPUT_BANDS](https://www.quantum-espresso.org/Doc/INPUT_BANDS.html)
- QE 源码（输出格式）：[gitlab.com/QEF/q-e](https://gitlab.com/QEF/q-e)：`PW/src/electrons.f90`、`print_ks_energies.f90`、`forces.f90`、`stress.f90`、`output_tau.f90`，`PP/src/bands.f90`、`dos.f90`、`partialdos.f90`、`average.f90`，`PHonon/PH/matdyn.f90`
- [SSSP](https://www.materialscloud.org/discover/sssp/table/efficiency)（赝势与推荐截断 JSON）、[SeeK-path](https://github.com/materialscloud-org/seekpath)（k 路径）
- 相关工具：[ASE](https://wiki.fysik.dtu.dk/ase/)（`ase.io.espresso` 读写 QE 输入输出）、[QEPlotter](https://github.com/shubics/QEPlotter)（能带/PDOS 作图）、[aiida-quantumespresso](https://github.com/aiidateam/aiida-quantumespresso)（大规模工作流）

[返回目录](../README.md)
