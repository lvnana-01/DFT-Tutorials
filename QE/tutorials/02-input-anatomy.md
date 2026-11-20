# 第2章 pw.x 输入文件结构与第一个 SCF（pw.x Input Anatomy & Your First SCF Calculation）

本章逐段讲解 pw.x 输入文件的两类组成部分：namelist（`&CONTROL`、`&SYSTEM`、`&ELECTRONS`、`&IONS`、`&CELL`）和 card（`ATOMIC_SPECIES`、`ATOMIC_POSITIONS`、`K_POINTS`、`CELL_PARAMETERS`、`HUBBARD`）。重点讲清 `ibrav` 与 `celldm(1)`（Bohr）/ `A`（Å）的关系、为什么推荐 `ibrav=0` + 显式单位、Ry/Bohr 单位制，以及如何从输出中 grep 总能量、收敛信息、HOMO/LUMO、费米能、力和应力。三个示例覆盖半导体 Si、金属 Al（smearing）和放在大盒子里的 H2O 分子（`K_POINTS gamma` + `assume_isolated='mt'`）。

> [!NOTE]
> 本章全部输入在示例目录 [02-input-anatomy](../examples/02-input-anatomy/README.md)。赝势使用 SSSP 1.3 PBE Efficiency（下载方法见 [第3章](03-pseudopotentials.md)），仓库不提交 UPF 文件。文中输出均为**示意输出，数值仅供参考**，不是真实计算结果。

---

## 1. 本章目标与适用场景

| 目标 | 你将学会 |
|---|---|
| 语法 | namelist 与 card 的区别、`/` 结束符、注释、字符串引号、大小写规则 |
| 结构 | `ibrav` + `celldm`/`A`、`ibrav=0` + `CELL_PARAMETERS`；四种原子坐标单位 |
| k 点 | `automatic`、`gamma`、显式列表、能带路径 `crystal_b`/`tpiba_b` |
| 运行 | `pw.x -i`、`mpirun`、`-nk` |
| 输出 | `!` 总能量行、SCF 历史、HOMO/LUMO 与费米能、力、应力、`JOB DONE`、`prefix.save/` |
| 单位 | Ry ↔ eV、Bohr ↔ Å、Ry/Bohr、kbar |
| 分子 | 大盒子 + `K_POINTS gamma` + `assume_isolated='mt'` |

从 VASP 转过来的读者，下表给出概念上的对应关系：

| VASP | pw.x | 说明 |
|---|---|---|
| INCAR | `&CONTROL`、`&SYSTEM`、`&ELECTRONS`、`&IONS`、`&CELL` | QE 把所有参数写在一个文件里，按 namelist 分组 |
| POSCAR 晶格 | `ibrav` + `celldm`/`A`，或 `ibrav=0` + `CELL_PARAMETERS` | |
| POSCAR 坐标（Direct/Cartesian） | `ATOMIC_POSITIONS crystal` / `angstrom` | |
| Selective dynamics `F F F` | 坐标后加 `0 0 0`（`if_pos`） | 见 [第5章](05-relax-eos.md) |
| KPOINTS | `K_POINTS` card | |
| POTCAR | `ATOMIC_SPECIES` 中的 UPF 文件名 + `pseudo_dir` | |
| `ENCUT`（eV） | `ecutwfc`（Ry）+ `ecutrho`（Ry） | 数值不能直接换算后套用，见 [第4章](04-convergence.md) |
| `ISMEAR`/`SIGMA`（eV） | `occupations`、`smearing`、`degauss`（Ry） | |
| `EDIFF`（eV） | `conv_thr`（Ry，广延量） | |
| OUTCAR 中的 `TOTEN` | 以 `!` 开头的 `total energy` 行 | 金属时是自由能 F = E − TS |
| WAVECAR/CHGCAR | `outdir/prefix.save/` | |

更多 VASP 文件的细节见 [VASP 第2章](../../VASP/tutorials/02-input-output-files.md)。

---

## 2. 原理速览：pw.x 需要知道什么

一次 SCF 计算，pw.x 需要知道五类信息：

1. **做什么**：`calculation`（`'scf'`、`'nscf'`、`'bands'`、`'relax'`、`'md'`、`'vc-relax'`、`'vc-md'`），以及文件放在哪里（`prefix`、`outdir`、`pseudo_dir`）；
2. **体系是什么**：晶格、原子种类与位置、电子数（由赝势的价电子数决定）、总电荷、自旋；
3. **离散化精度**：平面波截断 `ecutwfc`/`ecutrho`、k 点网格；
4. **如何占据电子**：`occupations`（绝缘体 `'fixed'`，金属 `'smearing'`，DOS 用 `'tetrahedra_opt'`）；
5. **如何收敛**：`conv_thr`、`mixing_beta`、`electron_maxstep`。

QE 的单位制：

| 物理量 | QE 单位 | 换算 |
|---|---|---|
| 能量 | Ry | 1 Ry = 13.605693 eV；1 Ha = 2 Ry |
| 长度 | Bohr（a.u.） | 1 Bohr = 0.529177 Å |
| 力 | Ry/Bohr | 1 Ry/Bohr ≈ 25.711 eV/Å |
| 应力/压强 | 输出同时给 Ry/Bohr³ 和 kbar | 1 kbar = 0.1 GPa |
| 时间（pw.x MD） | Rydberg 原子单位 | 1 a.u. = 4.8378×10⁻¹⁷ s（cp.x 用 Hartree 单位，是它的一半） |
| `degauss` | Ry | 0.01 Ry ≈ 0.136 eV |

> [!IMPORTANT]
> `celldm(1)` 的单位是 **Bohr**，`A` 的单位是 **Å**。最常见的错误是把 Å 为单位的晶格常数（例如 5.43）写进 `celldm(1)`，结果晶胞小了将近一半，原子挤在一起，SCF 报 `S matrix not positive definite` 或能量异常。

---

## 3. 输入文件的整体结构与语法

```text
&CONTROL          <- namelist: keyword = value pairs, ends with '/'
  ...
/
&SYSTEM
  ...
/
&ELECTRONS
  ...
/
&IONS             <- required only for relax/md/vc-relax/vc-md
/
&CELL             <- required only for vc-relax/vc-md
/
ATOMIC_SPECIES    <- cards: fixed-format blocks, no '/'
ATOMIC_POSITIONS {unit}
K_POINTS {type}
CELL_PARAMETERS {unit}   <- only if ibrav = 0
HUBBARD {projector}      <- only for DFT+U (QE >= 7.1)
```

语法规则：

| 规则 | 说明 |
|---|---|
| namelist 顺序 | 必须按 `&CONTROL` → `&SYSTEM` → `&ELECTRONS` → `&IONS` → `&CELL` 的顺序出现，每个以单独一行 `/` 结束 |
| 空 namelist | `&IONS` 在 relax 时即使全用默认值也要写出 `&IONS` 和 `/` 两行 |
| card | 在所有 namelist 之后；card 名后面的选项（如 `crystal`）可以写成 `ATOMIC_POSITIONS crystal` 或 `ATOMIC_POSITIONS {crystal}` |
| 注释 | namelist 中用 `!` 开始注释（Fortran 规则）；card 之间以 `#` 或 `!` 开头的整行会被跳过。card 的数据行内最好不要写注释，把说明放到 namelist 里更稳妥 |
| 字符串 | 必须加引号：`calculation = 'scf'`；单引号、双引号都可以，但必须是 ASCII 引号 |
| 逻辑值 | `.true.` / `.false.`（也可写 `.t.`、`.f.`） |
| 数组 | `celldm(1) = 10.26`、`starting_magnetization(1) = 0.5`，括号内的序号对应 `ATOMIC_SPECIES` 的顺序 |
| 大小写 | 关键词不区分大小写；**UPF 文件名区分大小写**（Linux 文件系统） |
| 浮点数 | `1.0d-8` 与 `1.0e-8` 都可以 |
| 换行 | 必须是 LF（Unix）。CRLF、Tab、全角符号、文件末尾缺少换行都可能导致读取错误 |

> [!TIP]
> 想在提交长作业之前检查输入，有两种办法。一是“干跑”：INPUT_PW 说明，把 `nstep` 设为 0 时，pw.x 初始化后立即停止，只打印摘要，可用于检查输入是否正确。二是用仓库脚本 `QE/scripts/qe_check_input.py` 做静态检查（不运行 pw.x），在仓库根目录运行 `python QE/scripts/qe_check_input.py -h` 查看用法。

---

## 4. `&CONTROL`：做什么、文件放在哪里

| 关键词 | 默认值 | 推荐值 | 说明 |
|---|---|---|---|
| `calculation` | `'scf'` | 按任务 | `'scf'`、`'nscf'`、`'bands'`、`'relax'`、`'md'`、`'vc-relax'`、`'vc-md'` |
| `prefix` | `'pwscf'` | 体系名，如 `'si'` | 输出文件前缀；nscf/bands/ph.x/pp.x 必须使用**相同的** `prefix` 和 `outdir` |
| `outdir` | 环境变量 `ESPRESSO_TMPDIR`，否则 `'./'` | `'./tmp'` 或 scratch 路径 | 存放 `prefix.save/`（电荷密度、波函数、XML）等文件 |
| `pseudo_dir` | 环境变量 `ESPRESSO_PSEUDO`，否则 `'$HOME/espresso/pseudo/'` | 显式写出 | UPF 文件所在目录 |
| `tprnfor` | `relax`/`md`/`vc-md` 时自动为 `.true.` | SCF 时设 `.true.` | 计算并打印原子受力 |
| `tstress` | `.false.`（`vc-relax`/`vc-md` 时自动为 `.true.`） | SCF 时设 `.true.` | 计算并打印应力张量 |
| `verbosity` | `'low'` | `'high'` 用于调试 | `'high'` 会打印更多信息（如每个 k 点的本征值、对称操作） |
| `restart_mode` | `'from_scratch'` | 默认 | `'restart'` 从中断处继续，见 [第18章](18-automation-parallel-troubleshooting.md) |
| `max_seconds` | 1.0d7（相当于无限制） | 略小于作业墙钟时间 | 到时间后干净地停止，便于续算 |
| `etot_conv_thr`、`forc_conv_thr` | 1.0d-4、1.0d-3（a.u.） | 见 [第5章](05-relax-eos.md) | 结构优化的能量与力判据 |
| `disk_io` | scf 为 `'low'`，其它计算为 `'medium'` | 默认 | 控制写盘量；可选 `'high'`、`'medium'`、`'low'`、`'nowf'`、`'minimal'`、`'none'`。`'none'` 什么都不写，后续 nscf/ph.x 无法读取，慎用 |
| `nstep` | scf/nscf/bands 为 1，其它为 50 | 默认 | 离子步数 |

> [!WARNING]
> `wf_collect` 和 `lkpoint_dir` 在当前 INPUT_PW 中标注为 OBSOLETE（不再实现），网上的旧输入经常带有它们。不要再写。

---

## 5. `&SYSTEM`：体系与离散化

### 5.1 必填与常用关键词

| 关键词 | 默认值 | 推荐值 | 说明 |
|---|---|---|---|
| `ibrav` | 必填 | `0`（配合 `CELL_PARAMETERS`）或标准值 | Bravais 晶格编号 |
| `celldm(1)`…`celldm(6)` | — | 与 `A,B,C,cosAB,cosAC,cosBC` 二选一 | `celldm(1)` = a，单位 **Bohr** |
| `A`、`B`、`C`、`cosAB`、`cosAC`、`cosBC` | — | 同上 | `A` = a，单位 **Å** |
| `nat` | 必填 | | 原胞中的原子数 |
| `ntyp` | 必填 | | 原子种类数（= `ATOMIC_SPECIES` 的行数） |
| `ecutwfc` | 必填 | 收敛测试决定 | 波函数截断（Ry） |
| `ecutrho` | `4 * ecutwfc` | NC：4×；US：8–12×；PAW：常为 4–8×，必须测试 | 电荷密度截断（Ry） |
| `occupations` | 绝缘体相当于 `'fixed'` | 金属 `'smearing'` | 也可选 `'tetrahedra'`、`'tetrahedra_lin'`、`'tetrahedra_opt'`、`'from_input'` |
| `smearing` | `'gaussian'` | 金属常用 `'mv'`（冷展宽） | `'gaussian'`、`'mp'`、`'mv'`（= `'cold'`）、`'fd'` |
| `degauss` | 0.0 Ry | 0.01–0.02 Ry（需测试） | 展宽宽度 |
| `nbnd` | 绝缘体为价带数（电子数/2）；金属多 20%（至少多 4 条） | nscf/bands/DOS 时显式设置 | 计算的能带数 |
| `nspin` | 1 | 磁性体系 2 | 非共线用 `noncolin=.true.`，不要写 `nspin=4`（见 [第7章](07-magnetism-soc.md)） |
| `tot_charge` | 0.0 | | +1 表示少一个电子；周期性计算会自动加补偿背景电荷 |
| `input_dft` | 从赝势读取 | 一般不写 | 强制指定泛函，会覆盖赝势中的泛函，谨慎使用 |
| `assume_isolated` | `'none'` | 分子用 `'mt'` | 孤立体系校正，见 §11 |
| `nosym` | `.false.` | | 关闭对称性 |

### 5.2 `ibrav`：常见 Bravais 晶格

下表摘自 INPUT_PW 的 `ibrav` 说明（a = `celldm(1)` 或 `A`）：

| `ibrav` | 晶格 | 晶格矢量 | 典型材料 |
|---|---|---|---|
| 0 | 自由 | 由 `CELL_PARAMETERS` 给出 | 任意结构、从 CIF 转换的结构 |
| 1 | 简单立方 sc | v1 = a(1,0,0), v2 = a(0,1,0), v3 = a(0,0,1) | 分子盒子、钙钛矿 |
| 2 | 面心立方 fcc | v1 = (a/2)(−1,0,1), v2 = (a/2)(0,1,1), v3 = (a/2)(−1,1,0) | Si、Al、Cu、NaCl、GaAs |
| 3 | 体心立方 bcc | v1 = (a/2)(1,1,1), v2 = (a/2)(−1,1,1), v3 = (a/2)(−1,−1,1) | Fe、W |
| 4 | 六方/三方 P | v1 = a(1,0,0), v2 = a(−1/2,√3/2,0), v3 = a(0,0,c/a)；`celldm(3)` = c/a | 石墨、MoS2、ZnO |
| 6 | 简单四方 | `celldm(3)` = c/a | TiO2 金红石 |
| 8 | 简单正交 | `celldm(2)` = b/a，`celldm(3)` = c/a | |

其它 `ibrav`（5、−5、7、9、−9、91、10–14 等）请直接查 INPUT_PW。

> [!NOTE]
> **QE 的 fcc 原胞矢量与 VASP/Materials Project 常见的取法不同**（例如 VASP 常用 (0,a/2,a/2)…），所以从别处拷来的 `crystal` 坐标、k 点路径不能直接套到 `ibrav=2` 上。最安全的做法是 `ibrav=0` + 显式 `CELL_PARAMETERS`，坐标与晶格来自同一个来源（同一个 CIF/POSCAR）。

### 5.3 `ibrav=0` + `CELL_PARAMETERS`（推荐写法）

INPUT_PW 对 `ibrav` 有两条重要提示：

1. `ibrav=0` 时，晶格矢量要给出**足够多的有效数字和正确的对称性**，否则对称性识别可能失败，导致对称化出现奇怪问题；
2. **不要**把 `celldm(1)` 或 `A` 当作 Bohr→Å 的换算因子来用。要么给真实的晶格常数，要么都不给，并在 `CELL_PARAMETERS` 和 `ATOMIC_POSITIONS` 中显式写单位。

```text
&SYSTEM
  ibrav = 0          ! no celldm(1), no A
  nat = 2, ntyp = 1
  ...
/
CELL_PARAMETERS angstrom
  -2.715000000   0.000000000   2.715000000
   0.000000000   2.715000000   2.715000000
  -2.715000000   2.715000000   0.000000000
ATOMIC_POSITIONS crystal
  Si  0.00  0.00  0.00
  Si  0.25  0.25  0.25
```

这与 `ibrav=2, celldm(1)=10.26`（a ≈ 5.43 Å）描述的是同一个晶胞（2.715 = 5.43/2）。

| `CELL_PARAMETERS` 单位 | 含义 |
|---|---|
| `angstrom` / `bohr` | 晶格矢量以 Å / Bohr 为单位，此时 alat = \|v1\| |
| `alat` | 以 `celldm(1)` 或 `A` 为单位 |
| 不写 | DEPRECATED，未来将不允许 |

`ibrav=0` 时 `CELL_PARAMETERS` 必须存在；`ibrav≠0` 时必须**不存在**。同时给出 `celldm` 和 `A` 会报错。

---

## 6. `&ELECTRONS`、`&IONS`、`&CELL`

### 6.1 `&ELECTRONS`

| 关键词 | 默认值 | 推荐值 | 说明 |
|---|---|---|---|
| `conv_thr` | 1.0d-6 | SCF 1e-8 ~ 1e-10；声子前 ≤ 1e-12 | 估计能量误差 < `conv_thr` 时收敛。**它是广延量**（随体系大小变化），大体系可适当放宽 |
| `electron_maxstep` | 100 | 难收敛体系 200–300 | SCF 最大迭代步数 |
| `mixing_beta` | 0.7 | 难收敛时 0.3、0.1 | 密度混合系数 |
| `mixing_mode` | `'plain'` | 表面、slab、非均匀体系 `'local-TF'` | `'TF'` 用于高度均匀的体系 |
| `mixing_ndim` | 8 | 内存紧张时 4 | 混合使用的历史步数 |
| `diagonalization` | `'david'` | 出问题时 `'cg'` 或 `'paro'` | 还有 `'rmm-davidson'`、`'rmm-paro'`；`'ppcg'` 在 7.5 起已不支持 |
| `startingwfc` | `'atomic+random'` | 默认 | 初始波函数 |
| `startingpot` | `'atomic'` | 默认 | `'file'` 从已有电荷密度开始 |

### 6.2 `&IONS` 与 `&CELL`

只有结构优化和 MD 才需要它们，详见 [第5章](05-relax-eos.md) 与 [第17章](17-md-cp.md)。最小写法：

```text
&IONS
  ion_dynamics = 'bfgs'      ! default for 'relax' and 'vc-relax'
/
&CELL
  cell_dynamics = 'bfgs'
  press         = 0.0        ! target pressure, kbar
  press_conv_thr = 0.5       ! kbar (default)
  cell_dofree   = 'all'      ! '2Dxy' for 2D materials
/
```

---

## 7. Card：原子、坐标与 k 点

### 7.1 `ATOMIC_SPECIES`

```text
ATOMIC_SPECIES
  Si  28.0855  Si.pbe-n-rrkjus_psl.1.0.0.UPF
```

每行三列：**标签**、**质量**（amu，C = 12）、**UPF 文件名**。

- 标签可以是元素符号（1–2 个字符），也可以是元素符号加数字或字母，如 `Fe1`、`Fe2`、`C_h`，**总长度不超过 3 个字符**。磁性计算中把同一元素拆成两种标签就靠这一点（[第7章](07-magnetism-soc.md)）。
- 质量只在 MD 和阻尼动力学优化中使用，但会写进数据文件，声子计算会用到它，所以要写对。
- 同一个标签出现两次会报 `two occurrences of the same atomic label`。

### 7.2 `ATOMIC_POSITIONS`

| 单位 | 含义 | 何时用 |
|---|---|---|
| `crystal` | 分数坐标（相对晶格矢量） | 晶体，最常用；与 VASP 的 Direct 相同 |
| `angstrom` | 笛卡尔坐标，Å | 分子、从 xyz 转换的结构 |
| `bohr` | 笛卡尔坐标，Bohr | |
| `alat` | 笛卡尔坐标，以 alat 为单位 | 旧输入常见；不写单位时默认就是它，但已 DEPRECATED |
| `crystal_sg` | 只给不等价原子 + Wyckoff 位置，需要 `space_group` | 高对称结构的简写 |

坐标后可以加三个整数 `if_pos`（0 或 1，默认 1），0 表示该方向的力被置零（在优化/MD 中固定）。下例中 O 原子被完全固定，见 [第5章](05-relax-eos.md)：

```text
ATOMIC_POSITIONS angstrom
  O   6.000  6.000  6.000   0 0 0
  H   6.757  6.586  6.000
```

坐标还可以写成简单的算术表达式（不能有空格，不能以 `+` 开头），例如 `C 1/3 1/2*3^(-1/2) 0`。XCrySDen 不支持这种写法。

> [!NOTE]
> 在 `calculation='nscf'` 或 `'bands'` 时，`ATOMIC_POSITIONS` 中的坐标会被**忽略**，程序使用前一次 SCF 保存的结构。

### 7.3 `K_POINTS`

| 类型 | 写法 | 用途 |
|---|---|---|
| `automatic` | 下一行 `nk1 nk2 nk3 sk1 sk2 sk3` | Monkhorst-Pack 网格；`sk` = 0 不偏移（含 Γ），= 1 偏移半格 |
| `gamma` | 后面不写任何内容 | 只用 Γ 点，波函数取实数，内存和时间约减半；分子、大超胞 |
| `tpiba` | 列表：k 点数，然后每行 `kx ky kz wk`（单位 2π/a） | 默认类型 |
| `crystal` | 同上，倒格矢分数坐标 | Wannier90 全网格（[第15章](15-wannier90.md)） |
| `tpiba_b` / `crystal_b` | 列表，每个点的“权重”是到下一个点之间插入的点数 | 能带路径（[第6章](06-band-dos-pdos.md)） |

```text
K_POINTS automatic
  8 8 8 0 0 0
```

> [!WARNING]
> `K_POINTS gamma` 和 `K_POINTS automatic 1 1 1 0 0 0` 都只用 Γ 点，但前者使用实数波函数的专用算法，后者使用一般 k 点的复数算法。两者结果在数值精度内一致，前者更快更省内存。但 ph.x **不能**读取用 gamma 技巧算出的数据（报错 `cannot start from pw.x data file using Gamma-point tricks`），之后要做声子的分子计算请用 `automatic 1 1 1 0 0 0`（[第12章](12-phonon-dfpt.md)）。

### 7.4 其它 card（简介）

| card | 用途 | 章节 |
|---|---|---|
| `CELL_PARAMETERS {alat\|bohr\|angstrom}` | `ibrav=0` 时的晶格矢量 | 本章 §5.3 |
| `HUBBARD {atomic\|ortho-atomic\|norm-atomic\|wf\|pseudo}` | DFT+U/U+V（7.1+，取代旧关键词） | [第8章](08-dftu-hp.md) |
| `OCCUPATIONS` | `occupations='from_input'` 时手动指定占据 | |
| `CONSTRAINTS` | 约束优化/MD | |
| `ATOMIC_VELOCITIES`、`ATOMIC_FORCES` | MD 初速度、外力 | [第17章](17-md-cp.md) |
| `ADDITIONAL_K_POINTS` | 杂化泛函中附加零权重 k 点画能带 | [第14章](14-hybrid-metagga.md) |

---

## 8. 第一个 SCF：金刚石 Si

### 8.1 输入文件

示例文件 [`si.scf.in`](../examples/02-input-anatomy/si.scf.in)：

```fortran
&CONTROL
  calculation = 'scf'
  prefix      = 'si'
  outdir      = './tmp'
  pseudo_dir  = './pseudo'
  tprnfor     = .true.
  tstress     = .true.
/
&SYSTEM
  ibrav     = 2
  celldm(1) = 10.26          ! Bohr (a = 5.43 Angstrom)
  nat       = 2
  ntyp      = 1
  ecutwfc   = 30.0           ! Ry, SSSP 1.3 efficiency hint for Si
  ecutrho   = 240.0          ! Ry, dual 8 (ultrasoft)
/
&ELECTRONS
  conv_thr    = 1.0d-10
  mixing_beta = 0.7
/
ATOMIC_SPECIES
  Si  28.0855  Si.pbe-n-rrkjus_psl.1.0.0.UPF
ATOMIC_POSITIONS crystal
  Si  0.00  0.00  0.00
  Si  0.25  0.25  0.25
K_POINTS automatic
  8 8 8 0 0 0
```

关键词说明：

| 关键词 | 值 | 为什么 |
|---|---|---|
| `ibrav`, `celldm(1)` | 2, 10.26 Bohr | fcc；10.26 Bohr ≈ 5.43 Å（实验晶格常数附近） |
| `ecutwfc`, `ecutrho` | 30, 240 Ry | SSSP 1.3 Efficiency 对 `Si.pbe-n-rrkjus_psl.1.0.0.UPF` 的推荐值；正式计算前按 [第4章](04-convergence.md) 测试 |
| `occupations` | 不写（相当于 `'fixed'`） | Si 是有带隙的绝缘体，8 个价电子，4 条占据带 |
| `conv_thr` | 1e-10 | 2 原子小体系，严格一些代价很小，力和应力也更干净 |
| `K_POINTS` | 8×8×8，不偏移 | 半导体的常用起点，需要测试 |

同一个结构的 `ibrav=0` 写法见示例 [`si.scf.ibrav0.in`](../examples/02-input-anatomy/si.scf.ibrav0.in)，两者应给出相同的总能量（数值精度内）。

### 8.2 运行

```bash
cd QE/examples/02-input-anatomy
mkdir -p pseudo        # put Si.pbe-n-rrkjus_psl.1.0.0.UPF here (chapter 3)

# serial
pw.x -i si.scf.in > si.scf.out

# MPI, 4 ranks, 2 k-point pools
mpirun -np 4 pw.x -nk 2 -i si.scf.in > si.scf.out
```

| 命令行参数 | 说明 |
|---|---|
| `-i`（同义：`-in`、`-inp`、`-input`） | 输入文件。并行运行时推荐用它，而不是 `< si.scf.in`（部分 MPI 不支持标准输入重定向） |
| `-nk`（同义：`-npool`、`-npools`） | k 点池数，最有效的并行方式 |
| `-nd`（同义：`-ndiag`、`-northo`） | 子空间对角化使用的进程数（需要 ScaLAPACK/ELPA） |
| `-nb`、`-nt`、`-ni` | 能带组、任务组、image（NEB）并行，见 [第18章](18-automation-parallel-troubleshooting.md) |

集群上的作业脚本模板见仓库中的 `QE/scripts/job_qe_slurm.sh`（在仓库根目录运行 `bash QE/scripts/job_qe_slurm.sh -h`，或直接阅读文件头部注释），调度系统用法见 [公共篇第2章](../../Common/tutorials/02-job-scheduler.md)。

### 8.3 读输出：先看头部摘要

```bash
grep -A12 "bravais-lattice index" si.scf.out
```

（示意输出，数值仅供参考）

```text
     bravais-lattice index     =            2
     lattice parameter (alat)  =      10.2600  a.u.
     unit-cell volume          =     xxx.xxxx (a.u.)^3
     number of atoms/cell      =            2
     number of atomic types    =            1
     number of electrons       =         8.00
     number of Kohn-Sham states=            4
     kinetic-energy cutoff     =      30.0000  Ry
     charge density cutoff     =     240.0000  Ry
     scf convergence threshold =      1.0E-10
     mixing beta               =       0.7000
     number of iterations used =            8  plain     mixing
     Exchange-correlation= PBE
```

这一段可以确认：晶格常数是否被正确理解（a.u. 即 Bohr）、电子数是否等于“价电子数 × 原子数”、截断能是否正确、泛函是否来自赝势。还可以检查：

```bash
grep "Sym. Ops" si.scf.out               # number of symmetry operations found
grep "number of k points" si.scf.out     # irreducible k-points actually used
grep "PseudoPot. #" -A3 si.scf.out       # which UPF files were read
```

### 8.4 读输出：SCF 过程与最终能量

```bash
grep -e "total energy" -e "estimated scf accuracy" si.scf.out
grep "^!" si.scf.out
grep "convergence has been achieved" si.scf.out
```

（示意输出，数值仅供参考）

```text
     total energy              =     -xx.xxxxxxxx Ry
     estimated scf accuracy    <       0.0xxxxxxx Ry
     ...
     total energy              =     -xx.xxxxxxxx Ry
     estimated scf accuracy    <          x.xE-11 Ry

!    total energy              =     -xx.xxxxxxxx Ry
     estimated scf accuracy    <          x.xE-11 Ry

     convergence has been achieved in   x iterations
```

| 输出行 | 含义 |
|---|---|
| `total energy = ...`（无 `!`） | 每次 SCF 迭代的能量 |
| `estimated scf accuracy < ...` | 估计的能量误差，小于 `conv_thr` 即收敛 |
| `! total energy = ...` | **最终总能量**，行首有感叹号，最适合 grep；金属（smearing）时是自由能 F = E − TS |
| `Harris-Foulkes estimate` | 另一种能量估计，收敛后应与总能量非常接近 |
| `The total energy is the sum of the following terms:` | 能量分解：one-electron、hartree、xc、ewald 等（金属时另有 `smearing contrib. (-TS)`） |
| `convergence has been achieved in N iterations` | SCF 成功 |
| `convergence NOT achieved after N iterations: stopping` | SCF 失败，见“常见问题”一节 |

> [!NOTE]
> 总能量的绝对值依赖于赝势，不同赝势之间、不同代码之间都**不能直接比较**。有意义的是同一套参数下的**能量差**（形成能、吸附能、能垒等）。

### 8.5 读输出：带隙、费米能

```bash
grep -e "highest occupied" -e "Fermi energy" si.scf.out
```

（示意输出，数值仅供参考）

```text
     highest occupied, lowest unoccupied level (ev):     x.xxxx    x.xxxx
```

| 情况 | 打印内容 |
|---|---|
| `occupations='fixed'` 且 `nbnd` > 占据带数 | `highest occupied, lowest unoccupied level (ev):` 两个数 |
| `occupations='fixed'` 且 `nbnd` = 占据带数（Si 默认 4 条） | 只有 `highest occupied level (ev):` |
| `occupations='smearing'` 或四面体 | `the Fermi energy is ... ev` |

要得到 LUMO，需要在 `&SYSTEM` 中设 `nbnd = 8` 之类的值。注意 SCF 网格上得到的 HOMO/LUMO 只是在这些 k 点上的值，Si 的导带底不在 8×8×8 网格的点上，真正的带隙要看能带图（[第6章](06-band-dos-pdos.md)）。

### 8.6 读输出：力、应力与结束标志

```bash
grep -A4 "Forces acting on atoms" si.scf.out
grep "Total force" si.scf.out
grep -A3 "total   stress" si.scf.out
grep "JOB DONE" si.scf.out
```

（示意输出，数值仅供参考）

```text
     Forces acting on atoms (cartesian axes, Ry/au):

     atom    1 type  1   force =     0.00000000    0.00000000    0.00000000
     atom    2 type  1   force =     0.00000000    0.00000000    0.00000000

     Total force =     0.000000     Total SCF correction =     0.000000

          total   stress  (Ry/bohr**3)                   (kbar)     P=      -x.xx
  -0.0000xxxx   0.00000000   0.00000000           -x.xx        0.00        0.00
  ...

   JOB DONE.
```

金刚石 Si 中每个原子所在位置都有 T<sub>d</sub> 点群对称性，力由对称性决定为 0（输出中是数值上的零）。`P=` 是压强（kbar），负值表示晶胞“想要收缩”，正值表示“想要膨胀”；在 SSSP 推荐截断下，压强可能还没有完全收敛，变胞优化前要测试（[第4章](04-convergence.md)、[第5章](05-relax-eos.md)）。

### 8.7 输出文件

| 文件 / 目录 | 内容 |
|---|---|
| `si.scf.out` | 标准输出（你重定向的文件） |
| `tmp/si.save/data-file-schema.xml` | XML 数据文件：结构、参数、本征值、费米能、总能量等，可用 qeschema、ASE、pymatgen 读取（[第18章](18-automation-parallel-troubleshooting.md)） |
| `tmp/si.save/charge-density.dat` | 电荷密度，nscf/bands/pp.x 会读取 |
| `tmp/si.save/wfc*.dat` | 波函数 |
| `tmp/si.save/*.UPF` | 用到的赝势副本 |
| `tmp/si.wfc*`、`tmp/si.mix*` 等 | 运行中的临时文件 |

```bash
grep -i "Writing" si.scf.out              # e.g. "Writing all to output data dir ./tmp/si.save/"
ls tmp/si.save/
```

也可以用仓库脚本 `QE/scripts/qe_parse_out.py` 一次性提取总能量、费米能/HOMO-LUMO、力、应力等（在仓库根目录运行 `python QE/scripts/qe_parse_out.py -h` 查看用法）。

---

## 9. 金属：fcc Al 与 smearing

金属的费米面会切过能带，`occupations='fixed'` 不再适用。对电子数为奇数的体系（例如每个原胞 3 个价电子的 Al），如果不设 `occupations`，pw.x 7.6 只打印一条提示 `Message from routine setup: the system is metallic, specify occupations`，**不会停止**，而是继续用固定占据计算，结果没有意义；对自旋极化体系则直接报错 `spin-polarized system, specify occupations`。所以金属一定要自己写 `occupations='smearing'`。

示例文件 [`al.scf.in`](../examples/02-input-anatomy/al.scf.in) 的关键部分：

```fortran
&SYSTEM
  ibrav     = 2
  celldm(1) = 7.65           ! Bohr (about 4.05 Angstrom)
  nat       = 1
  ntyp      = 1
  ecutwfc   = 30.0
  ecutrho   = 240.0
  occupations = 'smearing'
  smearing    = 'mv'         ! Marzari-Vanderbilt cold smearing
  degauss     = 0.02         ! Ry
/
...
K_POINTS automatic
  12 12 12 0 0 0
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `occupations` | `'smearing'` | 金属必须 |
| `smearing` | `'mv'` 或 `'mp'` | 冷展宽和 MP 展宽对总能量的展宽误差较小；`'gaussian'` 需要更小的 `degauss`；`'fd'` 对应物理的电子温度 |
| `degauss` | 0.01–0.02 Ry | 必须和 k 点网格一起测试（[第4章](04-convergence.md)） |
| k 网格 | 比同尺寸半导体更密 | 费米面附近需要足够的 k 点 |

输出中多出来的行：

```bash
grep -e "Fermi energy" -e "smearing contrib" -e "number of k points" al.scf.out
```

（示意输出，数值仅供参考）

```text
     number of k points=    xx  Marzari-Vanderbilt smearing, width (Ry)=  0.0200
     the Fermi energy is     x.xxxx ev
     smearing contrib. (-TS)   =      -0.0000xxxx Ry
```

> [!TIP]
> 金属时能量分解的标题会变成 `The total energy is F=E-TS. E is the sum of the following terms:`。`smearing contrib. (-TS)` 的绝对值应当远小于你关心的能量精度；如果它达到你关心的精度量级（例如 meV/atom），说明 `degauss` 偏大，应和 k 网格一起重新测试。

---

## 10. 运行流程中的“一致性”规则

后续章节会在 SCF 之后运行 nscf、bands、pp.x、ph.x 等，它们之间靠 `prefix` 和 `outdir` 找到数据：

```text
pw.x (scf, prefix='si', outdir='./tmp')
   -> writes ./tmp/si.save/
pw.x (nscf/bands, same prefix and outdir)
   -> reads ./tmp/si.save/charge-density.dat
dos.x / projwfc.x / bands.x / pp.x / ph.x (same prefix and outdir)
```

| 规则 | 违反时的现象 |
|---|---|
| `prefix`、`outdir` 必须相同 | 后续程序报找不到文件，或读到另一个体系的数据 |
| 同一目录下不同作业不要共用同一个 `prefix` + `outdir` | 互相覆盖，出现 `davcio` 读写错误（PW 用户手册 Troubleshooting） |
| nscf/bands 的 `ecutwfc`、`ecutrho`、赝势、泛函与 SCF 一致 | 结果无意义或报错 |
| nscf/bands 中显式设置 `nbnd` | 默认只算占据带，看不到导带 |

---

## 11. 分子：放进大盒子的 H2O

平面波代码总是周期性的。计算孤立分子时，要把分子放进足够大的盒子，并只用 Γ 点。

示例文件 [`h2o.scf.in`](../examples/02-input-anatomy/h2o.scf.in)：

```fortran
&CONTROL
  calculation = 'scf'
  prefix      = 'h2o'
  outdir      = './tmp'
  pseudo_dir  = './pseudo'
  tprnfor     = .true.
/
&SYSTEM
  ibrav   = 1
  A       = 12.0             ! Angstrom, cubic box
  nat     = 3
  ntyp    = 2
  ecutwfc = 60.0             ! max of SSSP hints: H 60 Ry, O 50 Ry
  ecutrho = 480.0            ! max of SSSP hints: H 480 Ry, O 400 Ry
  nbnd    = 8                ! 4 occupied + 4 empty -> HOMO and LUMO printed
  assume_isolated = 'mt'     ! Martyna-Tuckerman correction
/
&ELECTRONS
  conv_thr = 1.0d-9
/
ATOMIC_SPECIES
  O  15.999  O.pbe-n-kjpaw_psl.0.1.UPF
  H   1.008  H.pbe-rrkjus_psl.1.0.0.UPF
ATOMIC_POSITIONS angstrom
  O   6.000000   6.000000   6.000000
  H   6.757000   6.586000   6.000000
  H   5.243000   6.586000   6.000000
K_POINTS gamma
```

| 关键词 | 值 | 说明 |
|---|---|---|
| `ibrav=1` + `A` | 12 Å 立方盒子 | 分子周围留出足够真空；需要对盒子大小做收敛测试 |
| `ecutwfc`/`ecutrho` | 60/480 Ry | 多元素体系取**各元素推荐值中的最大值**（[第3章](03-pseudopotentials.md)） |
| `assume_isolated` | `'mt'` | Martyna–Tuckerman 校正，同时修正总能量和 SCF 势，去除周期镜像间的静电相互作用 |
| `nbnd` | 8 | H2O 有 8 个价电子、4 条占据轨道；多算 4 条空轨道，输出才会给出 LUMO |
| `K_POINTS gamma` | | 分子没有色散，Γ 点足够，并启用实数波函数算法（之后若要用 ph.x，改用 `automatic 1 1 1 0 0 0`，见 §7.3） |
| 坐标 | `angstrom` | 从 xyz 文件拷贝最方便；分子放在盒子中间 |

`assume_isolated` 的选项（INPUT_PW）：

| 值 | 说明 |
|---|---|
| `'none'` | 默认，普通周期性计算 |
| `'makov-payne'`（`'m-p'`、`'mp'`） | 只修正总能量，并估计真空能级；**仅限立方体系**（`ibrav`=1、2、3） |
| `'martyna-tuckerman'`（`'m-t'`、`'mt'`） | 修正总能量和 SCF 势；分子/团簇推荐 |
| `'esm'` | 有效屏蔽介质，用于带电或极化的 slab，见 [第9章](09-vdw-surface-workfunction.md) |
| `'2D'` | z 方向截断库仑相互作用，用于二维材料 |

> [!TIP]
> Martyna–Tuckerman 校正并不意味着盒子可以很小：电荷密度在盒子边界附近必须已经衰减到可以忽略。官方示例 `PW/examples/cluster_example` 用 16 Bohr 的立方盒子，并提供了 12、16、20、24 Bohr 的参考输出，用来检查结果对盒子大小的收敛。对你自己的分子也应该做同样的测试。

分子的 HOMO/LUMO：对 `occupations='fixed'` 的闭壳层分子，设 `nbnd` 大于占据轨道数（H2O 有 8 个价电子、4 个占据轨道；示例中设了 `nbnd = 8`），输出中就会有 `highest occupied, lowest unoccupied level`。几何优化（键长、键角）见 [第5章](05-relax-eos.md)。

---

## 参数选择与收敛建议

| 参数 | 起点 | 何时调整 |
|---|---|---|
| `ecutwfc`/`ecutrho` | SSSP 推荐值（多元素取最大值） | 正式计算前按 [第4章](04-convergence.md) 测试；需要应力时通常要提高 |
| `K_POINTS` | 半导体 8×8×8、金属 12×12×12（小原胞）；分子 `gamma` | 按 k 点密度和晶格常数换算，测试至 ΔE < 1 meV/atom |
| `occupations` | 绝缘体默认；金属 `'smearing'` + `'mv'` + 0.01–0.02 Ry | 不确定是否为金属时，先用 smearing 跑一次看带隙 |
| `conv_thr` | 1e-8（一般 SCF）；1e-10（小原胞、需要精确的力和应力）；≤ 1e-12（声子前） | 大体系可按原子数适当放宽 |
| `mixing_beta` | 0.7 | 不收敛时依次试 0.3、0.1 |
| 坐标/晶格单位 | `ibrav=0` + `CELL_PARAMETERS angstrom` + `ATOMIC_POSITIONS crystal` | 永远显式写单位 |
| `outdir` | 本地 scratch 或 `./tmp` | 不同作业用不同 `prefix` 或不同 `outdir` |

---

## 常见问题与报错

报错格式为 `Error in routine <name> (<code>): <message>`，下面列出 QE 源码中真实存在的消息文本：

| 报错 / 现象 | 原因 | 解决 |
|---|---|---|
| `Error in routine read_namelists ... bad line in namelist &system: "..." (error could be in the previous line)` | 关键词拼错、字符串缺引号、全角引号、上一行有问题 | 检查报错行**和上一行**；用 `file x.in`、`cat -A x.in` 查看隐藏字符 |
| `could not find namelist &electrons` | 缺少某个 namelist，或前一个 namelist 漏了 `/`、顺序写错 | 按 `&CONTROL`→`&SYSTEM`→`&ELECTRONS`→`&IONS`→`&CELL` 的顺序写全，每个以 `/` 结束 |
| 读取错误但信息很模糊（`Error while reading XXX namelist`） | CRLF 换行、Tab、非 ASCII 字符、文件末尾没有换行（PW 用户手册 5.0.0.4） | `sed -i 's/\r$//' x.in`；确保最后一行后有换行 |
| `Error in routine cell_base_init ... ibrav=0: must read cell parameters` | `ibrav=0` 但没给 `CELL_PARAMETERS` | 补上 `CELL_PARAMETERS angstrom` |
| `do not specify both celldm and a,b,c!` | 同时给了 `celldm(1)` 和 `A` | 二选一 |
| `species X in ATOMIC_POSITIONS is nonexistent` | 坐标中的标签不在 `ATOMIC_SPECIES` 中（大小写/拼写不一致） | 统一标签 |
| `two occurrences of the same atomic label` | `ATOMIC_SPECIES` 中标签重复 | 删除重复行或改用 `Fe1`/`Fe2` |
| `Error in routine readpp ... file ./pseudo/Si.UPF not found` | `pseudo_dir` 错误或文件名大小写不对 | `ls $pseudo_dir`；Linux 下 `.UPF` 与 `.upf` 不同 |
| `inconsistent DFT read from PP files` | 不同元素的赝势泛函不同（如 LDA 与 PBE 混用） | 换成同一泛函的赝势（[第3章](03-pseudopotentials.md)） |
| `DEPRECATED: no units specified in ATOMIC_POSITIONS card` | 坐标没写单位 | 显式写 `crystal`/`angstrom` |
| `Message from routine setup: the system is metallic, specify occupations`（只是提示，程序继续运行） | 电子数为奇数且未设 `occupations` | 设 `occupations='smearing'` 后重算；带着这条提示得到的结果不可用 |
| `spin-polarized system, specify occupations` | `nspin=2` 但未设 `occupations` | 设 smearing 或固定磁矩 |
| `ecutrho < 4*ecutwfc, are you sure?` | `ecutrho` 小于 4 倍 `ecutwfc` | 一般是笔误，NC 至少 4×，US/PAW 8–12× |
| `convergence NOT achieved after 100 iterations: stopping` | SCF 不收敛 | 降低 `mixing_beta`（0.3/0.1）、增大 `electron_maxstep`；金属用 smearing；表面用 `mixing_mode='local-TF'`；检查结构是否合理 |
| `S matrix not positive definite` 或 `problems computing cholesky`（cdiaghg） | 原子过近（常见于 Bohr/Å 混用）、赝势有 ghost 态、`ecutrho` 太小 | 先检查结构（用 VESTA/XCrySDen 看一眼）；再试 `diagonalization='cg'` 或提高 `ecutrho` |
| `too many bands are not converged`（c_bands） | 对角化失败 | `diagonalization='cg'`；检查结构 |
| 晶胞体积异常小、能量离谱 | 把 Å 写进了 `celldm(1)` | `celldm(1)` 用 Bohr，或改用 `A`（Å） |
| 并行运行停在开头 | 用了 `< x.in` 重定向 | 改用 `-i x.in` |

---

## 练习

1. 把 `si.scf.in` 改写成 `ibrav=0` + `CELL_PARAMETERS angstrom` + `ATOMIC_POSITIONS crystal`，与示例 `si.scf.ibrav0.in` 对照，然后（自己运行）比较两者的总能量和对称操作数。
2. 故意把 `celldm(1)` 写成 5.43，预测 pw.x 会出现什么现象，再查看头部摘要中的 `unit-cell volume` 来验证。
3. 在 `si.scf.in` 中设 `nbnd = 8`，找到输出中的 `highest occupied, lowest unoccupied level`，解释为什么这个差值不是 Si 的真实带隙。
4. 对 `al.scf.in`，分别用 `smearing='gaussian'` 和 `'mv'`、`degauss=0.01/0.02/0.04` 运行，比较 `!` 总能量与 `smearing contrib. (-TS)`。
5. 把 `h2o.scf.in` 的盒子从 10 Å 逐步增大到 16 Å，分别在 `assume_isolated='none'` 和 `'mt'` 下比较总能量随盒子大小的变化。
6. 用 `grep` 写一个单行命令，从一批 `*.out` 文件中提取文件名和 `!` 总能量；再与 `qe_parse_out.py` 的输出对比。

---

## 参考资料

**官方文档**

- [INPUT_PW（pw.x 全部关键词与 card 说明）](https://www.quantum-espresso.org/Doc/INPUT_PW.html)：`ibrav`、`celldm`、`ATOMIC_POSITIONS`、`K_POINTS`、`assume_isolated`、`occupations` 等
- [PWscf User's Guide](https://www.quantum-espresso.org/Doc/pw_user_guide/)：输入输出、性能、Troubleshooting
- [QE 总用户手册（并行与命令行参数）](https://www.quantum-espresso.org/Doc/user_guide/)
- 官方示例：QE 源码中的 `PW/examples/example01`（Si、Al、Cu、Ni）与 `PW/examples/cluster_example`（Martyna–Tuckerman，N、NH4+、H2O），在线浏览：[QEF/q-e PW/examples](https://github.com/QEF/q-e/tree/qe-7.6/PW/examples)
- [Materials Cloud QE 输入生成器](https://qeinputgenerator.materialscloud.io/)：上传结构自动生成带 SSSP 参数的 pw.x 输入，适合作为起点

**课程与教程（推荐，链接而非复制）**

- [Imperial College MSE404-MM Lab 2: Quantum Espresso Input and Output for Molecules](https://imperialcollegelondon.github.io/MSE404-MM/labs/lab02/readme/)：从分子出发讲 pw.x 输入输出，与本章 §11 互补
- [Imperial College MSE404-MM Lab 3](https://imperialcollegelondon.github.io/MSE404-MM/labs/lab03/readme/)：周期性晶体与收敛
- [Cornell PARADIM Computational Tutorials](https://www.paradim.org/toolbox/computation/tutorials)：“QE 输入与收敛参数”一讲讲解清楚（基于 QE 6.x 与 LDA 赝势）
- [Cornell PARADIM Summer School 2016](https://www.paradim.org/cu_ss_2016)：分子与体相的平衡结构、测试运行（讲义只链接，不复制）
- [pranabdas/espresso：SCF](https://pranabdas.github.io/espresso/hands-on/scf)：带 Jupyter 绘图的 SCF 入门
- [nanoHUB dftqe](https://nanohub.org/resources/dftqe)：在浏览器中直接运行 pw.x

**相关仓库（仅链接）**

- [QEF/q-e](https://github.com/QEF/q-e)：官方源码镜像，`PW/examples` 是最权威的输入示例
- [aiidateam/qe-tools](https://github.com/aiidateam/qe-tools)：纯 Python 解析 pw.x 输入输出
- [QEF/qeschema](https://github.com/QEF/qeschema)：读取 `data-file-schema.xml`
- [materialscloud-org/seekpath](https://github.com/materialscloud-org/seekpath)：标准原胞与高对称 k 路径，在线版可导出 QE 格式
- [yyyu200/DFTbook](https://github.com/yyyu200/DFTbook)：以 QE 为例的中文 DFT 入门书

**关键文献**

- P. Giannozzi et al., *J. Phys.: Condens. Matter* 21, 395502 (2009). DOI: [10.1088/0953-8984/21/39/395502](https://doi.org/10.1088/0953-8984/21/39/395502)
- P. Giannozzi et al., *J. Phys.: Condens. Matter* 29, 465901 (2017). DOI: [10.1088/1361-648X/aa8f79](https://doi.org/10.1088/1361-648X/aa8f79)
- G. J. Martyna, M. E. Tuckerman, *J. Chem. Phys.* 110, 2810 (1999). DOI: [10.1063/1.477923](https://doi.org/10.1063/1.477923)
- G. Makov, M. C. Payne, *Phys. Rev. B* 51, 4014 (1995). DOI: [10.1103/PhysRevB.51.4014](https://doi.org/10.1103/PhysRevB.51.4014)
- N. Marzari, D. Vanderbilt, A. De Vita, M. C. Payne, *Phys. Rev. Lett.* 82, 3296 (1999). DOI: [10.1103/PhysRevLett.82.3296](https://doi.org/10.1103/PhysRevLett.82.3296)
- H. J. Monkhorst, J. D. Pack, *Phys. Rev. B* 13, 5188 (1976). DOI: [10.1103/PhysRevB.13.5188](https://doi.org/10.1103/PhysRevB.13.5188)

---

[← 上一章](01-install.md) | [返回目录](../README.md) | [下一章 →](03-pseudopotentials.md)
