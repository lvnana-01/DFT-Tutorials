# 第9章 范德华修正、二维材料、表面与功函数（vdW Corrections, 2D Materials, Surfaces & Work Function）

本章讲 Quantum ESPRESSO 中层状材料与表面的常用设置：经验色散修正（`vdw_corr` = DFT-D2 / DFT-D3 / TS / MBD / XDM）与非局域泛函（`input_dft` = `vdw-df2-b86r`、`rvv10` 等），二维材料的真空层、`cell_dofree='2Dxy'` 与 `assume_isolated='2D'`，slab 建模与固定底层（`if_pos`），偶极修正（`tefield` + `dipfield`、`edir`、`emaxpos`、`eopreg`、`eamp`）与 ESM，最后用 pp.x（`plot_num=11`）+ average.x 求功函数并计算表面能（Quantum ESPRESSO vdW-D3 / 2D materials / slab dipole correction / work function tutorial）。

> [!NOTE]
> 全部输入在示例目录 [09-vdw-surface-workfunction](../examples/09-vdw-surface-workfunction/README.md)。关键词对照 QE 7.5 的 `INPUT_PW`、`INPUT_PP` 与源码（`set_vdw_corr.f90`、`funct.f90`、`average.f90`、`punch_plot.f90`）核实。本仓库**不运行任何计算**，示意输出只演示格式，数值仅供参考。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 为层状/分子/吸附体系选择色散修正，写对 `vdw_corr`、`dftd3_version`、`dftd3_threebody`，或用 `input_dft` 选非局域泛函；
- 搭建二维材料模型：真空厚度、k 网格（`N N 1`）、只优化面内晶格（`cell_dofree='2Dxy'`）、必要时 2D 库仑截断；
- 搭建金属 slab：层数、真空、固定底层、`mixing_mode='local-TF'`；
- 正确打开偶极修正，并把锯齿势的位置放在真空里；
- 用 pp.x 输出 V_bare + V_H，average.x 求平面平均与宏观平均，得到 Φ = V_vac − E_F；
- 用取向一致的体相参考计算表面能，并做层数收敛。

| 场景 | 推荐做法 | 本章小节 |
|---|---|---|
| 层状晶体（石墨、h-BN、MoS₂ 体相） | D3(BJ) 或 rVV10 / vdW-DF2-B86R，vc-relax | §3 |
| 分子吸附、物理吸附 | D3(BJ)；需要多体效应时 MBD | §3 |
| 单层/少层二维材料 | 真空 ≥ 15 Å，`2Dxy`，k `N N 1` | §4 |
| 带电/门控二维体系、2D 声子 | `assume_isolated='2D'` | §4.3 |
| 金属表面、吸附能 | slab + 固定底层 + 偶极修正 | §5、§6 |
| 功函数、表面偶极 | pp.x `plot_num=11` + average.x | §7 |
| 带电 slab、电极界面 | ESM（`assume_isolated='esm'`） | §5.4 |

### 与 VASP 的对应关系

| VASP | QE（pw.x / pp.x） | 说明 |
|---|---|---|
| `IVDW = 11` / `12` | `vdw_corr='grimme-d3'` + `dftd3_version = 3` / `4` | D3 零阻尼 / BJ 阻尼 |
| `IVDW = 1`（D2） | `vdw_corr='grimme-d2'` | |
| `IVDW = 2/20`（TS） | `vdw_corr='ts'` | |
| `IVDW = 202`（MBD） | `vdw_corr='mbd'` | |
| `GGA = MK` + `LUSE_VDW`（optB86b-vdW 等） | `input_dft='vdw-df-ob86'` 等 | QE 不需要外部 kernel 文件 |
| `ISIF = 3` + 固定 c（需补丁或 `LATTICE_CONSTRAINTS`） | `cell_dofree = '2Dxy'` | QE 原生支持 |
| `LDIPOL = .TRUE.`、`IDIPOL = 3`、`DIPOL` | `tefield`、`dipfield`、`edir = 3`、`emaxpos`、`eopreg` | QE 用锯齿势实现 |
| `LVHAR = .TRUE.` → LOCPOT | pp.x `plot_num = 11` | V_bare + V_H，不含 V_xc |
| `vaspkit` 平面平均 | average.x（`avg.dat`） | |
| Selective dynamics `F F F` | `ATOMIC_POSITIONS` 末尾 `0 0 0` | `if_pos` |

VASP 侧的讲解见 [VASP 第8章](../../VASP/tutorials/08-vdw-2d.md) 与 [VASP 第9章](../../VASP/tutorials/09-surface-workfunction.md)。

---

## 2. 原理速览

### 2.1 为什么需要色散修正

LDA/GGA 只依赖局部密度及其梯度，无法描述长程电子关联引起的色散（London）作用。后果：PBE 下石墨层间几乎不束缚，层状材料 c 轴严重偏大，物理吸附能偏小。两类补救：

1. **原子对加和修正**（DFT-D2/D3、TS、MBD、XDM）：$E=E_{\rm DFT}+E_{\rm disp}$，$E_{\rm disp}=-\sum_{i<j}\sum_{n=6,8} s_n\,C_n^{ij}/r_{ij}^n\,f_{\rm damp}(r_{ij})$。计算量几乎为零，只修正能量、力、应力，不改变电子结构。
2. **非局域关联泛函**（vdW-DF 系列、rVV10）：在交换关联泛函中加入非局域关联项 $E_c^{\rm nl}[n]$，自洽地影响密度。计算量略高于 GGA。QE 6.5 起 kernel 在运行时现算，不再需要 `vdW_kernel_table` 文件。

### 2.2 表面能

对上下两个表面相同的化学计量 slab：

$$
\gamma=\frac{E_{\rm slab}(N)-N\,\varepsilon_{\rm bulk}}{2A}
$$

N 为 slab 原子数，$\varepsilon_{\rm bulk}$ 为体相每原子能量，A 为表面胞面积。只弛豫一侧（底层固定）的 slab，两侧不等价，可拆为

$$
\gamma_{\rm rel}=\frac{E_{\rm slab}^{\rm rel}-N\varepsilon_{\rm bulk}}{A}-\gamma_{\rm unrel},\qquad
\gamma_{\rm unrel}=\frac{E_{\rm slab}^{\rm unrel}-N\varepsilon_{\rm bulk}}{2A}
$$

单位换算：1 Ry = 13.6057 eV；1 eV/Å² = 16.02 J/m²。

### 2.3 功函数

$$
\Phi=V_{\rm vac}-E_F
$$

$V_{\rm vac}$ 是真空中平坦区域的电子静电势能，$E_F$ 取自同一次 SCF 的 `the Fermi energy is`。V_xc 在真空中趋于 0，所以用 V_bare + V_H（pp.x `plot_num=11`）即可读出真空能级。

### 2.4 偶极修正

不对称 slab（单侧吸附、极性面）带净偶极 $p_z$。周期边界下这些偶极层周期排列，在真空中产生人为电场，真空势出现斜率。偶极修正（Neugebauer–Scheffler 1992；Bengtsson 1999）在真空中加一个锯齿势抵消该场：两侧真空各自平坦，中间出现台阶，台阶高度即两侧功函数之差

$$
\Delta\Phi=\frac{e\,p_z}{\varepsilon_0 A}
$$

QE 中这是"外加锯齿电场"功能（`tefield`）的一个特例：`dipfield = .true.` 让程序根据自洽偶极自动设定锯齿幅度，`eamp = 0` 表示不再额外加外场。

---

## 3. 范德华修正

### 3.1 关键词

`vdw_corr` 的合法取值（`Modules/set_vdw_corr.f90`，QE 7.5）：

| `vdw_corr` 取值（任一别名） | 方法 | 可选参数 |
|---|---|---|
| `'grimme-d2'`, `'Grimme-D2'`, `'DFT-D'`, `'dft-d'` | DFT-D2 | `london_s6`（默认 0.75，适合 PBE）、`london_rcut`（默认 200 a.u.）、`london_c6`、`london_rvdw` |
| `'grimme-d3'`, `'Grimme-D3'`, `'DFT-D3'`, `'dft-d3'` | DFT-D3 | `dftd3_version`、`dftd3_threebody` |
| `'TS'`, `'ts'`, `'ts-vdw'`, `'ts-vdW'`, `'tkatchenko-scheffler'` | Tkatchenko–Scheffler | `ts_vdw_econv_thr`（默认 1.D-6）、`ts_vdw_isolated` |
| `'MBD'`, `'mbd'`, `'many-body-dispersion'`, `'mbd_vdw'` | 多体色散 | `ts_vdw_isolated` |
| `'XDM'`, `'xdm'` | 交换空穴偶极矩模型 | `xdm_a1`、`xdm_a2`（PBE、B86bPBE、PW86PBE、BLYP 不必给） |
| `'none'`（默认） | 无 | |

`dftd3_version`（默认 3）：

| 值 | 含义 |
|---|---|
| 2 | 原始 Grimme-D2 参数化（此时三体项总是关闭） |
| 3 | D3 零阻尼（默认） |
| 4 | D3 BJ 阻尼 |
| 5 | D3M 零阻尼 |
| 6 | D3M BJ 阻尼 |

`dftd3_threebody` 默认 `.true.`（打开 ATM 三体项）。`INPUT_PW` 提示：并非所有泛函都有 D3 参数。

> [!IMPORTANT]
> TS 与 MBD 在 pw.x 中**只能配合模守恒（NC）赝势**：`PW/src/setup.f90`（QE 7.5）对超软/PAW 赝势直接报 `Tkatchenko-Scheffler not implemented with USPP/PAW`（`vdw_corr='mbd'` 内部同样打开 TS，受同一限制）。本章示例使用 SSSP 中的 PAW 碳赝势，因此 `run_vdw_compare.sh` 默认不跑 `ts`；要比较 TS，请换成 NC 赝势（例如 PseudoDojo 的 C 赝势，见 [第3章](03-pseudopotentials.md)）。D2/D3/XDM 与非局域泛函没有这一限制。

> [!WARNING]
> `vdw_corr` 拼错（例如写成 `'d3'`）在 QE 7.5 中**不会报错**：程序只打印 `WARNING: unknown vdw correction (vdw_corr): d3. No vdw correction used.` 然后照常运行。每次都要 grep 输出中的色散能量项确认生效。另外 `london`、`xdm`、`ts_vdw`、`mbd_vdw` 这类旧逻辑开关已标记为 OBSOLESCENT，会被翻译为对应的 `vdw_corr` 并打印提示。

非局域泛函通过 `input_dft` 选择（`Modules/funct.f90` 中的名字，大小写不敏感）：

| `input_dft` | 泛函 | 备注 |
|---|---|---|
| `'vdw-df'` | vdW-DF1（Dion 等 2004） | 层间距常偏大 |
| `'vdw-df2'` | vdW-DF2（Lee 等 2010） | |
| `'vdw-df-ob86'` | optB86b-vdW | 固体常用 |
| `'vdw-df-obk8'` | optB88-vdW | |
| `'vdw-df2-b86r'` | rev-vdW-DF2（Hamada 2014） | 层状材料表现好 |
| `'vdw-df-cx'` | vdW-DF-cx（Berland–Hyldgaard 2014） | |
| `'vdw-df3-opt1'`, `'vdw-df3-opt2'` | vdW-DF3 | 7.6 另增 vdW-DF3-mc |
| `'rvv10'` | rVV10（Sabatini 等 2013） | |

`INPUT_PW` 明确：非局域泛函**不在 `vdw_corr` 中指定**，而是用 `input_dft`。`input_dft` 会覆盖赝势中记录的泛函（"Use with care"）；严格的做法是使用与其交换部分一致的赝势，实际工作中用 PBE 赝势是常见近似，论文中应注明。

### 3.2 石墨：输入（`graphite.vdw.in`）

```fortran
&CONTROL
  calculation   = 'vc-relax'
  prefix        = 'graphite'
  outdir        = './tmp'
  pseudo_dir    = './pseudo'
  tprnfor       = .true.
  tstress       = .true.
  forc_conv_thr = 1.0d-4
  etot_conv_thr = 1.0d-6
/
&SYSTEM
  ibrav       = 4                ! hexagonal
  A           = 2.464            ! in-plane lattice constant (Angstrom)
  C           = 6.80             ! start slightly above the exp. 6.71 Angstrom
  nat         = 4
  ntyp        = 1
  ecutwfc     = 50.0
  ecutrho     = 400.0
  occupations = 'smearing'
  smearing    = 'mv'
  degauss     = 0.01
  vdw_corr    = 'grimme-d3'      ! @VDW@
  dftd3_version = 4              ! D3(BJ); ignored by other methods
/
&ELECTRONS
  conv_thr    = 1.0d-10
  mixing_beta = 0.5
/
&IONS
  ion_dynamics = 'bfgs'
/
&CELL
  cell_dynamics = 'bfgs'
  press_conv_thr = 0.5
/
ATOMIC_SPECIES
  C  12.011  C.pbe-n-kjpaw_psl.1.0.0.UPF
ATOMIC_POSITIONS crystal
  C  0.000000  0.000000  0.250000
  C  0.000000  0.000000  0.750000
  C  0.333333  0.666667  0.250000
  C  0.666667  0.333333  0.750000
K_POINTS automatic
  12 12 4 0 0 0
```

C 的赝势与截断取 SSSP 1.3 PBE Efficiency（`C.pbe-n-kjpaw_psl.1.0.0.UPF`，45 / 360 Ry），这里略微提高到 50 / 400 Ry，因为 vc-relax 需要更收敛的应力（Pulay 应力，见 [第5章](05-relax-eos.md)）。

### 3.3 运行与检查

```bash
cd QE/examples/09-vdw-surface-workfunction
bash run_vdw_compare.sh d3bj rvv10 none
cat vdw_summary.dat

grep 'DFT-D3 Dispersion' vdw_d3bj/graphite.out | tail -1
grep 'Dispersion Correction' vdw_d2/graphite.out | tail -1
grep -A4 'CELL_PARAMETERS' vdw_d3bj/graphite.out | tail -4
```

总能分解中的色散项（`PW/src/electrons.f90` 的格式）：D3 打印 `DFT-D3 Dispersion =`，D2 与 TS/MBD 打印 `Dispersion Correction =`，XDM 打印 `Dispersion XDM Correction =`。

（示意输出，数值仅供参考）

```text
method       c_final_A         E_final_Ry
none             7.xxx       -xx.xxxxxxxx
d3bj             6.7xx       -xx.xxxxxxxx
rvv10            6.7xx       -xx.xxxxxxxx
```

成功判据：

- 有色散修正的方法 c 落在实验值 6.71 Å 附近（层间距 3.35 Å），而纯 PBE 的 c 偏大很多甚至持续增大；
- `bfgs converged` 出现，末尾以新晶胞重做的 SCF 与最后一步压强差别不大（否则提高截断）；
- **不同方法之间只比较几何量**，不比较总能（它们是不同的能量泛函）。

> [!NOTE]
> `vdw_summary.dat` 中的 c 取自 `Begin final coordinates` 之后 `CELL_PARAMETERS` 的第三行。`PW/src/output_tau.f90` 按输入的单位打印该块：用 `ibrav=4` + `A` 输入时标题为 `CELL_PARAMETERS (alat= <alat in bohr>)`，数值以 alat 为单位，脚本会读出 alat 并换算成 Å（同时兼容 `bohr`、`angstrom` 标题）。

### 3.4 方法选择建议

| 体系 | 首选 | 备选 |
|---|---|---|
| 层状晶体结构（石墨、h-BN、TMD） | D3(BJ)、rVV10、vdW-DF2-B86R | optB86b-vdW |
| 分子在金属表面的物理吸附 | D3(BJ) | MBD（多体极化重要时，需 NC 赝势）、vdW-DF-cx |
| 分子晶体 | D3(BJ)、MBD（需 NC 赝势） | XDM |
| 与杂化泛函联用 | D3（参数化 HSE06/PBE0） | 非局域杂化（7.x 有 vdW-DF 杂化） |

---

## 4. 二维材料

### 4.1 建模要点

| 项目 | 建议 |
|---|---|
| 真空 | ≥ 15 Å；做收敛（15、20、25 Å 比较能带、能量） |
| k 网格 | `N N 1 0 0 0`，面外只取 1 |
| 晶胞优化 | `cell_dofree = '2Dxy'`：只动 x、y 分量，真空厚度不变 |
| 层间相互作用 | 多层/异质结必须加色散修正 |
| 偶极 | 非对称 2D（Janus MoSSe、单侧吸附）加偶极修正（§6） |
| 带电/门控 | `assume_isolated = '2D'`，或 `gate = .true.` |

`cell_dofree` 的相关取值（`INPUT_PW`）：`'2Dxy'` 只允许 x、y 分量变化；`'2Dshape'` 同上但保持 xy 面积不变；`'epitaxial_ab'` 固定 a、b 轴只动 c；`'ibrav'` 保持 ibrav 的晶格类型（可组合为 `'ibrav+2Dxy'`，对六角晶格可用）。`INPUT_PW` 提醒：轴不正交时某些选项会破坏对称性。

### 4.2 单层 MoS₂（`mos2_mono.vc-relax.in`）

```fortran
&SYSTEM
  ibrav       = 4
  A           = 3.19             ! in-plane lattice constant (Angstrom)
  C           = 20.0             ! about 17 Angstrom of vacuum
  nat         = 3
  ntyp        = 2
  ecutwfc     = 45.0
  ecutrho     = 360.0
  occupations = 'fixed'
  vdw_corr    = 'grimme-d3'
  dftd3_version = 4
/
&CELL
  cell_dynamics = 'bfgs'
  cell_dofree   = '2Dxy'         ! only x and y components of the cell move
  press_conv_thr = 0.5
/
ATOMIC_SPECIES
  Mo  95.95   Mo_ONCV_PBE-1.0.oncvpsp.upf
  S   32.06   s_pbe_v1.4.uspp.F.UPF
ATOMIC_POSITIONS angstrom
  Mo  0.000000  0.000000  10.000000
  S   1.595000  0.920873  11.560000
  S   1.595000  0.920873   8.440000
K_POINTS automatic
  12 12 1 0 0 0
```

赝势取 SSSP 1.3 PBE Efficiency：Mo 为 SG15 ONCV（35 / 140 Ry），S 为 GBRV 超软（35 / 280 Ry），这里统一提高到 45 / 360 Ry。单层 MoS₂ 本身不需要色散修正（只有一层），加上 D3 是为了与后续双层/异质结计算保持一致的能量面。

```bash
mpirun -np 4 pw.x -nk 4 -in mos2_mono.vc-relax.in > mos2_mono.vc-relax.out
grep -e 'bfgs converged' -e 'P=' mos2_mono.vc-relax.out | tail -3
grep -A4 'Begin final coordinates' mos2_mono.vc-relax.out
```

检查：最终 `CELL_PARAMETERS` 的第三行（c 轴）应与输入相同；面内应力分量收敛到 `press_conv_thr` 以内。在 2D 体系中报告的"压强"含有真空，数值本身没有物理意义，只看面内分量。

### 4.3 2D 库仑截断：`assume_isolated = '2D'`

`INPUT_PW`：对 xy 面内周期的结构截断 z 方向的库仑相互作用（Sohier, Calandra, Mauri 2017），总能、力、应力在 2D 框架下计算；在其基础上的线性响应（ph.x）也自动使用 2D 框架。要求：

- c 轴长度应大于 2D 材料厚度（含电子）的 2 倍；厚度可用对应层状体相的层间距估计，或"原子厚度 + 10 bohr"；程序另有 20 bohr 的下限；
- vc-relax 时只有面内应力有意义，应配合 `cell_dofree = '2Dxy'`；
- 与 ESM 一样，要求 c 轴沿 z、垂直于 xy 面。

示例 `bilayer_graphene.scf.in` 在 AB 双层石墨烯上同时使用 `assume_isolated='2D'` 与 D3(BJ)。中性、对称体系用不用 2D 截断差别不大；它的价值在于带电/门控体系与 2D 声子（避免 LO-TO 型的层间伪相互作用，见 [第12章](12-phonon-dfpt.md)）。

> [!TIP]
> `assume_isolated` 的其他选项：`'makov-payne'`（只限立方胞）、`'martyna-tuckerman'`（孤立分子，[第2章](02-input-anatomy.md)）、`'esm'`（§5.4）。注意 hp.x 不支持 `'2D'` 截断（[第8章](08-dftu-hp.md)）。

---

## 5. Slab 模型：Al(111)

### 5.1 构建

Al(111) 的 1×1 表面胞：面内 $a_s=a_0/\sqrt2$，层间距 $d=a_0/\sqrt3$。以 $a_0=4.04$ Å 为例：$a_s=2.8567$ Å，$d=2.3325$ Å，ABC 堆垛。示例取 6 层、c = 30 Å（真空约 18 Å），底两层固定。

建模工具：ASE `ase.build.fcc111('Al', size=(1,1,6), vacuum=..., a=...)` 或 pymatgen `SlabGenerator`，再用 `ase.io.write(..., format='espresso-in')` 输出 pw.x 输入，见 [Common 第4章](../../Common/tutorials/04-ase-pymatgen.md)；建模原则（终止面、极性、对称 slab）见 [VASP 第9章 §3](../../VASP/tutorials/09-surface-workfunction.md)。

> [!IMPORTANT]
> $a_0$ 必须来自**同一套参数**（泛函、赝势、截断、smearing）下的体相优化（[第5章](05-relax-eos.md)），否则 slab 存在人为应变。示例中的 4.04 Å 只是占位值。

### 5.2 输入（`al111_slab.relax.in`）

```fortran
&CONTROL
  calculation   = 'relax'
  prefix        = 'al111'
  outdir        = './tmp'
  pseudo_dir    = './pseudo'
  tprnfor       = .true.
  forc_conv_thr = 1.0d-4
  tefield       = .true.         ! required for the dipole correction
  dipfield      = .true.         ! Bengtsson dipole correction
/
&SYSTEM
  ibrav       = 0
  nat         = 6
  ntyp        = 1
  ecutwfc     = 40.0
  ecutrho     = 320.0
  occupations = 'smearing'
  smearing    = 'mv'
  degauss     = 0.02
  edir        = 3                ! correction along the third lattice vector (z)
  emaxpos     = 0.75             ! sawtooth maximum: inside the vacuum
  eopreg      = 0.05             ! width of the decreasing region (in vacuum)
  eamp        = 0.0              ! no external field, dipole correction only
/
&ELECTRONS
  conv_thr    = 1.0d-9
  mixing_beta = 0.3
  mixing_mode = 'local-TF'       ! helps for slabs (inhomogeneous systems)
/
&IONS
  ion_dynamics = 'bfgs'
/
ATOMIC_SPECIES
  Al  26.9815  Al.pbe-n-kjpaw_psl.1.0.0.UPF
CELL_PARAMETERS angstrom
  2.856711  0.000000   0.000000
  1.428356  2.473985   0.000000
  0.000000  0.000000  30.000000
ATOMIC_POSITIONS angstrom
  Al  0.000000  0.000000   2.000000  0 0 0
  Al  1.428356  0.824662   4.332495  0 0 0
  Al  2.856711  1.649323   6.664990
  Al  0.000000  0.000000   8.997485
  Al  1.428356  0.824662  11.329980
  Al  2.856711  1.649323  13.662475
K_POINTS automatic
  16 16 1 0 0 0
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `if_pos`（坐标后三个整数） | `0 0 0` 固定 | 默认 1（可动）；对应力分量乘以 0/1，用于结构优化与 MD 中固定原子或分量 |
| `mixing_mode` | `'local-TF'` | `INPUT_PW`：适用于高度不均匀体系 |
| `mixing_beta` | 0.2–0.3 | slab 金属的电荷振荡（charge sloshing）比体相严重 |
| `occupations`/`degauss` | `'smearing'` / 0.01–0.02 Ry | 金属 slab |
| K 网格 | 面内与体相相同密度，面外 1 | 本例 16×16×1 |

`ATOMIC_POSITIONS` 中固定原子的写法同样适用于吸附计算（固定底层、弛豫顶层与吸附物）。

### 5.3 运行与检查

```bash
mpirun -np 8 pw.x -nk 8 -in al111_slab.relax.in > al111_slab.relax.out
grep -e 'bfgs converged' -e 'Total force' al111_slab.relax.out | tail -3
grep -A12 'Begin final coordinates' al111_slab.relax.out
```

判据：`bfgs converged`；顶层相对体相层间距的变化在百分之几以内（Al(111) 表面弛豫很小）；固定原子的力不为零属正常，`Total force` 只统计可动原子。

### 5.4 ESM（简介）

`assume_isolated = 'esm'`（Otani–Sugino 2006）把超胞沿 z 嵌入有效屏蔽介质，适合带电或极化的 slab：

| `esm_bc` | 边界 |
|---|---|
| `'pbc'`（默认） | 普通周期（不启用 ESM） |
| `'bc1'` | 真空–slab–真空 |
| `'bc2'` | 金属–slab–金属（可加 `esm_efield`） |
| `'bc3'` | 真空–slab–金属 |

要求：c 轴沿 z 且垂直于 xy 面，**slab 居中于 z = 0**；vc-relax 时用 `cell_dofree='2Dxy'` 等不动 c 轴的设置；平面平均电荷与静电势写到 `prefix.esm1`。QE 源码 `PW/examples/ESM_example` 演示了 H₂O 分子与 Al(001) slab 的 ESM 计算，以及恒电势的 FCP、GCSCF 方法。

---

## 6. 偶极修正详解

### 6.1 关键词（`INPUT_PW`）

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `tefield` | `.true.` | 在离子势上加锯齿势；默认 `.false.` |
| `dipfield` | `.true.` | 仅当 `tefield=.true.` 时有效：加偶极修正（Bengtsson PRB 59, 12301）。**只能用于 slab，锯齿的不连续处必须落在真空中** |
| `edir` | 3 | 锯齿方向平行于倒格矢 bg(:,edir) |
| `emaxpos` | 真空中某点，默认 0.5 | 锯齿最大值位置（晶体坐标，0–1） |
| `eopreg` | 0.05–0.1，默认 0.1 | 锯齿下降区宽度（晶体坐标） |
| `eamp` | 0.0（只做偶极修正） | 外场幅度，Hartree 原子单位（1 a.u. = 51.422×10¹⁰ V/m），默认 0.001 |

`INPUT_PW` 对锯齿形状的描述：势在 (`emaxpos`+`eopreg`−1) 到 `emaxpos` 区间以斜率上升，然后在 `emaxpos` 到 `emaxpos`+`eopreg` 之间下降——**这个下降区必须完全在真空中，否则会产生非物理的力**。

> [!WARNING]
> `eamp` 默认是 0.001 a.u.（不是 0）。只想做偶极修正时必须显式写 `eamp = 0.0`，否则相当于额外加了约 0.05 V/Å（0.001 × 51.42 V/Å）的外电场。

### 6.2 如何放置 `emaxpos` / `eopreg`

示例中 slab 占据 z = 2.0–13.66 Å，即晶体坐标 0.067–0.455（c = 30 Å），真空为 0.455–1.067（跨越周期边界）。取 `emaxpos = 0.75`、`eopreg = 0.05`：下降区为 0.75–0.80（22.5–24.0 Å），距两侧表面都在 8 Å 以上。一般原则：

- 下降区放在**真空中部**，与两侧表面至少 3–5 Å；
- slab 不要跨越晶胞边界（把 slab 整体平移到胞内一侧更方便放置）；
- 改变真空厚度后要重新计算 `emaxpos`。

### 6.3 输出

偶极修正打开后，pw.x 每步打印（`PW/src/add_efield.f90`；`Elec. dipole`/`Ion. dipole` 两行只在 `verbosity='high'` 时出现，`E field amplitude` 一行只在 `eamp` ≠ 0 时出现）：

```bash
grep -A8 'Adding external electric field' al111_slab.scf.out | tail -9
```

（示意输出，数值仅供参考）

```text
     Adding external electric field

     Computed dipole along edir(3) : 
        Elec. dipole          x.xxxx Ry au,          x.xxxx Debye
        Ion. dipole           x.xxxx Ry au,          x.xxxx Debye
        Dipole                x.xxxx Ry au,          x.xxxx Debye
        Dipole field          x.xxxx Ry au, 
        Potential amp.        x.xxxx Ry
        Total length         xx.xxxx bohr
```

对称的 Al(111) slab，`Dipole` 应接近 0；单侧吸附后它显著非零，两侧真空能级出现台阶。

另一个相关选项 `gate = .true.`（`INPUT_PW`）：带电胞（`tot_charge` ≠ 0）中用一个位于 `zgate` 的带电平板代替均匀背景电荷，模拟场效应掺杂（Brumme 等 2014）；对 z 不对称的体系需同时开偶极修正。

---

## 7. 功函数：pp.x + average.x

### 7.1 流程

```text
pw.x scf (al111_slab.scf.in, tefield + dipfield)
   -> pp.x   plot_num = 11  (V_bare + V_H; dipole sawtooth added)  -> al111.v_bare_h
   -> average.x  (planar + macroscopic average along z)            -> avg.dat
   -> Phi = V_vac - E_F
```

`PP/src/punch_plot.f90` 中 `plot_num = 11` 在 `tefield` 与 `dipfield` 同时为真时会把偶极修正锯齿势加进去，因此输出的势在两侧真空中都是平的；`plot_num = 12` 则只输出锯齿势本身，可用于检查其位置。

### 7.2 pp.x 输入（`pp_v.in`）

```fortran
&INPUTPP
  prefix   = 'al111'
  outdir   = './tmp'
  filplot  = 'al111.v_bare_h'
  plot_num = 11
/
&PLOT
  iflag         = 3
  output_format = 6              ! Gaussian cube, for qe_cube_diff.py / VESTA
  fileout       = 'al111.v_bare_h.cube'
/
```

| 关键词 | 值 | 说明（`INPUT_PP`） |
|---|---|---|
| `plot_num` | 11 | V_bare + V_H；1 为总势（含 V_xc），真空中二者趋同，但 11 不受 V_xc 数值噪声影响 |
| `filplot` | 中间文件 | average.x 读取的就是这个文件 |
| `iflag` | 3 | 3D |
| `output_format` | 6 | cube（5 为 XCrySDen xsf） |

pp.x 输出的势单位为 Ry，量纲是能量（e·V，不是 V）。

### 7.3 average.x 输入（`average.in`）

average.x 没有 namelist，按行读取（`PP/src/average.f90`）：

```text
1                 ! nfile: number of files
al111.v_bare_h    ! file name (filplot of pp.x)
1.0               ! weight of this file
1000              ! npt: points for interpolation along idir
3                 ! idir: average over planes perpendicular to direction 3
4.4078            ! awin: window (bohr) of the macroscopic average
```

> [!NOTE]
> 上面每行的 `! ...` 注释只为说明；实际文件 `average.in` 中**不写注释**（文件名一行按整行读取，注释会被当成文件名的一部分）。宏观平均窗口 `awin` 取一个层间距（Al(111)：2.3325 Å = 4.4078 bohr），可以消除体内的层振荡。

average.x 必须**单进程**运行（多进程报 `run on a single processor`），输出 `avg.dat`，三列分别为 z（bohr）、平面平均 p(z)（Ry）、宏观平均 m(z)（Ry）。

### 7.4 运行与读取

```bash
mpirun -np 4 pw.x -nk 4 -in al111_slab.scf.in > al111_slab.scf.out
mpirun -np 4 pp.x -in pp_v.in > pp_v.out
average.x < average.in > average.out
grep 'the Fermi energy' al111_slab.scf.out
python ../../scripts/qe_workfunction.py -h        # helper: V_vac, Phi, plot
```

`run_workfunction.sh` 也给出一个最简的 awk 估计：在两侧真空中部各取一点，把 p(z)（Ry）乘 13.6057 换成 eV 再减 E_F。（示意输出，数值仅供参考）

```text
E_Fermi (eV) = x.xxxx
side 1 (frac z=0.60): V_vac =   x.xxxx eV   Phi =   4.xxxx eV
side 2 (frac z=0.93): V_vac =   x.xxxx eV   Phi =   4.xxxx eV
```

成功判据：

- p(z) 在两侧真空区各自平坦（起伏 < 0.01 eV）；若有斜率，检查偶极修正是否打开、`eamp` 是否为 0、下降区是否在真空中；
- 对称 slab 两侧 Φ 相同；不对称 slab 两侧差值与输出的偶极一致；
- Φ 随层数收敛（4、6、8 层差 < 0.05 eV）；Al(111) 实验功函数约 4.2–4.3 eV，PBE 结果通常在其附近。

---

## 8. 表面能计算

Al(111) 是对称 slab 的典型例子，但示例中底两层固定、只弛豫顶部，所以严格地应按 §2.2 的拆分公式计算。最简单的版本（未弛豫对称 slab）：

```bash
# E_slab (Ry) of the unrelaxed 6-layer slab: run al111_slab.scf.in with the
# initial coordinates; E_bulk per atom from the oriented 3-atom bulk cell
mpirun -np 4 pw.x -nk 4 -in al_bulk_oriented.scf.in > al_bulk_oriented.scf.out
grep '^!' al111_slab.scf.out al_bulk_oriented.scf.out
```

$$
\gamma=\frac{E_{\rm slab}-6\cdot E_{\rm bulk}^{(3)}/3}{2A},\qquad A=a_s^2\sin 60^\circ\approx 7.07\ \text{Å}^2
$$

取向一致的体相胞（`al_bulk_oriented.scf.in`：同样的面内晶胞与 16×16 面内 k 点，c = 3 层）比 fcc 原胞更好地抵消 k 点误差。另一种常用做法是用不同层数 N 的 slab 能量对 N 线性拟合，斜率即 $\varepsilon_{\rm bulk}$（避免体相参考不一致导致的 γ 随 N 发散）。

（示意输出，数值仅供参考）

| 层数 N | E_slab (Ry) | γ (J/m²) |
|---|---|---|
| 4 | −xx.xxxx | 0.xx |
| 6 | −xx.xxxx | 0.xx |
| 8 | −xx.xxxx | 0.xx |

单位：γ[J/m²] = γ[Ry/Å²] × 13.6057 × 16.02。

---

## 参数选择与收敛建议

| 项目 | 建议 |
|---|---|
| 色散方法 | 层状/吸附默认 D3(BJ)（`dftd3_version = 4`）；与非局域泛函交叉验证几何 |
| vc-relax | 截断能比 SSSP 推荐值高 10–20%；末尾重做 SCF 的压强与最后一步一致 |
| 2D 真空 | ≥ 15 Å，带电或做功函数时 ≥ 20 Å |
| 2D k 点 | `N N 1`，面内密度与体相一致 |
| slab 层数 | 金属 ≥ 5–7 层；以表面能、功函数、吸附能收敛为准 |
| 固定层 | 底部 1/3–1/2 固定在体相位置；或对称 slab 两侧同时弛豫 |
| 混合 | `mixing_mode='local-TF'`，`mixing_beta` 0.2–0.3 |
| 偶极修正 | 不对称体系必须开；`eamp=0`；下降区在真空中部 |
| average.x | `npt` ≥ FFT 网格数；`awin` 取一个层间距 |
| 能量比较 | 不同色散方法/泛函之间不比较总能；吸附能的各项必须用同一方法 |

---

## 常见问题与报错

| 报错 / 现象 | 原因 | 解决 |
|---|---|---|
| `WARNING: unknown vdw correction (vdw_corr): ... No vdw correction used.` | `vdw_corr` 拼错 | 用 §3.1 表中的合法取值；grep `DFT-D3 Dispersion` 确认 |
| `Tkatchenko-Scheffler not implemented with USPP/PAW` | `vdw_corr='ts'` 或 `'mbd'` 配了超软/PAW 赝势 | 换 NC 赝势，或改用 D3/XDM/非局域泛函 |
| `must choose a unique vdW correction!` | 同时打开了两种色散修正（如 `vdw_corr` 加旧开关 `london=.true.`） | 只保留一个 `vdw_corr` |
| 加了 D3 但输出没有 `DFT-D3 Dispersion` | 同上，或把 vdW-DF 写进了 `vdw_corr` | 非局域泛函用 `input_dft` |
| `XDM not implemented for noncolinear case` | XDM + 非共线 | 换 D3/TS 或改用共线 |
| 石墨 vc-relax 时 c 持续增大 | 没有色散修正 | 加 `vdw_corr` 或非局域泛函 |
| 2D vc-relax 后真空被压缩 | 没有设 `cell_dofree='2Dxy'` | 加 `cell_dofree='2Dxy'` |
| `dipole correction is not active if tefield = .false.` | 设了 `dipfield` 但没开 `tefield` | 同时设置 `tefield=.true.` |
| `gate cannot be used with tefield if dipole correction is not active` | `gate` + `tefield` 但缺 `dipfield` | 加 `dipfield=.true.` |
| 真空势有斜率 | 没开偶极修正，或 `eamp` 用了默认 0.001 | `tefield`+`dipfield`，`eamp = 0.0` |
| 开偶极修正后受力异常、SCF 振荡 | 锯齿下降区与原子或电子密度重叠 | 移动 `emaxpos`，减小 `eopreg`，增大真空 |
| `run on a single processor`（average.x） | 用 mpirun 启动了 average.x | 串行运行 |
| `convergence NOT achieved after 100 iterations: stopping`（slab） | 金属 slab 电荷振荡 | `mixing_mode='local-TF'`，减小 `mixing_beta`，增大 `electron_maxstep` |
| 表面能随层数单调漂移 | 体相参考与 slab 的 k 点/截断不一致 | 用取向一致的体相胞或线性拟合法 |
| 功函数随真空厚度变化 | 真空不足、平台不平 | 真空 ≥ 15–20 Å，检查 p(z) 平台 |

---

## 练习

1. 用 `run_vdw_compare.sh` 比较默认 6 种方法（另可在换成 NC 赝势后加 `ts`）下石墨的 c，与实验值 6.71 Å 列表对比；再固定 c 计算层间结合能曲线 E(d)（层间距 3.0–5.0 Å），注意只在同一方法内比较。
2. 对 D3 分别设 `dftd3_threebody = .true./.false.`，比较石墨 c 的差别。
3. 对单层 MoS₂ 测试真空 12、15、20、25 Å 对带隙（[第6章](06-band-dos-pdos.md) 流程）的影响。
4. Al(111)：比较 4、6、8 层 slab 的功函数与表面能，确定收敛层数；再把 `emaxpos` 移到距表面仅 1 Å 处，观察受力和功函数的变化。
5. 在 Al(111) 顶层放一个 Na 原子（2×2 超胞），比较开/关偶极修正时两侧真空的电势分布与功函数变化 ΔΦ。
6. 用 `python ../../scripts/qe_cube_diff.py -h` 查看用法，对 `al111.v_bare_h.cube` 做平面平均，与 average.x 的 `avg.dat` 比较。

---

## 参考资料

**官方文档**

- [INPUT_PW](https://www.quantum-espresso.org/Doc/INPUT_PW.html)：`vdw_corr`、`dftd3_version`、`dftd3_threebody`、`london_*`、`ts_vdw_*`、`xdm_*`、`input_dft`、`cell_dofree`、`assume_isolated`、`esm_*`、`tefield`、`dipfield`、`edir`、`emaxpos`、`eopreg`、`eamp`、`gate`、`mixing_mode`
- [INPUT_PP](https://www.quantum-espresso.org/Doc/INPUT_PP.html)：`plot_num`（11、12）、`iflag`、`output_format`
- average.x 的输入说明：源码 `PP/src/average.f90` 文件头注释（<https://gitlab.com/QEF/q-e>）
- QE 源码中的官方示例：`PW/examples/vdwDF_example`、`PW/examples/ESM_example`、`PW/examples/example02`（Al(001) slab 弛豫）、`test-suite/pw_vdw`、`test-suite/pw_dipole`
- [PWscf User's Guide](https://www.quantum-espresso.org/Doc/pw_user_guide/)
- 赝势：[SSSP 1.3 Efficiency](https://www.materialscloud.org/discover/sssp/table/efficiency)

**课程与教程（只链接）**

- [nguyen-group/QE-SSP](https://github.com/nguyen-group/QE-SSP)：《QE Course for Solid-State Physics》配套算例，含石墨烯、MoS₂、双层石墨烯，与本章 2D 部分高度相关
- [jkitchin/dft-book-espresso](https://github.com/jkitchin/dft-book-espresso)：CMU dft-book 的 QE 版，用 ASE 驱动表面、吸附与功函数计算
- [Imperial College MSE404-MM](https://imperialcollegelondon.github.io/MSE404-MM/)：QE 实验课，基础设置可作补充
- [Cornell PARADIM NSF 2D workshop 2018 数据](https://github.com/paradimdata/nsf_2d_workshop_2018)：PARADIM/NIST 的 2D 材料数据培训文件（非 QE 教程，作背景阅读）
- [pranabdas/espresso](https://pranabdas.github.io/espresso/)：入门到进阶的 QE 教程网站

**相关仓库（只链接）**

- [QEF/q-e](https://github.com/QEF/q-e)：官方镜像（开发在 GitLab）
- [pierrehirel/atomsk](https://github.com/pierrehirel/atomsk)：建超胞与格式转换，可输出 pw.x 输入
- [materialsproject/pymatgen](https://github.com/materialsproject/pymatgen)：`SlabGenerator` 切表面
- [aiidalab/aiidalab-qe](https://github.com/aiidalab/aiidalab-qe)：浏览器中建模与运行 QE 工作流

**关键文献**

- S. Grimme, *J. Comput. Chem.* 27, 1787 (2006)（DFT-D2）. DOI: [10.1002/jcc.20495](https://doi.org/10.1002/jcc.20495)
- S. Grimme et al., *J. Chem. Phys.* 132, 154104 (2010)（DFT-D3）. DOI: [10.1063/1.3382344](https://doi.org/10.1063/1.3382344)
- S. Grimme, S. Ehrlich, L. Goerigk, *J. Comput. Chem.* 32, 1456 (2011)（BJ 阻尼）. DOI: [10.1002/jcc.21759](https://doi.org/10.1002/jcc.21759)
- A. Tkatchenko, M. Scheffler, *Phys. Rev. Lett.* 102, 073005 (2009)（TS）. DOI: [10.1103/PhysRevLett.102.073005](https://doi.org/10.1103/PhysRevLett.102.073005)
- A. Ambrosetti et al., *J. Chem. Phys.* 140, 18A508 (2014)（MBD）. DOI: [10.1063/1.4865104](https://doi.org/10.1063/1.4865104)
- M. Dion et al., *Phys. Rev. Lett.* 92, 246401 (2004)（vdW-DF）. DOI: [10.1103/PhysRevLett.92.246401](https://doi.org/10.1103/PhysRevLett.92.246401)
- K. Lee et al., *Phys. Rev. B* 82, 081101 (2010)（vdW-DF2）. DOI: [10.1103/PhysRevB.82.081101](https://doi.org/10.1103/PhysRevB.82.081101)
- R. Sabatini, T. Gorni, S. de Gironcoli, *Phys. Rev. B* 87, 041108 (2013)（rVV10）. DOI: [10.1103/PhysRevB.87.041108](https://doi.org/10.1103/PhysRevB.87.041108)
- K. Berland, P. Hyldgaard, *Phys. Rev. B* 89, 035412 (2014)（vdW-DF-cx）. DOI: [10.1103/PhysRevB.89.035412](https://doi.org/10.1103/PhysRevB.89.035412)
- I. Hamada, *Phys. Rev. B* 89, 121103 (2014)（rev-vdW-DF2）. DOI: [10.1103/PhysRevB.89.121103](https://doi.org/10.1103/PhysRevB.89.121103)
- J. Neugebauer, M. Scheffler, *Phys. Rev. B* 46, 16067 (1992)（偶极修正）. DOI: [10.1103/PhysRevB.46.16067](https://doi.org/10.1103/PhysRevB.46.16067)
- L. Bengtsson, *Phys. Rev. B* 59, 12301 (1999)（偶极修正）. DOI: [10.1103/PhysRevB.59.12301](https://doi.org/10.1103/PhysRevB.59.12301)
- M. Otani, O. Sugino, *Phys. Rev. B* 73, 115407 (2006)（ESM）. DOI: [10.1103/PhysRevB.73.115407](https://doi.org/10.1103/PhysRevB.73.115407)
- T. Sohier, M. Calandra, F. Mauri, *Phys. Rev. B* 96, 075448 (2017)（2D 库仑截断）. DOI: [10.1103/PhysRevB.96.075448](https://doi.org/10.1103/PhysRevB.96.075448)
- T. Brumme, M. Calandra, F. Mauri, *Phys. Rev. B* 89, 245406 (2014)（gate）. DOI: [10.1103/PhysRevB.89.245406](https://doi.org/10.1103/PhysRevB.89.245406)

---

[← 上一章](08-dftu-hp.md) | [返回目录](../README.md) | [下一章 →](10-charge-elf-bader.md)
