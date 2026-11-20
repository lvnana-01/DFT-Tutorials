# 第4章 收敛测试：ecutwfc/ecutrho、k 点与 smearing（Convergence Tests in QE: Cutoffs, k-points & Smearing）

本章讲 Quantum ESPRESSO（pw.x）中最基础也最容易被省略的一步：收敛测试（convergence test）。我们分别扫描波函数截断能 `ecutwfc`、电荷密度截断能 `ecutrho`（dual）、Monkhorst-Pack k 点网格（`K_POINTS automatic`）以及金属的展宽方式 `smearing`/`degauss`，并用 `conv_thr` 控制 SCF 噪声。示例以半导体 Si（超软赝势）和金属 Al（PAW）对比，给出 bash 批量脚本、要 grep 的输出行和"收敛"的判据（ΔE < 1 meV/atom、力、应力）。

> [!NOTE]
> 本章全部输入与脚本在示例目录 [04-convergence](../examples/04-convergence/README.md)。赝势来自 SSSP 1.3 PBE Efficiency（下载方法见 [第3章](03-pseudopotentials.md)），仓库不提交 UPF 文件。文中所有"示意输出"只演示格式与趋势，**数值仅供参考，不是真实计算结果**。

---

## 1. 本章目标与适用场景

| 目标 | 你将学会 |
|---|---|
| 截断能 | 固定 dual 扫描 `ecutwfc`，再固定 `ecutwfc` 扫描 `ecutrho`；理解 NC 与 US/PAW 的差别 |
| k 点 | 用 `K_POINTS automatic N N N 0 0 0` 扫描网格，理解偏移（offset）与 k 点密度–晶格常数的反比关系 |
| smearing | 在 `gaussian`/`mv`/`mp`/`fd` 中选择，把 `degauss` 与 k 网格一起测试；绝缘体用 `occupations='fixed'` |
| SCF 精度 | 理解 `conv_thr` 是**广延量**，按原子数缩放；读懂 `estimated scf accuracy` |
| 判据 | 能量差 < 1 meV/atom、力与应力的收敛，比较能量差而不是绝对能 |
| 自动化 | bash 循环 + `sed` 替换占位符 + `grep '^!'` 汇总 |

适用场景：**每一个新体系、新赝势组合、新泛函**都要做。测试得到的参数会被后续所有章节沿用：结构优化（[第5章](05-relax-eos.md)）、能带与 DOS（[第6章](06-band-dos-pdos.md)）、磁性（[第7章](07-magnetism-soc.md)）、声子（[第12章](12-phonon-dfpt.md)）等。

与 VASP 用户的对应关系（概念层面，数值不能直接照搬）：

| VASP | QE（pw.x） | 说明 |
|---|---|---|
| `ENCUT`（eV） | `ecutwfc`（Ry） | 1 Ry = 13.6057 eV；两者的赝势不同，数值不可直接换算后套用 |
| （PREC 决定的 FFT 网格） | `ecutrho`（Ry） | QE 显式给出密度截断；默认 4×`ecutwfc`（NC 用默认即可），US 常需 8–12×，PAW 必须测试 |
| `KPOINTS` Gamma / Monkhorst-Pack | `K_POINTS automatic` + 偏移 `0 0 0` / `1 1 1` | 偏移为 0 时网格包含 Γ 点 |
| `ISMEAR=0` / `1` / `-1` | `smearing='gaussian'` / `'mp'` / `'fd'` | QE 另有 Marzari-Vanderbilt 冷展宽 `'mv'` |
| `ISMEAR=-5` | `occupations='tetrahedra'`（或 `'tetrahedra_opt'`） | 适合 DOS，不适合力/优化 |
| `SIGMA`（eV） | `degauss`（Ry） | 0.01 Ry ≈ 0.136 eV |
| `EDIFF`（eV） | `conv_thr`（Ry） | 两者都是总量（不是每原子） |

VASP 的同类测试见 [VASP 第3章](../../VASP/tutorials/03-convergence.md)。

---

## 2. 原理速览

### 2.1 两个截断能：`ecutwfc` 与 `ecutrho`

pw.x 用平面波展开 Kohn-Sham 轨道，只保留动能低于截断的平面波。Hartree 原子单位下平面波动能为 $\frac{1}{2}|\mathbf{k}+\mathbf{G}|^2$，pw.x 使用的 Rydberg 原子单位下为 $|\mathbf{k}+\mathbf{G}|^2$，所以条件就是 $|\mathbf{k}+\mathbf{G}|^2 < E_\text{cut}$（单位 Ry）：

- `ecutwfc`：波函数截断能（Ry），**必填**；
- `ecutrho`：电荷密度与势的截断能（Ry），默认 `4 * ecutwfc`。

为什么是 4 倍？密度 $n(\mathbf r)=\sum|\psi|^2$ 是波函数的平方，其傅里叶分量最高到波函数的 2 倍波矢，即 4 倍动能。官方 INPUT_PW 文档对 `ecutrho` 的说明可以概括为：

| 赝势类型 | `ecutrho` 建议 |
|---|---|
| 模守恒（NC，如 PseudoDojo、SG15） | 用默认 4×；可略降，但会给力与应力引入噪声 |
| 超软（US，如 `*rrkjus*`） | 常需 8–12× |
| PAW（如 `*kjpaw*`） | 很多情况下 4× 就够，但取决于增广电荷形状，**必须测试** |
| 梯度修正泛函 + 真空层、无 NLCC 的赝势 | 往往需要更高的 `ecutrho` |

我们把比值 `ecutrho/ecutwfc` 称为 **dual**。SSSP 表格给出的推荐值（如 Si 30/240 Ry）就是 dual = 8。

### 2.2 变分性与"比较能量差"

平面波基组越大，总能越低（对给定结构单调下降）。但我们真正关心的通常是**能量差**：形成能、吸附能、两个相的能量差、EOS 曲线的形状。结构相近的构型，截断误差大部分会相互抵消，所以能量差比绝对总能收敛得快。反过来，若两个构型差别很大（孤立原子 vs 体相），抵消就差，需要更高的截断。

> [!IMPORTANT]
> 只要两次计算的能量要相减，就必须使用**相同的赝势、`ecutwfc`、`ecutrho`、泛函**。不同截断下的总能不能互相比较。

### 2.3 k 点积分

布里渊区积分离散为加权求和。`K_POINTS automatic` 的语法是：

```text
K_POINTS automatic
  nk1 nk2 nk3  sk1 sk2 sk3
```

`nk1..3` 与 Monkhorst-Pack 网格含义相同；`sk1..3` 只能是 0（不偏移，网格包含 Γ）或 1（沿该方向偏移半个格点）。官方文档提示：四面体方法只适用于具有晶体完整对称性的网格，有些偏移网格可能不工作。

k 点密度与实空间晶格常数成**反比**：晶胞越大，倒格子越小，同样精度所需的 $N_i$ 越少。一个常用的经验是让 $N_i \times a_i$ 大致相同（例如对小原胞 $N_i a_i \approx 30$–$50$ Å 作为起点，具体以测试为准）。对表面/二维体系，真空方向取 1。

### 2.4 金属与 smearing

绝缘体被积函数光滑，k 点误差下降很快；金属在费米面处占据数从 1 跳到 0，直接求和收敛极慢。pw.x 的处理办法是 `occupations='smearing'`，再用 `smearing` 选择展宽函数、`degauss`（Ry，默认 0）给出宽度：

| `smearing` | 别名 | 特点 |
|---|---|---|
| `'gaussian'` | `'gauss'` | 默认值；简单，但总能对 `degauss` 依赖较强 |
| `'methfessel-paxton'` | `'m-p'`、`'mp'` | 一阶 MP，金属总能对展宽不敏感；占据数可为负 |
| `'marzari-vanderbilt'` | `'cold'`、`'m-v'`、`'mv'` | 冷展宽，金属结构优化与声子中常用 |
| `'fermi-dirac'` | `'f-d'`、`'fd'` | 物理意义是电子温度，收敛最慢 |

加了展宽之后，pw.x 最小化的是自由能 $F = E - TS$。输出中 `!    total energy` 那一行就是 $F$（输出里会写 `The total energy is F=E-TS`），同时打印：

```text
     smearing contrib. (-TS)   =      -0.000xxxxx Ry
     internal energy E=F+TS    =     -xx.xxxxxxxx Ry
```

（示意输出，数值仅供参考）

`-TS` 项的大小是判断 `degauss` 是否过大的直接指标。冷展宽和 MP 的设计目标就是让 $F$ 对展宽宽度的依赖变弱。

### 2.5 SCF 收敛阈值 `conv_thr`

`conv_thr`（默认 `1.D-6` Ry）是"估计能量误差 < conv_thr"的判据。官方文档特别强调它是**广延量**（像总能一样随体系大小增长）。所以：

- 2 原子的 Si 用 `1e-8`～`1e-10` 很便宜；
- 200 原子的超胞若仍用 `1e-10`，等于每原子要求高了 100 倍，浪费时间；
- 做收敛测试时 SCF 噪声必须远小于判据（1 meV/atom ≈ 7.3×10⁻⁵ Ry/atom），一般取 `conv_thr` ≤ 1e-8 × 原子数 量级即可。

每一步 SCF 都会打印：

```text
     estimated scf accuracy    <       0.0000xxxx Ry
```

（示意输出，数值仅供参考）

当它低于 `conv_thr` 时出现 `convergence has been achieved in   N iterations`。

### 2.6 "收敛"的标准

| 物理量 | 常用目标 | 备注 |
|---|---|---|
| 能量差（每原子） | < 1 meV/atom | 形成能、相对稳定性 |
| 力 | 变化 < 1e-3 Ry/Bohr（约 0.026 eV/Å）量级 | 与 `forc_conv_thr` 默认值同量级；声子要更严 |
| 应力 / 压强 | 变化 < 0.5 kbar 量级 | 与 `press_conv_thr` 默认值 0.5 kbar 对应 |
| 金属 `-TS` | 每原子 < 1 meV 量级 | 判断 `degauss` 是否过大 |
| 目标物理量 | 直接测 | 带隙、频率、吸附能、能垒 |

> [!TIP]
> "收敛"永远针对你关心的物理量。能量收敛了不代表应力收敛了——变胞优化和 EOS 要以**压强**的收敛为准（见第 4.4 节和 [第5章](05-relax-eos.md)）。

---

## 3. 测试策略总览

```text
1. 固定: 泛函、赝势、结构、occupations; conv_thr 足够小
2. 中等 k 网格 (Si: 8x8x8) 下, 固定 dual 扫描 ecutwfc     -> 选 ecutwfc
3. (US/PAW) 固定 ecutwfc, 扫描 ecutrho                    -> 选 ecutrho
4. 用选定截断扫描 k 网格                                  -> 选 k 网格
5. 金属: k 网格 x degauss 二维扫描                         -> 选 smearing/degauss/k
6. 用最终参数复核: 力、应力、(必要时) 目标物理量
```

几条原则：

- **一次只改一个参数**；
- 测试体系尽量接近生产体系（同元素、同赝势、同泛函、同自旋设置）；
- 从 SSSP 推荐截断出发：SSSP 的推荐值已经针对能量、力、应力、声子等做过系统测试（Prandini et al., 2018），但它是"对单质/氧化物的推荐"，你的体系和目标量仍需复核；
- 选值时留余量，向上取一档。

---

## 4. 实战一：Si 的截断能测试

### 4.1 模板输入

[`si_ecut_template.in`](../examples/04-convergence/si_ecut_template.in)：

```fortran
&CONTROL
! Si (diamond) template for the ecutwfc / ecutrho convergence scan.
! Placeholders: @ECUT@, @ECUTRHO@, @K@
  calculation = 'scf'
  prefix      = 'si'
  outdir      = './tmp'
  pseudo_dir  = './pseudo'
  tprnfor     = .true.
  tstress     = .true.
/
&SYSTEM
  ibrav     = 2
  celldm(1) = 10.26        ! lattice parameter a in Bohr (about 5.43 Angstrom)
  nat       = 2
  ntyp      = 1
  ecutwfc   = @ECUT@
  ecutrho   = @ECUTRHO@
  occupations = 'fixed'    ! Si is an insulator: no smearing needed
/
&ELECTRONS
  conv_thr    = 1.0d-10    ! 2 atoms; much tighter than the 1 meV/atom target
  mixing_beta = 0.7
/
ATOMIC_SPECIES
  Si  28.0855  Si.pbe-n-rrkjus_psl.1.0.0.UPF
ATOMIC_POSITIONS crystal
  Si  0.00  0.00  0.00
  Si  0.25  0.25  0.25
K_POINTS automatic
  @K@ 0 0 0
```

| 关键词 | 取值 | 说明 |
|---|---|---|
| `calculation` | `'scf'` | 单点自洽 |
| `tprnfor` / `tstress` | `.true.` | scf 默认不算力和应力；收敛测试要打开，观察应力收敛 |
| `ibrav` / `celldm(1)` | `2` / `10.26` | fcc 晶格，`celldm(1)` 单位为 **Bohr** |
| `ecutwfc` / `ecutrho` | 占位符 | 由脚本替换；dual 默认为 8（超软赝势） |
| `occupations` | `'fixed'` | 有带隙的绝缘体/半导体 |
| `conv_thr` | `1.0d-10` | 2 原子体系，SCF 噪声远低于 1 meV/atom |
| `K_POINTS automatic` | `8 8 8 0 0 0` | 截断测试时固定在中等网格 |

> [!NOTE]
> Fortran namelist 中 `!` 之后是注释，本仓库把说明性注释放在 `&CONTROL` 之后、`/` 之前。卡片之间以 `!` 或 `#` 开头的行会被跳过；卡片**数据行内部**只建议用第一列为 `#` 的整行注释或行尾 `!` 注释，避免把注释写成数据。

### 4.2 生成与运行

```bash
cd QE/examples/04-convergence
mkdir -p pseudo        # 放入 Si.pbe-n-rrkjus_psl.1.0.0.UPF

# 只生成输入（默认不运行）
bash run_conv.sh ecut

# 运行（4 个 MPI 进程）
RUN=1 PW_CMD="mpirun -np 4 pw.x" bash run_conv.sh ecut
```

[`run_conv.sh`](../examples/04-convergence/run_conv.sh) 的核心就是一个 `sed` 循环：

```bash
for e in 25 30 35 40 45 50 60 70 80; do
    d=conv_ecut/ecut_$e; mkdir -p "$d"
    ln -sf ../../pseudo "$d/pseudo"     # pseudo_dir = './pseudo' inside $d
    sed -e "s/@ECUT@/$e/" -e "s/@ECUTRHO@/$((e*8))/" \
        -e "s/@K@/8 8 8/" si_ecut_template.in > "$d/pw.in"
    (cd "$d" && mpirun -np 4 pw.x -in pw.in > pw.out)
done
```

仓库里还有通用版本 `QE/scripts/qe_conv_test.sh`（占位符相同，支持 ecut/kpt/smear 三种模式）。注意它的 dual 默认是 4，超软/PAW 赝势要显式加 `-d 8`；在本示例目录下用 `bash ../../scripts/qe_conv_test.sh -h` 查看用法。

### 4.3 要 grep 的输出

```bash
cd conv_ecut/ecut_40
grep -e 'kinetic-energy cutoff' -e 'charge density cutoff' pw.out   # 确认截断被正确读入
grep -e 'number of k points' pw.out                                 # 不可约 k 点数
grep '^!' pw.out                                                    # 最终总能 (Ry)
grep -e 'estimated scf accuracy' pw.out | tail -3                   # SCF 收敛过程
grep -e 'convergence has been achieved' pw.out
grep -A3 'total   stress' pw.out                                    # 应力张量与 P=
grep 'Total force' pw.out                                           # 总力
grep 'JOB DONE' pw.out
```

批量汇总：

```bash
for d in conv_ecut/ecut_*; do
    printf "%s %s %s\n" "${d##*_}" \
        "$(grep '^!' $d/pw.out | awk '{print $5}')" \
        "$(grep 'P=' $d/pw.out | awk '{print $NF}')"
done | sort -n
```

`run_conv.sh` 会自动生成 `conv_ecut/summary.dat`，以**最高截断**为参考计算 ΔE（meV/atom）：

```text
# label  E_total(Ry)  dE_vs_last(meV/atom)  P(kbar)  n_scf  status
ecut_25        -15.8xxxxxxx       xx.xxx      xx.xx      8 done
ecut_30        -15.8xxxxxxx        x.xxx      xx.xx      8 done
...
ecut_60        -15.8xxxxxxx        0.xxx      -x.xx      8 done
ecut_70        -15.8xxxxxxx        0.0xx      -x.xx      8 done
ecut_80        -15.8xxxxxxx        0.000      -x.xx      8 done
```

（示意输出，数值仅供参考）

你也可以在示例目录下用 `python ../../scripts/qe_parse_out.py -h` 查看解析脚本的用法，它能一次性提取总能、应力、力和 SCF 历史。

### 4.4 如何解读

1. **能量**：ΔE 随 `ecutwfc` 增大单调下降（变分性）。若出现明显非单调，检查是否混用了不同 `ecutrho`、或 `conv_thr` 太松；
2. **压强**：压强通常比能量收敛得慢。如果之后要做 vc-relax 或 EOS，以压强变化 < 0.5 kbar 为准选截断；
3. **力**：理想金刚石 Si 原子受力为零（对称性），要测力的收敛可把第二个原子挪到 `0.27 0.25 0.25` 再扫描一次；
4. **对照 SSSP**：Si.pbe-n-rrkjus_psl.1.0.0.UPF 在 SSSP 1.3 Efficiency 中的提示值是 30/240 Ry。你的测试结果若明显高于它，通常是因为判据更严（如应力）或目标量不同，属正常现象。

### 4.5 第二步：扫描 `ecutrho`（US/PAW）

对超软/PAW 赝势，选定 `ecutwfc` 后再固定它扫 dual：

```bash
for dual in 4 6 8 10 12; do
    e=40; rho=$((e*dual))
    d=conv_rho/dual_$dual; mkdir -p $d; ln -sf ../../pseudo $d/pseudo
    sed -e "s/@ECUT@/$e/" -e "s/@ECUTRHO@/$rho/" -e "s/@K@/8 8 8/" \
        si_ecut_template.in > $d/pw.in
done
```

判断要点：dual 太小时能量和**应力**会出现明显噪声；对 NC 赝势这一步通常可以省略（保持默认 4×）。多元素体系的 `ecutwfc`/`ecutrho` 要取各元素推荐值中的**最大值**（例如 SSSP 1.3 Efficiency 中 O 的提示值是 50/400 Ry，含 O 的氧化物就至少从这里起测）。

---

## 5. 实战二：Si 的 k 点测试

用上一步选定的截断（脚本默认 `ECUT_FIX=40`，`DUAL=8`）扫描 N×N×N：

```bash
RUN=1 PW_CMD="mpirun -np 4 pw.x" ECUT_FIX=40 bash run_conv.sh kpt
cat conv_kpt/summary.dat
```

每个子目录中的 K_POINTS 形如：

```text
K_POINTS automatic
  8 8 8 0 0 0
```

要点：

- 看 `number of k points=` 一行：对称性把 8×8×8 = 512 个点约化为几十个不可约点。不可约点数决定计算量；
- Si 是半导体，k 收敛通常较快；
- **偶数网格是否偏移**：`0 0 0` 包含 Γ，`1 1 1` 为偏移网格。对 fcc 原胞，偏移网格有时能用更少的不可约点达到同样精度；六角晶格在面内（a、b 方向）偏移会破坏六重对称，最多只沿 c 方向偏移（如 `N N M 0 0 1`）。最稳妥的做法是 Γ 中心（`0 0 0`），在需要时再对比偏移网格；
- 之后的 DOS 计算（[第6章](06-band-dos-pdos.md)）需要比 SCF 密得多的网格，那是另一回事，不要和 SCF 的收敛混淆。

（示意输出，数值仅供参考）

```text
# label  E_total(Ry)  dE_vs_last(meV/atom)  P(kbar)  n_scf  status
k_2            -15.7xxxxxxx      xxx.xxx      xx.xx      7 done
k_4            -15.8xxxxxxx        x.xxx      -x.xx      7 done
k_6            -15.8xxxxxxx        0.xxx      -x.xx      7 done
...
k_16           -15.8xxxxxxx        0.000      -x.xx      7 done
```

---

## 6. 实战三：金属 Al 的 k 点 × smearing 测试

### 6.1 为什么必须二维扫描

对金属，k 点误差与展宽误差耦合：

- `degauss` 大 → 被积函数平滑 → 少量 k 点就"收敛"，但收敛到的是**偏离真实值**的结果；
- `degauss` 小 → 结果更接近 σ→0 极限，但需要很密的 k 网格才收敛。

所以要在几个 `degauss` 下分别扫 k 网格，画出"E vs k"的一组曲线，选一个在合理 k 网格下已经收敛、且 `-TS` 足够小的组合。

### 6.2 模板

[`al_kpt_template.in`](../examples/04-convergence/al_kpt_template.in)：

```fortran
&CONTROL
! fcc Al template for the k-point x smearing convergence scan.
! Placeholders: @K@, @DEGAUSS@, @SMEAR@
  calculation = 'scf'
  prefix      = 'al'
  outdir      = './tmp'
  pseudo_dir  = './pseudo'
  tprnfor     = .true.
  tstress     = .true.
/
&SYSTEM
  ibrav     = 2
  celldm(1) = 7.65         ! a in Bohr (about 4.05 Angstrom)
  nat       = 1
  ntyp      = 1
  ecutwfc   = 40.0
  ecutrho   = 320.0
  occupations = 'smearing' ! Al is a metal: smearing is mandatory
  smearing    = '@SMEAR@'
  degauss     = @DEGAUSS@
/
&ELECTRONS
  conv_thr    = 1.0d-10
  mixing_beta = 0.7
/
ATOMIC_SPECIES
  Al  26.9815  Al.pbe-n-kjpaw_psl.1.0.0.UPF
ATOMIC_POSITIONS crystal
  Al  0.00  0.00  0.00
K_POINTS automatic
  @K@ 0 0 0
```

| 关键词 | 取值 | 说明 |
|---|---|---|
| `occupations` | `'smearing'` | 金属必须用；若不写（按 `'fixed'` 处理）而价电子数为奇数（如 Al 的 3 个），pw.x 会给出提示 `the system is metallic, specify occupations`；即使电子数为偶数，金属用 `'fixed'` 也会 SCF 难收敛、结果错误 |
| `smearing` | `'mv'`（默认扫描） | 冷展宽；也可设 `SMEAR=mp` 或 `SMEAR=gaussian` 对比 |
| `degauss` | 0.005–0.04 Ry | 单位是 **Ry**，不是 eV |
| `ecutwfc`/`ecutrho` | 40/320 Ry | 比 SSSP 提示的 30/240 留了余量；应先按第 4 节测过 |

### 6.3 运行与汇总

```bash
RUN=1 PW_CMD="mpirun -np 4 pw.x" bash run_conv.sh alsmear
# 生成 conv_al_mv_0.005/ conv_al_mv_0.01/ conv_al_mv_0.02/ conv_al_mv_0.04/
for d in conv_al_mv_*; do echo "== $d"; cat $d/summary.dat; done

# 每个计算的 -TS 项
grep 'smearing contrib' conv_al_mv_0.02/k_*/pw.out
grep 'the Fermi energy is' conv_al_mv_0.02/k_16/pw.out
```

把四组数据画在一张图上（横轴 k 网格 N，纵轴每原子能量），典型的定性特征是：

- 大 `degauss` 的曲线很早变平，但平台位置与小 `degauss` 的平台不同；
- 小 `degauss` 的曲线在稀疏网格时振荡明显，需要更大的 N 才平稳；
- 冷展宽/MP 的平台对 `degauss` 的依赖比 `gaussian` 小。

（示意输出，数值仅供参考）

```text
== conv_al_mv_0.02
# label  E_total(Ry)  dE_vs_last(meV/atom)  P(kbar)  n_scf  status
k_6            -39.xxxxxxxx       xx.xxx      xx.xx      9 done
k_8            -39.xxxxxxxx        x.xxx      xx.xx      9 done
...
k_24           -39.xxxxxxxx        0.000      xx.xx      9 done
```

### 6.4 如何选

1. 选一个 `degauss`，使 `smearing contrib. (-TS)` 每原子在 1 meV 量级以下，且与更小 `degauss` 的平台差 < 1 meV/atom；
2. 在这个 `degauss` 下选 k 网格（ΔE < 1 meV/atom，压强变化 < 0.5 kbar）；
3. 常见经验：金属结构优化用 `mv` 或 `mp`，`degauss` 0.01–0.02 Ry 量级；但具体以测试为准；
4. 做声子（[第12章](12-phonon-dfpt.md)）时 smearing 设置必须与 scf 一致；做 DOS 时改用四面体（[第6章](06-band-dos-pdos.md)）。

> [!WARNING]
> MP 展宽的占据数可以为负、也可以大于 1，不要用于半导体/绝缘体。对有带隙的体系，用 `occupations='fixed'`；如果不确定是否有带隙，先用小 `degauss` 的 `gaussian` 或 `mv` 跑一次，看 `highest occupied, lowest unoccupied level` 或 `the Fermi energy is` 的输出。

---

## 7. `conv_thr`、力与应力的配套设置

| 关键词 | 默认值（INPUT_PW） | 收敛测试建议 | 说明 |
|---|---|---|---|
| `conv_thr` | `1.D-6` | 小原胞 1e-8～1e-10 | 广延量，按原子数放宽或收紧 |
| `electron_maxstep` | 100 | 100–200 | 难收敛体系可增大，但先查原因 |
| `mixing_beta` | 0.7 | 0.7（Si/Al） | 磁性、表面、大体系可降到 0.3 或更低 |
| `tprnfor` | scf 下不计算 | `.true.` | relax/md 会自动打开 |
| `tstress` | `.false.` | `.true.` | vc-relax/vc-md 会自动打开 |
| `occupations` | 文档未列默认值；代码内部按 `'fixed'` 处理 | 绝缘体 `'fixed'`，金属 `'smearing'` | 金属务必显式写出；DOS 用 `'tetrahedra_opt'` |
| `degauss` | `0.D0` Ry | 测试决定 | Ry 单位 |

另外一条在 relax 中才会遇到的约束（见 [第5章](05-relax-eos.md)）：`conv_thr` 应比离子步判据 `etot_conv_thr` 小 2–3 个数量级，否则力带有 SCF 噪声，BFGS 会反复打印 `uphill step: resetting bfgs history` 并来回振荡。

---

## 8. 大体系与特殊体系的收敛建议

- **超胞/表面**：k 点按原胞结果换算（$N_i \propto 1/a_i$），真空方向取 1；截断能沿用小原胞测试值；
- **分子/孤立体系**：`K_POINTS gamma` 可启用 Γ 点专用算法（内存和时间约减半），真空盒子大小也要测（见 [第2章](02-input-anatomy.md)）；
- **磁性体系**：收敛测试要在自旋极化（`nspin=2`）下做，磁矩也是需要收敛的量（[第7章](07-magnetism-soc.md)）；
- **DFT+U、杂化泛函**：截断能一般沿用 PBE 测试值；杂化泛函的 q 点网格（`nqx1..3`）另需测试（[第14章](14-hybrid-metagga.md)）；
- **声子/介电**：力常数对截断和 k 点更敏感，常需比能量收敛更高一档；
- **多元素**：取各元素推荐截断的最大值，再整体测试。

---

## 参数选择与收敛建议

| 参数 | 起点 | 判据 | 典型陷阱 |
|---|---|---|---|
| `ecutwfc` | SSSP/PseudoDojo 推荐值 | ΔE < 1 meV/atom；压强变化 < 0.5 kbar | 用 VASP 的 eV 数值直接填 Ry |
| `ecutrho` | NC 4×，US/PAW 8× | 应力、力不再抖动 | US 赝势沿用默认 4× |
| k 网格 | $N_i a_i \approx 30$–$50$ Å 起测 | ΔE < 1 meV/atom | 超胞仍用原胞网格（浪费）或反之（不足） |
| `smearing` | 金属 `'mv'` 或 `'mp'` | 与更小展宽的平台差 < 1 meV/atom | 半导体用 MP |
| `degauss` | 0.01–0.02 Ry | `-TS` 每原子 ≲ 1 meV | 把 eV 当 Ry |
| `conv_thr` | 1e-8（小原胞测试 1e-10） | 远小于判据 | 大体系仍用极小值导致过慢 |

推荐记录方式：把每个体系的最终参数（赝势文件名、`ecutwfc`、`ecutrho`、k 网格、smearing、`degauss`）写进一个 README 或 JSON，后续章节统一引用。

---

## 常见问题与报错

| 报错 / 现象 | 原因 | 解决办法 |
|---|---|---|
| `convergence NOT achieved after 100 iterations: stopping` | 金属用了 `'fixed'`；`mixing_beta` 过大；磁性初猜不当 | 金属改 `occupations='smearing'`；降低 `mixing_beta`（0.7→0.3）；增大 `electron_maxstep`；表面体系可试 `mixing_mode='local-TF'` |
| `Message from routine setup: the system is metallic, specify occupations` | 未设置 `occupations`（按 `'fixed'` 处理），而价电子数为奇数 | 显式写 `occupations='smearing'` 和 `degauss` |
| relax 中反复出现 `uphill step: resetting bfgs history` | `conv_thr` 相对 `etot_conv_thr`/`forc_conv_thr` 太松，力有噪声 | 收紧 `conv_thr`（比 `etot_conv_thr` 小 2–3 个数量级），见 [第5章](05-relax-eos.md) |
| `Error in routine readpp ... file ... not found` | `pseudo_dir` 路径或 UPF 文件名大小写不对 | 检查 `pseudo_dir` 与 `ATOMIC_SPECIES` 中的文件名，或设置 `ESPRESSO_PSEUDO` |
| `S matrix not positive definite` / `problems computing cholesky` | 原子过近、单位混用（Bohr/Å）、`ecutrho` 太低 | 检查结构与单位；提高 `ecutrho`；试 `diagonalization='cg'` |
| `too many bands are not converged` | 对角化不稳定，常与结构问题相关 | 检查结构；`diagonalization='cg'` |
| 能量随 `ecutwfc` 非单调变化 | 同时改了 `ecutrho`，或 `conv_thr` 太松 | 固定 dual 或固定 `ecutrho`，收紧 `conv_thr` |
| 各 k 网格能量剧烈振荡（金属） | `degauss` 太小或用了 `'fixed'` | 用 smearing 并做二维扫描 |
| 各 k 网格能量随 N 奇偶交替 | 偏移网格与对称性不匹配，或偶/奇网格是否含 Γ 的差别 | 统一用 Γ 中心网格，分开比较奇数/偶数 |
| `bad line in namelist &system: "..." (error could be in the previous line)` | 漏 `/`、关键词拼错、字符串没加引号、CRLF 换行 | 在示例目录下用 `python ../../scripts/qe_check_input.py -h` 查看静态检查脚本用法；`dos2unix` 转换换行 |

---

## 练习

1. 用 `run_conv.sh ecut` 扫描 Si，分别以"ΔE < 1 meV/atom"和"压强变化 < 0.5 kbar"为判据选 `ecutwfc`，比较两个结果并解释。
2. 把 Si 的 dual 从 8 改成 4（`DUAL=4`），重复第 1 题，观察应力的变化，说明为何超软赝势需要更高的 `ecutrho`。
3. 对 Al 分别用 `SMEAR=gaussian`、`SMEAR=mv`、`SMEAR=mp` 做二维扫描，画出 E(N) 曲线族，比较平台位置对 `degauss` 的依赖。
4. 删掉 Al 输入中的 `occupations`/`smearing`/`degauss` 三行运行一次，记录 pw.x 的提示与 SCF 行为。
5. 对 Si 比较 `8 8 8 0 0 0` 和 `8 8 8 1 1 1` 两种网格的不可约 k 点数与能量差。
6. 选一个含 O 的氧化物（如 MgO，参考 [第3章](03-pseudopotentials.md) 的示例），先查 SSSP 中各元素的推荐截断，再设计完整的收敛测试方案。

---

## 参考资料

**官方文档**

- [INPUT_PW（pw.x 全部关键词）](https://www.quantum-espresso.org/Doc/INPUT_PW.html) — `ecutwfc`、`ecutrho`、`occupations`、`smearing`、`degauss`、`conv_thr`、`K_POINTS`
- [PWscf User's Guide](https://www.quantum-espresso.org/Doc/pw_user_guide/)
- [SSSP 1.3 Efficiency 表](https://www.materialscloud.org/discover/sssp/table/efficiency) — 推荐截断与赝势文件名
- [Materials Cloud QE 输入生成器](https://tools.materialscloud.org/qeinputgenerator/) — 自动给出 SSSP 截断与 k 网格，可作为起点

**课程与教程（推荐，链接而非复制）**

- [Imperial College MSE404-MM Lab 3: Convergence Testing](https://imperialcollegelondon.github.io/MSE404-MM/labs/lab03/readme/) — 用 QE 系统讲解截断能与 k 点收敛，讲解清楚、持续维护
- [Imperial College MSE404-MM Lab 5: Metals and the Density of States](https://imperialcollegelondon.github.io/MSE404-MM/labs/lab05/readme/) — 金属 smearing 的直观讨论
- [Cornell PARADIM Computational Tutorials](https://www.paradim.org/toolbox/computation/tutorials) — "输入与收敛参数"部分适合理解概念（基于 QE 6.x 与 LDA 赝势，下载/编译部分已过时）
- [PARADIM Hands-on tutorial（2021 暑校衍生，PDF）](https://indico.truba.gov.tr/event/47/contributions/427/attachments/110/269/Hands-onTutorial_DFT_BP.pdf) — 截断能/k 点收敛与晶格常数
- [pranabdas/espresso: Convergence](https://pranabdas.github.io/espresso/hands-on/convergence) — 带 Python 绘图的收敛测试示例
- [Ghent Computational Materials Physics](https://www.compmatphys.org) — 零基础概念课，解释"为什么要收敛"

**相关仓库（仅链接）**

- [QEF/q-e](https://github.com/QEF/q-e) — 官方源码镜像，`PW/examples` 是权威输入示例
- [aiidateam/aiida-pseudo](https://github.com/aiidateam/aiida-pseudo) — 一键安装 SSSP/PseudoDojo 并给出推荐截断
- [nguyen-group/QE-SSP](https://github.com/nguyen-group/QE-SSP) — 教材配套算例，含收敛测试

**关键文献**

- G. Prandini, A. Marrazzo, I. E. Castelli, N. Mounet, N. Marzari, *npj Comput. Mater.* 4, 72 (2018)（SSSP）. DOI: [10.1038/s41524-018-0127-2](https://doi.org/10.1038/s41524-018-0127-2)
- H. J. Monkhorst, J. D. Pack, *Phys. Rev. B* 13, 5188 (1976). DOI: [10.1103/PhysRevB.13.5188](https://doi.org/10.1103/PhysRevB.13.5188)
- M. Methfessel, A. T. Paxton, *Phys. Rev. B* 40, 3616 (1989). DOI: [10.1103/PhysRevB.40.3616](https://doi.org/10.1103/PhysRevB.40.3616)
- N. Marzari, D. Vanderbilt, A. De Vita, M. C. Payne, *Phys. Rev. Lett.* 82, 3296 (1999). DOI: [10.1103/PhysRevLett.82.3296](https://doi.org/10.1103/PhysRevLett.82.3296)
- K. Lejaeghere et al., *Science* 351, aad3000 (2016)（Δ 测试与代码可重复性）. DOI: [10.1126/science.aad3000](https://doi.org/10.1126/science.aad3000)
- E. Bosoni et al., *Nat. Rev. Phys.* 6, 45 (2024)（ACWF 验证）. DOI: [10.1038/s42254-023-00655-3](https://doi.org/10.1038/s42254-023-00655-3)
- P. Giannozzi et al., *J. Phys.: Condens. Matter* 29, 465901 (2017). DOI: [10.1088/1361-648X/aa8f79](https://doi.org/10.1088/1361-648X/aa8f79)

---

[← 上一章](03-pseudopotentials.md) | [返回目录](../README.md) | [下一章 →](05-relax-eos.md)
