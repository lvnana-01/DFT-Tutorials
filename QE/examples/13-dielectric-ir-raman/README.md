# 示例 13：GaAs 介电张量、Born 有效电荷、LO-TO、拉曼与 H2O 红外/拉曼

对应教程：[第13章 介电常数、Born 有效电荷、红外与拉曼](../../tutorials/13-dielectric-ir-raman.md)

## 目的

- 用 `ph.x`（`epsil = .true.`、`lraman = .true.`）在 Γ 点计算闪锌矿 GaAs 的高频介电张量 ε∞、Born 有效电荷 Z* 和拉曼张量；
- 用 `dynmat.x` 施加声学求和规则（ASR）、加上非解析项得到 LO 模，输出红外（IR）强度与拉曼活性；
- 用 `ph.x` 的 q 网格 + `q2r.x` + `matdyn.x` 得到包含 LO-TO 劈裂的声子色散；
- 用 `epsilon.x` 计算独立粒子近似（IPA）下的频率相关介电函数；
- 对 H2O 分子做 Γ 点简正模分析，得到 IR 强度和拉曼活性。

## 文件说明

| 文件 | 说明 |
|---|---|
| `gaas.scf.in` | GaAs 自洽，LDA（PseudoDojo PW 模守恒赝势），`conv_thr = 1e-12`，8×8×8 k 网格 |
| `gaas.ph_epsil.in` | Γ 点 DFPT：`epsil`、`trans`、`lraman`，写出 `gaas.dynG` |
| `gaas.dynmat.in` | `dynmat.x`：`asr='crystal'`，`q = (1,0,0)` 方向的 LO 非解析项，输出 IR 强度与拉曼活性 |
| `gaas.dynmat_eps0.in` | `dynmat.x`：不设 `q`，`lperm = .true.` 输出各 TO 模对静态介电常数 ε0 的贡献（与 LO 方向的计算分开，否则 dynmat.x 会警告结果混入 TO-LO 劈裂） |
| `gaas.ph_disp.in`, `gaas.q2r.in`, `gaas.matdyn.in` | 4×4×4 q 网格声子 → 实空间力常数（带 Z* 与 ε∞）→ X–Γ–L–U–X 色散 |
| `gaas.nscf_eps.in`, `gaas.eps.in` | 12×12×12 全网格非自洽 + 40 条能带 → `epsilon.x`（`calculation='eps'`） |
| `h2o.relax.in` | H2O 放在 12 Å 立方盒子中紧收敛弛豫（k 网格写成 `automatic 1 1 1 0 0 0`，**不能**用 `K_POINTS gamma`） |
| `h2o.ph_gamma.in`, `h2o.dynmat.in` | H2O 的 Γ 点 `ph.x`（`epsil` + `lraman`；开了 `lraman` 后 ph.x 内部的 `asr` 不起作用，所以不写）+ `dynmat.x`（`asr='zero-dim'`） |
| `job.slurm` | 依次运行四组步骤，可用 `STEPS` 变量选择 |

## 赝势（不随仓库提供）

拉曼张量（`lraman`）只支持**模守恒 + LDA**。本例使用 PseudoDojo v0.4 的 LDA（PW）标准精度 UPF 包：

```bash
mkdir -p pseudo && cd pseudo
wget https://www.pseudo-dojo.org/pseudos/nc-sr-04_pw_standard_upf.tgz
tar -xzf nc-sr-04_pw_standard_upf.tgz Ga.upf As.upf O.upf H.upf
grep -m1 'functional=' Ga.upf        # should show: SLA  PW   NOGX NOGC (LDA)
grep -m1 'pseudo_type=' Ga.upf       # should show: NC
cd ..
```

包内文件名就是 `Ga.upf`、`As.upf`、`O.upf`、`H.upf`。PseudoDojo 网页对每个元素给出 low/normal/high 截断提示（单位 Ha，乘 2 换成 Ry）；Ga（含 3d 半芯）normal 提示约 40 Ha，所以输入里先用 `ecutwfc = 80` Ry，请自行做收敛测试。

如果只需要 ε∞ 和 Z*（不算拉曼），可以改用 PBE 模守恒赝势或 SSSP 中的 US/PAW 赝势，并删去 `lraman`。

## 运行步骤

```bash
cd QE/examples/13-dielectric-ir-raman
sbatch job.slurm                                   # all steps
sbatch --export=ALL,STEPS="gaas_gamma" job.slurm   # only Gamma-point dielectric + Raman
```

手动运行 GaAs Γ 点部分：

```bash
mpirun -np 8 pw.x -nk 4 -in gaas.scf.in      > gaas.scf.out
mpirun -np 8 ph.x -nk 4 -in gaas.ph_epsil.in > gaas.ph_epsil.out
dynmat.x -in gaas.dynmat.in      > gaas.dynmat.out
dynmat.x -in gaas.dynmat_eps0.in > gaas.dynmat_eps0.out
```

## 预期结果 / 如何检查

```bash
grep -A4 "Dielectric constant in cartesian axis" gaas.ph_epsil.out
grep -A8 "Effective charges (d Force / dE) in cartesian axis with asr applied" gaas.ph_epsil.out
grep -A1 "Raman tensor" gaas.dynG | head
grep -A8 "# mode" gaas.dynmat.out
grep -A4 "with zone-center polar mode contributions" gaas.dynmat_eps0.out   # eps0 = eps_inf + polar modes
grep -c "JOB DONE" gaas.ph_epsil.out
```

判断标准（定性）：

- 立方晶体的 ε∞ 和每个原子的 Z* 都应是对角且三个对角元相等；Ga 与 As 的 Z* 大小相等、符号相反（施加 ASR 后求和为零）。
- `dynmat.x` 输出中前 3 个模接近 0 cm⁻¹（ASR 生效），其余 3 个光学模中有 2 个 TO、1 个 LO（沿 q 方向），LO 频率高于 TO。
- `matdyn.x` 得到的色散在 Γ 点附近，沿不同方向趋近 Γ 时光学支频率应出现 LO-TO 劈裂。
- PBE/LDA 的带隙偏小，ε∞ 通常高于实验值；具体数值请以自己的计算为准。
- H2O：`dynmat.x` 中 6 个零频（3 平动 + 3 转动）附近的模，加上 3 个内振动模（弯曲、对称伸缩、反对称伸缩），每个模都有 IR 和 Raman 列。

## 注意事项

- `epsil = .true.` 只能用于绝缘体、只能在 q = 0；scf 中不要设置 smearing，否则 ph.x 报 `no elec. field with metals`。
- `lraman` 需要模守恒赝势、LDA、非自旋极化、CPU 版本；使用 US/PAW 赝势时 ph.x 报 `The phonon code with US-PP and raman or elop not yet available`（PAW 对应 `paw and raman or elop`）。
- QE 7.6 的 release notes 说明：二阶响应拉曼计算从 7.4.1 起出现问题，到 7.6 修复。请使用 7.6 及以上版本做拉曼计算。
- H2O 的 scf 不能用 `K_POINTS gamma`，否则 ph.x 报 `cannot start from pw.x data file using Gamma-point tricks`。
- `gaas.nscf_eps.in` 会覆盖 `tmp/gaas.save` 中的波函数；如果之后还要跑 ph.x，请先重新执行 scf，或把 nscf 的 `outdir` 改成另一个目录（同时复制 scf 的 save 目录）。
- 本目录所有数值都需要你自行计算得到；教程中的输出片段只是格式示意。
