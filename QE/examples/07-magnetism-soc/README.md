# 示例 07：自旋极化、反铁磁、非共线与自旋轨道耦合（bcc Fe / fcc Ni / NiO / GaAs）

配套教程：[第7章 磁性、非共线与自旋轨道耦合](../../tutorials/07-magnetism-soc.md)

## 目的

- 用 `nspin = 2` + `starting_magnetization` 计算 bcc Fe、fcc Ni 的铁磁基态，并与非磁（NM）态比较能量；
- 用"同一元素拆成两个 species"的方法建立 NiO 的 AFM-II 构型；
- 用 `noncolin = .true.` + `lspinorb = .true.` + 全相对论（FR）赝势计算 GaAs 的 SOC 能带（价带顶劈裂）；
- 用 `lforcet = .true.` 的力定理流程估算 bcc Fe 的磁晶各向异性能（MAE）；
- 用无 SOC 的非共线计算检查"磁矩方向不影响能量"。

## 文件说明

| 文件 | 说明 |
|---|---|
| `fe.scf.in` | bcc Fe 铁磁 SCF（`nspin=2`，PAW，`mv` 展宽） |
| `fe.nm.scf.in` | bcc Fe 非磁参考（`nspin=1`），用于计算磁性稳定化能 |
| `ni.scf.in` | fcc Ni 铁磁 SCF（GBRV 超软赝势） |
| `nio_afm.scf.in` | NiO AFM-II，4 原子菱方胞，`Ni1`/`Ni2` 两个 species，初始磁化 ±0.5 |
| `gaas_soc.scf.in` | GaAs 非磁 SOC 自洽（`noncolin` + `lspinorb`，FR PAW 赝势） |
| `gaas_soc.bands.in` | GaAs SOC 能带 L–Γ–X（`tpiba_b`） |
| `gaas_soc.bands_pp.in` | `bands.x` 后处理，`lsigma(3)` 输出 ⟨σz⟩ |
| `fe_mae.lsda.in` | MAE 第 1 步：用 FR 超软赝势做共线 LSDA 自洽 |
| `fe_mae.nscf_001.in` / `_100.in` / `_110.in` | MAE 第 2 步：`lforcet` 非自洽 SOC，磁化分别沿 [001]、[100]、[110] |
| `fe_noncolin.scf.in` | 无 SOC 的非共线 Fe（磁化沿 x），用于自检 |
| `run_all.sh` | 依次运行各任务并 grep 关键结果；`bash run_all.sh mae` 只跑 MAE |
| `job_slurm.sh` | Slurm 包装脚本（`sbatch --export=ALL,TASK=mae job_slurm.sh`） |

## 赝势准备（不提交 UPF 文件）

```bash
mkdir -p pseudo && cd pseudo
# 1) SSSP 1.3 PBE Efficiency（https://www.materialscloud.org/discover/sssp/table/efficiency）
#    下载 SSSP_1.3.0_PBE_efficiency.tar.gz，解压后取：
#      Fe.pbe-spn-kjpaw_psl.0.2.1.UPF   (SSSP: 90 / 1080 Ry)
#      ni_pbe_v1.4.uspp.F.UPF           (SSSP: 45 / 360 Ry)
#      O.pbe-n-kjpaw_psl.0.1.UPF        (SSSP: 50 / 400 Ry)
# 2) pslibrary 1.0.0 全相对论赝势（QE 官方赝势页）：
wget https://pseudopotentials.quantum-espresso.org/upf_files/Ga.rel-pbe-dn-kjpaw_psl.1.0.0.UPF
wget https://pseudopotentials.quantum-espresso.org/upf_files/As.rel-pbe-n-kjpaw_psl.1.0.0.UPF
wget https://pseudopotentials.quantum-espresso.org/upf_files/Fe.rel-pbe-spn-rrkjus_psl.1.0.0.UPF
cd ..
grep -i 'relativistic\|has_so' pseudo/*.rel-*.UPF | head    # 应看到 relativistic="full" has_so="true"
```

赝势选择的背景见 [第3章](../../tutorials/03-pseudopotentials.md)。

## 运行步骤

```bash
cd QE/examples/07-magnetism-soc
bash run_all.sh fe      # bcc Fe FM + NM
bash run_all.sh ni      # fcc Ni
bash run_all.sh nio     # NiO AFM-II (plain PBE)
bash run_all.sh soc     # GaAs SOC scf + bands + bands.x
bash run_all.sh mae     # Fe force-theorem MAE (dense k mesh, the most expensive part)
bash run_all.sh nc      # noncollinear Fe without SOC

# custom launcher
PW_CMD="mpirun -np 8 pw.x -nk 4" bash run_all.sh fe

# on a Slurm cluster
sbatch --export=ALL,TASK=soc job_slurm.sh
```

## 预期结果 / 如何检查

```bash
grep -e 'total magnetization' -e 'absolute magnetization' fe.scf.out | tail -2
grep -A2 'Magnetic moment per site' fe.scf.out | tail -2
grep '^!' fe.scf.out fe.nm.scf.out
grep -A5 'Magnetic moment per site' nio_afm.scf.out | tail -5
grep 'Noncollinear calculation\|Non magnetic calculation' gaas_soc.scf.out
grep 'eband, Ef (eV)' fe_mae.nscf_*.out
```

- bcc Fe：总磁化接近文献中 PBE 的 2.2 μB/atom 左右，FM 能量明显低于 NM；fcc Ni 约 0.6 μB/atom（量级参考，以你的计算为准）。
- NiO：`total magnetization` ≈ 0，两个 Ni 位点磁矩大小相等、符号相反；纯 PBE 的带隙远小于实验值，是第8章 DFT+U 的动机。
- GaAs：输出应出现 `Non magnetic calculation with spin-orbit`；Γ 点价带顶由 3 重（不计自旋）劈裂为 Γ8（4 重）和 Γ7（2 重），劈裂量级约 0.3 eV（实验约 0.34 eV）。
- MAE：三个方向的 `eband` 之差就是力定理 MAE，bcc Fe 的 MAE 只有 μeV/atom 量级，**比 k 点误差还小是常态**，必须做 k 点收敛（见教程 §7）。

## 注意事项

- AFM 必须把同一元素拆成两个 species（`Ni1`、`Ni2`），否则无法给出相反的初始磁化；拆分后对称性降低、k 点数增加属正常。
- `noncolin = .true.` 时**不要**写 `nspin`；`lspinorb = .true.` 需要全相对论赝势（文件名中常带 `rel`）。
- `lforcet` 在 QE 7.5 源码中对 PAW 直接报错 `force theorem with PAW not implemented`，因此 MAE 示例使用超软 FR 赝势。
- 力定理比较的是**同一套 k 点、同一组 smearing 参数**下的 `eband`；nscf 输入设 `nosym = .true.`（与官方 `PP/examples/ForceTheorem_example` 一致），因此 24×24×24 网格展开为 13824 个 k 点，是本目录最贵的计算；三个方向的 nscf 都从同一份 LSDA 电荷出发，脚本为每个方向复制一份 `fe_mae.save`。
- 所有截断能、k 点、`degauss` 均为教学用起点，生产计算前请按 [第4章](../../tutorials/04-convergence.md) 做收敛测试。
