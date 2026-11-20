# 示例 09：范德华修正、二维材料、表面 slab 与功函数（graphite / MoS2 / Al(111)）

配套教程：[第9章 范德华修正、二维材料、表面与功函数](../../tutorials/09-vdw-surface-workfunction.md)

## 目的

- 用不同色散方案（无修正、DFT-D2、DFT-D3(zero/BJ)、rVV10、vdW-DF2-B86R；TS 需要模守恒赝势，见注意事项）弛豫石墨，比较层间距；
- 用 `cell_dofree = '2Dxy'` 只优化单层 MoS2 的面内晶格常数；
- 用 `assume_isolated = '2D'` 计算双层石墨烯（2D 库仑截断）；
- 建立 Al(111) 6 层 slab（底两层固定），加偶极修正，经 pp.x（`plot_num = 11`）→ average.x 得到平面平均静电势和功函数；
- 用取向相同的体相胞计算表面能。

## 文件说明

| 文件 | 说明 |
|---|---|
| `graphite.vdw.in` | AB 石墨 vc-relax，默认 `vdw_corr='grimme-d3'` + `dftd3_version=4`（D3-BJ）；`@VDW@` 行供脚本替换 |
| `run_vdw_compare.sh` | 依次运行 none/d2/d3bj/d3zero/rvv10/df2b86r（`ts` 需手动指定且要换 NC 赝势），汇总到 `vdw_summary.dat`（c 已换算为 Å） |
| `mos2_mono.vc-relax.in` | 单层 MoS2，`cell_dofree='2Dxy'`，约 17 Å 真空，k 12×12×1 |
| `bilayer_graphene.scf.in` | AB 双层石墨烯，`assume_isolated='2D'`，D3-BJ |
| `al111_slab.relax.in` | Al(111) 1×1、6 层、约 18 Å 真空，底两层 `0 0 0` 固定，偶极修正 |
| `al111_slab.scf.in` | 功函数单点（把弛豫后坐标粘贴进来） |
| `al_bulk_oriented.scf.in` | 与 slab 同面内胞、沿 [111] 取向的 3 原子体相胞（表面能参考） |
| `pp_v.in` | pp.x：`plot_num=11`（V_bare + V_H），输出 filplot 与 cube |
| `average.in` | average.x 输入：1 个文件、权重 1.0、1000 点、沿方向 3、宏观平均窗口 4.4078 bohr（= Al(111) 层间距） |
| `run_workfunction.sh` | relax → scf → pp.x → average.x → 粗略给出两侧真空能级与 Φ |
| `job_slurm.sh` | Slurm 包装（`WHAT=wf/vdw/mos2`） |

## 赝势准备（不提交 UPF 文件）

从 [SSSP 1.3 PBE Efficiency](https://www.materialscloud.org/discover/sssp/table/efficiency) 下载 `SSSP_1.3.0_PBE_efficiency.tar.gz`，解压后放入 `pseudo/`：

```bash
mkdir -p pseudo
#   C.pbe-n-kjpaw_psl.1.0.0.UPF     (SSSP: 45 / 360 Ry)
#   Mo_ONCV_PBE-1.0.oncvpsp.upf     (SSSP: 35 / 140 Ry)
#   s_pbe_v1.4.uspp.F.UPF           (SSSP: 35 / 280 Ry)
#   Al.pbe-n-kjpaw_psl.1.0.0.UPF    (SSSP: 30 / 240 Ry)
ls pseudo/
```

> rVV10 与 vdW-DF 系列是非局域泛函，`input_dft` 会覆盖赝势中的 PBE 标记；严格地说应配合与其交换部分一致的赝势，教学比较中使用 PBE 赝势是常见近似（见教程 §3.4）。

## 运行步骤

```bash
cd QE/examples/09-vdw-surface-workfunction

# vdW comparison on graphite
bash run_vdw_compare.sh                  # all methods
bash run_vdw_compare.sh d3bj rvv10       # selected

# MoS2 monolayer, in-plane relaxation only
mpirun -np 4 pw.x -nk 4 -in mos2_mono.vc-relax.in > mos2_mono.vc-relax.out

# bilayer graphene with 2D Coulomb cutoff
mpirun -np 4 pw.x -nk 4 -in bilayer_graphene.scf.in > bilayer_graphene.scf.out

# Al(111): relax -> scf -> pp.x -> average.x
bash run_workfunction.sh
SKIP_RELAX=1 bash run_workfunction.sh    # after pasting relaxed coordinates

# surface energy reference
mpirun -np 4 pw.x -nk 4 -in al_bulk_oriented.scf.in > al_bulk_oriented.scf.out

# work function with the repository helper
python ../../scripts/qe_workfunction.py -h
```

## 预期结果 / 如何检查

```bash
cat vdw_summary.dat
grep -e 'DFT-D3 Dispersion' -e 'Dispersion Correction' vdw_*/graphite.out | tail
grep -A3 'Begin final coordinates' mos2_mono.vc-relax.out
grep -A6 'Computed dipole along edir' al111_slab.scf.out | tail -7
grep 'the Fermi energy' al111_slab.scf.out
head avg.dat
```

- 石墨：无色散修正的 PBE 层间几乎不束缚（c 偏大很多或继续增大），D3/rVV10/vdW-DF2-B86R 应给出接近实验 c ≈ 6.71 Å（层间距约 3.35 Å）的结果。
- MoS2：`cell_dofree='2Dxy'` 下 c 轴保持 20 Å 不变，面内 a 在 3.18 Å 附近（PBE 量级参考）。
- Al(111)：`avg.dat` 第 2 列在两侧真空区都应平坦；对称 slab 两侧平台应基本相同（偶极近 0）。Al(111) 的实验功函数约 4.2–4.3 eV，PBE 结果通常在此附近（量级参考，以计算为准）。
- 表面能：γ = (E_slab − N·E_bulk/atom) / (2A)，A 为 slab 面积（本例 1×1 胞 A ≈ 7.07 Å²），注意单位换算（Ry → J 或 eV）。

## 注意事项

- `vdw_corr` 拼错**不会报错**，QE 7.5 只打印 `WARNING: unknown vdw correction (vdw_corr): ... No vdw correction used.` 并照常运行；务必 grep 输出中的色散能量项确认生效。
- `vdw_corr = 'ts'` / `'mbd'` 在 pw.x 中不支持超软/PAW 赝势（报 `Tkatchenko-Scheffler not implemented with USPP/PAW`）；本目录的 C 赝势是 PAW，所以脚本默认不跑 `ts`。
- 偶极修正要求 `tefield = .true.` **且** `dipfield = .true.`，`eamp = 0`；`emaxpos`/`eopreg` 定义的锯齿下降区必须落在真空中，否则会产生非物理受力。
- average.x 只能单进程运行（多进程报 `run on a single processor`）。
- `avg.dat` 的势能单位是 Ry、坐标是 bohr；pp.x 输出的"势"是电子势能（e·V），真空能级减费米能（eV）即功函数。
- Al(111) 的面内晶胞夹角为 60°，ABC 堆垛位置为分数坐标 A = (0,0)、B = (1/3,1/3)、C = (2/3,2/3)；换晶胞或换元素时请用 ASE `fcc111` 等工具重新生成坐标。
- 层数、真空厚度和 k 点都要收敛；固定层与表面能的关系见教程 §8。
