# 示例 08：DFT+U、DFT+U+V 与 hp.x 线性响应（NiO / LiCoO2）

配套教程：[第8章 DFT+U/U+V 与 hp.x 线性响应](../../tutorials/08-dftu-hp.md)

## 目的

- 用 QE ≥ 7.1 的 `HUBBARD` 卡对 AFM NiO 做 DFT+U（Dudarev，`ortho-atomic` 投影）；
- 按官方 `HP/examples/example02` 的思路，用三步法（假金属 SCF → 绝缘体 SCF → hp.x）从第一性原理计算 Ni-3d 的 U，并可选地迭代到自洽 U；
- 用 hp.x 计算 LiCoO2 的 DFT+U+V 参数，得到可直接粘贴到 pw.x 输入的 `HUBBARD.dat`。

## 文件说明

| 文件 | 说明 |
|---|---|
| `nio.dftu.scf.in` | NiO AFM-II DFT+U 单点（`U Ni1-3d 6.0` 为**演示值**） |
| `nio.lr_step1.in` | 线性响应第 1 步模板：小 smearing 的"假金属" SCF，占位符 `@U_IN@` |
| `nio.lr_step2.in` | 第 2 步模板：`occupations='fixed'` + 整数 `tot_magnetization` + 同样 `nbnd`，占位符 `@U_IN@`、`@NBND@`、`@TOTMAG@` |
| `nio.hp.in` | 第 3 步：hp.x 输入（`nq1=nq2=nq3=2`） |
| `run_hp.sh` | 自动执行 1→2→3，并可迭代（`NCYCLE`、`TOL`），结果记录在 `u_history.dat` |
| `licoo2.uv.scf.in` | LiCoO2 的 DFT+U+V 基态（U、V 初值 0.0001 eV） |
| `licoo2.hp.in` | LiCoO2 的 hp.x 输入（`num_neigh=6`、`lmin=2`） |
| `run_uv.sh` | LiCoO2：pw.x → hp.x → 打印 `HUBBARD.dat` |
| `job_slurm.sh` | Slurm 包装（`WHAT=hp/uv/dftu`） |

## 赝势准备（不提交 UPF 文件）

从 [SSSP 1.3 PBE Efficiency](https://www.materialscloud.org/discover/sssp/table/efficiency) 下载 `SSSP_1.3.0_PBE_efficiency.tar.gz`，解压后把以下文件放入 `pseudo/`：

```bash
mkdir -p pseudo
#   ni_pbe_v1.4.uspp.F.UPF       (SSSP: 45 / 360 Ry)
#   O.pbe-n-kjpaw_psl.0.1.UPF    (SSSP: 50 / 400 Ry)
#   Co_pbe_v1.2.uspp.F.UPF       (SSSP: 45 / 360 Ry)
#   li_pbe_v1.4.uspp.F.UPF       (SSSP: 40 / 320 Ry)
grep -c 'PP_CHI\|3D' pseudo/ni_pbe_v1.4.uspp.F.UPF    # the 3d orbital must be present
```

HUBBARD 卡中的 `Ni1-3d`、`O-2p` 等流形必须在赝势的原子波函数（`PP_CHI` 的 label，或旧格式 GBRV 文件头的 `Wavefunctions nl` 表）中存在。

## 运行步骤

```bash
cd QE/examples/08-dftu-hp

# A) DFT+U single point with a fixed (demo) U
mpirun -np 4 pw.x -nk 2 -in nio.dftu.scf.in > nio.dftu.scf.out

# B) linear-response U for NiO (one cycle from plain PBE)
bash run_hp.sh
# self-consistent U: iterate until |U_out - U_in| < 0.01 eV (max 5 cycles)
NCYCLE=5 TOL=0.01 bash run_hp.sh

# C) DFT+U+V for LiCoO2
bash run_uv.sh

# on Slurm
sbatch --export=ALL,NCYCLE=5 job_slurm.sh
sbatch --export=ALL,WHAT=uv job_slurm.sh
```

## 预期结果 / 如何检查

```bash
# DFT+U run
grep -e '^!' -e 'Hubbard energy' -e 'highest occupied' nio.dftu.scf.out
grep -A3 'Tr\[ns' nio.dftu.scf.out | tail -8
grep 'Atomic magnetic moment for atom' nio.dftu.scf.out | tail -2

# hp.x
grep -A6 'Hubbard U parameters' nio_lr.Hubbard_parameters.dat
cat u_history.dat
cat HUBBARD.dat                        # after run_uv.sh
```

- `nio.dftu.scf.out` 中出现 `Hubbard energy` 项和 `HUBBARD OCCUPATIONS` 块；两个 Ni 的 `Atomic magnetic moment` 大小相等、符号相反；带隙相对纯 PBE 明显增大。
- hp.x 的 U 表中 `Ni1` 与 `Ni2` 的 U 相同，`new_label` 合并为同一类型。官方 `HP/examples/example02` 的参考输出（PBEsol、pslibrary 超软赝势、2×2×2 k 与 q）给出约 7.9 eV——**换了赝势、泛函或网格，数值都会变**，只作量级参考。
- LiCoO2 的 `HUBBARD.dat` 形如 `V Co-3d Co-3d 1 1 <U>` 加 6 行 `V Co-3d O-2p 1 <j> <V>`（八面体配位的 6 个 O 近邻）。

## 注意事项

- **HUBBARD 卡是 7.1 起的唯一写法**；`lda_plus_u`、`Hubbard_U(i)`、`U_projection_type` 等旧关键词会直接报 `DFT+Hubbard input syntax has changed since v7.1`。
- hp.x 只支持 `atomic` 与 `ortho-atomic` 投影；不支持 Liechtenstein（`J`）形式、`J0`、杂化泛函、meta-GGA、tetrahedra、Gamma 技巧、约束磁化、`assume_isolated='2D'`；也不能直接读取 vc-relax 的结果（先用新结构重做 scf）。
- 官方 HP 示例 README 要求：在 `ATOMIC_POSITIONS` 中**先列 Hubbard 原子**（U 或 V 非零），再列其他原子，否则 hp.x 停止。
- 线性响应 U 与投影类型、赝势、泛函、结构绑定；`ortho-atomic` 的 U 不能拿到 `atomic` 计算中用，更不能拿到 VASP 中用。
- 示例中的 k 点（4×4×4）与 q 点（2×2×2）只是起点，U 对二者都需要收敛。
