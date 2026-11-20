# 示例 03：赝势下载、检查与 MgO 多元素截断（SSSP / PseudoDojo）

配套教程：[第3章 赝势选择：SSSP、PseudoDojo 与 pslibrary](../../tutorials/03-pseudopotentials.md)

## 目的

- 下载 SSSP 1.3.0（Efficiency/Precision，PBE/PBEsol），用官方 JSON 中的 md5 校验每个 UPF 文件，并读出推荐截断；
- 用一个只读脚本检查 UPF 头部（类型、相对论、NLCC、价电子数、泛函），发现“泛函混用”“全相对论文件”等问题；
- 以 MgO 为例演示多元素体系的截断取法（取各元素推荐值的最大值），并给出 SSSP（PAW）与 PseudoDojo（NC）两个版本的输入。

> [!WARNING]
> 本目录**不包含任何 UPF 文件**，也不要把 UPF 提交到仓库。脚本只是写好供你自己运行的，本仓库没有执行过它们，也没有运行过 pw.x。

## 文件说明

| 文件 | 说明 |
|---|---|
| `download_sssp.sh` | 从 Materials Cloud Archive（记录 `rcyfm-68h65`）下载 `SSSP_1.3.0_<FUNC>_<PROT>.tar.gz` 与 `.json`，解压到 `~/pseudo/<NAME>/`，校验 md5，打印常用元素的推荐截断；`LINK="Mg O"` 时在当前目录的 `pseudo/` 中建立符号链接 |
| `upf_info.sh` | 只读：打印 UPF 头部字段与文件内的建议截断，混用多种泛函时给出警告 |
| `mgo.scf.in` | MgO（岩盐，`ibrav=2`，`celldm(1)=7.96` Bohr），SSSP 1.3 PBE Efficiency，`ecutwfc=50`、`ecutrho=400` Ry |
| `mgo.scf.dojo.in` | 同一结构，PseudoDojo `nc-sr-05` PBE standard（`Mg.upf`、`O.upf`），`ecutwfc=84` Ry（42 Ha × 2），`ecutrho` 用默认的 4 倍 |

## 赝势准备

| 输入 | 需要的文件 | 来源 | 推荐截断 |
|---|---|---|---|
| `mgo.scf.in` | `Mg.pbe-n-kjpaw_psl.0.3.0.UPF` | SSSP 1.3.0 PBE efficiency | 30 / 240 Ry |
| `mgo.scf.in` | `O.pbe-n-kjpaw_psl.0.1.UPF` | SSSP 1.3.0 PBE efficiency | 50 / 400 Ry |
| `mgo.scf.dojo.in` | `Mg.upf`、`O.upf` | `nc-sr-05_pbe_standard_upf.tgz` | normal 提示 42 Ha（= 84 Ry） |

```bash
cd QE/examples/03-pseudopotentials

# SSSP (PAW/US): download + verify + symlink Mg and O into ./pseudo
LINK="Mg O" bash download_sssp.sh

# PseudoDojo (NC): download the UPF tarball and copy the two files
mkdir -p pseudo_dojo
wget https://www.pseudo-dojo.org/pseudos/nc-sr-05_pbe_standard_upf.tgz
tar -xzf nc-sr-05_pbe_standard_upf.tgz
cp nc-sr-05_pbe_standard_upf/Mg.upf nc-sr-05_pbe_standard_upf/O.upf pseudo_dojo/
```

## 运行步骤

```bash
# 1) inspect the headers (read-only)
bash upf_info.sh pseudo/* pseudo_dojo/*

# 2) run the two SCF calculations yourself
mpirun -np 4 pw.x -nk 2 -i mgo.scf.in      > mgo.scf.out
mpirun -np 4 pw.x -nk 2 -i mgo.scf.dojo.in > mgo.scf.dojo.out
```

## 预期结果 / 如何检查

```bash
grep -A3 "atomic species   valence" mgo.scf.out mgo.scf.dojo.out
grep "number of electrons" mgo.scf.out mgo.scf.dojo.out
grep -e "kinetic-energy cutoff" -e "charge density cutoff" mgo.scf.out mgo.scf.dojo.out
grep "Exchange-correlation" mgo.scf.out mgo.scf.dojo.out
grep "convergence has been achieved" mgo.scf.out mgo.scf.dojo.out
grep "JOB DONE" mgo.scf.out mgo.scf.dojo.out
```

| 检查项 | `mgo.scf.in` | `mgo.scf.dojo.in` |
|---|---|---|
| `valence` 列 | Mg 2.00，O 6.00（来自 `z_valence`） | Mg 10.00，O 6.00（`nc-sr-05` 的 `Mg.upf` 头部 `z_valence` 为 10，含 2s2p 半芯态；换表格后以 `upf_info.sh` 读出的值为准） |
| `number of electrons` | 8.00 | 16.00（与左列不同是正常的） |
| `charge density cutoff` | 400 Ry | 336 Ry（4 × 84） |
| `Exchange-correlation` | PBE | PBE |

（示意输出，数值仅供参考）

```text
mgo.scf.out:     number of electrons       =         8.00
mgo.scf.out:!    total energy              =     -xx.xxxxxxxx Ry
```

两个输入的总能量**不能互相比较**（赝势不同，参考能量不同）。要比较的是收敛后的物理量，例如用 [第5章](../../tutorials/05-relax-eos.md) 的方法求出的晶格常数和体模量。

## 注意事项

- SSSP JSON 中 `cutoff_wfc`/`cutoff_rho` 的单位是 Ry；PseudoDojo 网站的截断提示单位是 Ha，用于 QE 要乘以 2。
- 文件名区分大小写：例如 GBRV 文件 `ni_pbe_v1.4.uspp.F.UPF` 是小写开头；部分 SSSP 文件扩展名是 `.upf`。
- PseudoDojo 不同表格中的文件同名（都叫 `Mg.upf`），请分目录存放。
- 赝势的选择与截断要写进论文的方法部分；换赝势必须重新做收敛测试（[第4章](../../tutorials/04-convergence.md)）。
