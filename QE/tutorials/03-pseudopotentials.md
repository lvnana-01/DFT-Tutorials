# 第3章 赝势选择：SSSP、PseudoDojo 与 pslibrary（Pseudopotentials in QE: SSSP, PseudoDojo, pslibrary）

本章讲解 Quantum ESPRESSO 中赝势（pseudopotential）的选择与使用：模守恒（NC）、超软（US）与 PAW 三类赝势的区别，UPF 文件头部（`pseudo_type`、`z_valence`、`functional`、`core_correction`、`relativistic`）怎么读，以及三大来源怎么下载、怎么核对：SSSP 1.3（Efficiency/Precision，带推荐截断能）、PseudoDojo（ONCVPSP 模守恒，声子/Raman/GW 首选）和 pslibrary（`ld1.x` 生成的 US/PAW）。最后讨论多元素体系的截断能取法（`ecutrho` 的 dual）、泛函一致性、半芯态与非线性芯修正，并以 MgO 为例写出完整输入。

> [!NOTE]
> 本章示例在 [03-pseudopotentials](../examples/03-pseudopotentials/README.md)。仓库**不提交任何 UPF 文件**，只提供下载脚本、文件名和校验方法。文中输出均为**示意输出，数值仅供参考**。

---

## 1. 本章目标与适用场景

| 目标 | 你将学会 |
|---|---|
| 类型 | NC、US、PAW 的区别，以及它们对 `ecutwfc`/`ecutrho` 和可用功能的影响 |
| UPF | 用 `grep` 读出类型、价电子数、泛函、NLCC、是否全相对论、建议截断 |
| SSSP | 下载 1.3.0 Efficiency/Precision 压缩包与 JSON，读推荐截断，用 md5 校验 |
| PseudoDojo | 选择 `nc-sr-05`/`nc-fr-04` 等表格，理解截断“提示”（hint）的单位是 Ha |
| pslibrary / SG15 | 文件命名规则；用 `ld1.x` 自己生成 |
| 多元素 | 截断取各元素推荐值的最大值；泛函一致；价电子数核对 |
| 管理 | `pseudo_dir`、环境变量 `ESPRESSO_PSEUDO`、aiida-pseudo |

与 VASP 用户的对应关系：

| VASP | QE | 区别 |
|---|---|---|
| POTCAR（PAW，随许可证分发） | UPF 文件（NC/US/PAW 均可，免费下载） | QE 中每个元素可以来自不同的库，但泛函必须一致 |
| `potpaw_PBE/Fe_pv` 等后缀 | 文件名中的 `spn`、`dn`、`n` 等字段（见 §5.2） | 都表示半芯态 / 芯修正的选择 |
| POTCAR 中的 `ENMAX` | UPF 头部的 `Suggested minimum cutoff`、SSSP JSON 中的 `cutoff_wfc`/`cutoff_rho` | QE 要分别给出波函数和密度截断 |
| `ZVAL` | `z_valence` | |
| VASP 推荐 POTCAR 列表 | SSSP Efficiency/Precision 表格 | SSSP 是经过系统验证的组合 |

适用场景：**每个新体系开始前**。选定的赝势会决定之后所有章节的截断能（[第4章](04-convergence.md)）和可用功能（例如 Raman 需要 NC，见 [第13章](13-dielectric-ir-raman.md)；SOC 需要全相对论赝势，见 [第7章](07-magnetism-soc.md)）。

---

## 2. 原理速览：三类赝势

赝势用一个平滑的有效势替代原子核和芯电子，只保留价电子（有时包括部分半芯电子）。在截断半径 $r_c$ 以外，赝波函数与全电子波函数相同；在 $r_c$ 以内，赝波函数没有节点，因此只需要较少的平面波。

| 类型 | 原理要点 | `ecutwfc` | `ecutrho` | 典型用途与限制 |
|---|---|---|---|---|
| 模守恒 NC（norm-conserving；ONCV、Troullier–Martins 等） | 赝波函数在 $r_c$ 内的电荷与全电子相同（模守恒条件） | 较高 | 默认 4×`ecutwfc` 即可 | 功能最全：ph.x 的 Raman/电光张量（`lraman`、`elop`）只支持 NC；GW（Yambo）通常也用 NC；杂化泛函/meta-GGA 的部分功能只对 NC 实现 |
| 超软 US（ultrasoft，Vanderbilt） | 放弃模守恒，用增强电荷 $Q_{ij}(\mathbf r)$ 补偿 | 低 | 通常 8–12×`ecutwfc` | 结构、能带、声子都常用；部分高级功能受限（例如 ph.x 遇到 US/PAW + `lraman` 会报错） |
| PAW（projector augmented wave，Blöchl） | 在球内重建全电子波函数 | 低 | 依增强电荷形状而定，常为 4–8×，必须测试 | 精度接近全电子；可重建全电子密度（Bader 分析，见 [第10章](10-charge-elf-bader.md)） |

INPUT_PW 中 `ecutrho` 的说明是：NC 赝势应保持默认值（4×`ecutwfc`）；有 US 赝势时通常需要 8–12 倍；PAW 常常可以用 4 倍，但取决于增强电荷的形状，**必须测试**。此外，梯度修正泛函（尤其在含真空的体系中）可能需要更高的 `ecutrho`。

> [!IMPORTANT]
> 截断能是“赝势的属性”，不是“材料的属性”。同一种材料换一套赝势，收敛所需的 `ecutwfc` 可能相差一倍以上。所以“文献中用了 40 Ry”不能直接照搬，除非用的是同一个 UPF 文件。

---

## 3. UPF 文件：怎么读头部

QE 使用 UPF（Unified Pseudopotential Format）格式。新版 UPF v2 是类 XML 格式，所有关键信息都在 `<PP_HEADER .../>` 标签的属性里。

### 3.1 一条命令看全部关键信息

```bash
grep -o -E 'pseudo_type="[^"]*"|relativistic="[^"]*"|core_correction="[^"]*"|functional="[^"]*"|z_valence="[^"]*"|wfc_cutoff="[^"]*"|rho_cutoff="[^"]*"' Si.pbe-n-rrkjus_psl.1.0.0.UPF
grep -i "Suggested" Si.pbe-n-rrkjus_psl.1.0.0.UPF
```

以 pslibrary 的 `Si.pbe-n-rrkjus_psl.1.0.0.UPF` 为例，其头部包含以下信息（这是文件内容本身，不是计算结果）：

```text
pseudo_type="USPP"
relativistic="scalar"
core_correction="true"
functional="PBE"
z_valence="4.000000000000e0"
    Suggested minimum cutoff for wavefunctions:  44. Ry
    Suggested minimum cutoff for charge density: 175. Ry
```

> [!NOTE]
> 这个 UPF 文件自己“建议”的 `ecutwfc` 是 44 Ry，而 SSSP 1.3 Efficiency 对**同一个文件**给出的推荐值是 30/240 Ry。原因是两者的判据不同：ld1.x 的建议值来自原子测试，SSSP 的推荐值来自对晶体性质（相图能量、声子、应力、能带等）的系统收敛测试。**以 SSSP 推荐值为起点，再按 [第4章](04-convergence.md) 测试你关心的性质**。

### 3.2 头部属性含义

| 属性 | 取值示例 | 含义与用途 |
|---|---|---|
| `pseudo_type` | `NC`、`SL`、`US`/`USPP`、`PAW` | 赝势类型；决定 `ecutrho` 的 dual 以及可用功能 |
| `relativistic` | `scalar`、`full` | `full` 表示全相对论（FR），含自旋–轨道，`lspinorb=.true.` 时必需 |
| `has_so` | `true`/`T` | 是否包含自旋–轨道信息 |
| `core_correction` | `true`/`false`（或 `T`/`F`） | 是否有非线性芯修正（NLCC） |
| `functional` | `PBE`、`SLA PW PBX PBC`（等价于 PBE）、`SLA PZ NOGX NOGC`（LDA） | 生成赝势所用的交换关联泛函；pw.x 默认从这里读取泛函 |
| `z_valence` | `4.0`、`16.0` | 价电子数（含半芯态）。总电子数 = Σ<sub>原子</sub> z_valence − `tot_charge` |
| `wfc_cutoff`、`rho_cutoff` | Ry | 生成程序写入的建议截断（不是所有库都填写） |
| `is_ultrasoft`、`is_paw` | `true`/`false` | |
| `number_of_wfc` | 整数 | 原子波函数个数，PDOS（projwfc.x）和 DFT+U 投影会用到 |

### 3.3 批量查看一个目录

```bash
cd pseudo/
for f in *.UPF *.upf; do
  [ -f "$f" ] || continue
  printf "%-40s " "$f"
  grep -o -E 'pseudo_type="[^"]*"|z_valence="[^"]*"|functional="[^"]*"|relativistic="[^"]*"' "$f" \
    | head -n 4 | tr '\n' ' '
  echo
done
```

示例目录中的 [`upf_info.sh`](../examples/03-pseudopotentials/upf_info.sh) 是同样功能的脚本，并会提示“混用了不同泛函”“有全相对论文件”等情况。

### 3.4 pw.x 输出中确认实际读到的赝势

```bash
grep -A4 "PseudoPot. #" mgo.scf.out
grep -A4 "atomic species   valence    mass" mgo.scf.out
grep "Exchange-correlation" mgo.scf.out
```

`atomic species valence mass pseudopotential` 表中的 `valence` 列应与 `z_valence` 一致；`number of electrons` 应等于 Σ(原子数 × valence)。

---

## 4. SSSP：首选的“验证过的组合”

### 4.1 SSSP 是什么

SSSP（Standard Solid-State Pseudopotentials，Materials Cloud）从多个库（pslibrary、GBRV、SG15、PseudoDojo、Wentzcovitch 等）中为每个元素挑选赝势，并给出经过系统测试的推荐截断。它有两种协议：

| 协议 | 目标 | 适用 |
|---|---|---|
| **Efficiency** | 在可接受精度下截断尽量低 | 高通量、结构优化、大体系、日常计算 |
| **Precision** | 更严格的收敛判据，截断更高，部分元素换成更“硬”的赝势 | 需要高精度的声子、弹性、能量差；做方法对比 |

当前稳定版是 **SSSP 1.3.0**（2023-04 发布），提供 PBE 和 PBEsol 两套，各 103 个元素。sssp.materialscloud.org 上的 v2.0 仍是 alpha 版，本系列不使用。

从 1.3.0 PBE 的 JSON 文件统计（这是对官方元数据文件的统计，不是计算结果）：Efficiency 中 85 个元素的推荐 dual（`cutoff_rho/cutoff_wfc`）为 8，16 个为 4，2 个为 12（Fe、Mn）；dual 为 4 的 16 个元素全部是 NC 赝势（15 个来自 SG15，Hf 来自 PseudoDojo）。

### 4.2 下载

SSSP 存放在 Materials Cloud Archive。每种配置（版本 × 泛函 × 协议）有两个文件，命名为 `SSSP_{version}_{functional}_{protocol}.{tar.gz,json}`：

| 文件 | 内容 |
|---|---|
| `SSSP_1.3.0_PBE_efficiency.tar.gz` | 全部 UPF 文件（压缩包内没有子目录） |
| `SSSP_1.3.0_PBE_efficiency.json` | 每个元素的 `filename`、`md5`、`pseudopotential`（来源库）、`cutoff_wfc`、`cutoff_rho`（单位 Ry） |
| `SSSP_1.3.0_PBE_precision.{tar.gz,json}` | Precision 协议 |
| `SSSP_1.3.0_PBEsol_{efficiency,precision}.{tar.gz,json}` | PBEsol 版本 |

```bash
mkdir -p ~/pseudo/SSSP_1.3.0_PBE_efficiency && cd ~/pseudo/SSSP_1.3.0_PBE_efficiency
REC=https://archive.materialscloud.org/api/records/rcyfm-68h65/files
wget -O SSSP_1.3.0_PBE_efficiency.tar.gz "$REC/SSSP_1.3.0_PBE_efficiency.tar.gz/content"
wget -O SSSP_1.3.0_PBE_efficiency.json   "$REC/SSSP_1.3.0_PBE_efficiency.json/content"
tar -xzf SSSP_1.3.0_PBE_efficiency.tar.gz
ls | head
```

> [!NOTE]
> 上面的记录号 `rcyfm-68h65` 是截至 2026-11 的最新版本记录。如果将来 Materials Cloud 发布了新的记录版本，可以在 [SSSP 1.3 表格页面（legacy 站点）](https://legacy.materialscloud.org/discover/sssp/table/efficiency) 找到下载入口，或者用 aiida-pseudo 安装（§8.2）。注意 `www.materialscloud.org/discover/sssp/...` 现在会重定向到 v2.0 alpha 站点 sssp.materialscloud.org。示例目录中的 [`download_sssp.sh`](../examples/03-pseudopotentials/download_sssp.sh) 封装了下载、解压和 md5 校验（脚本由你自己运行）。

### 4.3 读推荐截断并校验

```bash
python3 - <<'EOF'
import json
d = json.load(open("SSSP_1.3.0_PBE_efficiency.json"))
for el in ["Mg", "O", "Si", "Al", "H"]:
    v = d[el]
    print(f"{el:3s} {v['filename']:40s} {v['pseudopotential']:12s} "
          f"wfc={v['cutoff_wfc']:5.1f} rho={v['cutoff_rho']:6.1f} Ry")
EOF
md5sum Mg.pbe-n-kjpaw_psl.0.3.0.UPF O.pbe-n-kjpaw_psl.0.1.UPF
```

本系列用到的元素在 SSSP 1.3.0 PBE 中的条目（直接摘自官方 JSON，单位 Ry）：

| 元素 | Efficiency 文件 | 来源 | wfc / rho | Precision 文件 | wfc / rho |
|---|---|---|---|---|---|
| H | `H.pbe-rrkjus_psl.1.0.0.UPF` | pslibrary 1.0.0 US | 60 / 480 | `H_ONCV_PBE-1.0.oncvpsp.upf` | 80 / 320 |
| C | `C.pbe-n-kjpaw_psl.1.0.0.UPF` | pslibrary 1.0.0 PAW | 45 / 360 | 同左 | 45 / 360 |
| O | `O.pbe-n-kjpaw_psl.0.1.UPF` | pslibrary 0.1 PAW | 50 / 400 | 同左 | 75 / 600 |
| Mg | `Mg.pbe-n-kjpaw_psl.0.3.0.UPF` | pslibrary 0.3.0 PAW | 30 / 240 | `mg_pbe_v1.4.uspp.F.UPF` | 45 / 360 |
| Al | `Al.pbe-n-kjpaw_psl.1.0.0.UPF` | pslibrary 1.0.0 PAW | 30 / 240 | 同左 | 30 / 240 |
| Si | `Si.pbe-n-rrkjus_psl.1.0.0.UPF` | pslibrary 1.0.0 US | 30 / 240 | 同左 | 30 / 240 |
| Fe | `Fe.pbe-spn-kjpaw_psl.0.2.1.UPF` | pslibrary 0.2.1 PAW | 90 / 1080 | 同左 | 90 / 1080 |
| Ni | `ni_pbe_v1.4.uspp.F.UPF` | GBRV 1.4 US | 45 / 360 | 同左 | 50 / 400 |
| Cu | `Cu.paw.z_11.ld1.psl.v1.0.0-low.upf` | pslibrary 1.0.0 PAW（low） | 90 / 720 | 同左 | 90 / 720 |
| Ga | `Ga.pbe-dn-kjpaw_psl.1.0.0.UPF` | pslibrary 1.0.0 PAW | 70 / 560 | 同左 | 90 / 720 |
| As | `As.pbe-n-rrkjus_psl.0.2.UPF` | pslibrary 0.2 US | 35 / 280 | `As.nc.z_15.oncvpsp3.dojo.v4-std.upf` | 75 / 300 |
| Mn | `mn_pbe_v1.5.uspp.F.UPF` | GBRV 1.5 US | 65 / 780 | 同左 | 90 / 1080 |

> [!WARNING]
> 注意文件名的**大小写和扩展名**：GBRV 的文件名全是小写开头（`ni_pbe_v1.4.uspp.F.UPF`），部分文件的扩展名是小写 `.upf`。Linux 文件系统区分大小写，`ATOMIC_SPECIES` 中必须逐字符一致，否则报 `readpp: file ... not found`。

### 4.4 Materials Cloud 工具

- [SSSP 1.3 Efficiency 表](https://legacy.materialscloud.org/discover/sssp/table/efficiency) / [SSSP 1.3 Precision 表](https://legacy.materialscloud.org/discover/sssp/table/precision)（legacy 站点）：点击元素可查看收敛测试图；
- [QE 输入生成器](https://qeinputgenerator.materialscloud.io/)：上传结构后自动填入 SSSP 文件名和截断。

---

## 5. pslibrary 与官方赝势页

### 5.1 来源

pslibrary 是 A. Dal Corso 维护的一套 `ld1.x` 输入文件集合，用来生成标量相对论和全相对论的 US 赝势与 PAW 数据集，覆盖到 Pu（Z = 94）。SSSP 中很多元素直接选用了 pslibrary 的文件。

获取方式有两种：

1. 从 [QE 官方赝势页](https://pseudopotentials.quantum-espresso.org/legacy_tables) 下载已生成的单个文件，例如：

   ```bash
   wget https://pseudopotentials.quantum-espresso.org/upf_files/Mg.pbe-n-kjpaw_psl.0.3.0.UPF
   wget https://pseudopotentials.quantum-espresso.org/upf_files/O.pbe-n-kjpaw_psl.0.1.UPF
   ```

2. 用自己的 QE 编译出的 `ld1.x` 重新生成（见 §5.3）。

### 5.2 文件命名规则

QE 赝势页的命名约定为 `元素.描述.UPF`，其中描述 = `[字段1-][字段2-]字段3-[字段4-]字段5[_字段6]`：

| 字段 | 取值 | 含义 |
|---|---|---|
| 字段1（可选） | `rel` | 全相对论（fully relativistic），用于 SOC |
| 字段2（可选） | `starNl`、`starhNl` | 含一个/半个芯空穴（XSpectra 等用） |
| 字段3 | `pz`、`vwn`、`pbe`、`blyp`、`pw91`、`tpss`、`coulomb` | 交换关联泛函 |
| 字段4（可选） | `s`、`p`、`d` 的组合，以及 `n`、`f` | `s`/`p`/`d`：该半芯态放入价层；`n`：非线性芯修正；`f`：f 电子留在芯中 |
| 字段5 | `rrkjus`、`kjpaw`、`bpaw`、`rrkj`、`van`、`mt`、`bhs`、`vbc`、`ae` | 生成方法：`rrkjus` = 超软，`kjpaw` = PAW，`rrkj`/`mt`/`bhs`/`vbc` = 模守恒 |
| 字段6（可选） | `psl.1.0.0` 等 | 版本或作者；`psl` 表示 pslibrary |

解读示例：

| 文件名 | 含义 |
|---|---|
| `Si.pbe-n-rrkjus_psl.1.0.0.UPF` | Si，PBE，有 NLCC，超软，pslibrary 1.0.0 |
| `Fe.pbe-spn-kjpaw_psl.0.2.1.UPF` | Fe，PBE，3s3p 半芯态 + NLCC，PAW，pslibrary 0.2.1（`z_valence` = 16） |
| `Ga.pbe-dn-kjpaw_psl.1.0.0.UPF` | Ga，PBE，3d 半芯态 + NLCC，PAW |
| `Pt.rel-pbe-n-rrkjus_psl.1.0.0.UPF` | Pt，全相对论（字段1 = `rel`），PBE，NLCC，超软，用于 SOC |

pslibrary 1.0.0 的说明文件还指出：`*_ps_high.job` 生成半芯态最多、截断半径最小、最精确但最“硬”的赝势；`*_ps_low.job` 生成为低截断（目标约 40 Ry）优化的赝势；文件名中的一个或两个字母 `l` 表示“低截断版本”。

### 5.3 用 ld1.x 生成（概要）

```bash
git clone https://github.com/dalcorso/pslibrary.git
cd pslibrary
# edit QE_path: set PWDIR to your QE installation (the directory containing bin/ld1.x)
vi QE_path
# generate everything for the default functionals (pz, pbe, rel-pz, rel-pbe) ...
. ./make_all_ps
# ... or a single functional directory, after setting element='Si' in make_ps
cd pbe && . ../make_ps
ls PSEUDOPOTENTIALS/
```

注意事项：

- pslibrary 主页强调，不要用低于 5.1.2 的 QE 版本生成这些赝势；
- `ld1.x` 只能生成 LDA/GGA 赝势，不能生成杂化泛函或 meta-GGA 赝势；
- 自己生成的赝势必须自己测试（如对比全电子的晶格常数），一般用户直接下载即可。

---

## 6. PseudoDojo：模守恒赝势的首选

### 6.1 特点

PseudoDojo（<https://www.pseudo-dojo.org>）提供用 ONCVPSP 生成的优化模守恒赝势，并对每个元素做了 Δ-gauge、GBRV 测试集、声子、ghost 态等系统检验。适用场景：

- **DFPT 拉曼**（ph.x 的 `lraman`/`elop` 遇到 US/PAW 会直接报错，只能用 NC，见 [第13章](13-dielectric-ir-raman.md)）；
- **GW/BSE**（Yambo 通常使用 NC，见 [第13章](13-dielectric-ir-raman.md) 的指引）；
- **SOC**：`nc-fr-04` 表提供全相对论版本；
- 杂化泛函、meta-GGA（部分功能只支持 NC，见 [第14章](14-hybrid-metagga.md)）。

### 6.2 表格与下载

网站上的表格类型（`TYP`）包括：

| 表格 | 含义 |
|---|---|
| `nc-sr-05` | NC，标量相对论，PseudoDojo 表格 v0.5（只有 PBE；网站说明 v0.5 更新了 Ba、Bi、I、Pb、Po、Rb、Rn、S、Te、Tl、Xe，LDA/PBEsol 请用 v0.4.1） |
| `nc-sr-04` | NC，标量相对论，PseudoDojo 表格 v0.4.1（PBE、PBEsol、LDA） |
| `nc-fr-04` | NC，**全相对论**，PseudoDojo 表格 v0.4（PBE、PBEsol；SOC 用） |
| `nc-sr-04-3plus` | 镧系 3+ 组态，f 电子冻结在芯中 |
| `paw`、`paw-sr-11` | JTH PAW（主要供 ABINIT 使用） |

表名中的 `04`/`05` 是 PseudoDojo 表格的版本号，不是 ONCVPSP 程序的版本。各表可选的泛函不同（见上表），精度都分 standard/stringent 两档。下载 QE 用的 UPF 格式压缩包：

```bash
mkdir -p ~/pseudo && cd ~/pseudo
wget https://www.pseudo-dojo.org/pseudos/nc-sr-05_pbe_standard_upf.tgz
wget https://www.pseudo-dojo.org/pseudos/nc-fr-04_pbe_standard_upf.tgz   # for SOC
tar -xzf nc-sr-05_pbe_standard_upf.tgz
ls nc-sr-05_pbe_standard_upf/ | head      # files are named Si.upf, Mg.upf, O.upf, ...
```

压缩包解压后得到 `nc-sr-05_pbe_standard_upf/` 目录，文件名就是 `元素.upf`（例如 `Mg.upf`、`O.upf`）。因为不同表格中的文件同名，**建议每个表格单独放一个目录**，并在输入中用不同的 `pseudo_dir`，避免混淆。

### 6.3 截断“提示”的单位是 Hartree

PseudoDojo 为每个元素给出 low / normal / high 三档截断提示（hints），**单位是 Ha**。用于 QE 时要乘以 2 换算成 Ry：

| 元素（nc-sr-05, PBE, standard） | low（Ha） | normal（Ha） | high（Ha） | normal 换算为 Ry |
|---|---|---|---|---|
| Si | 14 | 18 | 24 | 36 |
| Al | 16 | 20 | 26 | 40 |
| Mg | 38 | 42 | 48 | 84 |
| O | 36 | 42 | 48 | 84 |
| H | 32 | 36 | 42 | 72 |

（上表摘自 PseudoDojo 网站的 `nc-sr-05_pbe_standard.json` 元数据，不是计算结果。）

NC 赝势的 `ecutrho` 用默认的 4×`ecutwfc` 即可。网站说明 normal 提示适合高通量计算的初始猜测，high 提示之上结果不应再有明显变化；具体仍需按 [第4章](04-convergence.md) 测试。

---

## 7. 其它来源

| 库 | 类型 | 说明 |
|---|---|---|
| [SG15 ONCV](http://www.quantum-simulation.org/potentials/sg15_oncv/) | NC（ONCVPSP） | Schlipf–Gygi 优化的 ONCV 赝势；文件名如 `Si_ONCV_PBE-1.2.upf`，全相对论版如 `Si_ONCV_PBE_FR-1.1.upf`；压缩包 `sg15_oncv_upf_2020-02-06.tar.gz` |
| GBRV | US | 高通量用超软赝势，SSSP 中很多过渡金属来自它（文件名如 `ni_pbe_v1.4.uspp.F.UPF`） |
| QE 源码中的 `pseudo/` | 混合 | 只供官方 examples/test-suite 使用（如 `Si.bhs`、`C.UPF`），**不要用于研究** |
| QE 赝势页 legacy tables | 混合 | 旧的 `Si.pz-vbc.UPF` 等 LDA 赝势；PARADIM 等旧教程常用，新计算不推荐 |

> [!TIP]
> 旧教程（包括 Cornell PARADIM 2016/2018、MIT 3.320 等）常用 `Si.pz-vbc.UPF` 等 LDA 模守恒赝势。它们用来学习原理没问题，但做研究时请换成 SSSP 或 PseudoDojo，并相应地重新测试截断能。

---

## 8. 在输入中使用赝势

### 8.1 `pseudo_dir` 与 `ESPRESSO_PSEUDO`

| 方式 | 写法 | 说明 |
|---|---|---|
| 输入中显式写 | `pseudo_dir = './pseudo'` 或绝对路径 | **推荐**，输入自包含、可复现 |
| 环境变量 | `export ESPRESSO_PSEUDO=$HOME/pseudo/SSSP_1.3.0_PBE_efficiency` | 未写 `pseudo_dir` 时使用 |
| 默认值 | `$HOME/espresso/pseudo/` | 两者都没设置时 |

推荐的目录结构：

```text
~/pseudo/
  SSSP_1.3.0_PBE_efficiency/     <- all UPF + SSSP_1.3.0_PBE_efficiency.json
  SSSP_1.3.0_PBE_precision/
  nc-sr-05_pbe_standard_upf/     <- PseudoDojo scalar-relativistic
  nc-fr-04_pbe_standard_upf/     <- PseudoDojo fully relativistic
```

每个计算目录下建一个 `pseudo/` 子目录，用符号链接指向需要的文件，这样输入可以直接打包分享：

```bash
mkdir -p pseudo
ln -s ~/pseudo/SSSP_1.3.0_PBE_efficiency/Mg.pbe-n-kjpaw_psl.0.3.0.UPF pseudo/
ln -s ~/pseudo/SSSP_1.3.0_PBE_efficiency/O.pbe-n-kjpaw_psl.0.1.UPF   pseudo/
```

### 8.2 aiida-pseudo（可选）

[aiida-pseudo](https://aiida-pseudo.readthedocs.io/en/latest/) 可以一键安装 SSSP 和 PseudoDojo 到 AiiDA 数据库，并提供推荐截断：

```bash
pip install aiida-pseudo
aiida-pseudo install sssp -v 1.3 -x PBE -p efficiency
aiida-pseudo install pseudo-dojo -v 0.4 -x PBE -r SR -p standard -f upf
aiida-pseudo list
```

`install sssp` 的可选值：版本 `1.0`–`1.3`（默认 1.3），泛函 `PBE`/`PBEsol`，协议 `efficiency`/`precision`；`install pseudo-dojo` 的版本为 `0.4`、`0.5`、`1.0`，相对论 `SR`/`SR3plus`/`FR`，格式 `psp8`/`upf`/`psml`/`jthxml`（默认 `psp8`，QE 要选 `upf`）。它需要一个配置好的 AiiDA profile；`--download-only` 选项可以只下载压缩包到当前目录。自动化工作流见 [第18章](18-automation-parallel-troubleshooting.md)。

### 8.3 多元素体系：截断取最大值

一个计算只能有一组 `ecutwfc`/`ecutrho`，所以**取所有元素推荐值中的最大值**：

```text
ecutwfc = max_i(cutoff_wfc_i)
ecutrho = max_i(cutoff_rho_i)
```

以 MgO（SSSP 1.3 PBE Efficiency）为例：

| 元素 | 文件 | wfc | rho |
|---|---|---|---|
| Mg | `Mg.pbe-n-kjpaw_psl.0.3.0.UPF` | 30 | 240 |
| O | `O.pbe-n-kjpaw_psl.0.1.UPF` | 50 | 400 |
| **MgO 取值** | | **50** | **400** |

这只是起点。氧化物、含过渡金属的体系常常由“最硬”的元素决定截断；若该元素在体系中含量很少，也不能降低截断。仓库脚本 `QE/scripts/qe_inputgen.py` 可以根据 SSSP JSON 自动完成这一步（在仓库根目录运行 `python QE/scripts/qe_inputgen.py -h` 查看用法）。

### 8.4 示例：MgO（SSSP 1.3 Efficiency）

示例文件 [`mgo.scf.in`](../examples/03-pseudopotentials/mgo.scf.in)：

```fortran
&CONTROL
  calculation = 'scf'
  prefix      = 'mgo'
  outdir      = './tmp'
  pseudo_dir  = './pseudo'
  tprnfor     = .true.
  tstress     = .true.
/
&SYSTEM
  ibrav     = 2
  celldm(1) = 7.96           ! Bohr (a = 4.21 Angstrom, rock salt)
  nat       = 2
  ntyp      = 2
  ecutwfc   = 50.0           ! max(Mg 30, O 50)
  ecutrho   = 400.0          ! max(Mg 240, O 400)
/
&ELECTRONS
  conv_thr = 1.0d-10
/
ATOMIC_SPECIES
  Mg  24.305  Mg.pbe-n-kjpaw_psl.0.3.0.UPF
  O   15.999  O.pbe-n-kjpaw_psl.0.1.UPF
ATOMIC_POSITIONS crystal
  Mg  0.0  0.0  0.0
  O   0.5  0.5  0.5
K_POINTS automatic
  8 8 8 0 0 0
```

同一体系用 PseudoDojo NC 赝势的版本见 [`mgo.scf.dojo.in`](../examples/03-pseudopotentials/mgo.scf.dojo.in)：文件改为 `Mg.upf`、`O.upf`（来自 `nc-sr-05_pbe_standard_upf/`），`ecutwfc = 84`（normal 提示 42 Ha），`ecutrho` 用默认的 4 倍。两个输入得到的**总能量不能直接比较**（赝势不同），但晶格常数、体模量、带隙等物理量在收敛后应接近。

运行后检查：

```bash
grep -A4 "atomic species   valence" mgo.scf.out
grep "number of electrons" mgo.scf.out
grep "Exchange-correlation" mgo.scf.out
```

（示意输出，数值仅供参考）

```text
     atomic species   valence    mass     pseudopotential
     Mg               2.00    24.30500     Mg( 1.00)
     O                6.00    15.99900     O ( 1.00)
     number of electrons       =         8.00
```

`Mg.pbe-n-kjpaw_psl.0.3.0.UPF` 的 `z_valence` 是 2（只有 3s），`O.pbe-n-kjpaw_psl.0.1.UPF` 是 6，所以 MgO 原胞有 8 个价电子。PseudoDojo `nc-sr-05_pbe_standard_upf/Mg.upf` 的头部写的是 `z_valence="   10.00"`（2s2p 半芯态放入价层，且 `core_correction="F"`），O.upf 是 6，所以 `mgo.scf.dojo.in` 的原胞有 16 个价电子（以上都是文件头部的内容，不是计算结果）。换一个表格或版本后，请用 §3 的 grep 重新核对，这正是“核对价电子数”的意义。

---

## 9. 泛函一致性、半芯态与 NLCC

### 9.1 泛函一致性

- pw.x 默认从赝势读取泛函，并检查**所有赝势的泛函是否一致**；不一致时报错 `inconsistent DFT read from PP files`；
- `input_dft` 可以强制覆盖，但 PW 用户手册明确说“风险自负”。例如在 PBE 赝势上做 HSE、PBE0 是常见且可以接受的做法（[第14章](14-hybrid-metagga.md)），但用 LDA 赝势做 PBE 计算是不推荐的；
- PBEsol 计算请用 SSSP PBEsol 或 PseudoDojo PBEsol 表格中的赝势。

### 9.2 半芯态（semicore）

把外层芯电子（如过渡金属的 3s3p、Ga 的 3d）放入价层：

| 何时需要 | 代价 |
|---|---|
| 高压、短键、强离子性体系；过渡金属磁性；需要准确的 d 带位置 | 更多电子、更高截断、更慢 |
| 声子、弹性等对芯区敏感的性质 | |

文件名中的 `spn`、`dn`、`sp` 等字段，以及 `z_valence` 的数值，告诉你包含了哪些半芯态。例如 `Fe.pbe-spn-kjpaw_psl.0.2.1.UPF` 的 `z_valence` = 16（3s² 3p⁶ 3d⁶ 4s² 共 16 个电子），而只含 3d4s 的 Fe 赝势是 8。

### 9.3 非线性芯修正（NLCC）

交换关联能是密度的非线性函数，芯电荷与价电荷在空间上重叠时，只用价电荷计算 XC 会引入误差。NLCC 在 XC 计算中加入一个部分芯电荷。文件名中的 `n` 字段和头部的 `core_correction="true"` 表示有 NLCC。磁性体系和有半芯重叠的体系尤其需要它。

### 9.4 相对论

| 计算 | 需要的赝势 |
|---|---|
| 普通计算（含共线磁性） | 标量相对论（`relativistic="scalar"`） |
| `noncolin=.true.` 不含 SOC | 标量相对论即可 |
| `lspinorb=.true.`（SOC） | **全相对论**（`relativistic="full"`，文件名含 `rel-` 或 `_FR`），见 [第7章](07-magnetism-soc.md) |

---

## 参数选择与收敛建议

| 任务 | 推荐赝势 | 截断起点 |
|---|---|---|
| 日常结构、能带、DOS、表面吸附 | SSSP 1.3 Efficiency | JSON 中的最大 `cutoff_wfc`/`cutoff_rho` |
| 声子、弹性常数、精确能量差、方法对比 | SSSP 1.3 Precision | 同上，并做收敛测试 |
| Raman、GW/BSE、部分 meta-GGA/杂化泛函功能 | PseudoDojo NC（`nc-sr-05`）或 SG15 | normal 提示 × 2（Ha→Ry），`ecutrho` = 4× |
| SOC | PseudoDojo `nc-fr-04`、SG15 `_FR`、pslibrary `rel-` | 同上 |
| PBEsol 计算 | SSSP PBEsol 或 PseudoDojo PBEsol | |
| 杂化泛函 vc-relax（7.6.1 及以后） | NC（develop 分支 release notes 说明 7.6.1 起 USPP/PAW 下不再计算杂化泛函应力） | 见 [第14章](14-hybrid-metagga.md) |
| Bader 分析需要全电子密度 | PAW | 见 [第10章](10-charge-elf-bader.md) |

通用建议：

- **整个项目固定一套赝势**，写进论文方法部分（逐个列出文件名，最好附 md5）；
- 换赝势就要重新做截断收敛测试；
- `ecutrho/ecutwfc`：NC 4；US 8–12；PAW 4–8，以测试为准；
- 应力、声子对截断更敏感，比总能量需要更高的截断（[第4章](04-convergence.md)）。

---

## 常见问题与报错

| 报错 / 现象 | 原因 | 解决 |
|---|---|---|
| `Error in routine readpp ... file ./pseudo/Mg.UPF not found` | 文件名与 `ATOMIC_SPECIES` 不一致（大小写、`.UPF`/`.upf`），或 `pseudo_dir` 错误 | `ls -l $pseudo_dir`；从 SSSP JSON 的 `filename` 字段复制文件名 |
| `Error in routine readpp ... inconsistent DFT read from PP files` | 混用了不同泛函的赝势（例如 LDA 的 `Si.pz-vbc.UPF` 与 PBE 的 O） | 统一换成同一泛函的赝势（PW 用户手册 5.0.0.6） |
| `number of electrons` 与预期不符 | 某个赝势包含（或不包含）半芯态 | 查看 `z_valence`；改变赝势后重新确认 `tot_charge`、`nbnd`、`starting_magnetization` 等设置 |
| `ecutrho < 4*ecutwfc, are you sure?` | `ecutrho` 太小 | 至少 4×；US/PAW 按 SSSP 推荐（多为 8×） |
| `S matrix not positive definite` / cdiaghg 报错 | 结构有问题（原子过近）；赝势有 ghost 态；或 US 赝势给出非正的电荷密度，使 S 矩阵失去正定性（PW 用户手册 5.0.0.7） | 先检查结构；提高 `ecutrho`；换一个经过验证的赝势（SSSP/PseudoDojo）；试 `diagonalization='cg'` |
| 能量随 `ecutwfc` 振荡、不单调 | US/PAW 的 `ecutrho` 不够 | 固定 `ecutwfc`，单独扫描 `ecutrho`（[第4章](04-convergence.md)） |
| `Spin-orbit calculations require at least one spin-orbit pseudo`；或提示 `Not all pseudopotentials have spin-orbit data` | `lspinorb=.true.` 但用了标量相对论赝势（全部或部分） | 换用全相对论赝势（`rel-`、`_FR`、`nc-fr-04`） |
| ph.x 报 `The phonon code with US-PP and raman or elop not yet available`（PAW 时为 `The phonon code with paw and raman or elop is not yet available`） | 用 US/PAW 赝势设了 `lraman` 或 `elop` | Raman 需要 NC 赝势（[第13章](13-dielectric-ir-raman.md)） |
| PseudoDojo 截断按“Ry”用，结果没收敛 | PseudoDojo 提示的单位是 Ha | 乘以 2 |
| 解压 PseudoDojo 后文件互相覆盖 | 不同表格的文件同名（`Si.upf`） | 每个表格单独一个目录 |
| 下载的 UPF 是 HTML 网页 | 链接错误或被重定向到登录页 | `head -n 3 file.UPF` 应看到 `<UPF version=` 或 `<PP_INFO>`；与 SSSP JSON 中的 md5 比对 |

---

## 练习

1. 下载 `SSSP_1.3.0_PBE_efficiency.json`，写一个 Python 片段列出所有 dual = 4 的元素及其来源库，并说明它们大多属于哪类赝势。
2. 对 `Fe.pbe-spn-kjpaw_psl.0.2.1.UPF`，用 §3 的 grep 命令读出 `z_valence`，并根据文件名解释为什么是 16。
3. 用 SSSP Efficiency 和 Precision 两套赝势分别写出 GaAs 的 pw.x 输入（`ibrav=2`），各自应取多大的 `ecutwfc`/`ecutrho`？哪一套适合做 Raman？
4. 把 `mgo.scf.in` 的 O 换成 SG15 的 `O_ONCV_PBE-1.2.upf`，Mg 保持 PAW，这样混用是否合法？截断应如何选？（提示：泛函是否一致、NC 与 PAW 的 dual 取哪个。）
5. 用 `mgo.scf.in` 和 `mgo.scf.dojo.in`（自己运行）分别做晶格常数扫描（[第5章](05-relax-eos.md)），比较两套赝势得到的晶格常数和体模量，而不是总能量。

---

## 参考资料

**官方文档与数据**

- [INPUT_PW：`ecutrho`、`pseudo_dir`、`input_dft`](https://www.quantum-espresso.org/Doc/INPUT_PW.html)
- [QE 官方赝势页](https://pseudopotentials.quantum-espresso.org/) / [legacy tables](https://pseudopotentials.quantum-espresso.org/legacy_tables) / [命名约定](https://pseudopotentials.quantum-espresso.org/home/naming-convention) / [UPF 格式说明](https://pseudopotentials.quantum-espresso.org/home/unified-pseudopotential-format)
- [SSSP 1.3 Efficiency 表](https://legacy.materialscloud.org/discover/sssp/table/efficiency) / [SSSP 1.3 Precision 表](https://legacy.materialscloud.org/discover/sssp/table/precision) / [Materials Cloud Archive 记录](https://archive.materialscloud.org/records/rcyfm-68h65) / [SSSP v2.0 alpha 站点](https://sssp.materialscloud.org/)
- [PseudoDojo](https://www.pseudo-dojo.org/) / [PseudoDojo FAQ](https://www.pseudo-dojo.org/faq.html)
- [pslibrary 主页](https://dalcorso.github.io/pslibrary/)
- [SG15 ONCV 赝势](http://www.quantum-simulation.org/potentials/sg15_oncv/)
- [INPUT_LD1（ld1.x）](https://www.quantum-espresso.org/Doc/INPUT_LD1.html)
- [aiida-pseudo 文档](https://aiida-pseudo.readthedocs.io/en/latest/)

**课程与教程（推荐，链接而非复制）**

- [Imperial College MSE404-MM Lab 2](https://imperialcollegelondon.github.io/MSE404-MM/labs/lab02/readme/)：在第一个计算中解释赝势文件与截断
- [Cornell PARADIM Computational Tutorials](https://www.paradim.org/toolbox/computation/tutorials)：讲解输入与收敛参数（原教程使用 LDA `pz-vbc` 赝势，学原理可以，研究请换成 SSSP/PseudoDojo）
- [Advanced QE school 2023（EPFL THEOS / Materials Cloud）](https://github.com/materialscloud-org/QuantumESPRESSO-school-2023)：结合 AiiDA 与 aiida-pseudo 管理赝势
- [Ghent Computational Materials Physics](https://www.compmatphys.org)：赝势概念的入门讲解

**相关仓库（仅链接）**

- [dalcorso/pslibrary](https://github.com/dalcorso/pslibrary)：ld1.x 输入集
- [abinit/pseudo_dojo](https://github.com/abinit/pseudo_dojo)：PseudoDojo 生成与验证框架
- [oncvpsp/oncvpsp](https://github.com/oncvpsp/oncvpsp)：ONCV 赝势生成程序
- [aiidateam/aiida-pseudo](https://github.com/aiidateam/aiida-pseudo)：一键安装 SSSP/PseudoDojo 并给推荐截断

**关键文献**

- G. Prandini, A. Marrazzo, I. E. Castelli, N. Mounet, N. Marzari, *npj Comput. Mater.* 4, 72 (2018)（SSSP）. DOI: [10.1038/s41524-018-0127-2](https://doi.org/10.1038/s41524-018-0127-2)
- M. J. van Setten et al., *Comput. Phys. Commun.* 226, 39 (2018)（PseudoDojo）. DOI: [10.1016/j.cpc.2018.01.012](https://doi.org/10.1016/j.cpc.2018.01.012)
- A. Dal Corso, *Comput. Mater. Sci.* 95, 337 (2014)（pslibrary）. DOI: [10.1016/j.commatsci.2014.07.043](https://doi.org/10.1016/j.commatsci.2014.07.043)
- D. R. Hamann, *Phys. Rev. B* 88, 085117 (2013)（ONCV）. DOI: [10.1103/PhysRevB.88.085117](https://doi.org/10.1103/PhysRevB.88.085117)
- M. Schlipf, F. Gygi, *Comput. Phys. Commun.* 196, 36 (2015)（SG15）. DOI: [10.1016/j.cpc.2015.05.011](https://doi.org/10.1016/j.cpc.2015.05.011)
- D. Vanderbilt, *Phys. Rev. B* 41, 7892 (1990)（超软赝势）. DOI: [10.1103/PhysRevB.41.7892](https://doi.org/10.1103/PhysRevB.41.7892)
- P. E. Blöchl, *Phys. Rev. B* 50, 17953 (1994)（PAW）. DOI: [10.1103/PhysRevB.50.17953](https://doi.org/10.1103/PhysRevB.50.17953)
- N. Troullier, J. L. Martins, *Phys. Rev. B* 43, 1993 (1991). DOI: [10.1103/PhysRevB.43.1993](https://doi.org/10.1103/PhysRevB.43.1993)
- E. Bosoni et al., *Nat. Rev. Phys.* 6, 45 (2024)（全电子与赝势代码的系统验证）. DOI: [10.1038/s42254-023-00655-3](https://doi.org/10.1038/s42254-023-00655-3)

---

[← 上一章](02-input-anatomy.md) | [返回目录](../README.md) | [下一章 →](04-convergence.md)
