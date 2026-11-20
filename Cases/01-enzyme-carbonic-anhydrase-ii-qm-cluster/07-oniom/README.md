# 07-oniom：ONIOM(B3LYP:Amber) 替代方案模板

本目录只有一个**骨架模板** `g16_oniom_template.gjf`，没有附带完整蛋白的输入。原因是：一个能用的 ONIOM(QM:MM) 输入要先完成整套 MM 准备，包括加氢、Zn 位点参数化（例如 MCPB.py）、溶剂化、MM 平衡和抽取快照，每个原子行都要写上 Amber 原子类型、电荷、冻结标记和层标记。这些都要由你自己生成，不能手写。

| 占位符 | 应替换为 |
|---|---|
| `@CHARGE_MULT_REAL_LOW@` 等 | 二层 ONIOM 的三对电荷/多重度：real-low、model-high、model-low（见 Gaussian 第15章） |
| `@ONIOM_ATOM_LINES@` | `元素-Amber类型-电荷  冻结标记  x y z  层标记 [链接原子信息]` |
| `@CONNECTIVITY@` | `Geom=Connectivity` 要求的连接表 |
| `@MM_PARAMETERS_IF_ANY@` | Zn 位点等非标准参数（HrmStr1/HrmBnd1/AmbTrs 等），没有就删掉这一行 |

QM 区建议与簇模型一致：Zn、OH⁻、CO2、三个 His 的咪唑环（在 Cα–Cβ 处切断，加链接原子）、Thr199 侧链，可选 Glu106 侧链。MM 区为其余蛋白和水。

步骤说明见主教程 [Step 11](../README.md#step-11-进阶oniomb3lypamber-与-cp2k-qmmm-替代方案)。语法细节见 [Gaussian 第15章 ONIOM](../../../Gaussian/tutorials/15-oniom.md)，CP2K 方案见 [CP2K 第16章](../../../CP2K/tutorials/16-qmmm-xtb.md)。

> [!WARNING]
> 不要把簇模型的冻结 Cβ 坐标直接搬进 ONIOM。ONIOM 中的 MM 环境本身就起到约束作用，一般只冻结离 QM 区较远的 MM 原子（例如 QM 区 10–15 Å 以外）。
