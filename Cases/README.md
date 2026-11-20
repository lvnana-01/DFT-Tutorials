# 计算案例：从结构到机理（Case Studies）

“计算案例”板块把各软件板块中分散讲解的技能（结构建模、优化、过渡态、IRC、频率、溶剂模型、单点、后处理）串成完整的研究流程，每个案例都从一个真实体系出发，一步步走到机理结论。案例覆盖 enzyme mechanism、QM cluster model、transition state search、reaction energy profile 等主题，重点是“为什么这样设置”和“怎样检验结论是否可靠”。

> [!IMPORTANT]
> 本板块的案例只提供**工作流、输入模板和脚本**，**不提供计算结果**。仓库中没有附带任何输出文件、优化结构或能量；模板中的坐标用占位符表示，需要你用案例自带的脚本从 PDB 等原始数据生成，再自己完成计算。文中引用的数值要么来自原始实验数据（如 PDB 文件中的配位距离），要么是带 DOI 的文献值，要么明确标注为“示意输出，数值仅供参考”。

## 案例列表

| 编号 | 案例 | 体系与方法 | 涉及软件 | 难度 |
|---|---|---|---|---|
| 01 | [酶催化计算实例：hCA II 反应机理的 QM 簇模型研究](01-enzyme-carbonic-anhydrase-ii-qm-cluster/README.md) | 人碳酸酐酶 II（PDB 2CBA / 3D92），Zn(II)–OH⁻ 亲核进攻 CO2；簇模型截取、xTB 预优化、B3LYP-D3(BJ)/def2-SVP + CPCM(ε=4) 优化、TS、IRC、quasi-RRHO、def2-TZVP 单点、模型收敛性检验、ONIOM/QM/MM 替代方案 | Gaussian 16、ORCA 6、xTB、PDB2PQR/PROPKA、Multiwfn | 进阶 |

后续计划加入更多案例（如异相催化表面反应、电池材料离子迁移等），格式与本案例相同。

## 与软件板块的关系

案例不重复讲软件语法，遇到具体关键词时会链接回对应章节。建议先掌握下列章节，再做案例：

| 案例中的步骤 | 对应的软件章节 |
|---|---|
| 结构建模、截取簇模型 | [Common 第3章 结构建模](../Common/tutorials/03-structure-modeling.md) |
| 输入结构、冻结原子、ModRedundant | [Gaussian 第2章](../Gaussian/tutorials/02-input-anatomy.md)、[第5章 几何优化](../Gaussian/tutorials/05-opt.md) |
| 频率与热化学、quasi-RRHO | [Gaussian 第6章](../Gaussian/tutorials/06-freq-thermo.md) |
| 过渡态、IRC、柔性扫描 | [Gaussian 第7章](../Gaussian/tutorials/07-ts-irc.md)、[第8章](../Gaussian/tutorials/08-scan.md) |
| 隐式溶剂、自定义介电常数 | [Gaussian 第9章](../Gaussian/tutorials/09-solvation.md) |
| 电荷布居、NBO | [Gaussian 第12章](../Gaussian/tutorials/12-population-charges.md) |
| ONIOM 多层方法 | [Gaussian 第15章](../Gaussian/tutorials/15-oniom.md) |
| Multiwfn、IGMH | [Gaussian 第16章](../Gaussian/tutorials/16-postprocess-multiwfn.md) |
| 集群作业提交 | [Gaussian 第17章](../Gaussian/tutorials/17-hpc-jobs.md)、[Common 第2章 作业调度](../Common/tutorials/02-job-scheduler.md) |
| QM/MM、xTB/DFTB | [CP2K 第16章](../CP2K/tutorials/16-qmmm-xtb.md) |

## 每个案例的组织方式

```text
Cases/<NN>-<slug>/
├── README.md          完整教程：背景 -> 路线图 -> Step 1..N -> 常见问题 -> 练习 -> 参考资料
├── 00-structure/      从原始数据（PDB/CIF）生成模型的脚本
├── 01-.../ 02-.../    按计算步骤分目录的输入模板（.gjf / .inp / xcontrol 等）
├── 08-profile/        结果汇总与作图脚本 + 只含占位行的 CSV
└── jobs/              Slurm 提交脚本
```

使用约定：

- 输入文件为纯 ASCII、LF 换行，注释为英文；教程正文为简体中文，关键词保留英文；
- 模板中的 `@COORDS@`、`@FREEZE@` 等占位符由案例自带的 `fill_template.py` 填写，不要手工拼接坐标；
- PDB 等原始文件不随仓库分发，README 中给出下载命令；
- 结果汇总用的 CSV 只含占位行（`NaN`, `placeholder=yes`），作图脚本遇到占位行会拒绝运行，避免把占位数据误当结果。

## 推荐的外部资源

- Himo 等人关于 QM 簇模型方法的综述：F. Himo, *J. Am. Chem. Soc.* 2017, 139, 6780，[DOI: 10.1021/jacs.7b02671](https://doi.org/10.1021/jacs.7b02671)；P. E. M. Siegbahn, F. Himo, *WIREs Comput. Mol. Sci.* 2011, 1, 323，[DOI: 10.1002/wcms.13](https://doi.org/10.1002/wcms.13)
- [CCPBioSim/qmmm-workshop](https://github.com/CCPBioSim/qmmm-workshop)：酶 QM/MM 动手教程
- [natedey/RINRUS](https://github.com/natedey/RINRUS)：按残基相互作用网络系统地构建簇模型
- [ORCA 6.0 tutorials](https://www.faccts.de/docs/orca/6.0/tutorials/)：ORCA 官方教程，包括优化、TS 与溶剂模型

---

[返回仓库首页](../README.md)
