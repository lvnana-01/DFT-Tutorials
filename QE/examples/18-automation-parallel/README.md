# 示例 18：并行参数扫描、Slurm 模板、ASE 驱动与 XML 解析

配套教程：[第18章 自动化、并行性能与报错排查](../../tutorials/18-automation-parallel-troubleshooting.md)

> [!IMPORTANT]
> 本目录的脚本默认**只生成文件、不运行计算**：`ase_si_scf.py` 不加 `--run` 时只写输入，`scan_parallel.sh` 不设置 `SUBMIT=1` 时只生成作业目录。仓库作者没有运行这些计算；文中的输出示例都是"（示意输出，数值仅供参考）"。

## 目的

- 用同一个小算例（Si，8×8×8 k 网格）练习 pw.x 的并行参数 `-nk` 和 `OMP_NUM_THREADS` 的选择（`-nd` 只对能带数上百的大体系有意义，扫描脚本不包含它），并比较墙钟时间和内存估计；
- 读懂 Slurm 模板中每一行的作用：OpenMP 绑定、`ESPRESSO_TMPDIR`/`ESPRESSO_PSEUDO`，以及用 `prefix.EXIT` 实现超时前正常退出；
- 用 ASE（≥ 3.23 的 `EspressoProfile` 写法）生成或运行 pw.x 输入；
- 只用 Python 标准库解析 `data-file-schema.xml`，取出总能、费米能/HOMO-LUMO、力、应力。

## 文件说明

| 文件 | 说明 |
|---|---|
| `si.scf.in` | Si SCF（SSSP 超软赝势，`ecutwfc=30`/`ecutrho=240`，8×8×8 k），`tprnfor`、`tstress`、`max_seconds` |
| `job_qe_slurm.sh` | 带注释的 Slurm 模板：MPI + OpenMP，`-nk` 由环境变量 `NK` 控制，收到 `USR1` 时写 `prefix.EXIT` |
| `scan_parallel.sh` | 生成 `scan_nk*_omp*/` 目录和作业脚本；`SUBMIT=1` 时提交；`collect` 汇总墙钟时间和内存估计 |
| `ase_si_scf.py` | ASE 写 pw.x 输入（默认），或通过 `EspressoProfile` 运行（`--run`） |
| `parse_xml.py` | 标准库 `xml.etree` 解析 `data-file-schema.xml`，可输出 JSON |

仓库级的通用模板是 `QE/scripts/job_qe_slurm.sh` 和 `QE/scripts/job_qe_pbs.sh`（用法：`bash ../../scripts/job_qe_slurm.sh -h` 或直接阅读文件头注释）。

## 赝势准备（不提交 UPF 文件）

从 [SSSP 1.3 PBE Efficiency](https://www.materialscloud.org/discover/sssp/table/efficiency) 下载 `SSSP_1.3.0_PBE_efficiency.tar.gz`，解压后把 `Si.pbe-n-rrkjus_psl.1.0.0.UPF` 放到本目录的 `pseudo/`。详见 [第3章](../../tutorials/03-pseudopotentials.md)。

## 运行步骤

```bash
cd QE/examples/18-automation-parallel

# 1) one job with the annotated template (NK = number of k-point pools)
NK=4 sbatch job_qe_slurm.sh si.scf.in

# 2) parallel scan: prepare, submit, then collect
bash scan_parallel.sh si.scf.in                       # only creates scan_* directories
SUBMIT=1 NTOT=16 bash scan_parallel.sh si.scf.in      # submit all
bash scan_parallel.sh collect                         # after the jobs finish

# 3) ASE: write the input, inspect it, then (optionally) run
python ase_si_scf.py -h
python ase_si_scf.py --pseudo-dir ./pseudo
python ase_si_scf.py --pseudo-dir ./pseudo --run --command "mpirun -np 4 pw.x -nk 2"

# 4) parse the XML of a finished run
python parse_xml.py tmp/si.save/data-file-schema.xml --json si_result.json
```

## 预期结果 / 如何检查

`scan_parallel.sh collect` 的输出格式（示意输出，数值仅供参考）：

```text
run                      PWSCF_WALL  max_RAM_estimate
scan_nk1_omp1                 xx.xxs        xx.xx MB
scan_nk2_omp1                 xx.xxs        xx.xx MB
scan_nk4_omp1                 xx.xxs        xx.xx MB
...
```

检查要点：

- 每个输出都要有 `JOB DONE`，且各次运行的总能量（`grep '^!' scan_*/pw.out`）一致到 `conv_thr` 的量级；不同并行方式只应改变耗时，不应改变结果；
- 对这个只有 2 个原子、几十个不可约 k 点的小体系，`-nk` 增大通常更快；`-nd` 对它没有帮助（能带数太少）；
- `Estimated max dynamical RAM per process` 随 `-nk` 增大而**增大**（pool 并行不分摊波函数以外的内存），这是正常的。

`parse_xml.py` 的输出格式（示意输出，数值仅供参考）：

```text
scf_converged             : True
n_scf_steps               : x
etot_eV                   : -xxx.xxxxxxx
highestOccupiedLevel_eV   : x.xxxxxxx
lowestUnoccupiedLevel_eV  : x.xxxxxxx
pressure_kbar             : x.xxxxxx
max |F_i| (eV/A)          : x.xxxxe-xx
```

`lowestUnoccupiedLevel` 只有在计算中包含空带（设置了 `nbnd`）时才会出现在 XML 中；金属体系给出的是 `fermi_energy`。

## 注意事项

- 进程数必须能被 `-nk` 整除，否则报 `invalid number of pools, parent_nproc /= nproc_pool * npool`。
- 只有 Γ 点（`K_POINTS gamma`）时不能用 `-nk`。
- ASE 3.23 起旧的 `Espresso(command=...)` 写法会直接报错，必须改用 `EspressoProfile`；`input_data` 中的单位保持 QE 原生单位（Ry、Bohr），ASE 不做换算。
- `parse_xml.py` 只依赖标准库；更完整的 XML 读取可以使用官方的 [qeschema](https://github.com/QEF/qeschema) 或 [qe-tools](https://github.com/aiidateam/qe-tools)。
