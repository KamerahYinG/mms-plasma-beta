# 基于FGM、FPI仪器与等离子体β参数的入门级复现
## 概述
本项目是一项入门级复现工作，利用NASA的（MMS）在2015年10月16日观测到的昼侧磁层顶穿越事件，复现事件周边磁场等环境特征。
分析采用MMS1卫星的FGM与FPI的Level-2数据。

---
## 事件信息
时间区间：
UTC 2015-10-16 13:03:30 -- 13:08:00

---
## 数据来源
本项目使用的全部MMS数据为公开科学数据，由MMS科学数据中心与NASA空间物理数据平台（SPDF/CDAWeb）分发。
使用PySPEDAS完成CDF数据的检索、下载、缓存与读取。

---

标量温度按下式计算：
$$T = (T_\parallel + 2T_\perp)/3$$

---
## 主要结果
对比13:05:45 UTC附近剧烈过渡区两侧两个时段：

| 物理量 | A区间：13:05:00–13:05:35 | \(\boldsymbol{B}\)区间：13:05:50–13:06:20 |
| --- | ---: | ---: |
| $|\(\boldsymbol{B}\)|$ | 39.16 nT | 24.88 nT |
| $n_e$ | 0.408 cm⁻³ | 9.53 cm⁻³ |
| $n_i$ | 0.497 cm⁻³ | 10.70 cm⁻³ |
| $T_e$ | 128 eV | 26.0 eV |
| $T_i$ | 2466 eV | 323 eV |
| 热压强 | 0.208 nPa | 0.583 nPa |
| 磁压强 | 0.610 nPa | 0.246 nPa |
| 总β | 0.331 | 2.12 |

A区间特征：密度低、温度高、磁场强、β小于1。
\(\boldsymbol{B}\)区间特征：密度更高、温度更低、磁场更弱、β大于1。
上述特征符合从**磁层等离子体环境**过渡到**磁鞘等离子体环境**的物理图像。

---
## 项目目录结构
```
mms-plasma-beta/
├── README.md
├── AGENTS.md
├── src/
│   ├── load_fgm.py
│   ├── plot_fgm.py
│   ├── load_fpi.py
│   └── compute_beta.py
├── figures/
│   ├── mms1_fgm_overview.png
│   ├── mms1_fpi_overview.png
│   └── mms1_beta_overview.png
├── results/
└── data/
    └── spedas/
```

---
## 相关文献
Le Contel, O., et al. (2016).
*Whistler mode waves and Hall fields detected by MMS during a dayside magnetopause crossing.*
Geophysical Research Letters, 43, 5943–5952.
doi:10.1002/2016GL068968

\(\boldsymbol{B}\)urch, J. L., et al. (2016).
*Electron-scale measurements of magnetic reconnection in space.*
Science, 352(6290), aaf2939.
doi:10.1126/science.aaf2939

Zhao, C., et al. (2016).
*Force balance at the magnetopause determined with MMS: Application to flux transfer events.*
Geophysical Research Letters.
doi:10.1002/2016GL071568
