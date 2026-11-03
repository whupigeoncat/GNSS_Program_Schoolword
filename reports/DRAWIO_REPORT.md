# 流程图生成报告（DrawIO + Archify）

## Archify 交互版（HTML）

### 函数级数据流程（当前推荐）

| 文件 | 类型 | 来源依据 | 用途 | 状态 |
| --- | --- | --- | --- | --- |
| `figures/rtk-fnflow.json`（规格） | dataflow 图 | `Data_IO.cpp`、`Decoding.cpp`、`Coord_trans.cpp`、`Time_trans.cpp` | 函数级链路：fread()/recv() → FeedNovData() → PumpStream()（DecodeNovStream() 循环）→ ParseNovMsg() → crc32() → MsgId 分发 → decode_bestpos()(BLHToXYZ())/decode_gpsephem()/decode_bdsephem()/decode_rangeb_oem7() → GetSynGNSSObs()(PairEpoch()/GetDiffTime()) → RAWDAT → 下游链（FormSDEpochObs→…→lambda，未实现） | 已交付 |
| `figures/RTK_函数级数据流程.html` | dataflow 交互图 | 同上 | 每个处理步骤标注对应函数名；卡片列出 WaitReadable()、DecodeNovOem7Dat()、BLHToXYZ() 等 | 已交付 |

- 校验：showcase 9/9 通过，composition 0 错误 0 警告。
- 交付回执：spec sha256 `0344ED4A…`、artifact sha256 `9A450217…`（823 KB）。
- 浏览器证据：`visual-check` 四档分辨率全部通过（0 失败，scrollHeight 均等于视口高度）。
- 视觉复核：已查看 1440×900 浅色截图，函数链完整、无遮挡。
- 未画函数说明：`InitNovDecoder()`（当前未被调用）、`ReadRTKConfigInfo()`（占位未调用）不在数据通路上；SPP/SPV 相关函数属于另一条流水线。

### 数据全流程版（二进制 → RTK）

| 文件 | 类型 | 来源依据 | 用途 | 状态 |
| --- | --- | --- | --- | --- |
| `figures/rtk-dataflow-byte2rtk.json`（规格） | dataflow 图 | `Decoding.cpp`、`Data_IO.cpp` | 按"数据形态"从二进制字节讲起：原始字节 → 字节缓冲(Buff[40960]/Len) → 帧(0xAA4412/28B头/CRC32) → EPOCHOBS/PPRESULT/星历 → 配对(RAWDAT) → RTK | 已交付 |
| `figures/RTK_数据全流程_二进制到RTK.html` | dataflow 交互图 | 同上 | 节点标注数据载体与字段（Time、P/L/D、44B/条、L1/L2、B1/B3 等） | 已交付 |

- 校验：showcase 9/9 通过，composition 0 错误 0 警告。
- 交付回执：spec sha256 `955A88A5…`、artifact sha256 `8203CB2E…`（820 KB）。
- 浏览器证据：`visual-check` 四档分辨率全部通过（scrollHeight 均等于视口高度，0 失败 / 17 通过）。
- 视觉复核：已查看 1440×900 浅色截图，字节→帧→结构体→配对→RTK 数据链完整、无遮挡。

### 详细版

| 文件 | 类型 | 来源依据 | 用途 | 状态 |
| --- | --- | --- | --- | --- |
| `figures/rtk-dataflow-detail.json`（规格） | dataflow 图 | `main.cpp`、`Data_IO.cpp`、`Decoding.cpp`、`RTK_Structs.h` | 生成源（冻结的规格） | 交付用规格 |
| `figures/RTK_数据处理详细流程图.html` | dataflow 交互图 | 同上 | 5 阶段 12 节点；节点含变量注释（nSyncLoss / nBadCrc / nRange / nPos / nGeph / nBeph / BasPosValid 等） | 已交付 |

- 校验：`archify validate dataflow ... --quality showcase`，9/9 项检查通过，composition 0 错误 0 警告。
- 交付回执：spec sha256 `FD8D58C3…`、artifact sha256 `AEA93B36…`（820 KB）。
- 浏览器证据：`visual-check` 全部通过（1440×900 / 1600×1000 / 1920×1080 / 2048×1320，含深色主题）。
- 视觉复核：已查看 1440×900 浅色截图，帧解析链、解码分派、配对汇聚与下游虚线接口均正常。
- 侧车文件：`figures/RTK_数据处理详细流程图.visual-check.json/.html` 与 4 张 `.visual-check.*.png`。

### 变量跟踪版（workflow，两页）

| 文件 | 类型 | 来源依据 | 用途 | 状态 |
| --- | --- | --- | --- | --- |
| `figures/rtk-decode-trace.workflow.json` → `figures/RTK_变量跟踪_解码链.html` | workflow 图 | `Decoding.cpp`（`ParseNovMsg`、`FeedNovData`） | 逐步跟踪 Buff/Len/i/MsgId/MsgLen/frameLen/CRC 与统计计数；含 i++ 重试与循环取帧回退 | 已交付 |
| `figures/rtk-pairing-trace.workflow.json` → `figures/RTK_变量跟踪_配对输出.html` | workflow 图 | `Data_IO.cpp`（`PumpStream`、`GetSynGNSSObs`、`PairEpoch`） | 跟踪 dt、阈值、RovEpk/BasEpk、BasPos/BasPosValid、RAWDAT；含丢弃/追基站/EOF 终止分支 | 已交付 |

- 校验：两图 `archify validate workflow ... --quality showcase` 均通过，composition 0 错误 0 警告。
- 交付回执：解码链 spec `3E721201…` / artifact `AD861F65…`；配对输出 spec `CF379489…` / artifact `DC1B5A47…`。
- 浏览器证据：两图 `visual-check` 全分辨率通过。
- 视觉复核：已查看 1440×900 浅色截图；图例已本地化（流缓冲泵浦 / 帧解析步骤 / 解码产物 等）。

### 概览版

| 文件 | 类型 | 来源依据 | 用途 | 状态 |
| --- | --- | --- | --- | --- |
| `figures/rtk-dataflow-current.json`（规格） | dataflow 图 | `main.cpp`、`Data_IO.cpp`、`Decoding.cpp`、`RTK_Structs.h` | 生成源（冻结的规格） | 交付用规格 |
| `figures/RTK_数据处理流程图.html` | dataflow 交互图 | 同上 | 深/浅主题、缩放、聚焦、导出；共 5 阶段 11 节点 | 已交付 |

- 校验：`archify validate dataflow ... --quality showcase`，9/9 项检查通过，composition 0 错误 0 警告。
- 交付回执：spec sha256 `D128EAEF…930E`（5768 B）、artifact sha256 `F8FF2692…EDD9`（818 KB）。
- 浏览器证据：`visual-check` 全部通过（1440×900 / 1600×1000 / 1920×1080 / 2048×1320，含深色主题）。
- 视觉复核：已查看 1440×900 浅色截图，双流路径、配对回退循环、下游虚线接口与标注均正常。
- 侧车文件：`figures/RTK_数据处理流程图.visual-check.json/.html` 与 4 张 `.visual-check.*.png`。

## DrawIO 版本

### 图示清单

| 文件 | 类型 | 来源依据 | 用途 | 状态 |
| --- | --- | --- | --- | --- |
| `figures/fig_rtk_dataflow.drawio` / `.pdf` / `.png` | RTK 数据处理总流程图 | `main.cpp`、`Data_IO.cpp`、`Decoding.cpp`、`RTK_Structs.h` | 说明双流输入 → 独立缓冲 → 共用解码 → 时间配对 → RAWDAT → 下游解算 | 已导出 PDF / PNG |
| `figures/fig_rtk_decode.drawio` / `.pdf` / `.png` | ParseNovMsg 帧解码流程图 | `Decoding.cpp`（`ParseNovMsg`、`FeedNovData`） | 说明帧同步、长度校验、CRC32、MsgId 分发与失败重同步逻辑 | 已导出 PDF / PNG |

### 未生成图示及原因

- 时间配对单独详图：配对循环（`GetSynGNSSObs` / `PairEpoch` 的两个回退分支）已并入总流程图，不再单独出图。
- 数据型图表：不属于本阶段职责。

### 导出与自检记录

- 导出命令：`drawio --export --format pdf|png --crop`，4 次导出全部成功（exit 0）。
- 自检：节点无重叠；回退/循环箭头不穿过核心节点；同类节点样式、字号、配色统一；两图内容与当前代码一致（下游解算链以虚线标注“已声明、未实现”）。

### 给论文阶段的嵌入建议

- `fig_rtk_dataflow`：建议放入“RTK 数据处理流程”一节；caption 例：“图 X RTK 双流数据处理总体流程”。
- `fig_rtk_decode`：建议放入解码实现相关小节；caption 例：“图 X ParseNovMsg 帧解码流程”。
- 引用时保留图中“下游解算未实现 / 未接入主流程”的标注，避免与代码实现状态不符。
