# CBC 仪器与跨周期测量可比性审核

审核日期：2026-10-19。依据：CDC/NCHS 官方周期代码本及对应实验方法手册。

## 1. 审核结论与分析决定

**2005–2006 到 2007–2008 年发生了 Beckman Coulter MAXM → HMX 仪器转换，恰好位于本模块开发集和后期评价集之间。** 2007–2008 年官方代码本专门记载这一变化，并说明仪器基本原理没有变化；不能据此进一步声称测量数值完全等效。

本模块继续使用既定时间切分：2003–2004 年训练、2005–2006 年验证选型、2003–2006 年开发集重拟合、2007–2010 年后期评价。将此评价明确描述为跨调查周期且跨仪器型号的固定模型替代评价。**不追加同仪器模型分支，不调整时间切分，不根据分布图重新选择周期。** 本审核仅补充测量背景、描述性分布审计和结论边界。

限制到 2009–2010 年可以避开 CDC 后续记载的 2012 年 HMX → DxH 转换，但不能写成“2003–2010 年没有仪器转换”或“开发与评价使用完全相同的仪器”。

## 2. 四周期逐项证据

| 周期与文件 | 官方代码本 | 对应方法手册 | 本模块采用的判断 |
| --- | --- | --- | --- |
| 2003–2004，L25_C | 描述 MAXM，并说明相对前两年设备、方法及地点无变化 | 第 3 页介绍 MAXM；校准见第 7 页 | 当周期使用 MAXM；该说明不证明以后周期无变化 |
| 2005–2006，CBC_D | 描述 MAXM，并说明相对前两年设备、方法及地点无变化 | 第 3 页介绍 MAXM；校准见第 7 页 | 与 2003–2004 年型号证据一致，支持开发周期为 MAXM |
| 2007–2008，CBC_E | 专门说明由 2005–2006 年 MAXM 改为 2007–2008 年 HMX，基本原理无变化；同页通用段落仍称 MAXM | 第 3 页明确要求在 HMX 上完成 CBC；第 14–15 页说明校准 | 采用专门变更说明与当周期手册共同支持的 HMX；记录网页内部冲突 |
| 2009–2010，CBC_F | 通用段落仍称 MAXM；只说相对前两年方法与地点无变化，没有声称设备无变化 | 第 3 页明确要求使用 HMX；第 14–15 页说明校准 | 方法手册支持 HMX；不能仅凭通用网页单句认定为 MAXM |

页码为 PDF 文件的顺序页码，从 1 开始。网页名称不一致的原因没有被上述材料明确解释，本模块不将“模板残留”写成已得到 CDC 确认的事实。

### 官方原始链接

**2003–2004**

- 代码本：<https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/L25_C.htm>
- CBC 方法手册：<https://wwwn.cdc.gov/nchs/data/nhanes/public/2003/labmethods/l25_c_met_complete_blood_count.pdf>

**2005–2006**

- 代码本：<https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/CBC_D.htm>
- CBC 方法手册：<https://wwwn.cdc.gov/nchs/data/nhanes/public/2005/labmethods/cbc_d_met.pdf>

**2007–2008**

- 代码本：<https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2007/DataFiles/CBC_E.htm>
- CBC 方法手册：<https://wwwn.cdc.gov/nchs/data/nhanes/public/2007/labmethods/cbc_e_met.pdf>
- 关键位置：代码本的 “Description of Laboratory Methodology” 中专门说明仪器变化的段落。

**2009–2010**

- 代码本：<https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/CBC_F.htm>
- CBC 方法手册：<https://wwwn.cdc.gov/nchs/data/nhanes/public/2009/labmethods/CBC_F_met_HE.pdf>

## 3. 校准证据及其边界

2003、2005 年 MAXM 手册与 2007、2009 年 HMX 手册均描述 S-CAL 常规校准程序。Hb 属于被校准的项目；RDW 来源于红细胞体积分布，手册说明其不作为单独参数校准。**RDW 不单独校准不等于 RDW 不受仪器和测量条件影响。** 常规仪器校准也不等于已经完成不同型号之间的测量等效验证。

在本次查阅的官方代码本和方法手册中，没有找到可直接用于本模块 MAXM 时期与 HMX 时期 RDW/Hb 的经验证跨型号转换方程。故本模块不自行设定校正系数，不借用其他分析物或后期仪器的转换公式，也不使用后期评价集的测量分布拟合标准化或分布对齐规则。

仅凭这组观察性调查数据，不能把某个跨周期数值差异唯一分解为仪器因素、人群组成变化、真实健康趋势或其他测量条件变化。五条路线仍在同一评价人群中比较，但这并不能消除开发与评价之间的测量条件变化；M 与 S 对该变化的敏感性可能不同。因此，不把观察到的替换损失全部归因于代理遗漏了某种生物学信息。

## 4. 2005–2011 年早期趋势与后续仪器转换

2011–2012 年 CBC_G 的 “Analytic Notes” 指出，2005–2006 年以后 PLT、HCT、RBC、MCHC 与 WBC 出现跨总体、性别或年龄组的趋势，可能涉及方法问题；其中还讨论了 EDTA 抗凝剂变化。RDW 和 Hb 没有被列入该段对应的异常项目，但这不构成这两项指标完全稳定的证据。

- 2011–2012 年代码本：<https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/CBC_G.htm>

2013–2014 年 CBC_H 的特别分析说明另行记载 2012 年 HMX → DxH 转换以及 RDW 等指标的变化，并说明新鲜全血检测条件限制了利用历史保存样本补做方法比较的可能性。

- 2013–2014 年代码本：<https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/CBC_H.htm>

另一个文档边界是：2007、2009 年 CBC 方法手册的采血段落写 K3 EDTA，耗材清单则写 K2 EDTA。本模块记录这种不一致，不据此为每个周期或参与者指定未经确认的抗凝剂转换时间，也不推断 RDW/Hb 的具体校正幅度。

## 5. 描述性分布审计

使用本模块已经规定的分析人群，按四个调查周期分别展示 RDW 和 Hb 的原始单位分布；不读取或利用任何预测模型结果来决定审计形式。

应提供：

- 每周期分析人数及两项测量的有效人数；
- 使用该周期 MEC 权重的加权均值、加权 SD、P5、P25、P50、P75 和 P95；
- 每项测量一幅四周期加权经验累积分布图，使用同一原始单位横轴；
- 可直接用于 SI 的分布表、图及图注。

各周期的相对权重与其原始 WTMEC2YR 一致；同周期内乘以共同常数不会改变上述标准化加权描述值。开发集及后期评价集合并权重、模型拟合权重和主区间计算仍按主协议执行。

分布图只描述测量在实际分析人群中的变化。分布接近不作为仪器等效的通过标准，分布差异也不作为自动剔除周期、调整测量值或更换模型的依据。此处不增加按周期重新训练、无监督域适配、同仪器敏感性拟合或额外模型选型。

## 6. 英文 SI 段落

### Instrument documentation and cycle-specific measurement distributions

CDC documented a change from the Beckman Coulter MAXM in 2005–2006 to the HMX in 2007–2008, coinciding with our development–evaluation boundary. The 2007–2010 laboratory procedure manuals specify HMX instruments, although generic passages in the online codebooks retain the MAXM name. CDC described the basic measurement principles as unchanged; we did not identify a validated RDW or hemoglobin conversion equation between these instrument periods in the examined documentation. We retained the specified temporal split and report cycle-specific weighted distributions. The later-cycle evaluation therefore includes possible changes in measurement conditions and population composition. These descriptive distributions cannot establish instrument equivalence or identify the causes of any observed shift. We did not derive a correction from the evaluation data or add an alternative split based on these distributions.

### Interpretation of fixed-model substitution

The measured-input and substitution routes were evaluated in the same participants with the same fitted downstream model. However, this paired comparison does not isolate the contribution of changes in measurement conditions between development and evaluation. The observed substitution loss should therefore be interpreted as the performance difference under the specified later-cycle conditions, rather than attributed exclusively to unrecovered biological information.

正式 SI 排版时，按稿件引用体例为上述段落添加本文件列出的代码本与方法手册引用，并将分布图表编号替换为实际编号。本审核不包含模型结果。
