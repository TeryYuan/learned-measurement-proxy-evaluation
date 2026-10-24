# Q3 的信息边界、测量块决定与跨仪器评价

设计说明日期：2026-10-19。本文补充 V2 协议的解释边界，不改变两项测量、16 项 X、五条路线及原有时间切分。本文件的撰写和核查没有运行任何预测模型。

## 1. 用普通语言说清研究价值

代理 g 拟合以后，只根据一个人的基础信息 X 生成预测。因此，它不能补回这个人尚未被 X 解释的真实测量差异，也不能为这个人创造 X 之外的新信息。如果基础模型已经充分、最优地利用了 X，使用同一 X 生成的代理就不可能进一步降低期望 Brier loss。

本研究实际使用的是有限数据和预定候选模型。基础模型并不保证达到最优，代理生成的表示也可能改变这些有限模型利用 X 的方式。因此，真实测量是否产生收益、固定替换损失多大，以及为代理重新训练是否得到不同表现，都仍需实际评价。正的代理增量或正的保留率表示这些具体程序之间的预测差异，不能解释成代理增加了患者信息。

不把本模块的贡献写成新发现了“确定性变换不增加信息”。本模块提供的是预先规定的真实数据评价：在可信基础模型之上核查实测收益，量化同一已拟合模型的输入替换后果，并区分固定替换与重新训练。

## 2. 简短 Brier 推导

以下期望均针对同一个目标评价分布，包括协议所定义的分析域；所有拟合程序均视为已经固定。调查加权样本估计用于估计相应目标量。这里的期望风险恒等式不要求某一个有限评价样本的点估计严格满足相同不等式。

定义使用 X 的最优预测、实际基础预测和实际固定替换预测：

\[
b^*(X)=E[Y\mid X],\qquad b(X)=f_B(X),\qquad s(X)=f_M(X,g(X)).
\]

对任意只依赖 X 的预测 a，Brier 风险为

\[
L(a)=E[(Y-a(X))^2]
=E[\operatorname{Var}(Y\mid X)]+E[(a(X)-b^*(X))^2].
\]

因此，若基线为最优预测 b*，替换后相对该基线的剩余收益为

\[
H^*=L(b^*)-L(s)=-E[(s-b^*)^2]\leq0.
\]

这给出了确定性代理的信息边界。它不等同于实际有限基础模型的 H。实际比较满足

\[
H=L(b)-L(s)
=E[(b-b^*)^2]-E[(s-b^*)^2],
\]

所以 H 可正、可零、可负。为代理独立训练的 P 同样是 X 的函数，也服从相同的信息边界。

为说明与保留率的关系，进一步考虑真实测量模型也是目标分布中的最优预测：

\[
m^*(X,Z)=E[Y\mid X,Z].
\]

此时最优实测收益为

\[
G^*=L(b^*)-L(m^*)=E[(m^*-b^*)^2].
\]

对任何确定的 X 函数 s，包括固定输入替换所形成的函数，

\[
D^*=L(s)-L(m^*)=G^*+E[(s-b^*)^2].
\]

当 G*>0 时，

\[
R^*_{\mathrm{fixed}}=\frac{H^*}{G^*}
=-\frac{E[(s-b^*)^2]}{G^*}\leq0.
\]

若 G*=0，保留率没有正收益分母，不能借此宣称完整替代。若 Z 本身可完全由 X 确定，则它没有理想基线以上的新信息，也不能产生正的 G*。

实际 B 与 M 不是这里的理想预测，所以该理想推导不能用来预填真实实验的 G、D、H 或 R。实际有限算法的正保留、负保留乃至超过 100% 的点估计均按主协议及区间规则报告。

## 3. Q1 准确度不能单独解释 Q3

不使用“保留率就是 Z 中可由 X 预测部分的占比”这类表述。理想基础预测已经利用 X 所包含的风险信息，而整体测量 R²、MAE 或 RMSE 没有区分哪些预测误差与结局有关。

例如，一项测量的大部分变异可能容易预测，但模型的增量预测能力可能集中在较小的残余变异、极端值或某些非线性区间。代理对这些读数的收缩或偏差可能造成明显替换损失，即使整体预测误差不大。反过来，对模型很少使用的测量方向预测不准，也可能对下游损失影响有限。

因此应把 Q1 和 Q3 并列报告。替换表现的解释要结合实际测量收益、与结局有关的代理误差、固定 M 对输入的使用方式、模型近似和跨周期变化。不能仅根据面板 B 的准确度为 Q3 作充分归因，也不新增结果驱动的残差筛选或代理选型。

## 4. 本轮测量块决定：保留 RDW 和血红蛋白

Z 继续固定为 RDW（LBXRDW，%）和血红蛋白（LBXHGB，g/dL）。不加入白蛋白或肌酐，不同时运行两个测量块后选择实测收益更高的一组，不因 G 的结果更换面板。

白蛋白和肌酐具有合理的既往文献依据。PhenoAge 研究的测量组合包括这两项和 RDW，并考察了死亡结局；但它不能证明新增两项在本协议的 16 项 X、五年结局、分析人群和模型集合上一定有 Brier 增益。该文使用的 NHANES IV 年份与本模块存在重叠，不能称为完全独立评价人群的增量验证。[Liu et al., PLOS Medicine, 2018](https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1002718)。

保留两项的决定基于本模块的有限范围：新增测量可能提高实测收益，也增加完整测量入组条件、代理目标和解释来源，且不能解除确定性代理的信息限制。生化测量在这些周期还涉及另一套仪器和校准变化。对一个计划在一周内完成、可独立从稿件拆出的固定替代模块，保持小型预定测量块更利于完成清楚的程序比较。本决定不声称四项设计不科学；它是本次首次拟合前确定的范围选择。

### 未采用的生化候选：已核查的官方事实

下列内容用于记录为何扩块并非单纯增加两列数据，不能据此对本轮 CBC 数据执行任何生化校准。

| 周期 | 官方说明 | 若另行设计生化研究时的处理边界 |
| --- | --- | --- |
| 2003–2004 | 肌酐与参考方法比较后，CDC 表示无需校正。 | 使用对应公开值；见 [L40_C](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/L40_C.htm)。 |
| 2005–2006 | CDC 强烈建议将肌酐校正为 −0.016＋0.978×公开未校准值，单位为 mg/dL。 | 该校准仅用于规定周期；见 [BIOPRO_D，Analytic Notes](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/BIOPRO_D.htm)。 |
| 2007–2008 | 2007 年使用 LX20，2008 年改为 DxC800。官方文件说明 2007 年非标准化肌酐已经通过交叉实验的 Deming 回归调整。 | 不对已发布调整值重复应用旧校准；见 [BIOPRO_E，仪器及肌酐说明](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2007/DataFiles/BIOPRO_E.htm)。 |
| 2009–2010 | DxC800 的肌酐方法为 IDMS 标准化。 | 使用对应公开值并记录方法；见 [BIOPRO_F，Creatinine](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/BIOPRO_F.htm)。 |

白蛋白在 LX20 和 DxC800 上均采用溴甲酚紫端点方法，详见上述 BIOPRO_D、BIOPRO_E 和 BIOPRO_F。相同化学方法不能独立证明跨仪器完全无偏差。

## 5. 保留时间切分，明确跨 CBC 仪器评价

2003–2004 训练、2005–2006 验证、2003–2006 最终重拟合、2007–2010 后期评价的安排保持不变。

CDC 的 [CBC_E 文档](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2007/DataFiles/CBC_E.htm)明确说明：2005–2006 使用 Beckman Coulter MAXM，2007–2008 使用 HMX；实验室方法和地点未变，仪器基本原理未变。该页前面的通用描述仍提及 MAXM，因此应以专门说明周期仪器变更的段落记录这次变化。方法和原理未变不等同于跨仪器数值完全一致。

本模块据此称为“跨后期调查周期及已记录 CBC 仪器转换的固定替代评价”。五条路线在同一后期人群、同一结局及权重下评价，S 与 MEAN 仍查询同一个已拟合 M。因此，G、D、H 继续是该测量环境下实际预测程序之间的有效比较。

这些量不单独识别仪器转换的影响，也不构成相同测量协议条件下的纯替代效应。它们综合反映预测程序、后期人群以及测量环境。G 属于此情景下测量增强程序相对基础程序的预测收益，不能直接改称真实生物学信息量或测量的因果效应。

SI 按主协议呈现 RDW 和 Hb 各周期的加权分布摘要及分布图。观察到分布差异时，只描述其大小和位置；不能把差异全部归因于仪器，也不能因仪器原理相同便假定没有漂移。分布审计用于限定解释，不用于根据结果自行校正、重选周期、改变主面板或追加同仪器交叉验证。

## 6. 供英文 Discussion 使用的文字

以下段落不预填任何真实结果，可在结果完成后按上下文作语言衔接：

> After fitting, the proxy was a deterministic function of the baseline variables and therefore added no participant-level information beyond them. An optimal predictor using those variables alone could not be improved by fixed substitution under expected Brier loss. Our analysis instead compared prespecified prediction procedures fitted to finite data. The empirical quantities were the measured-input gain and the changes produced by substitution and separate proxy training. Their interpretation depends on outcome-relevant proxy error, model approximation, and transport across survey cycles and documented measurement changes. Overall proxy accuracy alone does not determine these quantities.

测量迁移限制可在相邻位置另列一段，并引用 CBC_E 官方文档：

> Evaluation used later NHANES cycles across a documented transition from the MAXM to the HMX CBC instrument. The comparisons therefore describe performance in that later-cycle measurement setting. They do not isolate instrument effects or establish substitution performance under an unchanged measurement protocol. Cycle-specific distributions provide a descriptive audit of the observed measurements; differences cannot be attributed solely to instrumentation.

上述文字需要与主协议的公开死亡文件结局定义、调查域和固定模型不确定性边界一并使用。真实结果是否建立有用测量前提，按实际输出报告。
