# 公开死亡链接文件、调查权重与分析前记录

版本：V2，2026-10-19。适用于 NHANES 2003–2010 RDW／血红蛋白代理替代模块。

本文件规定死亡时间限制的披露方式、调查权重的使用方法及分析前记录的措辞。它是分析协议的支持文件；文件日期不构成 OSF 登记回执，也不证明已经向 GitHub 或其他网站提交。本文件不包含模型结果。

## 1. 公开 LMF 能识别什么，不能识别什么

使用 2019 public-use Linked Mortality File（LMF）。CDC 说明，部分记录的随访时间或根本死因使用了合成替换值；生死状态 `MORTSTAT` 未受这一数据扰动。这里的 `MORTSTAT=0/1` 仍分别表示链接程序判定的 assumed alive／assumed deceased，不等于已经独立核实每个人的生死状态。

已核对 2019 公开数据字典列出的全部 13 个变量。字典没有记录级“随访时间被合成替换”标志、未经替换的原始随访时间或每条记录的替换幅度。因此，无法仅从公开文件识别或计数本分析中被替换随访时间的记录，也无法知道五年标签因替换而改变的准确人数。这个结论来自完整字段清单，不能写成使用了某个不存在的官方标记。

主结局继续称为 **public-use-file-defined five-year all-cause mortality**，使用从 MEC 查体起算的 `PERMTH_EXM`，而非 `PERMTH_INT`。五年内死亡者属于事件；不得因其随访不足 60 个月而删除。恰好在 60 个月死亡者属于五年事件，假定存活且随访恰好 60 个月者可以判为五年非事件。完整标签规则以主协议及整理脚本为准。

官方来源：

- [2019 public-use LMF 文件说明](https://www.cdc.gov/nchs/data/datalinkage/public-use-linked-mortality-file-description.pdf)
- [2019 public-use LMF 完整数据字典](https://www.cdc.gov/nchs/data/datalinkage/public-use-linked-mortality-files-data-dictionary.pdf)
- [CDC 死亡链接数据页面](https://www.cdc.gov/nchs/linked-data/mortality-files/index.html)

## 2. SI 的 60 个月边界描述

### 2.1 描述人群与分组

在满足年龄、查体及调查设计、死亡链接资格、两项 CBC 真实测量、MORTSTAT为0或1，以及PERMTH_EXM有限且不小于0的条件下，**先保留尚未因五年结局可判定性而排除的记录**，再描述公开随访时间。按四个调查周期分别列出，并汇总为开发池 2003–2006 和评价池 2007–2010。该表不使用模型预测。

对于公开 `PERMTH_EXM` 有效的记录，固定使用三个互不重叠的整数月份区间：54–59、恰好 60、61–66。54–66 月窗口是在模型拟合前选定的描述窗口，不是 CDC 提供的扰动范围。

各周期及汇总池分别报告：

| 描述项 | 固定窗口 |
| --- | --- |
| 公开死亡状态为1的记录 | 分别报告54–59月、恰60月、61–66月；另给54–66月合计 |
| 公开死亡状态为0的记录 | 报告54–66月合计 |
| 五年状态不可判定记录 | 报告54–66月合计；可能与存活记录重叠 |

以上首先报告未加权记录数，并给出该周期／池的描述人群分母。补充加权比例使用第3节的当前池权重，分母为同一前置条件人群。未知状态不是可相加的第三种生死状态。无效死亡状态／随访时间的排除已在原始人群流程中处理；不为边界表扩大人群或增加新的计数分支。

本版已完成上述描述，实际结果见diagnostics。54–66月存活与五年状态未知计数均为实际计算值；不能仅凭随访期较长的预期填零。

### 2.2 解释边界

该表描述公开文件中的时间分布，**不识别被合成替换的个体，不估计被替换人数，也不提供误分类人数的上界**。不得把窗口以外的记录视为未受扰动，或因边界附近人数较少而断言扰动可以忽略。

不通过比较 2015 与 2019 LMF 来推断哪些个体被扰动。官方说明指出，两版还涉及链接算法及随访的变化。也不为本项描述增加删除边界记录、搜索其他时间窗或反复重训的分析分支。

### 2.3 可用于 SI 的英文措辞

> The public-use linked mortality file does not provide a record-level indicator identifying follow-up times replaced with synthetic values. We therefore could not determine the number of affected records or the number of five-year labels changed by this procedure. We counted publicly recorded deaths at 54–59, exactly 60, and 61–66 months by survey cycle and analysis period, and reported 54–66-month totals for deaths, survivors and records with unascertainable five-year status. The latter descriptive counts may overlap. These counts describe the released data around the decision boundary; they do not identify perturbed records or bound the number of misclassified outcomes.

## 3. 调查权重按当前分析池定义

原始字段 `WTMEC2YR` 始终保留。多周期权重按当前分析池包含的两年调查周期数构造，不能给所有用途统一贴一个含义不清的 `weight` 标签。

| 当前用途 | 调查周期 | 合并周期数 | 调查权重 |
|---|---|---:|---|
| 候选模型训练 | 2003–2004 | 1 | `WTMEC2YR` |
| 下游模型选择验证 | 2005–2006 | 1 | `WTMEC2YR` |
| 开发池重拟合及其内部交叉拟合 | 2003–2006 | 2 | `WTMEC2YR / 2` |
| 最终后期评价 | 2007–2010 | 2 | `WTMEC2YR / 2` |
| 单周期分布描述 | 每个单独周期 | 1 | `WTMEC2YR` |
| 如展示全部周期合并描述 | 2003–2010 | 4 | `WTMEC2YR / 4` |

交叉验证子集继承其所属分析池的权重定义，随后按当次拟合的训练成员进行均值归一化。周期分组变量、权重及设计字段只用于切分、拟合加权和统计评价，不能作为 X 特征输入模型。

### 3.1 拟合、模型选择和最终评价均使用权重

- B、M、P 的所有候选拟合及最终拟合使用当前训练池的调查权重。
- 代理 g 的拟合和候选选择也使用权重。代理选择依据规定的测量预测误差，不依据死亡结局的预测收益或保留率。
- 为保持预设惩罚参数的尺度一致，每次模型拟合前，将当次训练成员的权重除以其均值，使训练权重均值为 1。交叉拟合的训练成员、预处理及权重归一化不得使用相应留出成员的信息。
- 下游模型选择采用验证池的加权 Brier loss；最终评价采用评价池的调查权重。加权平均的分母为参与该项计算的有效权重总和。
- 不叠加类别平衡权重，不实施 SMOTE，不根据结果重新标定调查权重。

同一池所有权重乘以同一个正常数，会在加权平均与“训练权重均值为 1”的归一化中约掉。因此，统一缩放错误不一定改变 Brier 点估计，但会造成权重含义和加权人数错误；各池仍须使用并记录上表定义。

官方依据：[NHANES Weighting Tutorial](https://wwwn.cdc.gov/nchs/nhanes/tutorials/weighting.aspx)。CDC 在该教程中明确给出 2003–2006 和 2007–2010 使用 `WTMEC2YR / 2`，以及 2003–2010 使用 `WTMEC2YR / 4` 的公式。

## 4. 分析域、方差与代表性边界

### 4.1 保留调查设计，再定义分析域

最终评价的完整设计框架保留 2007–2010 年所有正 MEC 权重且设计字段有效的查体参与者。在此框架内，用年龄、链接资格、两项 CBC 实测和五年结局可判定性定义分析域。进行年龄 40–79 岁的有限描述时，也在同一完整设计框架中定义相应域。

成对损失及其差异的 Taylor 线性化计算保留完整框架中的 PSU 和分层；域外记录的线性化贡献明确设为零。不得先按年龄、CBC 完整性或链接资格删成小表，再从该小表重建调查设计。若实现保留了四周期的完整框架，最终评价中不属于 2007–2010 的完整周期不贡献任何线性化量，且不能计入评价域的自由度。

自由度按当前分析域中实际有成员的 PSU 数减去相应分层数计算。每层 PSU 数使用数据中的实际值，不硬编码为 2。五条路线使用相同人群、权重和联合协方差；具体公式及 Fieller 集合见 `STATISTICS.md`。

官方依据：[NHANES Variance Estimation Tutorial](https://wwwn.cdc.gov/nchs/nhanes/tutorials/varianceestimation.aspx)，特别是 Analyzing Subgroups 和 Degrees of Freedom when Analyzing Subgroups 两节。

### 4.2 权重不会自动消除分析人群选择

本模块的目标应表述为：**调查加权、具备死亡链接资格且两项 CBC 实测可用的规定分析人群中的预测表现**。不得仅因使用 MEC 权重，就声称结果无条件代表全部美国 40 岁以上成人。

官方 LMF 文件说明指出，死亡链接不完整时，原调查权重的性质不确定；NHANES 未像 NHIS 一样附带链接资格调整权重。两项 CBC 的缺失也会造成额外的分析人群选择。本研究报告各步排除人数与缺失情况，并将其作为外推限制；不默认原权重已纠正上述选择。

本次不新增链接资格／CBC 缺失的逆概率加权分支。模型性能的调查设计区间仍以实际分析域、固定已拟合模型和公开数据为条件，不覆盖所有训练不确定性、链接误差或公开随访时间扰动。

官方依据：[2019 public-use LMF 文件说明](https://www.cdc.gov/nchs/data/datalinkage/public-use-linked-mortality-file-description.pdf) 的 Sample Weights 节。

## 5. 分析前记录与正式登记的区别

在任何模型拟合前，保留一份完整分析计划及明确日期；首次真实模型拟合记录实际开始时间。后续必要改动以简短偏离记录追加，说明修改理由和当时是否已经看过模型表现。这个模块的分析前记录只覆盖这个模块，不能追溯改变既有 DILI／BSEP／ALT 分析的时间性质。

已经完成的可行性工作应如实列出，包括官方来源和字段核查、人群整理、样本及事件人数、跨周期测量分布、公开随访时间边界描述。不得写成未曾接触数据。模型运行前的分布描述可以支持可行性核查，但不能用于事后隐藏或改换已经确定的测量块、时间窗或评价集。

GitHub 的版本提交可以提供内容追溯。Git 官方文档同时说明，作者与提交者日期可以设置，因此 Git 提交中的日期本身不等于独立、不可修改的登记时间证明。若尚未实际推送，不得声称协议已经公开；若没有登记记录，不得写成已经 preregistered。

如果作者选择正式登记，可采用 OSF 的 Secondary Data Preregistration 模板，并披露已经查看的描述信息。OSF 已提交登记的附件不能修改；只有实际提交、确认记录状态并取得链接后，稿件才加入该登记的地址和准确日期。不要求反复冻结、反复哈希或增加多阶段审批。

官方来源：

- [Git commit documentation](https://git-scm.com/docs/git-commit)
- [OSF Registrations & Preregistrations](https://help.osf.io/article/330-welcome-to-registrations)

### 可用于 Methods 的英文措辞

以下句子须与实际记录一致。跨周期分布和边界描述完成后、首次模型拟合开始前，可使用：

> The analysis plan for this additional module was specified before fitting prediction models or evaluating their performance. Prior feasibility work included source and variable checks, cohort construction, sample and event counts, distributions of the measured inputs across survey cycles, and descriptive counts of publicly released follow-up times around the 60-month boundary. Subsequent departures from the plan were documented with their timing and rationale.

如果最终没有发生偏离，最后一句可改为：

> No departures from the specified analysis plan were made.

上述措辞记录实际分析顺序，不宣称整个项目已经正式预注册。只有在对应事实发生后，才增加 GitHub 公开版本或 OSF 登记的引用。
