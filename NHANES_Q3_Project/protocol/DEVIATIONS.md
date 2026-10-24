# Deviations after the V2 analysis plan

The entries below retain the recorded chronology of changes after the local V2 plan. Earlier pre-model decisions are described in [CHANGES_V2.md](CHANGES_V2.md). Event times, scientific rules and the order of events are preserved.

## 2026-10-19T09:39:53+00:00 — 仅调整执行节奏

作者于该记录时间取消固定7个日历日时窗及由此触发的停止／移除模块指令，改为立即、连续完成全部必需产物。现行 analysis_config.json 的 scope.timebox 和 PROTOCOL_RECORD.json 的期限字段置空，并记为 continuous_until_required_deliverables_complete。旧 PREMODEL_ANALYSIS_PLAN.md、PROTOCOL_zh.md、CHANGES_V2.md 和原始启动目录完整保留，所含旧期限为历史计划。科学参数、候选、数据角色、折分和统计方法没有改变。

截至记录时首次研究拟合字段为 None；原包没有模型、预测或执行日志，新工作目录也未拟合研究模型。远端协议记录尚未完成，不能称作已预注册。该修订为执行节奏调整，不是研究设计修订；没有结果重跑。实际合并工作开始时间为 2026-10-19T09:37:19+00:00，此前已读取材料。

## 2026-10-19T09:47:25+00:00 — 作者明确调整远端记录前置条件

作者于该记录时间决定先按本地V2计划开展实验，完成实验后再公开材料。根据该作者决定，取消本轮首次真实数据拟合前必须完成GitHub远端记录的前置条件；保留本地V2原计划及真实执行日志，完成实验后再处理上传。该变化独立于先前“仅调整执行节奏”的记录，不将两者混写。此时 first_study_fit_utc=None，尚未拟合真实研究模型、尚未看任何研究性能。科学参数、16项X、两项Z、周期角色、候选、严格OOF、调查权重、Fieller及敏感性均未改变。

repository_url、protocol_commit_url、remote_recorded_utc维持null。不能把后续上传倒称模型前远端记录，不称正式预注册。此次仅改变记录／发布顺序，后续上传仍需实际指定目标；没有擅自创建仓库或上传数据。
