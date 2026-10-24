"""Record only the authorized cadence amendment, preserving the V2 snapshot."""
from runtime import ROOT, OUT, read_json, write_json, now, event

def main():
    rec = read_json(ROOT/'protocol/PROTOCOL_RECORD.json')
    if rec.get('cadence_amendment_utc'):
        print('Existing cadence amendment retained; no duplicate entry.')
        return
    stamp = now()
    # Merge into the new local execution root began at the observed tool time.
    rec['execution_started_utc'] = rec.get('execution_started_utc') or '2026-10-19T09:37:19+00:00'
    rec['execution_start_basis'] = 'Observed UTC immediately before merging the five supplied extracted bundles; initial reading preceded this record.'
    historical = {k: rec.get(k) for k in ['timebox_calendar_days','execution_deadline_utc']}
    rec.update(timebox_calendar_days=None, execution_deadline_utc=None,
               cadence_amendment_utc=stamp, execution_mode='continuous_until_required_deliverables_complete',
               historical_timebox=historical)
    cfg = read_json(ROOT/'protocol/analysis_config.json')
    old = cfg['scope']['timebox']
    cfg['scope']['timebox'] = dict(calendar_days=None, starts='recorded actual execution start',
        includes=old['includes'], outcome_independent=True, on_timeout=None,
        mode='continuous_until_required_deliverables_complete', amended_utc=stamp)
    write_json(ROOT/'protocol/analysis_config.json', cfg)
    write_json(ROOT/'protocol/PROTOCOL_RECORD.json', rec)
    fitted = rec.get('first_study_fit_utc')
    text = f'''\n## {stamp} — 仅调整执行节奏\n\n作者于该记录时间取消固定7个日历日时窗及由此触发的停止／移除模块指令，改为立即、连续完成全部必需产物。现行 analysis_config.json 的 scope.timebox 和 PROTOCOL_RECORD.json 的期限字段置空，并记为 continuous_until_required_deliverables_complete。旧 PREMODEL_ANALYSIS_PLAN.md、PROTOCOL_zh.md、CHANGES_V2.md 和原始启动目录完整保留，所含旧期限为历史计划。科学参数、候选、数据角色、折分和统计方法没有改变。\n\n截至记录时首次研究拟合字段为 {fitted!r}；原包没有模型、预测或执行日志，新工作目录也未拟合研究模型。远端协议记录尚未完成，不能称作已预注册。该修订为执行节奏调整，不是研究设计修订；没有结果重跑。实际合并工作开始时间为 {rec['execution_started_utc']}，此前已读取材料。\n'''
    with (ROOT/'protocol/DEVIATIONS.md').open('a', encoding='utf-8') as f:
        f.write(text)
    write_json(ROOT/'protocol/EXECUTION_CADENCE_AMENDMENT.json', dict(
        amended_utc=stamp, authority='author decision recorded on 2026-10-19',
        historical_timebox=old, current_timebox=cfg['scope']['timebox'],
        scientific_parameters_changed=False, first_study_fit_utc_at_amendment=fitted,
        prior_study_models_found=False))
    event('execution_initialized', sources='five existing extracted V2 bundles', cadence_amendment_utc=stamp, study_fits=0)
    print(text)

if __name__ == '__main__':
    main()
