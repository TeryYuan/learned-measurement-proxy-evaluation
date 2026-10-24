#!/usr/bin/env python3
"""Prespecified measurement distributions and public follow-up boundary counts.
No model fitting, predictor-outcome screening, reweighting, or recalibration.
"""
from pathlib import Path
from datetime import datetime, timezone
import json
# The report-only branch reuses the already saved SVG curve geometry and never
# reads the cohort or recomputes distribution/boundary statistics.
if __name__ == '__main__':
    import sys
    if '--report-only' in sys.argv or '--outputs-only' in sys.argv:
        from render_saved_figures import main as render_saved_main
        render_saved_main()
        raise SystemExit(0)

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
CYCLES = ['2003-2004', '2005-2006', '2007-2008', '2009-2010']
Z = {'rdw_pct': 'RDW (%)', 'hemoglobin_g_dl': 'Hemoglobin (g/dL)'}
COLORS = ['#355C7D', '#68A4A5', '#C17C36', '#9D4B59']


def weighted_quantiles(values, weights, probabilities):
    order = np.argsort(values, kind='stable')
    x, w = np.asarray(values)[order], np.asarray(weights)[order]
    cdf = np.cumsum(w) / np.sum(w)
    return x[np.searchsorted(cdf, probabilities, side='left')]


def main():
    out = ROOT / 'diagnostics'
    out.mkdir(exist_ok=True)
    cohort = pd.read_csv(ROOT / 'data/analysis_cohort.csv')
    frame = pd.read_csv(ROOT / 'data/design_frame.csv')
    rows = []
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.7))
    for ax, (variable, label), panel in zip(axes, Z.items(), ['A', 'B']):
        for cycle, color in zip(CYCLES, COLORS):
            c = cohort.loc[cohort.cycle.eq(cycle)]
            x = c[variable].to_numpy(float)
            w = c['WTMEC2YR'].to_numpy(float)  # Single-cycle descriptive weights.
            assert np.isfinite(x).all() and np.isfinite(w).all() and (w > 0).all()
            mu = np.average(x, weights=w)
            q = weighted_quantiles(x, w, [.05, .25, .5, .75, .95])
            rows.append(dict(cycle=cycle, measurement=variable, unit=label, n=len(c),
                             instrument='MAXM' if cycle in CYCLES[:2] else 'HMX',
                             weighted_mean=float(mu),
                             weighted_SD=float(np.sqrt(np.average((x-mu)**2, weights=w))),
                             min=float(x.min()), p05=float(q[0]), p25=float(q[1]),
                             median=float(q[2]), p75=float(q[3]), p95=float(q[4]), max=float(x.max())))
            order = np.argsort(x, kind='stable')
            xs, ws = x[order], w[order]
            ax.step(xs, np.cumsum(ws)/ws.sum(), where='post', lw=1.65, color=color,
                    label=f'{cycle} (n={len(c):,})')
        ax.set(xlabel=label, ylabel='Survey-weighted cumulative proportion', ylim=(0, 1.01))
        ax.set_title(panel, loc='left', weight='bold', fontsize=12)
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(axis='y', alpha=.18)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=2, frameon=False, fontsize=9,
               bbox_to_anchor=(.5, .06))
    fig.suptitle('CBC distributions in the prespecified analytic cohort', fontsize=13)
    fig.text(.5, .018, '2003-2006: MAXM | 2007-2010: HMX. Descriptive curves do not identify instrument effects.',
             ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .18, 1, .95))
    for extension in ['png', 'svg', 'pdf']:
        fig.savefig(out / f'cbc_distributions_by_cycle.{extension}', dpi=220)
    plt.close(fig)
    pd.DataFrame(rows).to_csv(out / 'cbc_distributions_by_cycle.csv', index=False)

    # Prerequisites preceding the five-year ascertainability rule. The full
    # frame already requires a MEC examination and valid positive design weight.
    t = frame['PERMTH_EXM']
    prerequisite = (frame['RIDAGEYR'].ge(40) & frame['ELIGSTAT'].eq(1)
                    & frame['MORTSTAT'].isin([0,1]) & np.isfinite(t) & t.ge(0))
    for variable in Z:
        prerequisite &= np.isfinite(frame[variable]) & frame[variable].gt(0)
    c = frame.loc[prerequisite].copy()
    scopes = [(cycle, c.cycle.eq(cycle)) for cycle in CYCLES]
    scopes += [('development', c.period.eq('development')), ('test', c.period.eq('test')),
               ('all_2003_2010', pd.Series(True, index=c.index))]
    boundary_rows = []
    for name, mask in scopes:
        s = c.loc[mask]
        t = s['PERMTH_EXM']
        died = s['MORTSTAT'].eq(1)
        alive = s['MORTSTAT'].eq(0)
        windows = {
            'death_54_59': died & t.between(54,59),
            'death_exact_60': died & t.eq(60),
            'death_61_66': died & t.between(61,66),
            'death_54_66_total': died & t.between(54,66),
            'alive_54_66': alive & t.between(54,66),
            'unknown_5y_54_66': s['y_5y'].isna() & t.between(54,66),
        }
        weights = s['WTMEC2YR'].to_numpy(float) / s.cycle.nunique()
        row = {'scope':name, 'prerequisite_n':len(s), 'analytic_n':int(s['y_5y'].notna().sum())}
        for key, selected in windows.items():
            row[key+'_n'] = int(selected.sum())
            row[key+'_weighted_fraction'] = float(weights[selected.to_numpy()].sum()/weights.sum())
        boundary_rows.append(row)
    pd.DataFrame(boundary_rows).to_csv(out / 'public_followup_boundary_counts.csv', index=False)
    report = {
        'generated_utc':datetime.now(timezone.utc).isoformat(),
        'status':'pre_model_descriptive_audit_completed_no_models_fitted',
        'cycle_n':{cycle:int(cohort.cycle.eq(cycle).sum()) for cycle in CYCLES},
        'measurement_distribution_population':'prespecified linked CBC-complete analytic cohort',
        'quantiles':'inverse weighted empirical CDF; no interpolation',
        'distribution_weight':'WTMEC2YR within each cycle',
        'boundary_population':'MEC, age40+, linkage eligible, positive observed CBC, valid vital status and nonnegative follow-up; before 5-year ascertainability exclusion',
        'descriptive_month_window':[54,66],
        'window_is_not_a_CDC_perturbation_bound':True,
        'perturbed_record_count':'not_identifiable_from_public_file',
        'misclassified_5y_record_count':'not_identifiable_from_public_file',
        'no_horizon_change_no_boundary_exclusion_no_model_fit':True,
    }
    (out/'pre_model_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    lines=['# 模型前描述性审计','',
           '本次只比较公开测量分布并计数随访边界记录；未训练模型、未计算预测性能。','',
           '## 周期分布','',
           '| 周期 | n | RDW 加权均值 | Hb 加权均值 | 仪器 |','| --- | ---: | ---: | ---: | --- |']
    table=pd.DataFrame(rows)
    for cycle in CYCLES:
        sub=table.loc[table.cycle.eq(cycle)].set_index('measurement')
        lines.append(f"| {cycle} | {int(sub.loc['rdw_pct','n'])} | {sub.loc['rdw_pct','weighted_mean']:.3f}% | {sub.loc['hemoglobin_g_dl','weighted_mean']:.3f} g/dL | {sub.loc['rdw_pct','instrument']} |")
    lines += ['', '完整加权SD、分位数和范围见配套CSV。分布不能识别仪器的因果效应，也不证明数值可互换；不据此进行测试集分布对齐或改变周期。','',
              '## 公开随访边界','',
              '| 分析池 | 54–59月死亡 | 恰60月死亡 | 61–66月死亡 | 54–66月存活 | 其中五年状态未知 |',
              '| --- | ---: | ---: | ---: | ---: | ---: |']
    for row in boundary_rows:
        if row['scope'] in ['development','test','all_2003_2010']:
            lines.append(f"| {row['scope']} | {row['death_54_59_n']} | {row['death_exact_60_n']} | {row['death_61_66_n']} | {row['alive_54_66_n']} | {row['unknown_5y_54_66_n']} |")
    lines += ['', '未知状态列是边界内记录的附加描述，不能与存活人数相加。54–66月是本次预定描述窗口，不是官方时间扰动范围。公开文件未标识被替换记录，因此这里的数目既不是被替换人数，也不是五年结局误分类人数或其上界。','',
              '官方依据见 protocol/CBC_INSTRUMENT_AUDIT.md 与 protocol/LMF_WEIGHT_AUDIT.md。']
    (out/'DESCRIPTIVE_AUDIT_zh.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'output_directory':str(out),'n':len(cohort),'model_fits':0,
                      'boundary_development':boundary_rows[-3], 'boundary_test':boundary_rows[-2]},ensure_ascii=False))


if __name__=='__main__':
    main()
