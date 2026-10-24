"""Render existing estimates only. No cohort, fitted object, cache or fit entry is read."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
from datetime import datetime, timezone

ROOT=Path(__file__).resolve().parents[1]
REFERENCE=ROOT/'reference' if (ROOT/'reference').is_dir() else ROOT
OUT=REFERENCE/'outputs'

def copy_tree(source,destination):
    # Preserve deep scientific attachments when a review package is unpacked
    # under a long Windows workspace path.
    def path(p):
        absolute=str(Path(p).resolve())
        return '\\\\?\\'+absolute if os.name=='nt' and not absolute.startswith('\\\\?\\') else absolute
    shutil.copytree(path(source),path(destination),dirs_exist_ok=True)
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def write(p,obj):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def signature(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(str(x).replace('|','/') for x in row)+' |' for row in rows])
def ci(x,d=6):return f'{x["estimate"]:.{d}f} ({x["lower"]:.{d}f} to {x["upper"]:.{d}f})'
def load_saved():
    s=read(OUT/'tables/fixed_substitution_summary.json')
    q=read(OUT/'tables/q1_measurement_performance.json')
    f=read(OUT/'figures/Figure6_source_data.json')
    # Check report sources agree. This does not re-evaluate any study prediction.
    assert f['primary']==s['primary_40plus'] and f['q1']==q,'Figure and saved table sources disagree.'
    assert (s['primary_40plus']['n'],s['primary_40plus']['events_5y'])==(7477,648)
    assert [x['measurement'] for x in q]==['rdw_pct','hemoglobin_g_dl']
    for key in ['primary_40plus','sensitivity_40_79']:
        d=s[key]
        assert set(d['Brier'])=={'B','M','S','MEAN','P'}
        for v in [*d['Brier'].values(),*d['contrasts'].values()]:
            assert all(math.isfinite(v[x]) for x in ['estimate','lower','upper'])
            assert v['lower']<=v['estimate']<=v['upper']
        c=d['contrasts']
        assert math.isclose(c['G_measured_gain']['estimate']-c['D_substitution_loss']['estimate'],c['H_remaining_gain']['estimate'],abs_tol=1e-14)
    return s,q,f

def draw_figure(source,dest):
    from render_saved_figures import draw_figure6
    draw_figure6(source,dest)


def merge_main(spec):
    text=(ROOT/'input/context/current_manuscript_patterns_en.md').read_text(encoding='utf-8')
    for op in spec['operations']:
        kind=op['kind'];new=op.get('new','')
        if kind=='append':text=text.rstrip()+'\n\n'+new+'\n';continue
        if kind=='table_append':
            lines=text.splitlines();start=next(i for i,l in enumerate(lines) if l.startswith('| '+op['header']+' |'))
            end=start
            while end<len(lines) and lines[end].startswith('|'):end+=1
            lines.insert(end,'| '+' | '.join(op['row'])+' |');text='\n'.join(lines)+'\n';continue
        old=op['old'];assert text.count(old)==1,(old[:80],text.count(old))
        value=new if kind=='replace' else (new+'\n\n'+old if kind=='before' else old+'\n\n'+new)
        text=text.replace(old,value,1)
    from methods_revision_texts import replace_table3
    return replace_table3(text,spec)

def generate(report_dir=None):
    dest=Path(report_dir).resolve() if report_dir else ROOT/'output/Patterns_Integrated'
    if dest==OUT or OUT in dest.parents:raise ValueError('Use a report destination outside the frozen outputs tree.')
    s,q,f=load_saved()
    sources=list((OUT/'tables').glob('*.json'))+list((OUT/'tables').glob('*.csv'))+[OUT/'figures/Figure6_source_data.json']
    before={str(p.relative_to(ROOT)):signature(p) for p in sources}
    log=OUT/'logs/execution.jsonl';log_before=signature(log) if log.exists() else None
    dest.mkdir(parents=True,exist_ok=True)
    if not (ROOT/'input/editorial/integration_texts.json').exists():
        # Public code exports intentionally omit unpublished manuscript templates.
        draw_figure(f,dest/'figures')
        (dest/'source_data').mkdir(exist_ok=True)
        for p in sources:shutil.copy2(p,dest/'source_data'/p.name)
        sections=[]
        for key in ['primary_40plus','sensitivity_40_79']:
            d=s[key]
            sections.extend(['## '+key,table(['Route or contrast','Estimate (95% CI)'],[[k,ci(v)] for k,v in {**d['Brier'],**d['contrasts']}.items()])])
        (dest/'results_tables.md').write_text('\n\n'.join(sections)+'\n',encoding='utf-8')
        assert before=={str(p.relative_to(ROOT)):signature(p) for p in sources}
        assert (signature(log) if log.exists() else None)==log_before
        write(dest/'author_only/report_reuse_record.json',dict(mode='outputs-only',source_statistics_unchanged=True,execution_log_unchanged=True,manuscript_templates='not included in public code export',source_sha256=before))
        print(f'Saved tables and Figure 6 generated: {dest}. Local unpublished manuscript templates are not present.')
        return dest
    for folder in ['figures','scientific_SI','tables','release_candidate','change_notes','author_only']:
        if not (ROOT/'input/editorial'/folder).exists():continue
        copy_tree(ROOT/'input/editorial'/folder,dest/folder)
    for name in ['author_items_to_complete.md','claim_to_evidence.tsv']:
        p=ROOT/'input/editorial'/name
        if p.exists():shutil.copy2(p,dest/name)
    draw_figure(f,dest/'figures')
    from render_main_figures import draw_main_figures
    draw_main_figures(dest/'figures')
    for p in sources:
        target=dest/'source_data'/p.name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
    spec=read(ROOT/'input/editorial/integration_texts.json')
    (dest/'manuscript_patterns_en.md').write_text(merge_main(spec),encoding='utf-8')
    from clinical_supplement import clinical_text, SI_REPLACEMENTS
    clinical=clinical_text(s,q)
    (dest/'clinical_module.md').write_text(clinical,encoding='utf-8')
    si=(ROOT/'input/context/current_supplementary_patterns_en.md').read_text(encoding='utf-8')
    oldintro=si.split('\n\n')[1]
    from methods_revision_texts import si_intro,coverage_text
    newintro=si_intro(ROOT)
    si=si.replace(oldintro,newintro,1).replace('## References',clinical+'\n\n## References',1).rstrip()+'\n\n'+spec['references']+'\n'
    for old,new in SI_REPLACEMENTS.items():si=si.replace(old,new,1)
    si=si.replace('## References',coverage_text(ROOT)+'\n\n## References',1)
    (dest/'supplementary_patterns_en.md').write_text(si,encoding='utf-8')
    legends='# Main figure legends\n\n'+'\n\n'.join('## Figure '+n+'\n\n'+t for n,t in spec['main_figure_legends'].items())+'\n'
    (dest/'figure_legends.md').write_text(legends,encoding='utf-8')
    # Keep the manuscript Figure 1 legend synchronized with the separate legend file.
    main=dest/'manuscript_patterns_en.md';mt=main.read_text(encoding='utf-8')
    from visual_editorial_texts import replace_main_legends
    mt=replace_main_legends(mt,spec)
    main.write_text(mt,encoding='utf-8')
    sub=dest/'submission_texts';sub.mkdir(exist_ok=True)
    for name,key in [('summary.txt','summary'),('the_bigger_picture.txt','bigger_picture'),('keywords.txt','keywords')]:
        (sub/name).write_text(spec[key]+'\n',encoding='utf-8')
    (sub/'highlights.txt').write_text('\n'.join(spec['highlights'])+'\n',encoding='utf-8')
    from render_saved_figures import draw_s12
    draw_s12(dest/'figures')
    for name in ['cbc_distributions_by_cycle.csv','public_followup_boundary_counts.csv']:shutil.copy2(ROOT/'diagnostics'/name,dest/'source_data'/name)
    for name in ['cohort_flow.json','feasibility_profile.json']:
        p=ROOT/'data_checks'/name
        if p.exists():shutil.copy2(p,dest/'source_data'/name)
    from submission_finalize import finalize
    finalize(dest)
    assert before=={str(p.relative_to(ROOT)):signature(p) for p in sources},'Report generation changed saved statistics.'
    assert (signature(log) if log.exists() else None)==log_before,'Report generation changed the fit/execution log.'
    write(dest/'author_only/report_reuse_record.json',dict(utc=datetime.now(timezone.utc).isoformat(),mode='outputs-only',source_statistics_unchanged=True,execution_log_unchanged=True,source_sha256=before,source='existing CSV/JSON; no predictions recomputed',destination=str(dest)))
    print(f'Report-only artifacts generated: {dest}. Saved statistics and execution log unchanged.')
    return dest

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--outputs-only','--report-only',action='store_true')
    p.add_argument('--report-dir')
    a=p.parse_args();generate(a.report_dir)
if __name__=='__main__':main()
