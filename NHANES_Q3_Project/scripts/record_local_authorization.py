"""One-time record of the author's decision to change the execution prerequisite."""
from runtime import ROOT, read_json, write_json, now, event

def main():
    rec=read_json(ROOT/'protocol/PROTOCOL_RECORD.json')
    if rec.get('author_authorized_local_execution'):
        print('Existing direct-author local-execution instruction retained.')
        return
    stamp=now()
    statement='目前没有Git仓库，先进行实验，后续做完所有实验再上传'
    authorization=dict(recorded_utc=stamp,instruction=statement,
        authority='author decision recorded on 2026-10-19',
        supersedes='remote protocol record before any study-data model fit',
        scope='execute complete V2 experiment locally now; later upload only to an author-designated repository',
        remote_record_completed=False,first_study_fit_utc_at_change=rec.get('first_study_fit_utc'))
    rec['author_authorized_local_execution']=authorization
    rec['status']='author_authorized_local_execution_before_later_repository_upload_no_study_fit_yet'
    write_json(ROOT/'protocol/PROTOCOL_RECORD.json',rec)
    cfg=read_json(ROOT/'protocol/analysis_config.json')
    cfg['protocol_record']['remote_record_required_before_first_study_fit']=False
    cfg['protocol_record']['current_remote_status']='not_recorded_author_directed_experiments_before_upload'
    cfg['protocol_record']['author_override']=authorization
    write_json(ROOT/'protocol/analysis_config.json',cfg)
    with (ROOT/'protocol/DEVIATIONS.md').open('a',encoding='utf-8') as f:
        f.write(f'''\n## {stamp} — 作者明确调整远端记录前置条件\n\n2026-10-19记录的作者决定为：“{statement}”。按最新明确指令，取消本轮首次真实数据拟合前必须完成GitHub远端记录的前置条件；保留本地V2原计划及真实执行日志，完成实验后再处理上传。该变化独立于先前“仅调整执行节奏”的记录，不将两者混写。此时 first_study_fit_utc={rec.get('first_study_fit_utc')!r}，尚未拟合真实研究模型、尚未看任何研究性能。科学参数、16项X、两项Z、周期角色、候选、严格OOF、调查权重、Fieller及敏感性均未改变。\n\nrepository_url、protocol_commit_url、remote_recorded_utc维持null。不能把后续上传倒称模型前远端记录，不称正式预注册。此次仅改变记录／发布顺序，后续上传仍需实际指定目标；没有擅自创建仓库或上传数据。\n''')
    with (ROOT/'EXECUTION_CURRENT.md').open('a',encoding='utf-8') as f:
        f.write(f'\n## {stamp} 最新作者覆盖指令\n\n作者明确要求先实验、后上传；上述远端前置限制已由该直接指令覆盖。程序允许读取已记录的作者本地执行授权后启动真实拟合。远端字段仍为空，不能声称预注册或模型前远端时间戳。其余科学与资源设置不变。\n')
    event('author_authorized_local_execution',**authorization)
    print(authorization)

if __name__=='__main__': main()
