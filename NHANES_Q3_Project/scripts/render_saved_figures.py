"""Print-size figures from existing aggregate estimates and saved vector curves.

No cohort, individual predictions, fitted models or statistical estimators are read.
"""
from pathlib import Path
import argparse,hashlib,json,re,shutil,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1]
REFERENCE=ROOT/'reference' if (ROOT/'reference').is_dir() else ROOT
WIDTH_MM=167
def setup():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':9.5,'axes.labelsize':9,'xtick.labelsize':8.5,'ytick.labelsize':8.5,'legend.fontsize':8.5,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42,'path.simplify':False})
    return plt
def save(fig,dest,stem):
    dest.mkdir(parents=True,exist_ok=True)
    # Fixed canvas width is the intended manuscript width; no tight-bbox resizing.
    for ext in ['pdf','svg','png']:fig.savefig(dest/(stem+'.'+ext),dpi=360)
    pts=[x.get_fontsize() for x in fig.findobj() if hasattr(x,'get_fontsize') and hasattr(x,'get_text') and x.get_text() and x.get_visible()]
    (dest/(stem+'_layout.json')).write_text(json.dumps({'canvas_width_mm':WIDTH_MM,'canvas_height_mm':fig.get_figheight()*25.4,'minimum_visible_text_pt':min(pts),'dpi':360,'statistics_recomputed':False},indent=2)+'\n',encoding='utf-8')
def draw_figure6(source,dest):
    plt=setup()
    from matplotlib.ticker import FormatStrFormatter,MaxNLocator
    r,q=source['primary'],source['q1']
    fig=plt.figure(figsize=(WIDTH_MM/25.4,6.0))
    a=fig.add_axes([.03,.56,.43,.37]);a.axis('off');a.set_title('A  Temporal design',loc='left',weight='bold',pad=10)
    boxes=[(.84,'2003–2004 · MAXM\nFit: n=2,955; deaths=370'),(.60,'2005–2006 · MAXM\nSelect: n=2,804; deaths=288'),(.36,'2003–2006 · MAXM\nDevelop: n=5,759; deaths=658'),(.12,'2007–2010 · HMX\nTest: n=7,477; deaths=648')]
    stage_labels=[]
    for y,t in boxes:
        stage_labels.append(a.text(.5,y,t,ha='center',va='center',transform=a.transAxes,fontsize=8.5,linespacing=1.4,bbox=dict(boxstyle='round,pad=.35',fc='#eef3f7',ec='#8aa6b9')))
    # Connect the actual box borders. The old short, freely placed arrows lost
    # their stems after endpoint shrinking at the final print width.
    for upper,lower in zip(stage_labels,stage_labels[1:]):
        a.annotate('',xy=lower.get_position(),xytext=upper.get_position(),xycoords='axes fraction',arrowprops=dict(arrowstyle='-|>',patchA=upper.get_bbox_patch(),patchB=lower.get_bbox_patch(),shrinkA=.6,shrinkB=.6,mutation_scale=9,color='#35546b',lw=1.2))
    for j,row in enumerate(q):
        ax=fig.add_axes([.65,.795 if j==0 else .575,.31,.115])
        ax.set_title(('B  Measurement prediction\n' if j==0 else '')+('RDW' if j==0 else 'Hb')+f' · R²={row["R2"]:.3f}',loc='left',weight='bold' if j==0 else 'normal',pad=9)
        for i,(key,col) in enumerate([('MAE','#C77722'),('development_mean_MAE','#68727B')]):
            v=row[key];ax.errorbar(v['estimate'],i,xerr=[[v['estimate']-v['lower']],[v['upper']-v['estimate']]],fmt='o' if i==0 else 's',markersize=4.5,color=col,capsize=3,lw=1.25)
        ax.set_yticks([0,1],['Proxy','Mean']);ax.invert_yaxis();ax.set_xlim(left=0);ax.margins(y=.65)
        ax.xaxis.set_major_locator(MaxNLocator(3));ax.set_xlabel('MAE (percentage points)' if j==0 else 'MAE (g/dL)');ax.grid(axis='x',alpha=.18)
    c=fig.add_axes([.21,.235,.265,.23])
    fig.text(.03,.485,'C  Five-route Brier losses',weight='bold',fontsize=9.5)
    for i,key in enumerate(['B','M','S','MEAN','P']):
        col=['#68727B','#216A98','#C77722','#68727B','#C77722'][i]
        v=r['Brier'][key];c.errorbar(v['estimate'],i,xerr=[[v['estimate']-v['lower']],[v['upper']-v['estimate']]],fmt=['s','o','o','s','D'][i],markersize=4.5,color=col,capsize=3,lw=1.25)
    c.set_yticks(range(5),['B · baseline','M · measured','S · fixed proxy','MEAN · fixed mean','P · retrained']);c.invert_yaxis();c.margins(y=.2)
    c.xaxis.set_major_locator(MaxNLocator(3));c.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    c.set_xlabel('Weighted Brier loss');c.grid(axis='x',alpha=.18)
    d=fig.add_axes([.735,.235,.23,.23])
    fig.text(.52,.485,'D  Paired Brier contrasts',weight='bold',fontsize=9.5)
    for i,key in enumerate(['G_measured_gain','D_substitution_loss','H_remaining_gain']):
        v=r['contrasts'][key];d.errorbar(v['estimate'],i,xerr=[[v['estimate']-v['lower']],[v['upper']-v['estimate']]],fmt=['o','s','D'][i],markersize=5,color=['#216A98','#C77722','#247C72'][i],capsize=3,lw=1.4)
    d.axvline(0,color='#666',linestyle='--',lw=.8);d.set_yticks(range(3),['G: B − M\nMeasured gain','D: S − M\nSubstitution loss','H: B − S\nRemaining gain']);d.invert_yaxis();d.margins(y=.4)
    d.set_xticks([0,.0015,.003]);d.xaxis.set_major_formatter(FormatStrFormatter('%.4f'));d.set_xlabel('Paired Brier difference');d.grid(axis='x',alpha=.18)
    rr=r['retention'];lo,hi=rr['intervals'][0]
    fig.text(.03,.10,'S and MEAN use the same fitted M.\nP is trained separately.',fontsize=8.5,linespacing=1.4)
    fig.text(.52,.10,f'R = H/G: {rr["point"]:.1%} retained\n95% Fieller: {lo:.1%} to {hi:.1%}',fontsize=8.5,linespacing=1.4,
             bbox=dict(boxstyle='round,pad=.4',fc='#F3F6F8',ec='#D8E0E5',lw=.7))
    fig.text(.03,.025,'Primary: age ≥40 years. Positive remaining gain not established.\nAge 40–79 sensitivity: see Table S30; no age-interaction claim.',fontsize=8.5,linespacing=1.4)
    save(fig,dest,'Figure6_NHANES_Q3');plt.close(fig)
    shutil.copy2(REFERENCE/'outputs/figures/Figure6_source_data.json',dest/'Figure6_source_data.json')
def saved_svg_curves(path):
    # Recover only the displayed polyline geometry, including the original SVG's
    # display simplification. No empirical distribution is recomputed.
    ns={'s':'http://www.w3.org/2000/svg'};parser=ET.XMLParser(target=ET.TreeBuilder(insert_comments=True))
    root=ET.parse(path,parser=parser).getroot();panels=[]
    colors=['#355c7d','#68a4a5','#c17c36','#9d4b59']
    for i in [1,2]:
        ax=root.find(f'.//s:g[@id="axes_{i}"]',ns)
        bg=ax.find('s:g/s:path',ns).get('d');nums=[float(x) for x in re.findall(r'-?\d+(?:\.\d+)?',bg)]
        left,bottom,right,_,_,top=nums[:6]
        ticks=[]
        for g in ax.findall('.//s:g',ns):
            if not g.get('id','').startswith('xtick_'):continue
            use=g.find('.//s:use',ns);label=next((n.text.strip() for n in g.iter() if n.tag is ET.Comment),None)
            if use is not None and label is not None:ticks.append((float(use.get('x')),float(label.replace('−','-'))))
        (p1,v1),(p2,v2)=ticks[:2];scale=(v2-v1)/(p2-p1)
        xlim=[v1+(left-p1)*scale,v1+(right-p1)*scale]
        curves=[]
        for color in colors:
            paths=[p for p in ax.findall('.//s:path',ns) if f'stroke: {color}' in p.get('style','') and p.get('clip-path')]
            assert len(paths)==1,(i,color,len(paths));d=paths[0].get('d')
            assert not re.search('[A-KN-Z]',d), 'Expected saved M/L polyline geometry.'
            xy=[float(x) for x in re.findall(r'-?\d+(?:\.\d+)?(?:e[-+]?\d+)?',d)]
            assert len(xy)%2==0
            vertices=[[(xy[j]-left)/(right-left),(bottom-xy[j+1])/(bottom-top)] for j in range(0,len(xy),2)]
            curves.append({'color':color,'vertices_axes_fraction':vertices,'original_svg_path':d})
        panels.append({'xlim':xlim,'ylim':[0,1.01],'curves':curves})
    return {'source':'diagnostics/cbc_distributions_by_cycle.svg','source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'role':'existing rendered curve vertices; not participant records and not newly calculated ECDF values','panels':panels}
def draw_s12(dest):
    plt=setup();source=saved_svg_curves(REFERENCE/'diagnostics/cbc_distributions_by_cycle.svg')
    fig,axes=plt.subplots(1,2,figsize=(WIDTH_MM/25.4,3.55))
    labels=['2003–2004 (n=2,955)','2005–2006 (n=2,804)','2007–2008 (n=3,648)','2009–2010 (n=3,829)']
    for ax,panel,title,xlabel in zip(axes,source['panels'],['A  RDW','B  Hemoglobin'],['RDW (%)','Hemoglobin (g/dL)']):
        for curve,label in zip(panel['curves'],labels):
            xy=curve['vertices_axes_fraction'];ax.plot([p[0] for p in xy],[p[1] for p in xy],transform=ax.transAxes,color=curve['color'],lw=1.1,label=label,clip_on=True)
        ax.set_xlim(panel['xlim']);ax.set_ylim(panel['ylim']);ax.set_title(title,loc='left',weight='bold');ax.set_xlabel(xlabel);ax.set_yticks([0,.25,.5,.75,1]);ax.grid(axis='y',alpha=.18)
    axes[0].set_ylabel('Survey-weighted\ncumulative proportion');axes[1].tick_params(labelleft=False)
    h,l=axes[0].get_legend_handles_labels();fig.legend(h,l,loc='lower center',ncol=2,frameon=False,bbox_to_anchor=(.52,.01),handlelength=2)
    fig.subplots_adjust(left=.12,right=.98,top=.89,bottom=.30,wspace=.20)
    save(fig,dest,'FigureS12_CBC_distributions');plt.close(fig)
    (dest/'FigureS12_curve_source_data.json').write_text(json.dumps(source,separators=(',',':'))+'\n',encoding='utf-8')
def main():
    p=argparse.ArgumentParser();p.add_argument('--report-only','--outputs-only',action='store_true');p.add_argument('--report-dir',type=Path,default=ROOT/'output/Patterns_Integrated');a=p.parse_args()
    source=json.loads((REFERENCE/'outputs/figures/Figure6_source_data.json').read_text(encoding='utf-8'))
    draw_figure6(source,a.report_dir/'figures');draw_s12(a.report_dir/'figures')
    print('Saved Figure 6 estimates and Figure S12 vector curves rendered at 167 mm; no statistical recomputation.')
if __name__=='__main__':main()
