"""Shared publication drawing primitives; no estimation or model access."""
from pathlib import Path
import copy, io, json, os, re, sys, xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
if (ROOT / '.runtime').exists():
    sys.path.append(str(ROOT / '.runtime'))
_dll = None
if os.name == 'nt':
    import shutil
    cairo_dir = os.environ.get('CAIRO_DLL_DIRECTORY')
    if not cairo_dir and shutil.which('pdftoppm'):
        cairo_dir = str(Path(shutil.which('pdftoppm')).parent)
    if cairo_dir:
        _dll = os.add_dll_directory(cairo_dir)
        os.environ['PATH'] = cairo_dir + os.pathsep + os.environ.get('PATH', '')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager, ft2font
from matplotlib.patches import Rectangle, FancyBboxPatch
import pymupdf as fitz
import cairosvg

W = 174 / 25.4 * 72
BLUE, ORANGE, GRAY, TEAL = '#216B8C', '#C97A26', '#68727B', '#27867A'
INK, MUTED, SEP = '#25313A', '#77818A', '#DDE3E8'
PALEBLUE, PALEORANGE = '#F3F8FA', '#FCF6EF'
LEDGER = []
STYLE_LEDGER = []

# Use the author's installed Arial. Never substitute another family silently and
# never redistribute a font binary or the private path-bearing configuration.
font_config = ROOT/'private/arial_fonts.local.json'
if font_config.exists():
    for face in json.loads(font_config.read_text(encoding='utf-8'))['faces'].values():
        font_manager.fontManager.addfont(face['path'])
for weight,style in [('normal','normal'),('bold','normal'),('normal','italic')]:
    path=font_manager.findfont(font_manager.FontProperties(family='Arial',weight=weight,style=style),fallback_to_default=False)
    if ft2font.FT2Font(path).family_name != 'Arial':
        raise RuntimeError('Actual installed Arial is required; run source/check_arial.py.')

plt.rcParams.update({
    'font.family': 'Arial', 'font.sans-serif': ['Arial'], 'font.size': 7.5,
    'axes.titlesize': 8, 'axes.titleweight': 'bold', 'axes.labelsize': 7.5,
    'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
    'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic',
    'mathtext.bf': 'Arial:bold', 'mathtext.sf': 'Arial', 'mathtext.fallback': None,
    'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': .65,
    'lines.linewidth': 1.2, 'pdf.fonttype': 42, 'svg.fonttype': 'none',
    'path.simplify': False, 'savefig.facecolor': 'white', 'text.color': INK,
    'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK,
    'axes.edgecolor': INK, 'axes.unicode_minus': True,
})

def canvas(height):
    fig = plt.figure(figsize=(W / 72, height / 72))
    fig._height = height
    fig._icons = []
    ax = fig.add_axes([0, 0, 1, 1], zorder=0)
    ax.set(xlim=(0, W), ylim=(height, 0)); ax.axis('off')
    return fig, ax

def txt(ax, x, y, text, size=7.5, color=INK, weight='normal', ha='left', va='center', **kw):
    return ax.text(x, y, text, fontsize=size, color=color, weight=weight,
                   ha=ha, va=va, linespacing=1.35, **kw)

def line(ax, x1, y1, x2, y2, color=SEP, lw=.7, ls='-'):
    return ax.plot([x1, x2], [y1, y2], color=color, lw=lw, ls=ls, clip_on=False)

def arrow(ax, start, end, color=GRAY, ls='-', lw=.85):
    ax.annotate('', xy=end, xytext=start,
                arrowprops=dict(arrowstyle='-|>', color=color, lw=lw, ls=ls,
                                mutation_scale=7, shrinkA=1, shrinkB=1))

def card(ax, x, y, w, h, label='', edge=SEP, fill='white', size=7.5, color=INK, weight='normal'):
    ax.add_patch(FancyBboxPatch((x,y), w,h, boxstyle='round,pad=0,rounding_size=2',
                              lw=.75, edgecolor=edge, facecolor=fill))
    if label: txt(ax,x+w/2,y+h/2,label,size,color,weight,ha='center')

def icon(fig, name, x, y, w, h, opacity=1):
    fig._icons.append(dict(name=name,x=x,y=y,w=w,h=h,opacity=opacity))

def panel(ax, letter, title, x, y):
    txt(ax,x,y,letter,9,weight='bold')
    txt(ax,x+14,y,title,8,weight='bold')

def chart(fig, x, y, w, h):
    ax=fig.add_axes([x/W,1-(y+h)/fig._height,w/W,h/fig._height],zorder=1)
    ax.tick_params(length=3, width=.65, pad=4)
    ax.set_axisbelow(True)
    return ax

def grid(ax): ax.grid(axis='x',color='#E8ECEF',lw=.55)

def interval(ax,row,y,keys,color=BLUE,marker='o',lw=1.2,size=4.4,tag='',open_marker=False):
    point,low,high=[float(row[k]) for k in keys]
    ax.plot([low,high],[y,y],color=color,lw=lw,zorder=2)
    for x in (low,high): ax.plot([x,x],[y-.045,y+.045],color=color,lw=.7,zorder=2)
    ax.plot(point,y,marker=marker,color=color,mfc='white' if open_marker else color,
            mew=.95,ms=size,zorder=3)
    LEDGER.append(dict(tag=tag,point=point,low=low,high=high))
    STYLE_LEDGER.append(dict(tag=tag,marker=marker,color=color,open_marker=open_marker))

def asset_svg(spec, prefix):
    source=ROOT/'assets/icons'/f"{spec['name']}.svg"
    root=ET.fromstring(source.read_text(encoding='utf-8'))
    # Prefix every ID/reference so separately licensed icons can coexist.
    mapping={e.attrib['id']:prefix+e.attrib['id'] for e in root.iter() if 'id' in e.attrib}
    classes={c:prefix+c for e in root.iter() for c in e.get('class','').split()}
    recolor={
        'cell_culture_dish':{'#d10b48':'#B98089','#ea4b6a':'#CAA1A8'},
        'blood_sample_tube':{'#47c9fc':'#88B4C2','#f15a24':'#AC7660','#f7755c':'#C69B8C','#c44426':'#936553'},
        'rat':{'#f5f5f6':'#D4D7DA','#ececed':'#C5CBD0','#d8d9da':'#ACB4BB',
               '#d3d4d5':'#A5AEB7','#d0d1d3':'#A5AEB7','#1a171b':'#25313A'},
    }.get(spec['name'],{})
    for e in root.iter():
        if e.tag.endswith('}style') and e.text:
            css=e.text
            for old,new in classes.items():
                css=re.sub(r'\.'+re.escape(old)+r'(?![\w-])','.'+new,css)
            for old,new in mapping.items(): css=css.replace(f'url(#{old})',f'url(#{new})')
            for old,new in recolor.items(): css=css.replace(old,new)
            e.text=css
        for k,v in list(e.attrib.items()):
            if k=='id': e.set(k,mapping[v]);continue
            if k=='class': e.set(k,' '.join(classes[c] for c in v.split()));continue
            for old,new in mapping.items():
                v=v.replace(f'url(#{old})',f'url(#{new})')
                if v==f'#{old}': v=f'#{new}'
            if 'currentColor' in v: v=v.replace('currentColor',GRAY)
            for old,new in recolor.items(): v=v.replace(old,new)
            if spec['name']=='rat' and k=='stroke-width' and e.get('stroke') not in ('none',None,'#fff'):
                v=re.sub(r'[\d.]+',lambda m:str(float(m.group())*1.65),v,count=1)
            e.set(k,v)
    root.set('x',str(spec['x']));root.set('y',str(spec['y']))
    root.set('width',str(spec['w']));root.set('height',str(spec['h']))
    root.set('preserveAspectRatio','xMidYMid meet')
    root.set('opacity',str(spec['opacity']))
    return root

def save(fig,name):
    dest=ROOT/'main_figures';dest.mkdir(exist_ok=True)
    (ROOT/'preview').mkdir(exist_ok=True);(ROOT/'validation').mkdir(exist_ok=True)
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    sizes=[];outside=[]
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_visible() and t.get_text():
            sizes.append(t.get_fontsize())
            bb=t.get_window_extent(renderer)
            if bb.x0 < -1 or bb.y0 < -1 or bb.x1 > fig.bbox.width+1 or bb.y1 > fig.bbox.height+1:
                outside.append(t.get_text())
    svg=io.BytesIO();pdf=io.BytesIO()
    fig.savefig(svg,format='svg');fig.savefig(pdf,format='pdf')
    ET.register_namespace('','http://www.w3.org/2000/svg')
    ET.register_namespace('xlink','http://www.w3.org/1999/xlink')
    root=ET.fromstring(svg.getvalue())
    doc=fitz.open(stream=pdf.getvalue(),filetype='pdf')
    for i,s in enumerate(fig._icons):
        root.append(asset_svg(s,f'icon{i}_'))
        raw=asset_svg(dict(s,x=0,y=0,opacity=1),f'pdficon{i}_')
        raw.set('width',f"{s['w']}pt");raw.set('height',f"{s['h']}pt")
        icon_pdf=cairosvg.svg2pdf(bytestring=ET.tostring(raw))
        with fitz.open(stream=icon_pdf,filetype='pdf') as idoc:
            doc[0].show_pdf_page(fitz.Rect(s['x'],s['y'],s['x']+s['w'],s['y']+s['h']),idoc,0)
    ET.ElementTree(root).write(dest/f'{name}.svg',encoding='utf-8',xml_declaration=True)
    doc.save(dest/f'{name}.pdf',garbage=4,deflate=True)
    doc[0].get_pixmap(dpi=360,alpha=False).save(dest/f'{name}.png')
    doc[0].get_pixmap(dpi=150,alpha=False).save(ROOT/'preview'/f'{name}_150dpi.png')
    from PIL import Image,ImageOps
    im=Image.open(ROOT/'preview'/f'{name}_150dpi.png')
    ImageOps.grayscale(im).save(ROOT/'preview'/f'{name}_gray.png')
    rec=dict(canvas_width_mm=174,canvas_height_mm=fig._height/72*25.4,
             minimum_visible_text_pt=min(sizes),text_outside_canvas=outside,
             icons=fig._icons,statistics_recomputed=False,png_dpi=360,
             svg_embedded_raster_images=len(root.findall('.//{http://www.w3.org/2000/svg}image')))
    (ROOT/'validation'/f'{name}_layout.json').write_text(json.dumps(rec,indent=2)+'\n',encoding='utf-8')
    (ROOT/'validation'/f'{name}_drawing_ledger.json').write_text(json.dumps(LEDGER,indent=2)+'\n',encoding='utf-8')
    (ROOT/'validation'/f'{name}_symbol_ledger.json').write_text(json.dumps(STYLE_LEDGER,indent=2)+'\n',encoding='utf-8')
    plt.close(fig);doc.close();LEDGER.clear();STYLE_LEDGER.clear()
    print(name, 'written;', 'outside:',outside)
