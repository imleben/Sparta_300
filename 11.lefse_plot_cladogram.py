#!/usr/bin/env python3

import os, sys, matplotlib, argparse, string
matplotlib.use('Agg')
from pylab import *
from lefse.lefse import *
import numpy as np

colors = ['r','g','b','m','c',[1.0,0.5,0.0],[0.0,1.0,0.0],[0.33,0.125,0.0],[0.75,0.75,0.75],'k']
dark_colors = [[0.4,0.0,0.0],[0.0,0.2,0.0],[0.0,0.0,0.4],'m','c',[1.0,0.5,0.0],[0.0,1.0,0.0],[0.33,0.125,0.0],[0.75,0.75,0.75],'k']

class CladeNode:
    def __init__(self, name, abundance, viz=True):
        self.id = name
        self.name = name.split(".")
        self.last_name = self.name[-1]
        self.abundance = abundance
        self.pos = (-1.0, -1.0)
        self.children = {}
        self.isleaf = True
        self.color = 'y'
        self.next_leaf = -1
        self.prev_leaf = -1
        self.viz = viz

    def __repr__(self):
        return self.last_name

    def add_child(self,node):
        self.isleaf = False
        self.children[node.__repr__()] = node

    def get_children(self):
        ck = sorted(self.children.keys())
        return [self.children[k] for k in ck]

    def get_color(self):
        return self.color

    def set_color(self,c):
        self.color = c

    def set_pos(self,pos):
        self.pos = pos

def read_params(args):
    parser = argparse.ArgumentParser(description='Cladoplot')
    parser.add_argument('input_file', metavar='INPUT_FILE', type=str, help="tab delimited input file")
    parser.add_argument('output_file', metavar='OUTPUT_FILE', type=str, help="the file for the output image")
    parser.add_argument('--clade_sep',dest="clade_sep", type=float, default=1.5)
    parser.add_argument('--max_lev',dest="max_lev", type=int, default=-1)
    parser.add_argument('--max_point_size',dest="max_point_size", type=float, default=6.0)
    parser.add_argument('--min_point_size',dest="min_point_size", type=float, default=1)
    parser.add_argument('--point_edge_width',dest="markeredgewidth", type=float, default=.25)
    parser.add_argument('--siblings_connector_width',dest="siblings_connector_width", type=float, default=2)
    parser.add_argument('--parents_connector_width',dest="parents_connector_width", type=float, default=0.75)
    parser.add_argument('--radial_start_lev',dest="radial_start_lev", type=int, default=1)
    parser.add_argument('--labeled_start_lev',dest="labeled_start_lev", type=int, default=1)
    parser.add_argument('--labeled_stop_lev',dest="labeled_stop_lev", type=int, default=5)
    parser.add_argument('--abrv_start_lev',dest="abrv_start_lev", type=int, default=3)
    parser.add_argument('--abrv_stop_lev',dest="abrv_stop_lev", type=int, default=5)
    parser.add_argument('--expand_void_lev',dest="expand_void_lev", type=int, default=1)
    parser.add_argument('--class_legend_vis',dest="class_legend_vis", type=int, default=1)
    parser.add_argument('--colored_connector',dest="colored_connectors", type=int, default=1)
    parser.add_argument('--alpha',dest="alpha", type=float, default=0.2)
    parser.add_argument('--title',dest="title", type=str, default="Cladogram")
    parser.add_argument('--sub_clade',dest="sub_clade", type=str, default="")
    parser.add_argument('--title_font_size',dest="title_font_size", type=str, default="14")
    parser.add_argument('--right_space_prop',dest="r_prop", type=float, default=0.1)
    parser.add_argument('--left_space_prop',dest="l_prop", type=float, default=0.1)
    parser.add_argument('--label_font_size',dest="label_font_size", type=str, default="8")
    parser.add_argument('--background_color',dest="back_color", type=str, choices=["k","w"], default="w", help="set the color of the background")
    parser.add_argument('--colored_labels',dest="col_lab", type=int, choices=[0,1], default=1, help="draw the label with class color (1) or in black (0)")
    parser.add_argument('--class_legend_font_size',dest="class_legend_font_size", type=str, default="6")
    parser.add_argument('--dpi',dest="dpi", type=int, default=300)
    parser.add_argument('--format', dest="format", choices=["png","svg","pdf"], default="pdf", type=str, help="the format for the output file")
    parser.add_argument('--all_feats', dest="all_feats", type=str, default="")
    args = parser.parse_args()
    return vars(args)

def cmp_names(la,lb):
    if len(la) != len(lb): return False
    for a,b in zip(la,lb):
        if a != b:
            return False
    return True

def build_tree(father,all_nodes,l,depth,viz):
    cc = [n for n in all_nodes if len(n.name) > len(father.name) and cmp_names(father.name,n.name[:len(father.name)])]
    children = [n for n in cc if len(n.name) == len(father.name)+1]
    if len(children) == 0 and l < depth -1:
        nc = CladeNode(father.id+"."+father.id.split(".")[-1],1.0,viz)
        father.add_child(nc)
        children.append(nc)
    for child in children:
        build_tree(child,cc,l+1,depth,viz)
        father.add_child(child)

def get_all_nodes(father):
    ret = [father]
    children = father.get_children()
    for c in children:
        ret += get_all_nodes(c)
    return ret

def read_data(input_file,params):
    with open(input_file, 'r') as inp:
        if params['sub_clade'] == "":
            rows = [line.strip().split()[:-1] for line in inp.readlines()
                    if params['max_lev'] < 1 or line.split()[0].count(".") < params['max_lev']]
        else:
            rows = [line.split(params['sub_clade']+".")[1].strip().split()[:-1]
                    for line in inp.readlines()
                    if (params['max_lev'] < 1 or line.split()[0].count(".") < params['max_lev'])
                    and line.startswith(params['sub_clade']+".")]
    all_names = [lin[0] for lin in rows]
    abundances = [float(v) for v in list(zip(*rows))[1] if float(v) >= 0.0]

    tree = {}
    tree['classes'] = list(set([v[2] for v in rows if len(v)>2]))
    tree['classes'].sort()

    all_nodes = [CladeNode("root."+row[0],float(row[1])) for row in rows]

    depth = max([len(n.name) for n in all_nodes])
    n2 = ["_".join(nn.name) for nn in all_nodes]
    for i, nn in enumerate(all_nodes):
        n = nn
        while "_".join(n.name[:-1]) not in n2 and len(n.name) > 1:
            n = CladeNode(".".join(n.name[:-1]), n.abundance)
            all_nodes.append(n)
            n2.append("_".join(n.name))

    cls2 = []
    if params['all_feats'] != "":
        cls2 = sorted(params['all_feats'].split(":"))
    for i,v in enumerate(rows):
        if len(v)>2:
            if len(cls2) > 0:
                idx = cls2.index(v[2])%len(colors)
                all_nodes[i].set_color(colors[idx])
            else:
                if v[2].count('rgbcol') > 0:
                    ccc = [float(tt) for tt in v[2].split('_')[1:]]
                    all_nodes[i].set_color(ccc)
                else:
                    idx = sorted(tree['classes']).index(v[2])%len(colors)
                    all_nodes[i].set_color(colors[idx])

    root = CladeNode("root",-1.0)
    root.set_pos((0.0,0.0))
    build_tree(root, all_nodes, 0, depth, params['expand_void_lev']==1)

    all_nodes = get_all_nodes(root)
    tree['root'] = root
    tree['max_abs'] = max(abundances) if abundances else 1.0
    tree['min_abs'] = min(abundances) if abundances else 0.0

    levs = []
    for i in range(depth):
        depthi = [n for n in all_nodes if len(n.name) == i+1]
        levs.append(len(depthi))
    tree['nlev'] = levs
    return tree

def add_all_pos(father,n,distn,seps,tsep,mlev,last_leaf=-1,nc=1):
    children = father.get_children()
    leaves = True if children[0].isleaf else False
    for i,child in enumerate(children):
        if leaves:
            n += 1.0
            men = 0.5 if len(children) == 1 else 0.0
            child.set_pos((n*distn-men*float(distn)+tsep,(len(father.name))/float(mlev-1)))
            if last_leaf != -1:
                child.prev_leaf = last_leaf
                last_leaf.next_leaf = child
            last_leaf = child
        else:
            ln = n
            ltsep = tsep
            n,tsep,last_leaf = add_all_pos(child,n,distn,seps,tsep,mlev,last_leaf,len(children))
            nn = (ln + n)*0.5*distn
            ssep = (ltsep + tsep)*0.5
            if n-ln == 1:
                ssep = ltsep
            child.set_pos((nn+ssep,(len(father.name))/float(mlev-1)))
    tsep += seps[len(father.name)-1] if (len(father.name)-1)<len(seps) else 0
    return n,tsep,last_leaf

def plot_points(father,params,pt_scale,ax):
    children = father.get_children()
    children.sort(key=lambda a: -int(a.get_color() == 'y')*a.abundance)
    x,r = father.pos
    for child in children:
        xc,rc = plot_points(child,params,pt_scale,ax)
    if not father.viz:
        return x,r
    ps = pt_scale[0] + father.abundance/pt_scale[1] + pt_scale[0]
    col = father.get_color()
    pw = params['markeredgewidth'] if col == 'y' else params['markeredgewidth']*3.0
    if x==0 and r==0:
        ax.plot(x,r, 'o', markersize=ps, color=col, markeredgewidth=0.01, markeredgecolor=params['fore_color'])
    else:
        ax.plot(x,r, 'o', markersize=ps, color=col, markeredgewidth=pw, markeredgecolor=params['fore_color'])
    return x,r

def plot_lines(father,params,depth,ax,xf):
    children = father.get_children()
    x,r = father.pos
    for i,child in enumerate(children):
        xc,rc = plot_lines(child,params,depth,ax,x)
        if i == 0:
            x_first, r_first = xc, rc
        if len(father.name) >= depth-params['radial_start_lev']:
            col = params['fore_color']
            lw = params['parents_connector_width']
            if not child.viz:
                continue
            if father.get_color() != 'y' and father.get_color() == child.get_color() and params['colored_connectors']:
                col = child.get_color()
                lw *=2.5
            if col != params['fore_color']:
                ax.plot([x,xc],[r,rc],"-",color=params['fore_color'],lw=lw*1.5)
            ax.plot([x,xc],[r,rc],"-",color=col,lw=lw)
    if not father.viz or (len(children) == 1 and not children[0].viz):
        return x,r
    if len(father.name) < depth-params['radial_start_lev']:
        col = params['fore_color']
        lw = params['parents_connector_width']
        if father.get_color() != 'y':
            same_flag = True
            for c in children:
                if c.get_color() != father.get_color() or not params['colored_connectors']:
                    same_flag = False
                    break
            if same_flag:
                col = father.get_color()
                lw *=2.5
        if not (x==0 and r==0):
            xx = xc if len(children) > 0 else x
            if len(children) == 0:
                rc = r
            xt = x if len(children) > 1 else xx
            if col != params['fore_color']:
                ax.plot([x,xt],[r,rc],"-",color=params['fore_color'],lw=lw*1.5)
            ax.plot([x,xt],[r,rc],"-",color=col,lw=lw)
    if len(children) > 0 and 1 < len(father.name) < depth-params['radial_start_lev']:
        xs = arange(x_first, xc, 0.01)
        ys = [rc for _ in xs]
        ax.plot(xs, ys, "-", color=col, lw=params['siblings_connector_width'], markeredgecolor=params['fore_color'])
    return x,r

def uniqueid():
    # 축약 기호를 a, b, ..., z, a0, a1, ... 식으로 생성
    for ch in string.ascii_lowercase:
        yield ch
    for ch in string.ascii_lowercase:
        for i in range(10):
            yield ch + str(i)
    i = 0
    while True:
        yield str(i)
        i += 1

def plot_names(father,params,depth,ax,u_i,seps):
    children = father.get_children()
    # 레벨에서 root(=1) 빼서 d__=1, ... g__=6
    l = len(father.name) - 1

    if len(children)==0:
        if father.prev_leaf == -1 or father.next_leaf == -1:
            fr_0, fr_1 = father.pos[0], father.pos[0]
        else:
            fr_0, fr_1 = (father.pos[0]+father.prev_leaf.pos[0])*0.5, (father.pos[0]+father.next_leaf.pos[0])*0.5
    for i, child in enumerate(children):
        fr,to = plot_names(child, params, depth, ax, u_i, seps)
        if i == 0:
            fr_0 = fr
        fr_1 = to

    if father.get_color() != 'y' and params['labeled_start_lev'] < l <= params['labeled_stop_lev']+1:
        col = father.get_color()
        dd = params['labeled_stop_lev'] - params['labeled_start_lev'] + 1
        de = depth - 1
        dim = 1.0/float(de)
        perc_ext = 0.65 if dim>0.1 else 1.0
        clto = (de-l+1)*dim+dim*(dd+1-(l-dd-1))*perc_ext
        clto = (de-l+1)*dim+dim*(dd-(l-params['labeled_start_lev'])+1)*perc_ext
        des = float(180.0*(fr_0+fr_1)/np.pi)*0.5-90
        lab = ""
        txt = father.last_name
        if params['abrv_start_lev']< l <= params['abrv_stop_lev']+1:
            ide = next(u_i)
            lab = str(ide)+": "+father.last_name
            txt = str(ide)

        # 실제 bar
        ax.bar(fr_0, clto, width=fr_1-fr_0,
               bottom=float(l-1)/float(de),
               alpha=params['alpha'], color=col, edgecolor=col)
        # 레전드 핸들(dummy bar)
        ax.bar(0.0, 0.0, width=0.0, bottom=0.0,
               alpha=1.0, color=col, edgecolor=params['fore_color'], label=lab)

        if l <= params['abrv_stop_lev']+1:
            if not params['col_lab']:
                col = params['fore_color']
            else:
                if col not in colors:
                    col = params['fore_color']
                else:
                    col = dark_colors[colors.index(col)%len(dark_colors)]
            ax.text((fr_0+fr_1)*0.5,
                    clto+float(l-1)/float(de)-dim*perc_ext/2.0,
                    txt, size=params['label_font_size'],
                    rotation=des, ha="center", va="center", color=col)
    return fr_0, fr_1

################# 여기서부터 핵심 수정 ####################

def draw_tree(out_file, tree, params):
    nlev = tree['nlev']
    pt_scale = (params['min_point_size'],
                max(1.0,((tree['max_abs']-tree['min_abs']))
                    /(params['max_point_size']-params['min_point_size'])))
    depth = len(nlev)
    sep = (2.0*np.pi)/float(nlev[-1]) if nlev[-1] else 0.1
    seps = [params['clade_sep']*sep/float(depth-i+1) for i in range(1,len(tree['nlev'])+1)]
    totseps = sum([s*nlev[i] for i,s in enumerate(seps[:-1])]) if len(nlev)>1 else 0
    clade_sep_err = (totseps>np.pi)
    while totseps>np.pi:
        params['clade_sep']*=0.75
        seps=[params['clade_sep']*sep/(float(depth-i+1)*0.25)
              for i in range(1,len(tree['nlev'])+1)]
        totseps=sum([s*nlev[i] for i,s in enumerate(seps[:-1])])
    if clade_sep_err:
        print("clade_sep parameter too large, lowered to", params['clade_sep'])

    # (A) Figure 생성 (넉넉히)
    fig = plt.figure(edgecolor=params['back_color'],facecolor=params['back_color'],
                     figsize=(12,8))

    # (B) 왼쪽(클라도그램), 오른쪽(레전드)
    ax_cladogram = fig.add_subplot(121, polar=True, frame_on=False,
                                   facecolor=params['back_color'])
    ax_legend = fig.add_subplot(122)
    ax_legend.axis("off")

    plt.subplots_adjust(left=params['l_prop'], right=0.98, wspace=0.05)

    ax_cladogram.grid(False)
    ax_cladogram.set_xticks([])
    ax_cladogram.set_yticks([])

    # (C) 클라도그램 그리기
    ds = (2.0*np.pi - totseps)/float(nlev[-1]) if nlev[-1] else 0.1
    add_all_pos(tree['root'], 0.0, ds, seps, 0.0, depth)
    plot_lines(tree['root'], params, depth, ax_cladogram, 0)
    plot_points(tree['root'], params, pt_scale, ax_cladogram)
    plot_names(tree['root'], params, depth, ax_cladogram, uniqueid(), seps)

    ax_cladogram.set_title(params['title'], size=params['title_font_size'],
                           color=params['fore_color'])

    # (D) 왼쪽 상단에 클래스 레전드(DQ, MW, TC)
    #   - tree['classes']에 저장된 클래스 중 DQ, MW, TC만 추출
    left_classes = ["DQ","MW","TC"]
    left_handles = []
    left_labels = []
    # 기존 코드에서 classes를 만들지만, 여기선 직접 DQ,MW,TC만 bar를 만들어 추가
    for i, c in enumerate(tree['classes']):
        if c in left_classes:
            # class 색상
            idx = i % len(colors)
            dummy_bar = ax_cladogram.bar(0,0,width=0,bottom=0,
                                         color=colors[idx], label=c)
            left_handles.append(dummy_bar)
            left_labels.append(c)

    if params['class_legend_vis'] and len(left_labels)>0:
        leg_left = ax_cladogram.legend(left_handles, left_labels,
                                       loc='upper left',
                                       bbox_to_anchor=(0.05,1.05),
                                       frameon=False,
                                       prop={'size': params['class_legend_font_size']},
                                       ncol=1)
        def get_col_attr(x):
            return hasattr(x, 'set_color') and not hasattr(x, 'set_facecolor')
        if leg_left is not None:
            for o in leg_left.findobj(get_col_attr):
                o.set_color(params['fore_color'])
        ax_cladogram.add_artist(leg_left)

    # (E) 오른쪽 레전드(박테리아 축약 라벨 + 그 외 클래스)
    #   1) 먼저 ax_cladogram에서 handles/labels 가져옴
    all_h, all_l = ax_cladogram.get_legend_handles_labels()

    #   2) 필터링: 빈 라벨, _nolegend 제외
    filtered = [(h,l) for (h,l) in zip(all_h, all_l)
                if l not in("", "_nolegend")]

    #   3) left_classes(DQ,MW,TC)는 이미 왼쪽 레전드로 표시했으므로 제외
    #      + 중복 제거
    final_right = []
    used_labels = set(left_labels)  # 왼쪽 레전드 라벨
    for (h,l) in filtered:
        # 만약 'DQ'... 등이랑 겹치지 않는 라벨이면 오른쪽에 추가
        if l not in used_labels:
            final_right.append((h,l))

    if len(final_right)>0:
        rh, rl = zip(*final_right)
        # 열 개수(ncol=2 -> 필요시 3,4)
        leg_right = ax_legend.legend(rh, rl, loc='upper left',
                                     bbox_to_anchor=(0,1),
                                     frameon=False,
                                     prop={'size': float(params['label_font_size'])},
                                     ncol=2)
        def get_col_attr(x):
            return hasattr(x, 'set_color') and not hasattr(x, 'set_facecolor')
        if leg_right is not None:
            for o in leg_right.findobj(get_col_attr):
                o.set_color(params['fore_color'])

    # (F) 저장
    plt.savefig(out_file, format=params['format'],
                facecolor=params['back_color'],
                edgecolor=params['fore_color'],
                dpi=params['dpi'],
                bbox_inches='tight')
    plt.close()

def plot_cladogram():
    params = read_params(sys.argv)
    params['fore_color'] = 'w' if params['back_color'] == 'k' else 'k'
    clad_tree = read_data(params['input_file'], params)
    draw_tree(params['output_file'], clad_tree, params)

if __name__ == '__main__':
    plot_cladogram()
