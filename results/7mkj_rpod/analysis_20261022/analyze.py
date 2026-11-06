"""Read-only source audit and descriptive analysis of the completed priority scan."""
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics as st

OUT=Path(__file__).resolve().parent
ROOT=OUT.parent
AA='ACDEFGHIKLMNPQRSTVWY'

def read(path):
    with path.open(newline='') as f: return list(csv.DictReader(f,delimiter='\t'))

def write(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)

rows=read(ROOT/'selected_completed_predictions.tsv')
ids={r['Result_Id']:r for r in rows}
assert len(ids)==len(rows)==665
for r in rows:
    assert r['Chain']=='L' and r['Status']=='complete'
    r['site']=int(r['Mutation'][1:-1]); r['ddg']=float(r['DeltaDeltaG_kcal_mol'])
    assert math.isfinite(r['ddg']) and math.isfinite(float(r['PROVEAN_full_sequence']))
sites=sorted({r['site'] for r in rows})
assert sites==list(range(375,391))+list(range(435,454))
for site in sites:
    items=[r for r in rows if r['site']==site]
    assert len(items)==19 and len({r['Mutation'][0] for r in items})==1
    assert {r['Mutation'][-1] for r in items}==set(AA)-{items[0]['Mutation'][0]}

plan=json.loads((ROOT/'selection_plan.json').read_text())
features=[]; seen=set()
for item in plan['batches']:
    base=ROOT/'batches'/f"batch_{item['batch']:04d}"
    assert json.loads((base/'status.json').read_text())['stage']=='complete'
    raw=read(base/'results.tsv'); ft=read(base/'features.tsv')
    assert len(raw)==len(ft)==19
    for r in raw:
        rid=r['Result_Id']; assert rid not in seen; seen.add(rid)
        assert ids[rid]['Mutation']==r['Mutation_PDB']
        assert ids[rid]['ddg']==float(r['SunDDG'])
        assert ids[rid]['Model_Confidence']==r['Confidence']
    features.extend(ft)
assert seen==set(ids)=={r['Result_Id'] for r in features}
write(OUT/'all_features.tsv',features)

def stats(group):
    v=[r['ddg'] for r in group]
    return {'count':len(v),'min':min(v),'median':st.median(v),'mean':st.mean(v),'max':max(v),
            'negative':sum(x<0 for x in v),'zero':sum(x==0 for x in v),
            'ge_1_5':sum(x>=1.5 for x in v),'interface_yes':sum(r['Model_Interface']=='yes' for r in group),
            'near_truncation':sum(r['Near_truncated_terminus_10A']=='1' for r in group)}

site_summary=[]
for site in sites:
    group=[r for r in rows if r['site']==site]; v=[r['ddg'] for r in group]
    site_summary.append({'Site':site,'WT':group[0]['Mutation'][0],'N':19,'Median_ddG':st.median(v),
        'Mean_ddG':round(st.mean(v),4),'Minimum_ddG':min(v),'Maximum_ddG':max(v),
        'Count_ge_1_5':sum(x>=1.5 for x in v),
        'Largest_effect':','.join(r['Mutation'] for r in group if r['ddg']==max(v)),
        'Smallest_effect':','.join(r['Mutation'] for r in group if r['ddg']==min(v)),
        'Partner_distance_A':float(group[0]['Partner_min_heavy_distance_A'])})
write(OUT/'site_summary.tsv',site_summary)
larger=sorted([r for r in rows if r['ddg']>=1.5],key=lambda r:(-r['ddg'],r['Mutation']))
smaller=sorted([r for r in rows if r['ddg']<0.5 and float(r['PROVEAN_full_sequence'])>-2.5],key=lambda r:r['ddg'])
summary={'verification':'ANALYZED; native calculations not rerun',
    'all':stats(rows),'375-390':stats([r for r in rows if r['site']<400]),
    '435-453':stats([r for r in rows if r['site']>400]),
    'larger_effects':larger,'small_effect_and_provean_above_minus_2_5':smaller,
    'provean_le_minus_2_5':sum(float(r['PROVEAN_full_sequence'])<=-2.5 for r in rows),
    'site_summary':site_summary,'source_sha256':{
        name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in
        ('selected_completed_predictions.tsv','selection_plan.json','site_annotations.tsv','batch_status.json')}}
(OUT/'analysis.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')

def table(data,fields):
    return '\n'.join(['| '+' | '.join(fields)+' |','| '+' | '.join('---' for _ in fields)+' |']+
                     ['| '+' | '.join(str(r[k]) for k in fields)+' |' for r in data])

report=f'''# RpoD 两个重点区段预测解读

## Material Passport

验证状态：ANALYZED。分析日期：2026-10-22。分析单位为每个独立单点替换，共665项；不是665个实验重复。未重新执行物理计算；逐项核对了35个原生批次的结果与特征，检查665个请求全部覆盖、无重复、每个位点19种非WT替换。

## 结果位置与完成范围

本次结果：[selected_completed_predictions.tsv](../selected_completed_predictions.tsv)。35个位点，L链375–390与435–453，每个位点19种替换。北京时间2026-10-22 01:34:47完成。`selected_predictions.tsv`为同范围的状态表。

`completed_predictions.tsv`另含3个早期试算，共668行；`predictions.tsv`包含全量8949项请求，未完成行的空白值不是0。不能把本次完成理解成全量8949项已完成。

新入口：`/home/zek/mutabind2-tool/projects/7mkj_rpod/`。它当前仍是指向原任务目录的兼容链接。

## 主要结论

665个预测全部为正，范围+0.11至+2.35 kcal/mol，中位数+0.71。此次未发现模型预测的结合增强突变。按现有打分程序的1.5 kcal/mol幅度阈值，24项（{24/665:.1%}）属于较大影响候选；另外641项低于该阈值，不能据此认定它们完全无影响。

| 区段 | 位点/突变数 | ΔΔG中位数 | 范围 | ≥1.5 kcal/mol |
|---|---:|---:|---:|---:|
| 375–390 | 16 / 304 | +0.70 | +0.11～+1.73 | 5（{5/304:.1%}） |
| 435–453 | 19 / 361 | +0.71 | +0.13～+2.35 | 19（{19/361:.1%}） |

两个区段的中位数很接近；435–453的较大影响候选更多，主要集中在少数位点。它不是整个区段都更敏感的证据，没有进行区段间显著性检验。

## 怎样理解分数

ΔΔG = ΔG结合(突变体) − ΔG结合(WT)。正值预测结合减弱；负值预测结合增强；绝对值反映模型估计变化幅度。本次问题是L链RpoD与I+J伙伴的结合，不直接预测整体折叠稳定性、转录活性、启动子识别、细胞生长或抗性。

`Model_Confidence` 的 High/Low 在代码中由 `interface` yes/no直接赋值，未给出每项预测的概率或置信区间。665项中284项为High、381项为Low。1.5是程序分类阈值，不是p值，也不是实验显著性标准；1.49与1.51的差别不应过度解读。

PROVEAN来自完整613 aa序列；≤−2.5是常用的序列影响筛选界限。本次642/665低于或等于该界限，表明这些替换常被该保守性指标判为有潜在功能影响。PROVEAN已经是MutaBind2模型输入，不能将二者一致当成独立工具的重复验证。PROVEAN也不提供结合增强或减弱的方向。

## 优先关注的位点

I450最突出：19种替换中12种≥1.5，位点中位数+1.62，最小也为I450L的+0.94；I450V为+1.10。I450K最高，为+2.35。L384有4种较大影响替换，I452有3种，P453有2种；V380、Q446、T449各1种。全部24项的界面标签为yes。

较大影响候选如下（用于选择实验验证对象，不能直接称为功能缺失突变）：

{table(larger,['Mutation','DeltaDeltaG_kcal_mol','PROVEAN_full_sequence','Model_Confidence'])}

若想验证结合减弱，可优先考虑I450K、T449R、L384E和I452W，覆盖不同位点。I450A可作为同一热点的另一种替换，便于判断是否存在位点效应。实际实验还需检查表达量、可溶性和折叠，避免将蛋白失稳误解为界面作用。

## 较小结合影响的对照候选

仅作为探索性筛选，使用ΔΔG<0.5且PROVEAN>−2.5得到以下8项；界限在查看结果后用于整理候选，不代表统计检验或独立验证：

{table(smaller,['Mutation','DeltaDeltaG_kcal_mol','PROVEAN_full_sequence','Model_Confidence'])}

V380I（+0.18、PROVEAN −0.572、界面yes）可优先作为较小影响对照，与WT和较大影响候选一起验证。它不是预测的增强突变。E381M虽为全表最低ΔΔG（+0.11），PROVEAN为−6.153，不宜只按最低ΔΔG称为最佳或最安全突变。

## 结构和模型限制

本次使用7个断开的连续结构片段及中性ACE/ACP、CT2封端，没有补建缺失部分。保留原随机森林与力场，复用固定WT基线。小分子和金属未纳入。预测属于现有结构条件下的近似筛选，不是完整蛋白的精确自由能。

这35个位点均未被标记为距任意截断端≤10 Å，也没有位于截断端本身；这是有利的局部几何条件，但不能消除缺失片段的远程影响或结构状态差异。原有两个端基敏感性对照差异0.14和0.42 kcal/mol，不能推广成所有突变的误差条。

所有分数为正是本次模型的输出，不能推出真实蛋白已达到最优结合状态，也不能排除模型偏移。原始ddVdw、ddPb和FoldX特征均包含正负值；预测代码直接调用随机森林，没有用绝对值把负预测改为正预测。独立模型打分复核见`rescore_audit.json`。

## 图与可复核数据

[热图](heatmap.png)：两段使用相同0–2.5 kcal/mol色标；灰色×是WT自身替换、未计算，不代表0。空心圈标出≥1.5；两段之间未展示的391–434不在本次范围内。

[35个位点汇总](site_summary.tsv)，以及机器可读统计`analysis.json`。原始表不改写，变换仅包括指定区段筛选、按位点汇总及排序；未删除离群值、插值、平滑或进行显著性检验。

## 解释偏差检查（11/11）

1. 辛普森悖论：分别报告两个区段，不将汇总排序当作分层因果结论。
2. 生态谬误：区段均值不能替代具体突变结果。
3. 选择偏差：只覆盖预先指定35个位点，不外推为全RpoD。
4. 碰撞变量偏差：未作条件回归或因果推断。
5. 基率忽视：没有已校准预测概率；High不视为正确率。
6. 回归均值：极端候选需独立验证，不宣称其排名稳定。
7. 幸存者偏差：665/665覆盖，无缺失预测被排除。
8. 多重搜索：披露全部665项及24项阈值筛选，不声称p值或显著性。
9. 分析路径：结构和封端选择明确；候选阈值整理属于探索性。
10. 相关与因果：计算预测不等于实验证明结合/功能改变。
11. 反向因果：未作观察性因果推断，也不以保守性解释唯一分子机制。
'''
(OUT/'结果解读.md').write_text(report)
print(json.dumps({k:summary[k] for k in ('all','375-390','435-453')},indent=2))

if __name__=='__main__' and '--plot' in __import__('sys').argv:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.colors import Normalize
    from matplotlib.lines import Line2D
    style={'font.family':'Noto Sans CJK JP','font.size':11,'axes.unicode_minus':False,
           'pdf.fonttype':42,'svg.fonttype':'none'}
    with plt.rc_context(style):
        fig,axes=plt.subplots(1,2,figsize=(16,8.8),gridspec_kw={'width_ratios':[16,19]},layout='constrained')
        for ax,(lo,hi) in zip(axes,[(375,390),(435,453)]):
            ss=list(range(lo,hi+1)); mat=np.full((20,len(ss)),np.nan)
            for r in rows:
                if lo<=r['site']<=hi: mat[AA.index(r['Mutation'][-1]),r['site']-lo]=r['ddg']
            image=ax.imshow(mat,cmap=plt.get_cmap('YlOrBr').with_extremes(bad='#dddddd'),
                            norm=Normalize(0,2.5),interpolation='nearest',aspect='auto')
            yy,xx=np.where(np.isnan(mat)); ax.scatter(xx,yy,marker='x',c='#555555',s=20,linewidths=.9)
            yy,xx=np.where(mat>=1.5);ax.scatter(xx,yy,s=65,facecolors='none',edgecolors='white',linewidths=1.2)
            labels=[next(s['WT'] for s in site_summary if s['Site']==n)+str(n) for n in ss]
            ax.set_xticks(range(len(ss)),labels,rotation=90,fontsize=10)
            ax.set_yticks(range(20),list(AA));ax.set_xlabel('原始 L 链位点（WT + 编号）')
            ax.set_title(f'{lo}–{hi}：{len(ss)} 位点，{len(ss)*19} 个独立替换',pad=12)
            ax.set_xticks(np.arange(-.5,len(ss),1),minor=True); ax.set_yticks(np.arange(-.5,20,1),minor=True)
            ax.grid(which='minor',color='white',linewidth=.5);ax.tick_params(which='minor',length=0)
        axes[0].set_ylabel('突变后的氨基酸（单字母）')
        cb=fig.colorbar(image,ax=axes,fraction=.025,pad=.02)
        cb.set_label('预测结合 ΔΔG（kcal/mol）；正值 = 结合减弱')
        fig.suptitle('RpoD 重点区段：665 个预测均为正值，较大影响集中于 I450 等位点\n'
                     '灰色 ×：WT 自身替换未计算；空心圈：ΔΔG ≥ 1.5（筛选阈值，非统计显著性）',fontsize=15)
        fig.savefig(OUT/'heatmap.png',dpi=180,facecolor='white')
        fig.savefig(OUT/'heatmap.pdf',facecolor='white')
        plt.close(fig)
