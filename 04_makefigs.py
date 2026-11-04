"""Normalisation, PCA and gene-level figures.
Reads the gene tables from 03_makecounts.py. Computes DESeq2-style
median-of-ratios size factors, then writes figures and relabelled tables."""
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
H=os.path.expanduser("~"); OUT=H+"/Desktop/RNAseq_for_meeting"; os.makedirs(OUT,exist_ok=True)
# old core name, new name, treatment, library prep run, original bench tube
S=[("3818_Sample1","Vehicle_1","Vehicle","original","V1"),
("3818_Sample3","Vehicle_2","Vehicle","original","V3"),
("3818_Sample4","Vehicle_3","Vehicle","original","V4"),
("3818_Sample5","Vehicle_4","Vehicle","original","V5"),
("3818_Sample8","Vehicle_5","Vehicle","original","V8"),
("3818_Sample9","SB273005_1","SB273005","original","300-1"),
("3818_Sample12","SB273005_2","SB273005","original","300-4"),
("3818_Sample13","SB273005_3","SB273005","original","300-5"),
("3818_1_300_T3","SB273005_4","SB273005","replacement","1_300_T3"),
("3818_5_T3","SB273005_5","SB273005","replacement","5_T3"),
("3818_6_T3","SB273005_6","SB273005","replacement","6_T3")]
OLD=[a for a,b,c,d,e in S]; NEW=[b for a,b,c,d,e in S]
RN=dict(zip(OLD,NEW)); G={b:c for a,b,c,d,e in S}; B={b:d for a,b,c,d,e in S}
V,D="#2a78d6","#eb6834"; INK,MU,GR="#1a1a18","#5f5e57","#e1e0d9"
plt.rcParams.update({"font.size":9,"axes.edgecolor":"#c3c2b7","text.color":INK,
"axes.labelcolor":INK,"xtick.color":MU,"ytick.color":MU,"figure.facecolor":"white",
"axes.spines.top":False,"axes.spines.right":False})
co=lambda s: V if G[s]=="Vehicle" else D
mk=lambda s: "o" if B[s]=="original" else "^"
cf=pd.read_csv(H+"/postn_rnaseq/gene_counts.csv"); tf=pd.read_csv(H+"/postn_rnaseq/gene_tpm.csv")
sym=cf.set_index("gene_id")["symbol"]
c=cf.set_index("gene_id").drop(columns=["symbol"])[OLD].rename(columns=RN)
# keep genes with at least 10 counts in at least 3 samples
c=c.loc[(c>=10).sum(axis=1)>=3]; print("genes after filtering:",len(c))
# median-of-ratios size factors, as used by DESeq2
lg=np.log(c.replace(0,np.nan)); sf=np.exp((lg.sub(lg.mean(axis=1),axis=0)).median())
rM=c.sum()/1e6; nM=(c/sf).sum()/1e6
print("RAW        %.1f-%.1f M   spread %.1f%%"%(rM.min(),rM.max(),100*rM.std()/rM.mean()))
print("NORMALISED %.1f-%.1f M   spread %.1f%%"%(nM.min(),nM.max(),100*nM.std()/nM.mean()))
k=pd.DataFrame(S,columns=["name_from_core","sample_name","treatment","library_prep","bench_tube"])
k["dose"]=np.where(k.treatment=="Vehicle","0 (DMSO)","300 nM")
k=k[["sample_name","treatment","dose","library_prep","bench_tube","name_from_core"]]
k.to_csv(OUT+"/0_SAMPLE_KEY.csv",index=False); print(); print(k.to_string(index=False))
t=tf.set_index(["gene_id","symbol"])[OLD].rename(columns=RN)
t.insert(0,"mean_Vehicle",t[[s for s in NEW if G[s]=="Vehicle"]].mean(axis=1).round(2))
t.insert(1,"mean_SB273005",t[[s for s in NEW if G[s]=="SB273005"]].mean(axis=1).round(2))
t.sort_values("mean_Vehicle",ascending=False).round(2).to_csv(OUT+"/1_expression_TPM.csv")
nc=(c/sf).round(1); nc.insert(0,"gene_symbol",sym.reindex(nc.index).values)
nc.to_csv(OUT+"/2_counts_normalised.csv")
pd.DataFrame({"sample":NEW,"treatment":[G[s] for s in NEW],
"raw_millions":[round(rM[s],2) for s in NEW],"size_factor":[round(sf[s],3) for s in NEW],
"normalised_millions":[round(nM[s],2) for s in NEW]}).to_csv(OUT+"/6_normalisation_table.csv",index=False)
x=np.arange(len(NEW))
f,ax=plt.subplots(1,2,figsize=(12.5,4.5),sharey=True)
for a,vals,ti in [(ax[0],rM,"BEFORE normalising"),(ax[1],nM,"AFTER normalising")]:
    a.grid(color=GR,lw=.7,axis="y",zorder=0); a.set_axisbelow(True)
    a.bar(x,[vals[s] for s in NEW],color=[co(s) for s in NEW],width=.68,zorder=3)
    a.axhline(np.mean([vals[s] for s in NEW]),color=INK,lw=1.2,ls="--",zorder=4)
    a.set_xticks(x); a.set_xticklabels(NEW,rotation=45,ha="right",fontsize=8.5)
    a.set_title("%s  -  spread %.1f%%"%(ti,100*vals.std()/vals.mean()),loc="left")
ax[0].set_ylabel("total counts, millions")
ax[0].legend(handles=[plt.Rectangle((0,0),1,1,fc=V),plt.Rectangle((0,0),1,1,fc=D)],
labels=["Vehicle","SB273005 300 nM"],frameon=False,fontsize=8.5,labelcolor=INK)
plt.tight_layout(); plt.savefig(OUT+"/3_figure_normalisation.png",dpi=200); plt.close()
ln=np.log2(c/sf+1); vv=ln.var(axis=1).sort_values(ascending=False)
Z=ln.loc[vv.index[:2000]].T.values; Z=Z-Z.mean(axis=0)
U,SS,_=np.linalg.svd(Z,full_matrices=False); p=U*SS; ev=SS**2/(SS**2).sum()*100
f,a=plt.subplots(figsize=(7,5.5)); a.grid(color=GR,lw=.7,zorder=0); a.set_axisbelow(True)
for i,s in enumerate(NEW):
    a.scatter(p[i,0],p[i,1],s=150,zorder=3,c=co(s),marker=mk(s),edgecolors="white",linewidths=1.8)
    a.annotate(s,(p[i,0],p[i,1]),fontsize=7.5,color=MU,xytext=(0,-15),
               textcoords="offset points",ha="center")
a.set_xlabel("PC1 - %.0f%% of variance"%ev[0]); a.set_ylabel("PC2 - %.0f%%"%ev[1])
a.set_title("Samples group by library prep, not treatment",loc="left")
a.legend(handles=[Line2D([],[],marker="o",ls="",mfc=V,mec="w",mew=1.4,ms=9,label="Vehicle"),
Line2D([],[],marker="o",ls="",mfc=D,mec="w",mew=1.4,ms=9,label="SB273005"),
Line2D([],[],marker="o",ls="",mfc=MU,mec="w",mew=1.4,ms=9,label="original prep"),
Line2D([],[],marker="^",ls="",mfc=MU,mec="w",mew=1.4,ms=9,label="replacement prep")],
frameon=False,fontsize=7.5,labelcolor=INK)
plt.tight_layout(); plt.savefig(OUT+"/4_figure_PCA.png",dpi=200); plt.close()
print("\nPC1 %.1f%%  PC2 %.1f%%"%(ev[0],ev[1]))
GN=["Postn","Itgav","Itgb1","Itgb3","Itgb5","Fn1","Col1a1","Runx2"]
rg=np.random.default_rng(0); f,axs=plt.subplots(2,4,figsize=(13,6.2))
for a,g in zip(axs.ravel(),GN):
    a.grid(color=GR,lw=.7,axis="y",zorder=0); a.set_axisbelow(True)
    r=t[t.index.get_level_values("symbol")==g]
    if len(r)==0: a.axis("off"); continue
    vl=r.iloc[0]
    for xi,gr,cc in [(0,"Vehicle",V),(1,"SB273005",D)]:
        ss=[s for s in NEW if G[s]==gr]
        for s in ss:
            a.scatter(xi+rg.uniform(-.1,.1),vl[s],s=55,c=cc,zorder=3,marker=mk(s),
                      edgecolors="white",linewidths=1.3)
        m=np.mean([vl[s] for s in ss]); a.plot([xi-.24,xi+.24],[m,m],color=INK,lw=1.8,zorder=4)
    a.set_title(g,loc="left",fontsize=10); a.set_xticks([0,1])
    a.set_xticklabels(["Vehicle","SB273005"],fontsize=8.5,color=INK)
    a.set_xlim(-.5,1.5); a.set_ylabel("TPM"); a.set_ylim(bottom=0)
f.suptitle("Genes of interest - each point is one sample, bar is group mean. Triangles = replacement prep.",
x=.006,ha="left",fontsize=10,fontweight="bold")
plt.tight_layout(rect=[0,0,1,.95]); plt.savefig(OUT+"/5_figure_genes.png",dpi=200); plt.close()
print("\nWritten to",OUT)
