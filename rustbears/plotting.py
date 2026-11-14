"""Matplotlib ancestral-range and stochastic-map summaries."""
import numpy as np


def plot_fit(fit,filepath,plot_what="pie",show_corners=True,title=None):
    import matplotlib.pyplot as plt
    from matplotlib.offsetbox import AnnotationBbox, DrawingArea
    from matplotlib.patches import Wedge

    if plot_what not in ("pie","text"):
        raise ValueError("plot_what must be 'pie' or 'text'")
    model=fit.model
    t=model.tree
    p=fit.evaluation
    y=np.zeros(len(t.labels))
    y[:t.tip_count]=np.arange(t.tip_count)[::-1]
    for node,_ in t.edges[::2]:
        y[node]=y[list(t.children[node])].mean()
    fig,ax=plt.subplots(figsize=(9,max(4,min(t.tip_count*.28,20))))
    try:
        colors=plt.get_cmap("turbo")(np.linspace(.04,.96,len(p.ranges_list)))

        def annotate(x,y,probs,small=False):
            if not np.isfinite(probs).all() or probs.sum()<=0:
                return
            if plot_what=="text":
                ax.annotate(p.ranges_list[np.argmax(probs)],(x,y),xytext=(1,3),
                            textcoords="offset points",fontsize=6 if small else 8)
            else:
                r=4 if small else 7
                art=DrawingArea(2*r,2*r)
                angle=0
                for probability,color in zip(probs,colors):
                    if probability>0:
                        art.add_artist(Wedge((r,r),r,angle,angle+360*probability,facecolor=color))
                        angle+=360*probability
                ax.add_artist(AnnotationBbox(art,(x,y),frameon=False,pad=0))

        for node,child in t.edges:
            a,b=t.ages[node],t.ages[child]
            ax.plot([a,a,b],[y[node],y[child],y[child]],color="0.2",lw=.7)
            if show_corners:
                annotate(a,y[child],p.ancestral_probs_branch_bottom[child],True)
        for node in range(len(t.labels)):
            if node<t.tip_count:
                ax.annotate(t.labels[node],(t.ages[node],y[node]),xytext=(8,0),
                            textcoords="offset points",fontsize=7,va="center")
            else:
                annotate(t.ages[node],y[node],p.ancestral_probs_branch_top[node])
        for age in model.times:
            if age<t.ages[t.root]:
                ax.axvline(age,ls="--",color="0.6",lw=.8)
        ax.set_xlim(t.ages[t.root]*1.04,-t.ages[t.root]*.22)
        ax.set_xlabel("Time before present")
        ax.set_yticks([])
        ax.set_title(title or f"{model.model.name} ancestral ranges")
        fig.savefig(filepath,bbox_inches="tight")
    finally:
        plt.close(fig)
    return filepath


def plot_bsm(result,filepath):
    import matplotlib.pyplot as plt

    counts=result.summary_counts["counts_list"]
    fig,axes=plt.subplots(2,4,figsize=(11,6))
    try:
        for ax,event in zip(axes.flat,counts.columns):
            ax.hist(counts[event],bins="auto",color="#216c83",edgecolor="white")
            ax.set_title(f"{event} events")
            ax.set_xlabel("Count per map")
            ax.set_ylabel("Maps")
        axes.flat[-1].set_visible(False)
        fig.tight_layout()
        fig.savefig(filepath)
    finally:
        plt.close(fig)
    return filepath
