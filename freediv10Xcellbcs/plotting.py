import matplotlib.pyplot as plt


def knee_plot(counts, thresh=None, good_label='Cells'):
    counts.sort(reverse=True)
    if thresh is None:
        thresh = 50
    threshold_idx = next(i for i, v in enumerate(counts) if v < thresh) 
    fig, ax = plt.subplots()
    ax.plot(range(1, threshold_idx+1), counts[:threshold_idx], linewidth=2, label=good_label)
    ax.plot(range(threshold_idx+1, len(counts)+1), counts[threshold_idx:], color='grey', label='Background')
    ax.set_xlabel('Observed barcode')
    ax.set_ylabel('Count')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend()
    return fig, ax 

def umi_len_plot(umi_len_cntr):
    fs = 12
    x = list(umi_len_cntr.keys())
    y = [umi_len_cntr[xx] for xx in x]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), gridspec_kw=dict(wspace=0.3))
    axes[0].stem(x, y, basefmt='C0:')
    axes[1].scatter(x, y)
    for ax in axes:
        ax.set_ylabel('Count', fontsize=fs)
        ax.set_xlabel('UMI length', fontsize=fs)
        ax.set_title(sample)
    ax = axes[1]
    ax.set_yscale('log')
    ylim = ax.get_ylim()
    ax.plot([10]*2, ylim, 'C0:', zorder=-1, linewidth=1)
    ax.set_ylim(ylim)
    return fig, ax
