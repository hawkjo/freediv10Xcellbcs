import matplotlib.pyplot as plt


def knee_plot(counts, thresh=None, good_label='Cells'):
    counts.sort(reverse=True)
    if thresh is None:
        thresh = 50
    threshold_idx = next(i for i, v in counts if v < thresh) 
    fig, ax = plt.subplots()
    ax.plot(range(1, threshold_idx+1), counts[:threshold_idx], linewidth=2, label=good_label)
    ax.plot(range(threshold_idx+1, len(counts)+1), counts[threshold_idx:], color='grey', label='Background')
    ax.set_xlabel('Observed barcode')
    ax.set_ylabel('Count')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend()
    return fig, ax 
