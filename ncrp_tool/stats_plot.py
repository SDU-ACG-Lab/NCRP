# stats_plot.py

from collections import Counter


def edge_overlap_histogram(adj, bin_width=100, out_path=None):
    """
    对已构建好的邻接表 adj（对称、无自环、每对仅保留 best overlap），
    统计“唯一无向边”的 overlap 直方图。仅计 u < v，避免双计。
    """
    max_ov = 0
    for u, nbrs in adj.items():
        for v, ov in nbrs.items():
            if u < v and ov > max_ov:
                max_ov = ov

    if max_ov <= 0:
        print("Edge overlap histogram: no edges to count.")
        return

    nbins = max(1, (max_ov + bin_width - 1) // bin_width)
    counts = [0] * nbins

    for u, nbrs in adj.items():
        for v, ov in nbrs.items():
            if u < v:
                idx = ov // bin_width
                if idx >= nbins:
                    idx = nbins - 1
                counts[idx] += 1

    print(f"Edge overlap histogram (unique undirected, bin={bin_width}bp):")
    print("bin_left\tbin_right\tcount")
    for i, c in enumerate(counts):
        left = i * bin_width
        right = (i + 1) * bin_width - 1
        print(f"{left}\t{right}\t{c}")

    if out_path:
        with open(out_path, "w", encoding="utf-8") as fout:
            fout.write("bin_left\tbin_right\tcount\n")
            for i, c in enumerate(counts):
                left = i * bin_width
                right = (i + 1) * bin_width - 1
                fout.write(f"{left}\t{right}\t{c}\n")
        print(f"Edge histogram written to: {out_path}")


def _save_bar_png(xs, heights, width, out_png, title, xlabel, ylabel="count", dpi=200, logy=False):
    import matplotlib.pyplot as plt

    plt.figure()
    plt.bar(xs, heights, width=width, align="edge")
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if logy:
        plt.yscale("log")
    plt.grid(True, axis="y")
    plt.tight_layout()
    plt.savefig(out_png, dpi=dpi)
    plt.close()
    print(f"Figure saved: {out_png}")


def compute_edge_overlap_hist(adj, bin_width=100):
    """返回 (bin_left_list, counts, bin_width)；只计唯一无向边（u < v）。"""
    max_ov = 0
    for u, nbrs in adj.items():
        for v, ov in nbrs.items():
            if u < v and ov > max_ov:
                max_ov = ov
    if max_ov <= 0:
        return [], [], bin_width

    nbins = max(1, (max_ov + bin_width - 1) // bin_width)
    counts = [0] * nbins
    for u, nbrs in adj.items():
        for v, ov in nbrs.items():
            if u < v:
                idx = ov // bin_width
                if idx >= nbins:
                    idx = nbins - 1
                counts[idx] += 1
    bin_left = [i * bin_width for i in range(nbins)]
    return bin_left, counts, bin_width


def plot_edge_overlap_histogram(adj, bin_width, out_png, dpi=200, logy=False):
    bin_left, counts, bw = compute_edge_overlap_hist(adj, bin_width=bin_width)
    if not counts:
        print("No edges to plot for overlap histogram.")
        return
    _save_bar_png(
        xs=bin_left,
        heights=counts,
        width=bw,
        out_png=out_png,
        title=f"Edge overlap histogram (bin={bw} bp)",
        xlabel="overlap length (bp)",
        dpi=dpi,
        logy=logy,
    )


def plot_degree_histogram_from_adj(adj, out_png, dpi=200, logy=False):
    """度分布（基于邻接表度数）。"""
    degs = [len(nbrs) for nbrs in adj.values()]
    if not degs:
        print("No nodes to plot for degree histogram.")
        return
    c = Counter(degs)
    xs = sorted(c.keys())
    ys = [c[d] for d in xs]
    _save_bar_png(
        xs=xs,
        heights=ys,
        width=1.0,
        out_png=out_png,
        title="Degree histogram",
        xlabel="degree",
        dpi=dpi,
        logy=logy,
    )


def plot_weight_histogram_from_graph(graph, bins=50, out_png="weights_hist.png", dpi=200, logy=False):
    """归一化权重直方图（w=overlap/Lmax）。"""
    import matplotlib.pyplot as plt

    vals = []
    for _, (_, wts) in graph.items():
        vals.extend(wts.tolist())
    if not vals:
        print("No weights to plot.")
        return
    plt.figure()
    plt.hist(vals, bins=bins)
    plt.title("Edge weight histogram")
    plt.xlabel("weight (normalized overlap)")
    plt.ylabel("count")
    if logy:
        plt.yscale("log")
    plt.grid(True, axis="y")
    plt.tight_layout()
    plt.savefig(out_png, dpi=dpi)
    plt.close()
    print(f"Figure saved: {out_png}")

