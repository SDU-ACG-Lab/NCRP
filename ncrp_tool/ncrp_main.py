# main.py

import argparse
import sys
import time
import resource

from idmap import IdMap
from kraken_io import load_labels_and_order
from graph_io import parse_paf_to_adj, parse_edges_to_adj, finalize_graph
from stats_plot import (
    edge_overlap_histogram,
    plot_edge_overlap_histogram,
    plot_degree_histogram_from_adj,
    plot_weight_histogram_from_graph,
)
from refine import refine_labels_linear
from propagate import label_propagation_layered


def main():
    ap = argparse.ArgumentParser(
        description="NCRP label refinement & propagation (PAF/edge-list 兼容，带边统计与直方图)"
    )
    ap.add_argument("--paf", help="minimap2 PAF")
    ap.add_argument("--edges", help="简化边文件：read1 read2 overlap")
    ap.add_argument("--kraken", required=True, help="Kraken2 per-read 输出")
    ap.add_argument("--min-overlap", type=int, default=130, help="overlap 阈值（bp）小于此值的边将删除")
    ap.add_argument("--max-iter", type=int, default=20, help="（保留参数）")
    ap.add_argument("--hist-bin", type=int, default=100, help="overlap 直方图分箱宽度（bp）")
    ap.add_argument("--hist-out", default=None, help="直方图 TSV 输出路径，可选")
    ap.add_argument("--output", default="final_labels.tsv", help="最终标签输出文件")

    # 图像输出选项
    ap.add_argument("--plot-overlap-png", default=None, help="输出 overlap 直方图 PNG 路径")
    ap.add_argument("--plot-degree-png", default=None, help="输出度分布直方图 PNG 路径")
    ap.add_argument(
        "--plot-weight-png",
        default=None,
        help="输出归一化权重直方图 PNG 路径（finalize 后）",
    )
    ap.add_argument("--dpi", type=int, default=300, help="PNG DPI")
    ap.add_argument("--logy", action="store_true", help="y 轴使用对数坐标")
    args = ap.parse_args()

    rid_map = IdMap()   # read id
    tax_map = IdMap()   # taxid（字符串）也压缩为 int

    # 读 Kraken（先建 rid 映射与初始标签；order 保存输出顺序）
    print("Loading Kraken2 labels ...")
    init_labels_int, order = load_labels_and_order(args.kraken, rid_map, tax_map)
    print(f"Initial labeled reads: {len(init_labels_int)} | All reads in Kraken: {len(order)}")

    # 读图（边读边累加 best overlap），并统计阈值筛选前后数量
    if args.edges:
        print(f"Parsing edge list: {args.edges}")
        adj, lengths, estats = parse_edges_to_adj(args.edges, rid_map, min_overlap=args.min_overlap)
    elif args.paf:
        print(f"Parsing PAF: {args.paf}")
        adj, lengths, estats = parse_paf_to_adj(args.paf, rid_map, min_overlap=args.min_overlap)
    else:
        raise ValueError("必须提供 --paf 或 --edges")

    # 打印边统计信息
    print(
        "Edge stats | "
        f"raw records: {estats['raw']}  |  "
        f"removed (<{args.min_overlap}): {estats['too_short']}  |  "
        f"self-loops skipped: {estats['self_loop']}  |  "
        f"passed threshold (records): {estats['kept_records']}  |  "
        f"unique undirected edges (kept): {estats['unique_edges']}"
    )

    # 文本版直方图（唯一无向边，已去自环）
    edge_overlap_histogram(adj, bin_width=args.hist_bin, out_path=args.hist_out)

    # PNG 版直方图
    if args.plot_overlap_png:
        plot_edge_overlap_histogram(
            adj,
            bin_width=args.hist_bin,
            out_png=args.plot_overlap_png,
            dpi=args.dpi,
            logy=args.logy,
        )
    if args.plot_degree_png:
        plot_degree_histogram_from_adj(
            adj,
            out_png=args.plot_degree_png,
            dpi=args.dpi,
            logy=args.logy,
        )

    print(f"Nodes(with degree>0): {len(adj)}")

    # 固化邻接并做权重归一化
    graph, Lmax = finalize_graph(adj, lengths, default_len=1000)
    print(f"Graph fixed. Lmax={Lmax}")

    # 权重直方图（可选）
    if args.plot_weight_png:
        plot_weight_histogram_from_graph(
            graph,
            bins=50,
            out_png=args.plot_weight_png,
            dpi=args.dpi,
            logy=args.logy,
        )

    # refine
    print("Refining labels ...")
    labels_after_refine = refine_labels_linear(graph, init_labels_int)

    # 传播（分层）
    print("Running label propagation ...")
    final_labels_int = label_propagation_layered(graph, labels_after_refine)

    # 输出：严格按 Kraken 出现顺序；无法标注输出 U
    print(f"Writing: {args.output}")
    with open(args.output, "w", encoding="utf-8") as fout:
        for rid in order:
            lab_int = final_labels_int.get(rid, labels_after_refine.get(rid))
            if lab_int is None:
                fout.write(f"U\t{rid_map.get_str(rid)}\t0\n")
            else:
                fout.write(f"C\t{rid_map.get_str(rid)}\t{tax_map.get_str(lab_int)}\n")
    print("Done.")


def _unit_kb(ru_maxrss):
    # Linux: ru_maxrss 单位=KB；macOS: 单位=Bytes
    return ru_maxrss if sys.platform.startswith("linux") else ru_maxrss / 1024.0


if __name__ == "__main__":
    t0 = time.perf_counter()
    main()
    t1 = time.perf_counter()
    ru_self = _unit_kb(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    ru_child = _unit_kb(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    peak_kb = max(ru_self, ru_child)
    print(f"[PROF] Elapsed: {t1 - t0:.2f}s  |  Peak RSS: {peak_kb/1024:.2f} MB")

