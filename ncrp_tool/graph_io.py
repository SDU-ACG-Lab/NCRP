# graph_io.py

from collections import defaultdict
from array import array
from idmap import IdMap


def parse_paf_to_adj(paf_file, rid_map: IdMap, min_overlap=130):
    """
    读取 PAF，构建对称邻接表 adj[u][v] = best_overlap。
    返回：adj, lengths, stats（包含 raw/too_short/self_loop/kept_records/unique_edges）
    """
    adj = defaultdict(dict)
    lengths = {}
    stats = dict(raw=0, too_short=0, self_loop=0, kept_records=0, unique_edges=0)

    with open(paf_file, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line or line[0] == "#" or line == "\n":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue

            q, t = fields[0], fields[5]
            try:
                qlen = int(fields[1])
                tlen = int(fields[6])
            except Exception:
                continue

            # overlap 估计：优先用 col[9]（nmatch），否则回退到 qEnd - qStart
            ov = None
            try:
                ov = int(fields[9])
            except Exception:
                try:
                    ov = max(0, int(fields[3]) - int(fields[2]))
                except Exception:
                    ov = None
            if ov is None:
                continue

            stats["raw"] += 1

            if ov < min_overlap:
                stats["too_short"] += 1
                continue

            uq = rid_map.get_int(q)
            vt = rid_map.get_int(t)
            if uq == vt:
                stats["self_loop"] += 1
                continue

            stats["kept_records"] += 1

            # 对称存储 best overlap（取最大）
            prev = adj[uq].get(vt)
            if prev is None or ov > prev:
                adj[uq][vt] = ov
                adj[vt][uq] = ov

            lengths[uq] = qlen
            lengths[vt] = tlen

    # 统计唯一无向边（对称存一半）
    undirected = sum(len(nbrs) for nbrs in adj.values()) // 2
    stats["unique_edges"] = undirected
    return adj, lengths, stats


def parse_edges_to_adj(edge_file, rid_map: IdMap, min_overlap=130):
    """
    读取简单边文件：r1 r2 overlap（整数），构建对称邻接表。
    """
    adj = defaultdict(dict)
    lengths = {}  # edge-list 无长度，留空即可
    stats = dict(raw=0, too_short=0, self_loop=0, kept_records=0, unique_edges=0)

    with open(edge_file, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            parts = line.split()
            if len(parts) != 3:
                continue
            r1, r2, ov_s = parts
            try:
                ov = int(ov_s)
            except Exception:
                continue

            stats["raw"] += 1

            if ov < min_overlap:
                stats["too_short"] += 1
                continue

            u = rid_map.get_int(r1)
            v = rid_map.get_int(r2)
            if u == v:
                stats["self_loop"] += 1
                continue

            stats["kept_records"] += 1

            prev = adj[u].get(v)
            if prev is None or ov > prev:
                adj[u][v] = ov
                adj[v][u] = ov

    undirected = sum(len(nbrs) for nbrs in adj.values()) // 2
    stats["unique_edges"] = undirected
    return adj, lengths, stats


def finalize_graph(adj, lengths, default_len=1000):
    """
    把 dict(dict) 变为紧凑结构：
      graph[u] = (array('I', nei_ids), array('f', weights))
    """
    Lmax = max(lengths.values()) if lengths else default_len
    graph = {}
    for u, nbrs in adj.items():
        if not nbrs:
            continue
        neis = array("I")
        wts = array("f")
        for v, ov in nbrs.items():
            neis.append(v)
            wts.append(ov / Lmax)
        graph[u] = (neis, wts)
    return graph, Lmax

