# propagate.py

from collections import deque, defaultdict


def label_propagation_layered(graph, init_labels_int):
    """
    分层传播：只从上一层（dist=d-1）的已标注邻居聚合权重并赋给当前层。
    """
    labels = dict(init_labels_int)
    if not labels:
        return labels

    INF = 10**15
    N = max(graph.keys(), default=0) + 1
    dist = [INF] * N
    q = deque()

    for s in labels.keys():
        if s < N:
            dist[s] = 0
            q.append(s)

    # 多源 BFS 计算层次
    while q:
        u = q.popleft()
        neis, _ = graph.get(u, ((), ()))
        du = dist[u]
        for j in range(len(neis)):
            v = neis[j]
            if v >= N:
                continue
            if dist[v] == INF:
                dist[v] = du + 1
                q.append(v)

    # 按层赋值
    buckets = defaultdict(list)
    maxd = 0
    for v in graph.keys():
        if v < N and dist[v] != INF and dist[v] > 0:
            buckets[dist[v]].append(v)
            if dist[v] > maxd:
                maxd = dist[v]

    for d in range(1, maxd + 1):
        for v in buckets.get(d, []):
            neis, wts = graph.get(v, ((), ()))
            acc = {}
            for j in range(len(neis)):
                u = neis[j]
                if u < N and dist[u] == d - 1:
                    lab_u = labels.get(u)
                    if lab_u is not None:
                        w = wts[j]
                        acc[lab_u] = acc.get(lab_u, 0.0) + w
            if acc:
                # 稳定 tie-break：权重大，其次 taxid 小
                best_lab = max(acc.items(), key=lambda kv: (kv[1], -kv[0]))[0]
                labels[v] = best_lab

    return labels

