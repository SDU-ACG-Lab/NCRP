# refine.py

from collections import deque


def refine_labels_linear(graph, init_labels_int):
    """
    无权多源 BFS，找“最近相遇”的对手标签集合；
    若集合包含多个或不含自身标签，则删为歧义。
    """
    seeds = list(init_labels_int.keys())
    if not seeds:
        return dict(init_labels_int)

    N = max(graph.keys(), default=0) + 1
    owner = [-1] * N
    dist = [-1] * N

    INF = 10**15
    bestD = {s: INF for s in seeds}
    bestLabSet = {s: set() for s in seeds}

    q = deque()
    for s in seeds:
        owner[s] = s
        dist[s] = 0
        q.append(s)

    while q:
        u = q.popleft()
        own_u = owner[u]
        du = dist[u]
        neis, _ = graph.get(u, ((), ()))
        for j in range(len(neis)):
            v = neis[j]
            if v >= N:
                continue
            if owner[v] == -1:
                owner[v] = own_u
                dist[v] = du + 1
                q.append(v)
            else:
                own_v = owner[v]
                if own_v != own_u:
                    cand = du + dist[v] + 1
                    for s, t in ((own_u, own_v), (own_v, own_u)):
                        if cand < bestD[s]:
                            bestD[s] = cand
                            bestLabSet[s] = {init_labels_int[t]}
                        elif cand == bestD[s]:
                            bestLabSet[s].add(init_labels_int[t])

    refined = dict(init_labels_int)
    removed = 0
    for s, self_lab in init_labels_int.items():
        labs = bestLabSet[s]
        if labs and (len(labs) > 1 or self_lab not in labs):
            refined.pop(s, None)
            removed += 1

    print(f"Refined labels (linear): removed {removed} ambiguous nodes.")
    return refined

