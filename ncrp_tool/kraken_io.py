# kraken_io.py

from idmap import IdMap


def load_labels_and_order(kraken_file, rid_map: IdMap, tax_map: IdMap):
    """
    读取 Kraken2 per-read 输出。

    返回：
      labels: {rid_int: taxid_int}
      order:  [rid_int]  （按 Kraken 出现顺序，含 U/C 全部）
    """
    labels = {}
    order = []
    seen = set()

    with open(kraken_file, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            status, rid_s, tax_s = parts[0], parts[1], parts[2]

            # 无论 U/C 都注册，确保可输出
            rid = rid_map.get_int(rid_s)
            if rid not in seen:
                order.append(rid)
                seen.add(rid)

            if status == "C" and tax_s != "0":
                labels[rid] = tax_map.get_int(tax_s)

    return labels, order

