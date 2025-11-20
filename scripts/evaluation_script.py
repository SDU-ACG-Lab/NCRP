import argparse, csv
from collections import defaultdict
from taxonomy import Taxonomy

# ---------- helpers ----------
def get_node_safe(tx, taxid):
    try:
        return tx.node(str(taxid))
    except Exception:
        return None

def lift_to_rank(tx, taxid, rank):
    """
    把 taxid 提升/上卷到指定 rank。
    返回：该 rank 的 taxid（字符串），或 None（不存在/未知）。
    """
    n = get_node_safe(tx, taxid)
    if n is None:
        return None
    if n.rank == rank:
        return str(n.id)
    try:
        p = tx.parent(str(taxid), at_rank=rank)
        return str(p.id) if p else None
    except Exception:
        return None

def is_ancestor(tx, ancestor_taxid, descendant_taxid):
    """ancestor_taxid 是否为 descendant_taxid 的祖先（含自身）。"""
    cur = get_node_safe(tx, descendant_taxid)
    anc = str(ancestor_taxid)
    while cur:
        if str(cur.id) == anc:
            return True
        cur = tx.parent(cur.id)
    return False

# ---------- main eval ----------
def evaluate(kraken_file, read_mapping_file, taxonomy_dir, eval_rank="genus"):
    tx = Taxonomy.from_ncbi(dump_dir=taxonomy_dir)

    # 1) ground truth -> eval_rank
    read2truth = {}
    lookup_fail_truth = 0
    with open(read_mapping_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rid = row["#anonymous_read_id"]
            tid = row["tax_id"]
            # 未知/未识别
            n = get_node_safe(tx, tid)
            if (n is None) or ("unidentified" in (n.name or "").lower()):
                read2truth[rid] = "0"     # 约定 0 表示“无标签”
                continue
            # 升到评估层级
            t_rank = lift_to_rank(tx, tid, eval_rank)
            if t_rank is None:
                # 该 reads 的标签在 eval_rank 以上（上卷不到）
                lookup_fail_truth += 1
                continue
            read2truth[rid] = t_rank

    total_reads = len(read2truth)

    # 2) parse kraken output
    tp = fp = fn = vp = 0
    lookup_fail_pred = 0
    with open(kraken_file) as f:
        for line in f:
            if not line.strip(): 
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            status, rid, pred_taxid = parts[0], parts[1], parts[2]
            if rid not in read2truth:
                continue
            truth = read2truth[rid]

            if status == "U":
                # 预测未分类
                if truth == "0":
                    tp += 1           # 真也未标（当作正确未检出，保持与原脚本口径）
                else:
                    fn += 1
                continue

            if status != "C":
                continue

            # 预测已分类：先把预测 taxid 归一到评估层级
            pred_rank = lift_to_rank(tx, pred_taxid, eval_rank)
            if pred_rank is None:
                # 预测的 taxid 不在谱系里 或 在评估层级以上：算一次查找失败并记 FP
                lookup_fail_pred += 1
                if truth != "0":
                    fp += 1
                continue

            if truth == "0":
                # 真值无标签，但模型给了分类 -> 记 FP
                fp += 1
            elif pred_rank == truth:
                tp += 1
            elif is_ancestor(tx, pred_rank, truth):
                vp += 1
            else:
                fp += 1

    # 3) metrics
    denom = tp + vp + fp + fn
    sen = tp / denom if denom else 0.0
    pre = tp / (tp + fp) if (tp + fp) else 0.0
    f1  = 2 * sen * pre / (sen + pre) if (sen + pre) else 0.0

    print(f"\nTotal Reads in Truth:   {total_reads}")
    print(f"Evaluation Rank:         {eval_rank}")
    print("Kraken2 Evaluation:")
    print(f"TP: {tp}, VP: {vp}, FP: {fp}, FN: {fn} "
          f"(Taxonomy lookup failed: truth={lookup_fail_truth}, pred={lookup_fail_pred})")
    print(f"Sensitivity (SEN):       {sen:.4f}")
    print(f"Precision (PRE):         {pre:.4f}")
    print(f"F1-score:                {f1:.4f}")

def main():
    ap = argparse.ArgumentParser("Evaluate Kraken2 classification")
    ap.add_argument("--reads", required=True, help="reads_mapping.tsv (ground truth)")
    ap.add_argument("--kraken", required=True, help="Kraken2 output (default format)")
    ap.add_argument("--taxonomy", required=True, help="NCBI taxonomy dump dir")
    ap.add_argument("--rank", default="genus",
                    choices=["species","genus","family","order","class","phylum"])
    args = ap.parse_args()
    evaluate(args.kraken, args.reads, args.taxonomy, eval_rank=args.rank)

if __name__ == "__main__":
    main()
