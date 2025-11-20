# idmap.py

class IdMap:
    """双向 ID 映射：字符串 <-> 连续整数"""
    __slots__ = ("to_int", "to_str")

    def __init__(self):
        self.to_int = {}
        self.to_str = []

    def get_int(self, s):
        """把字符串映射为整数 ID；第一次遇到会分配新 ID。"""
        m = self.to_int.get(s)
        if m is None:
            m = len(self.to_str)
            self.to_int[s] = m
            self.to_str.append(s)
        return m

    def get_str(self, idx):
        """根据整数 ID 取回原始字符串。"""
        return self.to_str[idx]

    def __len__(self):
        return len(self.to_str)

