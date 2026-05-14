"""
core/id_map.py – Bidirectional string↔integer ID compressor.

Keeps insertion order; IDs are assigned sequentially starting from 0.
"""

from __future__ import annotations
from typing import List


class IdMap:
    """
    Compress arbitrary string identifiers to compact integers.

    Attributes:
        to_int: Maps string → integer ID.
        to_str: Maps integer ID → original string (index = ID).
    """

    __slots__ = ("to_int", "to_str")

    def __init__(self) -> None:
        self.to_int: dict[str, int] = {}
        self.to_str: List[str] = []

    def get_int(self, s: str) -> int:
        """Return existing ID for *s*, or register and return a new one."""
        existing = self.to_int.get(s)
        if existing is not None:
            return existing
        new_id = len(self.to_str)
        self.to_int[s] = new_id
        self.to_str.append(s)
        return new_id

    def lookup(self, i: int) -> str:
        """Return string for integer ID *i*. Raises IndexError if out of range."""
        return self.to_str[i]

    def __len__(self) -> int:
        return len(self.to_str)

    def __contains__(self, s: str) -> bool:
        return s in self.to_int

    def __repr__(self) -> str:
        return f"IdMap(size={len(self)})"
