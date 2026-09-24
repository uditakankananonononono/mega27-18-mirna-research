"""Duplex CNN: separate conv towers for UTR window and mature miRNA, plus an
explicit pairwise complementarity map (outer product of one-hots through the
Watson-Crick/wobble kernel) processed by a 2D CNN."""
from __future__ import annotations
import torch
import torch.nn as nn

# pairing energy proxy: A-U, U-A, G-C, C-G = 1 ; G-U, U-G wobble = 0.5
_PAIR = torch.tensor([[0, 0, 0, 1],
                      [0, 0, 1, 0],
                      [0, 1, 0, .5],
                      [1, 0, .5, 0]], dtype=torch.float32)


def pair_map(u: torch.Tensor, m: torch.Tensor) -> torch.Tensor:
    """u: (B,4,L_u) m: (B,4,L_m) -> (B,1,L_m,L_u) complementarity."""
    P = _PAIR.to(u.device)
    return torch.einsum("bim,ij,bju->bmu", m, P, u).unsqueeze(1)


class DuplexCNN(nn.Module):
    def __init__(self, n_site_types: int = 4, width: int = 32, n_feat: int = 0):
        super().__init__()
        self.utr = nn.Sequential(nn.Conv1d(4, width, 7, padding=3), nn.ReLU(),
                                 nn.Conv1d(width, width, 5, padding=2), nn.ReLU(),
                                 nn.AdaptiveMaxPool1d(1))
        self.mir = nn.Sequential(nn.Conv1d(4, width, 5, padding=2), nn.ReLU(),
                                 nn.AdaptiveMaxPool1d(1))
        self.pair = nn.Sequential(nn.Conv2d(1, width, (7, 7), padding=3), nn.ReLU(),
                                  nn.MaxPool2d(2),
                                  nn.Conv2d(width, width, 3, padding=1), nn.ReLU(),
                                  nn.AdaptiveMaxPool2d(1))
        self.st = nn.Embedding(n_site_types, 8)
        self.head = nn.Sequential(nn.Linear(3 * width + 8 + n_feat, 64), nn.ReLU(), nn.Linear(64, 1))

    def forward(self, u, m, st, f=None):
        parts = [self.utr(u).flatten(1), self.mir(m).flatten(1),
                       self.pair(pair_map(u, m)).flatten(1), self.st(st)]
        if f is not None:
            parts.append(f)
        h = torch.cat(parts, 1)
        return self.head(h).squeeze(1)
