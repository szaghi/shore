"""Spacing1D: a single edge's spacing law and parameters."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from shore.exceptions import ParameterError
from shore.volume.spacing import (
    build_k_steps,
    tanh_one_sided,
    tanh_two_sided,
    uniform,
)


@dataclass
class Spacing1D:
    """One-dimensional spacing law for a block edge.

    Parameters
    ----------
    law :
        Spacing algorithm: ``'uniform'``, ``'tanh'``, ``'tanh2'``,
        or ``'geometric'``.
    beta :
        Clustering strength for tanh/tanh2 laws.
    ds :
        First-layer thickness for geometric law.
    growth :
        Geometric growth ratio (geometric law only).
    total_thickness :
        Total extent required for tanh/tanh2 laws.
    cluster_at :
        Asymmetric clustering target for tanh2 (longitude axis use case).
    """

    law: Literal["uniform", "tanh", "tanh2", "geometric"]
    beta: float = 3.0
    ds: float | None = None
    growth: float = 1.15
    total_thickness: float | None = None
    cluster_at: float | None = None

    def build(self, n: int, start: float = 0.0, stop: float = 1.0) -> np.ndarray:
        """Return n coordinates on [start, stop] according to the spacing law.

        For ``'geometric'`` this returns step sizes (n-1,), not coordinates.
        For all other laws it returns coordinates (n,).

        Parameters
        ----------
        n :
            Number of points (or layers for geometric).
        start, stop :
            Interval endpoints.

        Returns
        -------
        numpy.ndarray
            Shape ``(n,)`` for uniform/tanh/tanh2; ``(n-1,)`` for geometric.
        """
        if self.law == "uniform":
            return uniform(n, start, stop)
        if self.law == "tanh":
            return tanh_one_sided(n, start, stop, self.beta)
        if self.law == "tanh2":
            return tanh_two_sided(n, start, stop, self.beta, self.cluster_at)
        if self.law == "geometric":
            if self.ds is None:
                raise ParameterError("Spacing1D with law='geometric' requires ds.")
            return build_k_steps(n, self.ds, self.growth)
        raise ParameterError(f"Unknown spacing law '{self.law}'.")

    def build_steps(self, nk: int) -> np.ndarray:
        """Return (nk-1,) wall-normal step sizes.

        Convenience wrapper for the k-direction; delegates to
        ``volume.spacing.build_k_steps`` for all laws.

        Parameters
        ----------
        nk :
            Total number of layers including k=0.
        """
        if self.law == "geometric":
            if self.ds is None:
                raise ParameterError("Spacing1D with law='geometric' requires ds.")
            return build_k_steps(nk, self.ds, self.growth, spacing="geometric")
        if self.total_thickness is None:
            raise ParameterError(
                f"Spacing1D with law='{self.law}' requires total_thickness."
            )
        return build_k_steps(
            nk,
            ds=self.ds or 0.0,
            growth=self.growth,
            spacing=self.law,
            beta=self.beta,
            total_thickness=self.total_thickness,
        )
