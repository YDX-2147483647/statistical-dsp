from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from json import dumps
from pathlib import Path

from numpy import random


@dataclass
class Simulation:
    base: float
    real_value: float
    max_n: int

    gains: list[float]
    data: list[float]
    estimators: list[float]
    variances: list[float]


def simulate(
    *,
    base: float,
    real_value: float,
    max_n: int,
    rng: random.Generator,
) -> Simulation:
    """Simulate recursive least squares.

    The n-th noise variance is `base ** n`, where `n` = 0, 1, …, `max_n`.
    """
    gains = deque([0.0])  # gain[0] is unused and fixed to zero
    data = deque([rng.normal(real_value, 1.0)])
    estimators = deque(data)
    variances = deque([1.0])

    for n in range(1, max_n + 1):
        noise_variance = base**n
        data.append(rng.normal(real_value, noise_variance**0.5))

        gains.append(1 / (1 + noise_variance / variances[-1]))

        estimators.append(estimators[-1] + gains[-1] * (data[-1] - estimators[-1]))
        variances.append((1 - gains[-1]) * variances[-1])

    return Simulation(
        base=base,
        real_value=real_value,
        max_n=max_n,
        gains=list(gains),
        data=list(data),
        estimators=list(estimators),
        variances=list(variances),
    )


if __name__ == "__main__":
    rng = random.default_rng(seed=42)
    simulations = [
        simulate(base=0.95, real_value=10, max_n=100, rng=rng),
        simulate(base=1.00, real_value=10, max_n=100, rng=rng),
        simulate(base=1.05, real_value=10, max_n=100, rng=rng),
    ]
    Path(__file__).with_name("recursive-data.json").write_text(
        dumps([s.__dict__ for s in simulations]),
        encoding="utf-8",
    )
