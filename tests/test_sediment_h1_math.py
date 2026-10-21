#!/usr/bin/env python3
"""Independent mathematical H1 donor-budget invariant checks; not CoLM tests."""
from math import expm1, isfinite

def h1_extract(exchange_stock, total_stock, gross_susp, gross_bed, dt):
    assert min(exchange_stock, total_stock, gross_susp, gross_bed, dt) >= 0.0
    assert exchange_stock <= total_stock
    demand = (gross_susp + gross_bed) * dt
    if total_stock == 0 or demand == 0:
        return (0.0, 0.0)
    taken = exchange_stock * (-expm1(-demand / total_stock))
    return (taken * (gross_susp * dt / demand),
            taken * (gross_bed * dt / demand))

def run():
    for se, sb, es, bed, dt in (
        (6.0, 60.0, 1e5, 1e3, 3600.0),
        (6.0, 60.0, 1.0, 2.0, 1e-4),
        (0.0, 60.0, 1.0, 2.0, 3600.0),
        (0.0, 0.0, 1.0, 2.0, 3600.0),
        (6.0, 60.0, 0.0, 0.0, 3600.0),
    ):
        susp, transport = h1_extract(se, sb, es, bed, dt)
        assert all(isfinite(x) for x in (susp, transport))
        assert susp >= 0 and transport >= 0
        assert susp + transport <= se * (1.0 + 1e-12)
        if es > 0 and bed > 0 and susp > 0:
            assert abs(susp / transport - es / bed) < 1e-10 * es / bed
    # Under fixed positive total stock, a single analytic step equals
    # repeated consumption of a *depleting* exchange stock, without refill.
    se, sb, rate, T = 6.0, 60.0, 3.0, 3600.0
    ref = se * (-expm1(-rate*T/sb))
    for n in (1, 2, 4, 8, 16):
        remaining = se
        for _ in range(n):
            remaining *= 1.0 + expm1(-rate*T/(n*sb))
        assert abs((se-remaining)-ref) < 1e-11
    print("PASS: independent H1 algebra, bounds, proportional allocation, exact integration")

if __name__ == "__main__":
    run()
