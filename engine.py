from database import build_default_registry
import numpy as np


def separator(title=""):
    if title:
        print(f"\n{'═'*60}")
        print(f"  {title}")
        print(f"{'═'*60}")
    else:
        print(f"{'─'*60}")


def main():
    registry = build_default_registry()

    # ----------------------------------------------------------------
    # OVERVIEW
    # ----------------------------------------------------------------
    separator("PHASE BEHAVIOR ONTOLOGY ENGINE")
    stats = registry.stats()
    print(f"  Phases:      {stats['total_phases']}")
    print(f"  Transitions: {stats['total_transitions']}")
    print(f"  Materials:   {stats['total_materials']}")
    print(f"  Categories:  {', '.join(sorted(stats['categories']))}")
    print(f"  Year range:  {stats['year_range'][0]} – {stats['year_range'][1]}")
    print(f"  Topological: {stats['topological']}")
    print(f"  Non-equil.:  {stats['non_equilibrium']}")

    # ----------------------------------------------------------------
    # ALL PHASES — summary table
    # ----------------------------------------------------------------
    separator("ALL PHASES")
    print(registry.all().summary())

    # ----------------------------------------------------------------
    # COMPOUND QUERIES
    # ----------------------------------------------------------------
    separator("COMPOUND QUERIES")

    print("\n  Topological + equilibrium:")
    print(f"  → {registry.all().topological().equilibrium().names()}")

    print("\n  Breaks symmetry + 3D:")
    print(f"  → {registry.all().breaks_symmetry().dimensionality(3).names()}")

    print("\n  Has edge states:")
    print(f"  → {registry.all().has_edge_states().names()}")

    print("\n  Has anyons:")
    print(f"  → {registry.all().has_anyons().names()}")

    print("\n  Berry phase = π:")
    print(f"  → {registry.all().berry_phase(np.pi).names()}")

    print("\n  Discovered after 1990:")
    print(f"  → {registry.all().discovered_after(1990).names()}")

    print("\n  Quantum computing applications:")
    print(f"  → {registry.all().has_application('quantum computing').names()}")

    print("\n  Topological + anyons + 2D:")
    print(f"  → {registry.all().topological().has_anyons().dimensionality(2).names()}")

    # ----------------------------------------------------------------
    # PHASE TRANSITION GRAPH
    # ----------------------------------------------------------------
    separator("PHASE TRANSITION GRAPH")
    for t in registry.transitions:
        order_str = "continuous" if t.continuous else "first-order"
        print(f"\n  {t.phase_from}  →  {t.phase_to}")
        print(f"    Control: {t.control_parameter} at {t.critical_value}")
        print(f"    Type: {order_str}, universality: {t.universality_class or 'unknown'}")
        print(f"    Mechanism: {t.mechanism}")

    print("\n  Neighbors of 'Mott Insulator':")
    print(f"  → {registry.neighbors('Mott Insulator')}")

    # ----------------------------------------------------------------
    # MATERIAL SURVEY
    # ----------------------------------------------------------------
    separator("MATERIALS SURVEY")
    print(f"  {'Formula':<20} {'Phase':<30} {'T range (K)':<20} {'Cost $/g'}")
    print(f"  {'─'*18} {'─'*28} {'─'*18} {'─'*10}")
    for phase in registry.phases:
        for mat in phase.materials:
            t_lo, t_hi = mat.temperature_range
            t_hi_str = f"{t_hi:.2e}" if t_hi < 1 else str(int(t_hi))
            cost = f"${mat.cost_per_gram}" if mat.cost_per_gram else "N/A"
            print(f"  {mat.formula:<20} {phase.name:<30} {t_lo}-{t_hi_str:<17} {cost}")

    # ----------------------------------------------------------------
    # OPEN QUESTIONS COUNT
    # ----------------------------------------------------------------
    separator("OPEN QUESTIONS BY PHASE")
    for phase in registry.phases:
        n = len(phase.open_questions)
        bar = "■" * n
        print(f"  {phase.name:<35} {bar} ({n})")

    # ----------------------------------------------------------------
    # EXPORT
    # ----------------------------------------------------------------
    separator("EXPORT")
    registry.to_json("/tmp/phases.json")
    print("  JSON saved to /tmp/phases.json")

    separator()


if __name__ == "__main__":
    main()
