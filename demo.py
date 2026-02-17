#!/usr/bin/env python3
"""
demo.py — Detailed phase walkthrough
Shows the full richness of one phase entry: Bi₂Se₃ Topological Insulator
"""

from database import build_default_registry
import numpy as np


def print_section(title):
    print(f"\n{'─'*60}")
    print(f"  {title}")
    print(f"{'─'*60}")


def main():
    registry = build_default_registry()

    # Pull the Topological Insulator entry
    ti = registry.query(name="Topological Insulator")[0]

    print("=" * 60)
    print(f"  PHASE DETAIL: {ti.name}")
    print(f"  Category: {ti.category}  |  Dimensionality: {ti.dimensionality}D")
    print(f"  Discovered: {ti.discovery_year} by {ti.discovered_by}")
    print("=" * 60)

    # --- Symmetry ---
    print_section("SYMMETRY")
    if ti.symmetry_breaking is None:
        print("  Breaks NO symmetry — purely topological phase")
    else:
        sb = ti.symmetry_breaking
        print(f"  Order parameter: {sb.order_parameter.name} ({sb.order_parameter.symbol})")
        print(f"  Universality class: {sb.universality_class}")

    # --- Topology ---
    print_section("TOPOLOGY")
    topo = ti.topology
    print(f"  Has topological order: {topo.has_topological_order} ({topo.topological_order_type})")
    print(f"  Bulk gap: {topo.bulk_gap}")
    print(f"  Edge states: {topo.edge_states}")
    print(f"  Dispersion: {topo.edge_state_dispersion}")
    print(f"  Field theory: {topo.topological_field_theory}")
    print(f"\n  Topological Invariants:")
    for inv in topo.invariants_3d:
        print(f"    • {inv.name} ({inv.symbol}) = {inv.value}")
        print(f"      Meaning: {inv.physical_meaning}")
        print(f"      Math: {inv.k_theory}")
        print(f"      Consequences:")
        for c in inv.observable_consequence:
            print(f"        - {c}")

    # --- Quantum ---
    print_section("QUANTUM PROPERTIES")
    q = ti.quantum
    print(f"  Berry phase: {q.berry_phase/np.pi:.2f}π")
    print(f"  Berry curvature: {q.berry_curvature}")
    if q.coherence:
        print(f"  Decoherence time: {q.coherence.decoherence_time*1e12:.1f} ps")
        print(f"  Decoherence mechanisms:")
        for m in q.coherence.decoherence_mechanism:
            print(f"    - {m}")
    if q.entanglement:
        print(f"  Entanglement: {q.entanglement.entanglement_entropy_formula}")
        print(f"  Long-range entanglement: {q.entanglement.long_range_entanglement}")

    # --- Experimental ---
    print_section("EXPERIMENTAL SIGNATURES")
    exp = ti.experimental
    print("  Primary signatures:")
    for sig in exp.primary_signatures:
        print(f"    • {sig}")
    print("\n  Spectroscopy:")
    for spec in exp.spectroscopy:
        print(f"    [{spec.technique}] ({spec.energy_range[0]} to {spec.energy_range[1]} eV)")
        for feat in spec.characteristic_features:
            print(f"      - {feat}")
    if exp.transport:
        print(f"\n  Transport:")
        print(f"    Conductivity: {exp.transport.conductivity}")
        print(f"    Hall: {exp.transport.hall_coefficient}")
    print(f"\n  Required techniques:")
    for tech in exp.required_techniques:
        print(f"    • {tech.name}")
        print(f"      Measures: {tech.measures}")
        print(f"      Setup: {tech.typical_setup}")
    print(f"\n  Sample purity required: {exp.sample_purity_required}")
    print(f"  Temperature control: {exp.temperature_control[0]}-{exp.temperature_control[1]} K")

    # --- Materials ---
    print_section("MATERIALS")
    print(f"  Prototype: {ti.prototype_material}")
    for mat in ti.materials:
        print(f"\n  • {mat.formula}")
        print(f"    Structure: {mat.crystal_structure} ({mat.space_group})")
        print(f"    Temperature range: {mat.temperature_range[0]}-{mat.temperature_range[1]} K")
        if mat.band_gap:
            print(f"    Band gap: {mat.band_gap} eV")
        if mat.lattice_constants:
            print(f"    Lattice constants: {mat.lattice_constants} Å")
        print(f"    Synthesis: {', '.join(mat.synthesis_method)}")
        if mat.cost_per_gram:
            print(f"    Cost: ${mat.cost_per_gram}/g")
        if mat.discovery_paper:
            print(f"    Discovery: {mat.discovery_paper}")

    # --- References ---
    print_section("KEY REFERENCES")
    print("  Theory:")
    for p in ti.key_theoretical_papers:
        print(f"    • {p}")
    print("  Experiment:")
    for p in ti.key_experimental_papers:
        print(f"    • {p}")
    print("  Reviews:")
    for p in ti.review_papers:
        print(f"    • {p}")
    print("  Textbooks:")
    for p in ti.textbook_references:
        print(f"    • {p}")

    # --- Applications ---
    print_section("APPLICATIONS")
    print("  Current:")
    for app in ti.technological_applications:
        print(f"    • {app}")
    print("  Potential:")
    for app in ti.potential_applications:
        print(f"    • {app}")

    # --- Open questions ---
    print_section("OPEN QUESTIONS")
    for q in ti.open_questions:
        print(f"  ? {q}")
    print("\n  Controversies:")
    for c in ti.controversies:
        print(f"  ⚠ {c}")

    print("\n" + "=" * 60)


if __name__ == "__main__":
    main()
