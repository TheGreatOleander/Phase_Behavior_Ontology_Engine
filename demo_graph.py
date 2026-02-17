"""
demo_graph.py — Graph engine demo with real phase names
"""
from database import build_default_registry
from graph_engine import PhaseGraph

registry = build_default_registry()
graph = PhaseGraph(registry)

print(graph)

print("\nShortest paths:")
path_tests = [
    ("Ferromagnet", "Paramagnet"),
    ("Normal Metal", "Topological Superconductor"),
    ("Spin Glass", "BCS Superconductor"),
    ("Dirac Semimetal", "Weyl Semimetal"),
    ("Topological Insulator", "Floquet Topological Insulator"),
    ("Wigner Crystal", "Fractional Quantum Hall State"),
    ("Time Crystal", "Floquet Topological Insulator"),
    ("Quantum Spin Liquid", "Mott Insulator"),
]
for a, b in path_tests:
    try:
        path = graph.shortest_path(a, b)
        print(f"  {a}  →  {b}")
        print(f"    {' → '.join(path)}")
    except ValueError as e:
        print(f"  {a} → {b}: {e}")

print("\nCentrality (top 8):")
for node, score in sorted(graph.centrality().items(), key=lambda x: -x[1])[:8]:
    print(f"  {node:<40} {score:.3f}")

print("\nStrongly connected components:")
for comp in graph.strongly_connected_components():
    if len(comp) > 1:
        print(f"  {sorted(comp)}")

print("\nReachable from Normal Metal:")
print(f"  {sorted(graph.reachable_from('Normal Metal'))}")
