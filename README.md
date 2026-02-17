# Phase Behavior Ontology Engine

A computational "periodic table of behaviors" for phases of matter.

Phases are encoded by symmetry properties, topological invariants, quantum
characteristics, experimental signatures, material realizations, and how they
connect to each other via phase transitions.

## Install

```bash
git clone <repo>
cd phase-ontology

python -m venv venv && source venv/bin/activate   # optional but recommended
pip install -r requirements.txt
cp _env.example .env                               # add MP_API_KEY if you have one
```

## Docker

```bash
docker build -t phase-ontology .
docker run -p 5000:5000 --env-file .env phase-ontology
```

Or with Compose:

```bash
docker-compose up
```

## Files

| File | Purpose |
|------|---------|
| `phase.py` | Core schema — all dataclasses, PhaseQuery, PhaseRegistry |
| `database.py` | All phase entries and the transition graph |
| `discovery.py` | Phase space coordinate system, gap finder, theorem checker |
| `materials.py` | Material suggestion engine |
| `graph_engine.py` | NetworkX phase transition graph with path-finding |
| `engine.py` | Query, explore, and export the registry |
| `demo.py` | Deep walkthrough of one phase (Bi₂Se₃ Topological Insulator) |
| `demo_graph.py` | Graph engine demo — shortest paths, centrality |
| `app.py` | Flask REST API |
| `pipeline.py` | Materials Project API integration (requires MP_API_KEY) |

## Run

```bash
python engine.py        # Full query demo + export
python demo.py          # Deep dive into one phase
python demo_graph.py    # Graph traversal demo
python discovery.py     # Gap finder report
python materials.py     # Material suggestion report
flask run               # Start API server (port 5000)
```

## Tests

```bash
pytest test_*.py -v
```

## Phases in Registry (21)

| Phase | Category | Dim | Year |
|-------|----------|-----|------|
| Crystal | Symmetry-Broken | 3D | 1912 |
| Ferromagnet | Symmetry-Broken | 3D | 1895 |
| BCS Superconductor | Symmetry-Broken | 3D | 1911 |
| Bose-Einstein Condensate | Symmetry-Broken | 3D | 1995 |
| Nematic Liquid Crystal | Symmetry-Broken | 3D | 1888 |
| Paramagnet | Symmetry-Broken | 3D | 1895 |
| Normal Metal | Symmetry-Broken | 3D | 1900 |
| Normal Liquid | Symmetry-Broken | 3D | 1687 |
| Normal Bose Gas | Symmetry-Broken | 3D | 1924 |
| Mott Insulator | Strongly-Correlated | 3D | 1937 |
| Spin Glass | Strongly-Correlated | 3D | 1972 |
| Wigner Crystal | Strongly-Correlated | 2D | 1934 |
| Topological Insulator | Topological | 3D | 2007 |
| Weyl Semimetal | Topological | 3D | 2015 |
| Dirac Semimetal | Topological | 3D | 2014 |
| Topological Superconductor | Topological | 1D | 2012 |
| Integer Quantum Hall State | Topological | 2D | 1980 |
| Fractional Quantum Hall State | Topological | 2D | 1982 |
| Quantum Spin Liquid | Topological | 2D | 1973 |
| Floquet Topological Insulator | Topological | 2D | 2013 |
| Time Crystal | Non-Equilibrium | 1D | 2016 |

## Phase Transitions (17)

- Crystal → Normal Liquid (melting, first-order)
- Ferromagnet ↔ Paramagnet (Curie point, 3D Heisenberg)
- BCS Superconductor → Normal Metal (T_c, mean-field)
- BCS Superconductor ↔ Topological Superconductor (Zeeman + SOC, class D)
- BCS Superconductor ↔ Mott Insulator (doping-driven, cuprates)
- BEC → Normal Bose Gas (λ-point, 3D XY)
- Normal Metal ↔ Topological Insulator (SOC band inversion)
- Normal Metal ↔ Spin Glass (frustration + disorder)
- Normal Metal → Topological Superconductor (proximity + field)
- Normal Liquid ↔ Nematic Liquid Crystal (isotropic-nematic, weakly first-order)
- Integer QHE ↔ Fractional QHE (interactions, plateau transition)
- Dirac Semimetal ↔ Topological Insulator (gap opening)
- Dirac Semimetal ↔ Weyl Semimetal (TRS/inversion breaking)
- Wigner Crystal ↔ FQHE (filling factor competition)
- Topological Insulator ↔ Floquet TI (drive amplitude)
- Quantum Spin Liquid → Mott Insulator (pressure/frustration)
- Time Crystal → Floquet TI (disorder tuning)

## Query System

Queries are chainable:

```python
from database import build_default_registry
import numpy as np

registry = build_default_registry()

registry.all().topological().equilibrium().names()
registry.all().has_anyons().names()
registry.all().berry_phase(np.pi).names()
registry.all().discovered_after(1990).names()
registry.all().has_application("quantum computing").names()
registry.all().topological().has_anyons().dimensionality(2).names()
registry.all().dimensionality(1).names()
registry.all().non_equilibrium().names()

registry.neighbors("Mott Insulator")
registry.stats()
registry.to_json("phases.json")
```

## Graph Engine

```python
from graph_engine import PhaseGraph

graph = PhaseGraph(registry)
graph.shortest_path("Spin Glass", "BCS Superconductor")
# → ['Spin Glass', 'Normal Metal', 'BCS Superconductor']

graph.reachable_from("Normal Metal")
graph.centrality()
graph.betweenness()
```

## Discovery Engine

```python
from discovery import find_gaps, KNOWN_PHASE_COORDINATES, THEORETICAL_PREDICTIONS

gaps = find_gaps(KNOWN_PHASE_COORDINATES, THEORETICAL_PREDICTIONS)
# Scans ~5700 coordinates in phase space
# Returns: known / forbidden (by theorem) / predicted / candidate gaps
```

## Extend

Add new phases in `database.py` and a coordinate entry in
`discovery.py:KNOWN_PHASE_COORDINATES`. Copy any existing entry as a template.

Good next candidates: Anderson Insulator, Kitaev Honeycomb, ν=5/2 Moore-Read,
Chern Insulator (zero-field), Haldane Phase, Polar Metal.
