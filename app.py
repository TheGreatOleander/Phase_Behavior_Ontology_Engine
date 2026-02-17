"""
app.py — Flask backend for the Phase Ontology Engine
"""

import os
import numpy as np
from flask import Flask, jsonify, request, abort
from flask_cors import CORS
from dotenv import load_dotenv

load_dotenv()

from database import build_default_registry
from discovery import (
    PhaseCoordinate, TopoType, SymBreaking,
    DynClass, AnyonType,
    KNOWN_PHASE_COORDINATES,
    THEORETICAL_PREDICTIONS,
    find_gaps
)
from materials import extract_requirements, suggest_materials

app = Flask(__name__)
CORS(app)

registry = build_default_registry()
_gaps_cache = find_gaps(KNOWN_PHASE_COORDINATES, THEORETICAL_PREDICTIONS)


# ============================================================
# Helpers
# ============================================================

def phase_to_dict(phase):
    return {
        "name": phase.name,
        "category": phase.category,
        "dimensionality": phase.dimensionality,
        "topological": phase.topology.has_topological_order,
        "breaks_symmetry": phase.symmetry_breaking is not None,
        "equilibrium": phase.dynamics.equilibrium,
        "prototype_material": phase.prototype_material,
        "applications": phase.technological_applications,
    }


def coord_from_dict(d):
    return PhaseCoordinate(
        dimensionality=int(d.get("dimensionality", 3)),
        symmetry_breaking=SymBreaking(d.get("symmetry_breaking", "none")),
        topological=TopoType(d.get("topological", "none")),
        dynamics=DynClass(d.get("dynamics", "equilibrium")),
        anyons=AnyonType(d.get("anyons", "none")),
        edge_states=bool(d.get("edge_states", False)),
        bulk_gap=bool(d.get("bulk_gap", False)),
        long_range_entanglement=bool(d.get("long_range_entanglement", False)),
    )


# ============================================================
# API Routes
# ============================================================

@app.route("/api/phases")
def get_phases():
    q = registry.all()

    if request.args.get("topological") == "true":
        q = q.topological()

    if request.args.get("equilibrium") == "true":
        q = q.equilibrium()

    if request.args.get("breaks_symmetry") == "true":
        q = q.breaks_symmetry()

    if d := request.args.get("dimensionality"):
        q = q.dimensionality(int(d))

    phases = list(q)

    return jsonify({
        "count": len(phases),
        "phases": [phase_to_dict(p) for p in phases]
    })


@app.route("/api/phases/<name>")
def get_phase(name):
    results = registry.query(name=name)
    phase = results.first()

    if not phase:
        abort(404, f"Phase '{name}' not found")

    return jsonify(phase_to_dict(phase))


@app.route("/api/stats")
def get_stats():
    stats = registry.stats()
    stats["categories"] = sorted(stats["categories"])
    return jsonify(stats)


@app.route("/api/transitions")
def get_transitions():
    return jsonify({
        "count": len(registry.transitions),
        "transitions": [
            {
                "from": t.phase_from,
                "to": t.phase_to,
                "order": t.order,
                "continuous": t.continuous,
                "control_parameter": t.control_parameter,
                "critical_value": (
                    None if t.critical_value == float("inf")
                    else t.critical_value
                ),
                "universality_class": t.universality_class,
            }
            for t in registry.transitions
        ]
    })


@app.route("/api/gaps")
def get_gaps():
    return jsonify({
        "total": len(_gaps_cache),
        "gaps": [
            {
                "coordinate": str(g.coordinate),
                "status": g.status,
                "interest": round(g.interest_score, 3),
            }
            for g in _gaps_cache
        ]
    })


@app.route("/api/materials", methods=["POST"])
def get_materials():
    data = request.get_json()

    if not data:
        abort(400, "JSON body required")

    coord = coord_from_dict(data)
    reqs = extract_requirements(coord)
    candidates = suggest_materials(coord, top_n=5)

    return jsonify({
        "coordinate": str(coord),
        "requirements": [r.name for r in reqs],
        "candidates": [
            {
                "name": c.name,
                "score": round(c.overall_score, 3),
                "rationale": c.rationale,
            }
            for c in candidates
        ]
    })


# ============================================================
# Boot
# ============================================================

if __name__ == "__main__":
    port = int(os.getenv("PORT", 8000))
    debug = os.getenv("FLASK_DEBUG", "true").lower() == "true"
    app.run(host="0.0.0.0", port=port, debug=debug)
