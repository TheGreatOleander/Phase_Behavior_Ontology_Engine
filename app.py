"""
app.py — Flask backend for the Phase Ontology Engine

Endpoints:
  GET  /                        → Web UI
  GET  /api/phases              → All phases (with optional filters)
  GET  /api/phases/<name>       → Single phase detail
  GET  /api/stats               → Registry statistics
  GET  /api/gaps                → Discovery engine results
  POST /api/materials           → Material suggestions for a coordinate
  GET  /api/transitions         → Phase transition graph
  POST /api/pipeline            → Live MP API query (if key set)
"""

import os
import json
import numpy as np
from flask import Flask, jsonify, request, render_template, send_from_directory, abort
from flask_cors import CORS
from dotenv import load_dotenv

load_dotenv()

# Import our engine
from database import build_default_registry
from discovery import (
    PhaseCoordinate, TopoType, SymBreaking, DynClass, AnyonType,
    KNOWN_PHASE_COORDINATES, THEORETICAL_PREDICTIONS,
    find_gaps
)
from materials import extract_requirements, suggest_materials, MATERIAL_FAMILIES

app = Flask(__name__)
CORS(app)  # Allow cross-origin requests from any frontend origin
registry = build_default_registry()
# Pre-compute gap analysis once at startup (expensive combinatorial scan)
_gaps_cache = find_gaps(KNOWN_PHASE_COORDINATES, THEORETICAL_PREDICTIONS)


# ============================================================================
# HELPERS
# ============================================================================

def phase_to_dict(phase):
    """Serialize a Phase object to JSON-safe dict."""
    sb = phase.symmetry_breaking
    return {
        "name": phase.name,
        "category": phase.category,
        "dimensionality": phase.dimensionality,
        "discovery_year": phase.discovery_year,
        "theoretical_prediction_year": phase.theoretical_prediction_year,
        "experimental_confirmation_year": phase.experimental_confirmation_year,
        "discovered_by": phase.discovered_by,
        "topology": {
            "has_topological_order": phase.topology.has_topological_order,
            "topological_order_type": phase.topology.topological_order_type,
            "edge_states": phase.topology.edge_states,
            "bulk_gap": phase.topology.bulk_gap,
            "anyonic_excitations": phase.topology.anyonic_excitations,
            "anyon_types": phase.topology.anyon_types,
            "long_range_entanglement": phase.topology.long_range_entanglement,
            "topological_entanglement_entropy": phase.topology.topological_entanglement_entropy,
            "topological_field_theory": phase.topology.topological_field_theory,
            "edge_state_dispersion": phase.topology.edge_state_dispersion,
            "invariants": [
                {
                    "name": inv.name,
                    "symbol": inv.symbol,
                    "value": inv.value,
                    "meaning": inv.physical_meaning,
                    "k_theory": inv.k_theory,
                    "consequences": inv.observable_consequence,
                }
                for inv in (
                    phase.topology.invariants_0d +
                    phase.topology.invariants_1d +
                    phase.topology.invariants_2d +
                    phase.topology.invariants_3d
                )
            ],
        },
        "symmetry_breaking": {
            "order_parameter": sb.order_parameter.name,
            "symbol": sb.order_parameter.symbol,
            "physical_meaning": sb.order_parameter.physical_meaning,
            "units": sb.order_parameter.units,
            "universality_class": sb.universality_class,
            "critical_temperature": sb.critical_temperature,
            "num_goldstone_modes": sb.num_goldstone_modes,
            "goldstone_dispersion": sb.goldstone_dispersion,
            "topological_defects": sb.topological_defects,
            "measurement_techniques": sb.order_parameter.measurement_technique,
        } if sb else None,
        "quantum": {
            "berry_phase": phase.quantum.berry_phase,
            "berry_phase_over_pi": (
                round(phase.quantum.berry_phase / np.pi, 2)
                if phase.quantum.berry_phase is not None else None
            ),
            "berry_curvature": phase.quantum.berry_curvature,
            "coherence": {
                "decoherence_time": phase.quantum.coherence.decoherence_time,
                "mechanisms": phase.quantum.coherence.decoherence_mechanism,
            } if phase.quantum.coherence else None,
            "entanglement": {
                "area_law": phase.quantum.entanglement.area_law,
                "long_range": phase.quantum.entanglement.long_range_entanglement,
                "formula": phase.quantum.entanglement.entanglement_entropy_formula,
                "tensor_network": phase.quantum.entanglement.tensor_network_structure,
            } if phase.quantum.entanglement else None,
        },
        "dynamics": {
            "equilibrium": phase.dynamics.equilibrium,
            "driven": phase.dynamics.driven is not None,
            "floquet": phase.dynamics.driven.floquet_system if phase.dynamics.driven else False,
            "mbl": phase.dynamics.time_evolution.many_body_localized,
            "ergodic": phase.dynamics.time_evolution.ergodic,
        },
        "experimental": {
            "primary_signatures": phase.experimental.primary_signatures,
            "sample_purity": phase.experimental.sample_purity_required,
            "temperature_control": list(phase.experimental.temperature_control)
                if phase.experimental.temperature_control else None,
            "spectroscopy": [
                {
                    "technique": s.technique,
                    "energy_range": list(s.energy_range),
                    "features": s.characteristic_features,
                }
                for s in phase.experimental.spectroscopy
            ],
            "transport": {
                "conductivity": phase.experimental.transport.conductivity,
                "resistivity": phase.experimental.transport.resistivity,
                "hall": phase.experimental.transport.hall_coefficient,
            } if phase.experimental.transport else None,
            "required_techniques": [
                {
                    "name": t.name,
                    "measures": t.measures,
                    "setup": t.typical_setup,
                }
                for t in phase.experimental.required_techniques
            ],
        },
        "materials": [
            {
                "formula": m.formula,
                "crystal_structure": m.crystal_structure,
                "space_group": m.space_group,
                "temperature_range": list(m.temperature_range),
                "band_gap_eV": m.band_gap,
                "synthesis": m.synthesis_method,
                "cost_per_gram_usd": m.cost_per_gram,
                "commercially_available": m.commercially_available,
                "discovery_paper": m.discovery_paper,
            }
            for m in phase.materials
        ],
        "prototype_material": phase.prototype_material,
        "neighboring_phases": phase.neighboring_phases,
        "key_theoretical_papers": phase.key_theoretical_papers,
        "key_experimental_papers": phase.key_experimental_papers,
        "review_papers": phase.review_papers,
        "textbook_references": phase.textbook_references,
        "technological_applications": phase.technological_applications,
        "potential_applications": phase.potential_applications,
        "open_questions": phase.open_questions,
        "controversies": phase.controversies,
    }


def coord_from_dict(d):
    """Parse a PhaseCoordinate from a JSON dict."""
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


# ============================================================================
# ROUTES — Frontend
# ============================================================================

@app.route("/")
def index():
    # render_template looks in templates/; fall back to project root for
    # deployments that serve index.html directly from the project directory.
    try:
        return render_template("index.html")
    except Exception:
        import os
        return send_from_directory(
            os.path.dirname(os.path.abspath(__file__)), "index.html"
        )


# ============================================================================
# ROUTES — API
# ============================================================================

@app.route("/api/phases")
def get_phases():
    """Return all phases, with optional query filters."""
    q = registry.all()

    # Apply filters from query params
    if request.args.get("topological") == "true":
        q = q.topological()
    elif request.args.get("topological") == "false":
        q = q.non_topological()

    if request.args.get("equilibrium") == "true":
        q = q.equilibrium()
    elif request.args.get("equilibrium") == "false":
        q = q.non_equilibrium()

    if request.args.get("breaks_symmetry") == "true":
        q = q.breaks_symmetry()

    if request.args.get("has_anyons") == "true":
        q = q.has_anyons()

    if request.args.get("has_edge_states") == "true":
        q = q.has_edge_states()

    if d := request.args.get("dimensionality"):
        q = q.dimensionality(int(d))

    if cat := request.args.get("category"):
        q = q.category(cat)

    if year := request.args.get("discovered_after"):
        q = q.discovered_after(int(year))

    if app_kw := request.args.get("application"):
        q = q.has_application(app_kw)

    if request.args.get("berry_phase_pi") == "true":
        q = q.berry_phase(np.pi)

    phases = list(q)
    return jsonify({
        "count": len(phases),
        "phases": [phase_to_dict(p) for p in phases]
    })


@app.route("/api/phases/<name>")
def get_phase(name):
    """Return full detail for a single phase by name."""
    results = registry.query(name=name)
    phase = results.first()
    if not phase:
        abort(404, f"Phase '{name}' not found")
    return jsonify(phase_to_dict(phase))


@app.route("/api/stats")
def get_stats():
    """Registry statistics."""
    stats = registry.stats()
    # Make JSON serializable
    stats["categories"] = sorted(stats["categories"])
    return jsonify(stats)


@app.route("/api/transitions")
def get_transitions():
    """Phase transition graph."""
    return jsonify({
        "count": len(registry.transitions),
        "transitions": [
            {
                "from": t.phase_from,
                "to": t.phase_to,
                "order": t.order,
                "continuous": t.continuous,
                "control_parameter": t.control_parameter,
                "critical_value": t.critical_value if t.critical_value != float('inf') else None,
                "universality_class": t.universality_class,
                "mechanism": t.mechanism,
                "critical_exponents": t.critical_exponents,
                "experimental_signatures": t.experimental_signatures,
            }
            for t in registry.transitions
        ]
    })


@app.route("/api/gaps")
def get_gaps():
    """Run discovery engine and return gap analysis."""
    all_gaps = _gaps_cache

    known     = [g for g in all_gaps if g.status == "known"]
    forbidden = [g for g in all_gaps if g.status == "forbidden"]
    predicted = [g for g in all_gaps if g.status == "predicted"]
    candidate = [g for g in all_gaps if g.status == "candidate"]

    # Top candidates sorted by interest
    top_candidates = sorted(candidate, key=lambda x: x.interest_score, reverse=True)[:20]
    top_predicted  = sorted(predicted, key=lambda x: x.prediction.confidence, reverse=True)

    # Forbidden zone summary
    theorem_counts = {}
    for g in forbidden:
        name = g.forbidden_by.name
        theorem_counts[name] = theorem_counts.get(name, 0) + 1

    return jsonify({
        "summary": {
            "total_evaluated": len(all_gaps),
            "known": len(known),
            "forbidden": len(forbidden),
            "predicted": len(predicted),
            "candidate_gaps": len(candidate),
        },
        "predictions": [
            {
                "name": g.prediction.name,
                "coordinate": str(g.coordinate),
                "predicted_by": g.prediction.predicted_by,
                "year": g.prediction.prediction_year,
                "confidence": g.prediction.confidence,
                "mechanism": g.prediction.mechanism,
                "status": g.prediction.experimental_status,
                "candidate_materials": g.prediction.candidate_materials,
                "paper": g.prediction.prediction_paper,
                "novelty": g.novelty_score,
                "feasibility": g.feasibility_score,
                "interest": g.interest_score,
            }
            for g in top_predicted
        ],
        "top_candidates": [
            {
                "coordinate": str(g.coordinate),
                "novelty": round(g.novelty_score, 3),
                "feasibility": round(g.feasibility_score, 3),
                "interest": round(g.interest_score, 3),
                "reasoning": g.reasoning,
            }
            for g in top_candidates
        ],
        "forbidden_zones": [
            {
                "theorem": name,
                "count": count,
                "statement": next(
                    g.forbidden_by.statement for g in forbidden
                    if g.forbidden_by.name == name
                ),
                "reference": next(
                    g.forbidden_by.reference for g in forbidden
                    if g.forbidden_by.name == name
                ),
            }
            for name, count in sorted(theorem_counts.items(), key=lambda x: -x[1])
        ]
    })


@app.route("/api/materials", methods=["POST"])
def get_materials():
    """Suggest materials for a given phase coordinate."""
    data = request.get_json()
    if not data:
        abort(400, "JSON body required with phase coordinate fields")

    try:
        coord = coord_from_dict(data)
    except (ValueError, KeyError) as e:
        abort(400, f"Invalid coordinate: {e}")

    reqs = extract_requirements(coord)
    candidates = suggest_materials(coord, top_n=8)

    return jsonify({
        "coordinate": str(coord),
        "requirements": [
            {
                "name": r.name,
                "description": r.description,
                "difficulty": r.difficulty,
                "implies_elements": r.implies_elements,
                "implies_conditions": r.implies_conditions,
                "implies_techniques": r.implies_techniques,
            }
            for r in reqs
        ],
        "candidates": [
            {
                "name": c.name,
                "examples": c.formula_or_class,
                "rationale": c.rationale,
                "requirements_satisfied": c.requirements_satisfied,
                "requirements_missing": c.requirements_missing,
                "feasibility": round(c.feasibility_score, 3),
                "overall_score": round(c.overall_score, 3),
                "synthesis_difficulty": round(c.synthesis_difficulty, 3),
                "detection_difficulty": round(c.detection_difficulty, 3),
                "experiment": c.suggested_experiment,
                "timeline": c.estimated_timeline,
                "warning": c.warning,
            }
            for c in candidates
        ]
    })


@app.route("/api/neighbors/<name>")
def get_neighbors(name):
    """Return neighboring phases in the transition graph."""
    neighbors = registry.neighbors(name)
    transitions = []
    for neighbor in neighbors:
        t = registry.get_transition(name, neighbor)
        if t:
            transitions.append({
                "neighbor": neighbor,
                "control_parameter": t.control_parameter,
                "critical_value": t.critical_value if t.critical_value != float('inf') else None,
                "mechanism": t.mechanism,
                "continuous": t.continuous,
            })
    return jsonify({"phase": name, "neighbors": transitions})


@app.route("/api/pipeline", methods=["POST"])
def run_pipeline():
    """
    Live Materials Project API query.
    Requires MP_API_KEY in .env or passed in request body.
    """
    data = request.get_json() or {}
    api_key = data.get("api_key") or os.getenv("MP_API_KEY")

    if not api_key:
        return jsonify({
            "error": "No API key provided",
            "hint": "Set MP_API_KEY in .env or pass api_key in request body",
            "get_key": "https://materialsproject.org/api"
        }), 400

    try:
        coord = coord_from_dict(data.get("coordinate", {}))
    except Exception as e:
        return jsonify({"error": f"Invalid coordinate: {e}"}), 400

    try:
        from pipeline import query_materials_project, score_material, enrich_with_topology
        reqs = extract_requirements(coord)
        materials = query_materials_project(reqs, coord)
        for mat in materials:
            score_material(mat, coord)
        materials.sort(key=lambda x: x.coordinate_match_score, reverse=True)

        return jsonify({
            "coordinate": str(coord),
            "total_found": len(materials),
            "materials": [
                {
                    "mp_id": m.mp_id,
                    "formula": m.formula_pretty,
                    "crystal_system": m.crystal_system,
                    "space_group": m.space_group,
                    "band_gap_eV": m.band_gap,
                    "is_stable": m.is_stable,
                    "energy_above_hull": m.energy_above_hull,
                    "elements": m.elements,
                    "topological_class": m.topological_class,
                    "z2_invariant": m.z2_invariant,
                    "score": round(m.coordinate_match_score, 3),
                    "reasoning": m.reasoning,
                }
                for m in materials[:20]
            ]
        })

    except ImportError:
        return jsonify({
            "error": "mp-api not installed",
            "fix": "pip install mp-api pymatgen"
        }), 500
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/enum_options")
def enum_options():
    """Return all valid enum values for the coordinate builder UI."""
    return jsonify({
        "dimensionality": [1, 2, 3],
        "symmetry_breaking": [e.value for e in SymBreaking],
        "topological": [e.value for e in TopoType],
        "dynamics": [e.value for e in DynClass],
        "anyons": [e.value for e in AnyonType],
    })


if __name__ == "__main__":
    port = int(os.getenv("PORT", 5000))
    debug = os.getenv("FLASK_DEBUG", "false").lower() == "true"
    app.run(host="0.0.0.0", port=port, debug=debug)
