"""3D molecular visualization helpers based on py3Dmol."""

from __future__ import annotations

from pathlib import Path
from tempfile import NamedTemporaryFile

import mdtraj as md
import py3Dmol


def _validate_pdb_path(pdb_path: Path) -> Path:
    """Validate the public PDB-path boundary."""

    path = Path(pdb_path)
    if not path.exists():
        raise FileNotFoundError(f"pdb file not found: {path}")
    if path.suffix.lower() != ".pdb":
        raise ValueError("pdb_path must point to a .pdb file")
    return path


def _validate_frame_index(frame_index: int, n_frames: int) -> int:
    """Validate a trajectory frame index."""

    index = int(frame_index)
    if index < 0 or index >= n_frames:
        raise ValueError("frame_index must be within [0, N_frames)")
    return index


def _trajectory_frame_to_pdb_string(frame: md.Trajectory) -> str:
    """Serialize a single MDTraj frame to PDB text for py3Dmol."""

    with NamedTemporaryFile(suffix=".pdb", delete=True) as handle:
        frame.save_pdb(handle.name)
        return Path(handle.name).read_text(encoding="utf-8")


def _extract_chain_ids(pdb_text: str) -> list[str]:
    """Extract unique chain identifiers from PDB ATOM/HETATM records in order of appearance."""

    seen: dict[str, None] = {}
    for line in pdb_text.splitlines():
        record = line[0:6].strip()
        if record in {"ATOM", "HETATM"} and len(line) > 21:
            chain_id = line[21]
            if chain_id.strip() and chain_id not in seen:
                seen[chain_id] = None
    return list(seen)


def _apply_base_chain_style(view: py3Dmol.view, style: str, pdb_text: str = "") -> None:
    """Apply the blueprint chain styling to the complex widget."""

    if style not in {"cartoon", "stick", "sphere", "line"}:
        raise ValueError("style must be one of: cartoon, stick, sphere, line")

    chains = _extract_chain_ids(pdb_text) if pdb_text else []
    chain_colors = ["blue", "red", "green", "orange", "cyan", "magenta"]

    if len(chains) >= 2:
        for i, chain_id in enumerate(chains):
            color = chain_colors[i] if i < len(chain_colors) else "gray"
            view.setStyle({"chain": chain_id}, {style: {"color": color}})
    else:
        view.setStyle({"chain": "A"}, {style: {"color": "blue"}})
        view.setStyle({"chain": "B"}, {style: {"color": "red"}})
        view.setStyle({"not": {"chain": ["A", "B"]}}, {style: {"color": "gray"}})


def _highlight_interface_residues(view: py3Dmol.view, residue_ids: list[int] | None) -> None:
    """Highlight requested interface residues as green spheres."""

    if residue_ids is None or len(residue_ids) == 0:
        return
    residue_list = [int(residue_id) for residue_id in residue_ids]
    view.addStyle({"resi": residue_list}, {"sphere": {"color": "green", "radius": 0.9}})


def _highlight_active_site_triad(view: py3Dmol.view) -> None:
    """Highlight His/Asp/Ser residues as the active-site triad representation when present."""

    view.addStyle({"resn": ["HIS", "ASP", "SER"]}, {"stick": {"colorscheme": "yellowCarbon", "radius": 0.18}})
    view.addStyle({"resn": ["HIS", "ASP", "SER"]}, {"sphere": {"colorscheme": "yellowCarbon", "radius": 0.35}})


def render_complex(
    pdb_path: Path,
    highlight_interface_residues: list[int] | None = None,
    style: str = "cartoon",
) -> py3Dmol.view:
    """Render a prepared SPINK7-KLK5 complex as an interactive py3Dmol widget."""

    validated_path = _validate_pdb_path(pdb_path)
    pdb_text = validated_path.read_text(encoding="utf-8")

    view = py3Dmol.view(width=800, height=600)
    view.addModel(pdb_text, "pdb")
    _apply_base_chain_style(view, style, pdb_text)
    _highlight_interface_residues(view, highlight_interface_residues)
    _highlight_active_site_triad(view)
    view.zoomTo()
    return view


def render_trajectory_frame(
    trajectory: md.Trajectory,
    frame_index: int,
    color_by: str = "chain",
) -> py3Dmol.view:
    """Render a single MDTraj frame as an interactive py3Dmol widget."""

    if trajectory.n_frames < 1:
        raise ValueError("trajectory must contain at least one frame")

    index = _validate_frame_index(frame_index, trajectory.n_frames)
    frame = trajectory[index]
    pdb_text = _trajectory_frame_to_pdb_string(frame)

    view = py3Dmol.view(width=800, height=600)
    view.addModel(pdb_text, "pdb")

    if color_by == "chain":
        _apply_base_chain_style(view, "cartoon")
    elif color_by == "bfactor":
        view.setStyle({}, {"cartoon": {"colorscheme": {"prop": "b", "gradient": "roygb", "min": 0, "max": 100}}})
    elif color_by == "rmsf":
        view.setStyle({}, {"cartoon": {"colorscheme": "spectrum"}})
    else:
        raise ValueError("color_by must be one of: chain, bfactor, rmsf")

    view.zoomTo()
    return view