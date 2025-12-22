"""LaTeX/TikZ export functions for tiling visualizations."""

from typing import Dict, Tuple


def export_tiling_to_latex(tiling, filename: str = "tiling_visual.tex", compile_pdf: bool = True) -> str:
    """Export a tiling as a TikZ picture to `filename`.
    
    If `compile_pdf` is True, attempt to run `pdflatex` (from PATH) to
    produce a PDF next to the .tex file. Returns path to .tex file.
    """
    dim_x, dim_y = tiling._cached_properties['dimensions']

    # label assignment consistent with pretty_print
    obs_labels = {ob: chr(ord('A') + i % 26) + (str(i//26) if i//26 else '') 
                  for i, ob in enumerate(tiling._obstructions)}
    all_reqs = [r for req_list in tiling._requirements for r in req_list]
    req_labels = {rq: chr(ord('a') + i % 26) + (str(i//26) if i//26 else '') 
                  for i, rq in enumerate(all_reqs)}

    # Calculate scale factor based on number of obstructions and requirements
    total_constraints = len(tiling._obstructions) + len(all_reqs)
    scale_factor = 1.0 + max(0, (total_constraints - 3) * 0.15)

    tex_lines = []
    tex_lines.append("\\documentclass{standalone}")
    tex_lines.append("\\usepackage{tikz}")
    tex_lines.append("\\begin{document}")
    tex_lines.append(f"\\begin{{tikzpicture}}[scale={scale_factor}]")

    # draw rounded cell boxes
    tex_lines.append("  % draw rounded cell boxes")
    for x in range(dim_x):
        for y in range(dim_y):
            tex_lines.append(f"  \\draw[rounded corners=3pt] ({x},{y}) rectangle ({x+1},{y+1});")

    # helper: cell-local coordinate
    def cell_coord(cell, ux, uy):
        return (cell[0] + ux, cell[1] + uy)

    # draw obstructions: place endpoints inside each cell and draw connections
    for idx, ob in enumerate(tiling._obstructions):
        col = "red"
        lab = obs_labels[ob]
        
        # Determine vertical levels for each unique chord ID (bottom to top: 0 below 1 below 2, etc)
        unique_chords = sorted(set(ob._chord_dict.keys()))
        chord_to_level = {}
        
        # Distribute obstruction base levels across the vertical space
        num_obs = len(tiling._obstructions)
        if num_obs == 1:
            obs_base_uy = 0.5
        else:
            # Spread obstructions evenly from 0.1 to 0.9
            obs_base_uy = 0.1 + 0.8 * (idx / max(1, num_obs - 1))
        
        if len(unique_chords) == 1:
            # Single chord: use the obstruction's base level
            chord_to_level[unique_chords[0]] = obs_base_uy
        else:
            # Multiple chords: distribute around the obstruction's base level
            num_chords = len(unique_chords)
            for i, chord_id in enumerate(unique_chords):
                # Spread chords evenly: chord 0 at lowest, chord N-1 at highest
                uy = obs_base_uy - 0.1 + (i / max(1, num_chords - 1)) * 0.2 if num_chords > 1 else obs_base_uy
                chord_to_level[chord_id] = max(0.05, min(0.95, uy))
        
        # Place endpoints left-to-right in the order they appear in the pattern
        # Each endpoint's horizontal position is based on its position in the overall pattern
        pts: Dict[int, Tuple[float, float]] = {}
        num_endpoints = len(ob._pos)
        
        for endpoint_idx in range(num_endpoints):
            cell = ob._pos[endpoint_idx]
            
            # Find which chord uses this endpoint to determine vertical level
            uy = 0.5  # default
            for chord_id, (i1, i2) in ob._chord_dict.items():
                if endpoint_idx == i1 or endpoint_idx == i2:
                    uy = chord_to_level[chord_id]
                    break
            
            # Horizontal position: left-to-right based on endpoint index
            ux = 0.2 + 0.6 * (endpoint_idx / max(1, num_endpoints - 1)) if num_endpoints > 1 else 0.5
            
            pts[endpoint_idx] = cell_coord(cell, ux, uy)
        
        # create named nodes for all endpoints
        for i, (x, y) in pts.items():
            tex_lines.append(f"  \\node[circle, fill={col}, inner sep=0.03cm] (obs{id(ob)}_pt{i}) at ({x},{y}) {{}};")

        # draw polyline connecting all endpoints using node references
        visited = set()
        polyline_points = []
        for idx, chord_id in enumerate(ob._patt):
            if chord_id not in visited:
                visited.add(chord_id)
                i1, i2 = ob._chord_dict[chord_id]
                polyline_points.append(i1)
                polyline_points.append(i2)
        if len(polyline_points) >= 2:
            path = " -- ".join(f"(obs{id(ob)}_pt{i})" for i in polyline_points)
            tex_lines.append(f"  \\draw[{col}, line width=1.2pt] {path};")

    # draw requirements similarly but in blue dashed style
    for idx, rq in enumerate(all_reqs):
        col = "blue"
        lab = req_labels[rq]
        
        # Determine vertical levels for each unique chord ID (bottom to top: 0 below 1 below 2, etc)
        unique_chords = sorted(set(rq._chord_dict.keys()))
        chord_to_level = {}
        
        # Distribute requirement base levels across the vertical space
        num_reqs = len(all_reqs)
        if num_reqs == 1:
            req_base_uy = 0.5
        else:
            # Spread requirements evenly from 0.1 to 0.9
            req_base_uy = 0.1 + 0.8 * (idx / max(1, num_reqs - 1))
        
        if len(unique_chords) == 1:
            # Single chord: use the requirement's base level
            chord_to_level[unique_chords[0]] = req_base_uy
        else:
            # Multiple chords: distribute around the requirement's base level
            num_chords = len(unique_chords)
            for i, chord_id in enumerate(unique_chords):
                # Spread chords evenly: chord 0 at lowest, chord N-1 at highest
                uy = req_base_uy - 0.1 + (i / max(1, num_chords - 1)) * 0.2 if num_chords > 1 else req_base_uy
                chord_to_level[chord_id] = max(0.05, min(0.95, uy))
        
        # Place endpoints left-to-right in the order they appear in the pattern
        # Each endpoint's horizontal position is based on its position in the overall pattern
        pts: Dict[int, Tuple[float, float]] = {}
        num_endpoints = len(rq._pos)
        
        for endpoint_idx in range(num_endpoints):
            cell = rq._pos[endpoint_idx]
            
            # Find which chord uses this endpoint to determine vertical level
            uy = 0.5  # default
            for chord_id, (i1, i2) in rq._chord_dict.items():
                if endpoint_idx == i1 or endpoint_idx == i2:
                    uy = chord_to_level[chord_id]
                    break
            
            # Horizontal position: left-to-right based on endpoint index
            ux = 0.2 + 0.6 * (endpoint_idx / max(1, num_endpoints - 1)) if num_endpoints > 1 else 0.5
            
            pts[endpoint_idx] = cell_coord(cell, ux, uy)
        
        # create named nodes for all endpoints
        for i, (x, y) in pts.items():
            tex_lines.append(f"  \\node[circle, fill={col}, inner sep=0.025cm] (req{id(rq)}_pt{i}) at ({x},{y}) {{}};")

        # draw polyline connecting all endpoints using node references
        visited = set()
        polyline_points = []
        for idx_pat, chord_id in enumerate(rq._patt):
            if chord_id not in visited:
                visited.add(chord_id)
                i1, i2 = rq._chord_dict[chord_id]
                polyline_points.append(i1)
                polyline_points.append(i2)
        if len(polyline_points) >= 2:
            path = " -- ".join(f"(req{id(rq)}_pt{i})" for i in polyline_points)
            tex_lines.append(f"  \\draw[{col}, dashed, line width=1.2pt] {path};")

    tex_lines.append("\\end{tikzpicture}")
    tex_lines.append("\\end{document}")

    with open(filename, "w") as f:
        f.write("\n".join(tex_lines))

    # attempt to compile
    if compile_pdf:
        import shutil, subprocess, os
        pdflatex = shutil.which("pdflatex")
        if pdflatex:
            try:
                subprocess.run([pdflatex, "-interaction=nonstopmode", filename], 
                             cwd=os.path.dirname(os.path.abspath(filename)) or ".", 
                             check=True, stdout=subprocess.DEVNULL)
                # open the PDF on macOS
                pdf_file = filename.replace(".tex", ".pdf")
                if os.path.exists(pdf_file):
                    subprocess.run(["open", pdf_file], check=False)
            except Exception:
                pass

    return filename
