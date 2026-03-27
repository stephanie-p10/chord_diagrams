"""LaTeX/TikZ export functions for tiling visualizations."""

from typing import Dict, Tuple, List, Set


def export_tiling_to_latex(tiling, filename: str = "tiling_visual.tex", compile_pdf: bool = True) -> str:
    """Export a tiling as a TikZ picture to `filename`.
    
    If `compile_pdf` is True, attempt to run `pdflatex` (from PATH) to
    produce a PDF next to the .tex file. Returns path to .tex file.
    """
    dim_x, dim_y = tiling._cached_properties['dimensions']
    
    # Get chord_row_cells and determine which rows should show black chords
    chord_row_cells = tiling.chord_row_cells
    chord_rows = set(cell[1] for cell in chord_row_cells)
    
    # Group chord_row_cells by row to draw black chords
    chord_cells_by_row = {}
    for cell in chord_row_cells:
        row = cell[1]
        if row not in chord_cells_by_row:
            chord_cells_by_row[row] = []
        chord_cells_by_row[row].append(cell)
    
    # Filter out obstructions and requirements that are fully in chord_rows
    # (keep those that are only partially in chord_rows)
    filtered_obs = [ob for ob in tiling._obstructions 
                    if not all(cell[1] in chord_rows for cell in ob._pos)]
    filtered_reqs = [r for req_list in tiling._requirements 
                     for r in req_list 
                     if not all(cell[1] in chord_rows for cell in r._pos)]

    # Calculate scale factor based on number of obstructions and requirements
    # More adaptive: account for density of obstructions per cell
    total_constraints = len(filtered_obs) + len(filtered_reqs)
    
    # Calculate maximum obstructions and requirements in any single cell
    cell_ob_count = {}
    for ob in filtered_obs:
        for cell in set(ob._pos):
            cell_ob_count[cell] = cell_ob_count.get(cell, 0) + 1
    max_obs_per_cell = max(cell_ob_count.values()) if cell_ob_count else 1
    
    cell_req_count = {}
    for rq in filtered_reqs:
        for cell in set(rq._pos):
            cell_req_count[cell] = cell_req_count.get(cell, 0) + 1
    max_req_per_cell = max(cell_req_count.values()) if cell_req_count else 1
    
    # Maximum constraints per cell (obstructions + requirements)
    max_constraints_per_cell = max(max_obs_per_cell, max_req_per_cell)
    
    # Scale factor increases proportionally with max constraints per cell
    # Base scale gives minimum size, then add proportional scaling
    base_scale = 1.5
    constraint_scale = max(0, (total_constraints - 3) * 0.25)
    # Proportional scaling: each additional constraint per cell adds significant scale
    density_scale = max(0, (max_constraints_per_cell - 1) * 0.4)
    scale_factor = base_scale + constraint_scale + density_scale

    tex_lines = []
    tex_lines.append("\\documentclass{standalone}")
    tex_lines.append("\\usepackage{tikz}")
    tex_lines.append("\\begin{document}")
    tex_lines.append(f"\\begin{{tikzpicture}}[scale={scale_factor}]")

    # helper: cell-local coordinate
    def cell_coord(cell, ux, uy):
        return (cell[0] + ux, cell[1] + uy)

    ###########################################################################
    # 1. Build initial endpoint layouts for all obstructions/requirements
    ###########################################################################

    Shape = Dict[str, object]  # lightweight structure per obstruction/requirement
    shapes: List[Shape] = []

    # --- Obstructions (red, solid) ---
    for ob in filtered_obs:
        col = "red"

        # Vertical levels per chord within this obstruction (local, cell-relative)
        unique_chords = sorted(set(ob._chord_dict.keys()))
        chord_to_level: Dict[int, float] = {}
        num_chords = len(unique_chords)
        if num_chords == 0:
            continue
        if num_chords == 1:
            chord_to_level[unique_chords[0]] = 0.5
        else:
            # Spread chords in a wider vertical band around the cell midpoint for better separation
            band_center = 0.5
            band_half_height = 0.20  # results in [0.30, 0.70] - increased from 0.08
            for i, chord_id in enumerate(unique_chords):
                frac = i / max(1, num_chords - 1)
                uy = band_center - band_half_height + 2 * band_half_height * frac
                chord_to_level[chord_id] = uy

        # Endpoint positions (global coordinates) before packing
        pts: Dict[int, Tuple[float, float]] = {}
        num_endpoints = len(ob._pos)
        for endpoint_idx in range(num_endpoints):
            cell = ob._pos[endpoint_idx]

            # Vertical position from chord level
            uy = 0.5
            for chord_id, (i1, i2) in ob._chord_dict.items():
                if endpoint_idx == i1 or endpoint_idx == i2:
                    uy = chord_to_level[chord_id]
                    break

            # Horizontal position: left-to-right based on endpoint index
            # Increased range from 0.15-0.85 (was 0.2-0.8) for better spacing between endpoints
            ux = 0.15 + 0.70 * (endpoint_idx / max(1, num_endpoints - 1)) if num_endpoints > 1 else 0.5

            pts[endpoint_idx] = cell_coord(cell, ux, uy)

        # If this obstruction spans multiple rows, gently compress it vertically
        rows = [cell[1] for cell in ob._pos]
        min_row, max_row = min(rows), max(rows)
        if max_row - min_row == 1:
            center_y = (min_row + max_row + 1) / 2.0  # midpoint between cell centers
            condense_factor = 0.4
            for endpoint_idx, (x, y) in pts.items():
                y = center_y + (y - center_y) * condense_factor
                pts[endpoint_idx] = (x, y)

        # Build edges as in the original pretty-print logic
        edges: List[Tuple[int, int]] = []
        visited: Set[int] = set()
        for chord_id in ob._patt:
            if chord_id not in visited:
                visited.add(chord_id)
                i1, i2 = ob._chord_dict[chord_id]
                edges.append((i1, i2))

        # Compute initial bounding box
        xs = [x for (_, (x, _)) in pts.items()]
        ys = [y for (_, (_, y)) in pts.items()]
        if not xs or not ys:
            continue
        margin_x = 0.05
        margin_y = 0.05
        bbox = {
            "x_min": min(xs) - margin_x,
            "x_max": max(xs) + margin_x,
            "y_min": min(ys) - margin_y,
            "y_max": max(ys) + margin_y,
        }

        shapes.append(
            {
                "kind": "obstruction",
                "obj": ob,
                "color": col,
                "dashed": False,
                "node_radius": "0.03cm",
                "pts": pts,
                "edges": edges,
                "bbox": bbox,
                "y_offset": 0.0,
            }
        )

    # --- Requirements (blue, dashed) ---
    for rq in filtered_reqs:
        col = "blue"

        unique_chords = sorted(set(rq._chord_dict.keys()))
        chord_to_level: Dict[int, float] = {}
        num_chords = len(unique_chords)
        if num_chords == 0:
            continue
        if num_chords == 1:
            chord_to_level[unique_chords[0]] = 0.5
        else:
            band_center = 0.5
            band_half_height = 0.20  # results in [0.30, 0.70] - increased from 0.08 for better separation
            for i, chord_id in enumerate(unique_chords):
                frac = i / max(1, num_chords - 1)
                uy = band_center - band_half_height + 2 * band_half_height * frac
                chord_to_level[chord_id] = uy

        pts: Dict[int, Tuple[float, float]] = {}
        num_endpoints = len(rq._pos)
        for endpoint_idx in range(num_endpoints):
            cell = rq._pos[endpoint_idx]

            uy = 0.5
            for chord_id, (i1, i2) in rq._chord_dict.items():
                if endpoint_idx == i1 or endpoint_idx == i2:
                    uy = chord_to_level[chord_id]
                    break

            # Increased range from 0.15-0.85 (was 0.2-0.8) for better spacing between endpoints
            ux = 0.15 + 0.70 * (endpoint_idx / max(1, num_endpoints - 1)) if num_endpoints > 1 else 0.5

            pts[endpoint_idx] = cell_coord(cell, ux, uy)

        # If this requirement spans multiple rows, gently compress it vertically
        rows = [cell[1] for cell in rq._pos]
        min_row, max_row = min(rows), max(rows)
        if max_row - min_row == 1:
            center_y = (min_row + max_row + 1) / 2.0
            condense_factor = 0.4
            for endpoint_idx, (x, y) in pts.items():
                y = center_y + (y - center_y) * condense_factor
                pts[endpoint_idx] = (x, y)

        edges: List[Tuple[int, int]] = []
        visited: Set[int] = set()
        for chord_id in rq._patt:
            if chord_id not in visited:
                visited.add(chord_id)
                i1, i2 = rq._chord_dict[chord_id]
                edges.append((i1, i2))

        xs = [x for (_, (x, _)) in pts.items()]
        ys = [y for (_, (_, y)) in pts.items()]
        if not xs or not ys:
            continue
        margin_x = 0.05
        margin_y = 0.05
        bbox = {
            "x_min": min(xs) - margin_x,
            "x_max": max(xs) + margin_x,
            "y_min": min(ys) - margin_y,
            "y_max": max(ys) + margin_y,
        }

        shapes.append(
            {
                "kind": "requirement",
                "obj": rq,
                "color": col,
                "dashed": True,
                "node_radius": "0.025cm",
                "pts": pts,
                "edges": edges,
                "bbox": bbox,
                "y_offset": 0.0,
            }
        )

    ###########################################################################
    # 2. Iteratively separate rectangles vertically until none intersect
    ###########################################################################

    def rects_intersect(a, b, pad: float = 0.02) -> bool:
        """Return True if two rectangles (with padding) intersect."""
        ax1, ay1, ax2, ay2 = a["x_min"], a["y_min"], a["x_max"], a["y_max"]
        bx1, by1, bx2, by2 = b["x_min"], b["y_min"], b["x_max"], b["y_max"]

        ax1, ax2 = ax1 - pad, ax2 + pad
        ay1, ay2 = ay1 - pad, ay2 + pad
        bx1, bx2 = bx1 - pad, bx2 + pad
        by1, by2 = by1 - pad, by2 + pad

        return not (ax2 <= bx1 or bx2 <= ax1 or ay2 <= by1 or by2 <= ay1)

    max_iters = 80
    padding = 0.05  # vertical gap between rectangles (requirements vs obstructions too)

    for _ in range(max_iters):
        changed = False
        for i in range(len(shapes)):
            for j in range(i + 1, len(shapes)):
                si = shapes[i]
                sj = shapes[j]

                bi = si["bbox"]
                bj = sj["bbox"]

                # Effective bboxes with current y-offsets
                bi_eff = {
                    "x_min": bi["x_min"],
                    "x_max": bi["x_max"],
                    "y_min": bi["y_min"] + si["y_offset"],
                    "y_max": bi["y_max"] + si["y_offset"],
                }
                bj_eff = {
                    "x_min": bj["x_min"],
                    "x_max": bj["x_max"],
                    "y_min": bj["y_min"] + sj["y_offset"],
                    "y_max": bj["y_max"] + sj["y_offset"],
                }

                if not rects_intersect(bi_eff, bj_eff):
                    continue

                # Resolve intersection by moving the "later" shape upward
                # Compute minimal upward shift for sj so that it sits above si
                shift = (bi_eff["y_max"] - bj_eff["y_min"]) + padding
                if shift <= 0:
                    # If sj is actually above, move si up instead
                    shift = (bj_eff["y_max"] - bi_eff["y_min"]) + padding
                    si["y_offset"] += shift
                else:
                    sj["y_offset"] += shift

                changed = True
        if not changed:
            break

    ###########################################################################
    # 2b. Adjust vertical positions per row so points stay in their boxes
    ###########################################################################

    # Collect all endpoint y-positions per original row (after packing)
    # We also keep track of kind so we can separate obstructions/requirements vertically.
    row_points: Dict[int, List[Tuple[dict, int, float, str]]] = {}
    for shape in shapes:
        obj = shape["obj"]
        y_off = shape["y_offset"]
        pts = shape["pts"]
        for idx, (x, y) in pts.items():
            cell_row = obj._pos[idx][1]
            y_total = y + y_off
            row_points.setdefault(cell_row, []).append(
                (shape, idx, y_total, shape["kind"])
            )

    # For each row, linearly remap y so all endpoints for that row lie within it.
    # If both obstructions and requirements are present in the row, give them
    # slightly different vertical bands to reduce overlap.
    base_margin = 0.08
    for row, entries in row_points.items():
        has_obs = any(kind == "obstruction" for (_, _, _, kind) in entries)
        has_req = any(kind == "requirement" for (_, _, _, kind) in entries)

        def remap_band_grouped(local_entries, band_min, band_max):
            """
            Remap a list of (shape, idx, y_total) entries into [band_min, band_max],
            but keep different shapes in distinct vertical sub‑bands. This avoids
            different obstructions (or requirements) sitting exactly on top of
            each other in rows where several live in the same cells.
            """
            if not local_entries:
                return

            # Group endpoints by their parent shape object (use the obj field, not the dict itself)
            groups: Dict[object, List[Tuple[dict, int, float]]] = {}
            for s, i, y in local_entries:
                obj_key = s["obj"]  # Use the GriddedChord object as the key
                groups.setdefault(obj_key, []).append((s, i, y))

            num_groups = len(groups)
            if num_groups == 0:
                return

            # Add padding between groups for better visual separation
            # Use 30% of total height as padding, distributed between groups
            padding_per_group = 0.30 * (band_max - band_min) / max(1, num_groups - 1) if num_groups > 1 else 0.0
            total_padding = padding_per_group * (num_groups - 1)
            available_height = (band_max - band_min) - total_padding
            sub_height = available_height / num_groups

            # Add horizontal offset to separate different obstructions/requirements
            # when they share cells - spread them across the cell width
            horizontal_offset_range = 0.15  # Maximum horizontal offset per group
            for g_idx, (_, group_entries) in enumerate(groups.items()):
                # Add padding before each group (except the first)
                padding_offset = padding_per_group * g_idx
                g_band_min = band_min + g_idx * sub_height + padding_offset
                g_band_max = g_band_min + sub_height

                # Horizontal offset: distribute groups across the cell width
                # Center the offsets around 0 (so they range from -range/2 to +range/2)
                if num_groups > 1:
                    h_offset = -horizontal_offset_range/2 + (g_idx / (num_groups - 1)) * horizontal_offset_range
                else:
                    h_offset = 0.0

                ys_local = [y for (_, _, y) in group_entries]
                min_y = min(ys_local)
                max_y = max(ys_local)

                if max_y - min_y < 1e-6:
                    # Collapse this entire group to the center of its sub‑band
                    center = (g_band_min + g_band_max) / 2.0
                    for shape, idx, _ in group_entries:
                        x, _ = shape["pts"][idx]
                        shape["pts"][idx] = (x + h_offset, center)
                    continue

                scale = (g_band_max - g_band_min) / (max_y - min_y)
                for shape, idx, y_total in group_entries:
                    x, _ = shape["pts"][idx]
                    new_y = g_band_min + (y_total - min_y) * scale
                    shape["pts"][idx] = (x + h_offset, new_y)

        if has_obs and has_req:
            # Obstructions: lower band, requirements: upper band.
            # Within each band, different shapes get their own sub‑bands so they
            # don't visually overlap when they share cells/rows.
            obs_entries = [(s, i, y) for (s, i, y, k) in entries if k == "obstruction"]
            req_entries = [(s, i, y) for (s, i, y, k) in entries if k == "requirement"]

            # Choose two non-overlapping bands within the row:
            # obstructions strictly below requirements, with a clear vertical gap
            # between them so that their endpoints never overlap visually.
            #
            # We keep a slightly smaller band for obstructions, a smaller band
            # for requirements, and a larger gap in the middle.
            obs_min = row + base_margin
            obs_max = row + 0.40
            req_min = row + 0.60
            req_max = row + 1 - base_margin

            remap_band_grouped(obs_entries, obs_min, obs_max)
            remap_band_grouped(req_entries, req_min, req_max)
        else:
            # Only one kind present. If this row also has a black chord, keep
            # that chord centered and push the single-kind endpoints either
            # below or above the center band to avoid overlap. As above, keep
            # different shapes in the row in distinct vertical sub‑bands.
            ys = [y_total for (_, _, y_total, _) in entries]
            min_y = min(ys)
            max_y = max(ys)

            has_chord = row in chord_cells_by_row

            if has_chord and has_obs and not has_req:
                # Obstructions only and a black chord on this row:
                # put the obstruction entirely *above* the black chord so
                # diagonals to higher rows stay above the chord and do not
                # cross its center line.
                target_min = row + 0.55
                target_max = row + 1 - base_margin
            elif has_chord and has_req and not has_obs:
                # Requirements only: keep them in the upper half of the row.
                target_min = row + 0.55
                target_max = row + 1 - base_margin
            else:
                # No black chord on this row: use the full band.
                target_min = row + base_margin
                target_max = row + 1 - base_margin

            # Reuse the grouped remapping so multiple obstructions/requirements
            # in the same row do not overlap.
            flat_entries = [(s, i, y_total) for (s, i, y_total, _) in entries]
            remap_band_grouped(flat_entries, target_min, target_max)

    # After remapping, pts contain absolute y-coordinates; offsets can be zeroed
    for shape in shapes:
        shape["y_offset"] = 0.0

    # Draw the rounded cell boxes as exact squares (1x1 in TikZ units)
    tex_lines.append("  % draw rounded cell boxes")
    for x in range(dim_x):
        for y in range(dim_y):
            tex_lines.append(f"  \\draw[rounded corners=3pt] ({x},{y}) rectangle ({x+1},{y+1});")

    # Draw black chords for chord_row_cells *before* red/blue chords so that
    # we can treat their vertical position as reserved and keep other chords
    # away from them in the remapping logic above.
    tex_lines.append("  % draw black chords for chord_row_cells")
    for row, cells in chord_cells_by_row.items():
        # Keep black chord centered in the cell; obstructions/requirements
        # are packed away from this center band in the remapping above.
        chord_y = row + 0.5

        if len(cells) >= 2:
            # Sort cells by column for consistent left-to-right ordering
            sorted_cells = sorted(cells, key=lambda c: c[0])
            # Draw a black chord connecting the leftmost and rightmost cells in this row
            left_cell = sorted_cells[0]
            right_cell = sorted_cells[-1]

            # Place endpoints at the chosen height for this row
            left_x, left_y = left_cell[0] + 0.5, chord_y
            right_x, right_y = right_cell[0] + 0.5, chord_y

            # Create nodes for endpoints
            tex_lines.append(f"  \\node[circle, fill=black, inner sep=0.04cm] (chord_row{row}_left) at ({left_x},{left_y}) {{}};")
            tex_lines.append(f"  \\node[circle, fill=black, inner sep=0.04cm] (chord_row{row}_right) at ({right_x},{right_y}) {{}};")

            # Draw the chord
            tex_lines.append(f"  \\draw[black, line width=1.5pt] (chord_row{row}_left) -- (chord_row{row}_right);")
        elif len(cells) == 1:
            # Single cell: draw a chord within the cell
            cell = cells[0]
            left_x, cell_y = cell[0] + 0.3, chord_y
            right_x = cell[0] + 0.7

            tex_lines.append(f"  \\node[circle, fill=black, inner sep=0.04cm] (chord_row{row}_left) at ({left_x},{cell_y}) {{}};")
            tex_lines.append(f"  \\node[circle, fill=black, inner sep=0.04cm] (chord_row{row}_right) at ({right_x},{cell_y}) {{}};")
            tex_lines.append(f"  \\draw[black, line width=1.5pt] (chord_row{row}_left) -- (chord_row{row}_right);")

    ###########################################################################
    # 3. Emit TikZ nodes and edges using the packed layouts
    ###########################################################################

    for shape in shapes:
        kind = shape["kind"]
        obj = shape["obj"]
        color = shape["color"]
        dashed = shape["dashed"]
        node_radius = shape["node_radius"]
        pts = shape["pts"]
        edges = shape["edges"]
        y_off = shape["y_offset"]

        # Create nodes
        for i, (x, y) in pts.items():
            y_draw = y + y_off
            if kind == "obstruction":
                tex_lines.append(
                    f"  \\node[circle, fill={color}, inner sep={node_radius}] "
                    f"(obs{id(obj)}_pt{i}) at ({x},{y_draw}) {{}};"
                )
            else:
                tex_lines.append(
                    f"  \\node[circle, fill={color}, inner sep={node_radius}] "
                    f"(req{id(obj)}_pt{i}) at ({x},{y_draw}) {{}};"
                )

        # Draw edges as polylines
        if not edges:
            continue

        # Reconstruct the polyline order from the pattern (as before)
        visited: Set[int] = set()
        polyline_points: List[int] = []
        for chord_id in obj._patt:
            if chord_id not in visited:
                visited.add(chord_id)
                i1, i2 = obj._chord_dict[chord_id]
                polyline_points.append(i1)
                polyline_points.append(i2)

        if len(polyline_points) >= 2:
            if kind == "obstruction":
                path = " -- ".join(f"(obs{id(obj)}_pt{i})" for i in polyline_points)
                tex_lines.append(f"  \\draw[{color}, line width=1.2pt] {path};")
            else:
                path = " -- ".join(f"(req{id(obj)}_pt{i})" for i in polyline_points)
                tex_lines.append(f"  \\draw[{color}, dashed, line width=1.2pt] {path};")

    tex_lines.append("\\end{tikzpicture}")
    tex_lines.append("\\end{document}")

    with open(filename, "w") as f:
        f.write("\n".join(tex_lines))

    # attempt to compile
    if compile_pdf:
        import os
        import shutil
        import subprocess
        import tempfile
        from pathlib import Path

        pdflatex = shutil.which("pdflatex")
        if pdflatex:
            tex_path = Path(filename).resolve()
            tex_dir = tex_path.parent
            jobname = tex_path.stem
            try:
                # Compile into a temp directory so .aux/.log/etc. are not created
                # next to the requested .tex output.
                with tempfile.TemporaryDirectory(prefix="tiling_latex_") as tmpdir:
                    tmpdir_path = Path(tmpdir)
                    subprocess.run(
                        [
                            pdflatex,
                            "-interaction=nonstopmode",
                            "-halt-on-error",
                            f"-output-directory={str(tmpdir_path)}",
                            str(tex_path),
                        ],
                        cwd=str(tex_dir) or ".",
                        check=True,
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL,
                    )

                    built_pdf = tmpdir_path / f"{jobname}.pdf"
                    if built_pdf.exists():
                        out_pdf = tex_dir / f"{jobname}.pdf"
                        shutil.copy2(built_pdf, out_pdf)

                        # open the PDF on macOS
                        if os.name == "posix":
                            subprocess.run(["open", str(out_pdf)], check=False)
            except Exception:
                # Exporting the .tex file is still valuable even if compilation fails.
                pass
            finally:
                # Best-effort cleanup of common LaTeX sidecar files next to the .tex.
                # This also cleans up artifacts from older runs where we compiled in-place.
                sidecars = [
                    f"{jobname}.aux",
                    f"{jobname}.log",
                    f"{jobname}.out",
                    f"{jobname}.toc",
                    f"{jobname}.fls",
                    f"{jobname}.fdb_latexmk",
                    f"{jobname}.synctex.gz",
                ]
                for name in sidecars:
                    p = tex_dir / name
                    try:
                        if p.exists():
                            p.unlink()
                    except Exception:
                        pass

    return filename


