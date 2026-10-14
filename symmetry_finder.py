#!/usr/bin/env python3
"""
Point group symmetry finder for molecules from XYZ files.
Uses only numpy as external dependency.
"""
import numpy as np
import sys

# Tolerance for position matching (Angstrom)
TOL = 1e-3


def read_xyz(path):
    """Read XYZ file: returns (atom_symbols: list[str], coords: np.ndarray Nx3)"""
    with open(path, 'r') as f:
        lines = f.read().strip().splitlines()
    n_atoms = int(lines[0].strip())
    # Skip comment line (line 1)
    data_lines = lines[2:2 + n_atoms]
    symbols = []
    coords = []
    for line in data_lines:
        parts = line.split()
        symbols.append(parts[0])
        coords.append([float(x) for x in parts[1:4]])
    return symbols, np.array(coords)


def center_geom(coords):
    return np.mean(coords, axis=0)


def rotation_matrix(axis, angle_deg):
    """Return 3x3 rotation matrix for rotation by angle_deg around axis (unit vector)."""
    theta = np.radians(angle_deg)
    ux, uy, uz = axis / np.linalg.norm(axis)
    cos_t, sin_t = np.cos(theta), np.sin(theta)
    # Rodrigues' rotation formula as matrix
    R = np.array([
        [cos_t + ux*ux*(1-cos_t), ux*uy*(1-cos_t) - uz*sin_t, ux*uz*(1-cos_t) + uy*sin_t],
        [uy*ux*(1-cos_t) + uz*sin_t, cos_t + uy*uy*(1-cos_t), uy*uz*(1-cos_t) - ux*sin_t],
        [uz*ux*(1-cos_t) - uy*sin_t, uz*uy*(1-cos_t) + ux*sin_t, cos_t + uz*uz*(1-cos_t)]
    ])
    return R


def all_atom_types_match(symbols):
    # Helper: check if all symbols identical (simplifies matching)
    return len(set(symbols)) == 1


def find_symmetry_elements(symbols, coords):
    center = center_geom(coords)
    centered = coords - center
    n = len(centered)

    # Detect inversion center: every atom i has matching atom j of same type at -r
    has_inversion = True
    for i in range(n):
        matched = False
        for j in range(n):
            if i == j:
                # Atom at center required for self-match
                if np.linalg.norm(centered[i]) > TOL:
                    matched = False
                    break
                else:
                    matched = True
                    break
            if symbols[i] != symbols[j]:
                continue
            if np.linalg.norm(centered[i] + centered[j]) < TOL:
                matched = True
                break
        if not matched:
            has_inversion = False
            break

    # Test rotation axes. Try common angles around x, y, z and some diagonals.
    # For simplicity, test principal axes first.
    axes_to_test = {
        'C2_z': (np.array([0,0,1]), 180),
        'C2_y': (np.array([0,1,0]), 180),
        'C2_x': (np.array([1,0,0]), 180),
        'C3_z': (np.array([0,0,1]), 120),
        'C4_z': (np.array([0,0,1]), 90),
        'C5_z': (np.array([0,0,1]), 72),
        'C6_z': (np.array([0,0,1]), 60),
        'C3_y': (np.array([0,1,0]), 120),
        'C4_y': (np.array([0,1,0]), 90),
        'S2': (np.array([0,0,1]), 180),  # improper = rotation + reflection through plane perpendicular
    }
    # Note: full S_n detection requires reflection component, handled separately

    rotation_axes = []
    for name, (axis, angle) in axes_to_test.items():
        R = rotation_matrix(axis, angle)
        rotated = (R @ centered.T).T
        # Try to match rotated positions to original (with same symbol mapping)
        matched_all = True
        # For each original atom, find unique match in rotated
        used = [False] * n
        for i in range(n):
            found = False
            for j in range(n):
                if used[j] or symbols[i] != symbols[j]:
                    continue
                if np.linalg.norm(rotated[i] - centered[j]) < TOL:
                    used[j] = True
                    found = True
                    break
            if not found:
                matched_all = False
                break
        if matched_all:
            # Extract axis direction and order (angle/360)
            order = int(round(360 / angle)) if angle != 0 else 0
            rotation_axes.append((name, axis, order))

    # Deduplicate by axis direction/order
    unique_axes = {}
    for name, axis, order in rotation_axes:
        # Normalize axis direction
        axis_n = tuple(np.round(axis / np.linalg.norm(axis), decimals=3))
        key = (order, axis_n)
        if key not in unique_axes:
            unique_axes[key] = name

    # Detect mirror planes: reflection through plane with normal n (unit) passing through center.
    # Reflection of point r: r - 2*(r·n)*n = reflection if plane through origin.
    plane_normals = [
        np.array([0,0,1]), np.array([0,1,0]), np.array([1,0,0]),
    ]
    # Add diagonals
    plane_normals += [
        np.array([1,1,0]), np.array([1,-1,0]), np.array([1,0,1]), np.array([1,0,-1]),
        np.array([0,1,1]), np.array([0,1,-1])
    ]
    has_sigma_v = False
    has_sigma_h = False
    has_sigma_d = False
    # We'll just report presence of mirror planes
    mirror_planes = []
    for norm in plane_normals:
        norm = norm / np.linalg.norm(norm)
        reflected = centered - 2 * np.dot(centered, norm)[:, np.newaxis] * norm
        matched_all = True
        used = [False] * n
        for i in range(n):
            found = False
            for j in range(n):
                if used[j] or symbols[i] != symbols[j]:
                    continue
                if np.linalg.norm(reflected[i] - centered[j]) < TOL:
                    used[j] = True
                    found = True
                    break
            if not found:
                matched_all = False
                break
        if matched_all:
            mirror_planes.append(tuple(np.round(norm, 3)))

    # Deduplicate normals (and opposite directions represent same plane for reflection)
    unique_planes = []
    seen = set()
    for p in mirror_planes:
        # Represent plane by absolute direction and sorted tuple
        abs_p = tuple(np.round(np.abs(np.array(p)), 3))
        if abs_p not in seen:
            seen.add(abs_p)
            unique_planes.append(p)

    # Determine point group heuristically based on detected elements
    orders = sorted([k[0] for k in unique_axes.keys()], reverse=True)
    max_order = orders[0] if orders else 0
    num_planes = len(unique_planes)

    # Classification rules (simplified for common molecular point groups)
    point_group = "C1"
    if max_order == 0 and not unique_planes and not has_inversion:
        point_group = "C1"
    elif max_order == 0 and has_inversion:
        point_group = "Ci"
    elif max_order == 0 and num_planes >= 1:
        point_group = "Cs"
    elif max_order == 2:
        if has_inversion:
            point_group = "C2h"
        elif num_planes >= 1:
            # Distinguish C2v vs C2h (C2h already handled)
            # Check if plane perpendicular to axis exists (sigma_h) -> C2h handled
            point_group = "C2v"
        else:
            point_group = "C2"
    elif max_order == 3:
        if has_inversion:
            point_group = "S6"
        elif num_planes >= 3:
            point_group = "D3h"
        elif num_planes >= 1:
            point_group = "C3v"
        else:
            point_group = "C3"
    elif max_order == 4:
        if has_inversion:
            point_group = "D4h"
        elif num_planes >= 4:
            point_group = "D4h"  # approximate
        elif num_planes >= 2:
            point_group = "D4"
        elif num_planes >= 1:
            point_group = "C4v"
        else:
            point_group = "C4"
    elif max_order == 5:
        if num_planes >= 5:
            point_group = "D5h"
        elif num_planes >= 1:
            point_group = "C5v"
        else:
            point_group = "C5"
    elif max_order == 6:
        if num_planes >= 6:
            point_group = "D6h"
        elif num_planes >= 1:
            point_group = "C6v"
        else:
            point_group = "C6"

    # Linear molecules (all atoms collinear) have special groups Cinfv / Dinfh
    coords_flat = centered
    # Check collinearity
    if n > 1:
        # Check if all points lie on a line through center
        # Compute cross products with first non-zero vector
        v_ref = None
        for i in range(n):
            if np.linalg.norm(coords_flat[i]) > TOL:
                v_ref = coords_flat[i]
                break
        linear = True
        if v_ref is not None:
            for i in range(n):
                if np.linalg.norm(np.cross(v_ref, coords_flat[i])) > TOL:
                    linear = False
                    break
        else:
            linear = False  # all at center (degenerate)
        if linear:
            point_group = "Dinfh" if has_inversion else "Cinfv"

    return {
        'point_group': point_group,
        'inversion': has_inversion,
        'rotation_axes': unique_axes,
        'mirror_planes': unique_planes,
        'linear': linear if 'linear' in locals() else False,
    }


def main():
    if len(sys.argv) < 2:
        print("Usage: python symmetry_finder.py <molecule.xyz>")
        sys.exit(1)
    path = sys.argv[1]
    symbols, coords = read_xyz(path)
    result = find_symmetry_elements(symbols, coords)
    print(f"Point group: {result['point_group']}")
    print(f"Inversion center: {'Yes' if result['inversion'] else 'No'}")
    print(f"Rotation axes (order, axis): {[(k, v) for k,v in result['rotation_axes'].items()]}")
    print(f"Mirror planes (normals): {result['mirror_planes']}")
    if result.get('linear'):
        print("Molecule is linear.")


if __name__ == '__main__':
    main()
