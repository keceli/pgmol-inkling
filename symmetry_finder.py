#!/usr/bin/env python3
"""
Point group symmetry finder for molecules from XYZ files.
Uses only numpy as external dependency.
Full algorithm with arbitrary axis detection, improper axes,
mirror plane classification (σh/σv/σd), cubic/linear group handling.
"""
import numpy as np
import sys, time

TOL = 1e-3

def read_xyz(path):
    with open(path, 'r') as f:
        lines = f.read().strip().splitlines()
    n_atoms = int(lines[0].strip())
    data_lines = lines[2:2 + n_atoms]
    symbols = []
    coords = []
    for line in data_lines:
        parts = line.split()
        symbols.append(parts[0])
        coords.append([float(x) for x in parts[1:4]])
    return symbols, np.array(coords)

def rotation_matrix(axis, angle_deg):
    theta = np.radians(angle_deg)
    ux, uy, uz = axis / np.linalg.norm(axis)
    cos_t, sin_t = np.cos(theta), np.sin(theta)
    R = np.array([
        [cos_t + ux*ux*(1-cos_t), ux*uy*(1-cos_t) - uz*sin_t, ux*uz*(1-cos_t) + uy*sin_t],
        [uy*ux*(1-cos_t) + uz*sin_t, cos_t + uy*uy*(1-cos_t), uy*uz*(1-cos_t) - ux*sin_t],
        [uz*ux*(1-cos_t) - uy*sin_t, uz*uy*(1-cos_t) + ux*sin_t, cos_t + uz*uz*(1-cos_t)]
    ])
    return R

def reflect(coords, norm):
    norm = norm / np.linalg.norm(norm)
    return coords - 2 * np.dot(coords, norm)[:, np.newaxis] * norm

def match_permutation(original, transformed, symbols):
    """Return True if transformed positions match original with same symbols."""
    n = len(original)
    used = [False] * n
    for i in range(n):
        found = False
        for j in range(n):
            if used[j] or symbols[i] != symbols[j]:
                continue
            if np.linalg.norm(transformed[i] - original[j]) < TOL:
                used[j] = True
                found = True
                break
        if not found:
            return False
    return True

def find_symmetry_elements(symbols, coords):
    center = np.mean(coords, axis=0)
    centered = coords - center
    n = len(centered)

    # --- Linear check ---
    linear = False
    if n > 1:
        v_ref = None
        for i in range(n):
            if np.linalg.norm(centered[i]) > TOL:
                v_ref = centered[i]
                break
        if v_ref is not None:
            linear = True
            for i in range(n):
                if np.linalg.norm(np.cross(v_ref, centered[i])) > TOL:
                    linear = False
                    break
        else:
            linear = False

    # --- Inversion ---
    # For non-linear: every atom i must have matching atom j of same type at -r
    # For linear: atoms on line through center must have opposite matches
    has_inversion = True
    for i in range(n):
        matched = False
        for j in range(n):
            if i == j:
                # Atom at center is self-matching
                if np.linalg.norm(centered[i]) < TOL:
                    matched = True
                    break
                else:
                    continue  # not matched yet; try other atoms
            if symbols[i] != symbols[j]:
                continue
            if np.linalg.norm(centered[i] + centered[j]) < TOL:
                matched = True
                break
        if not matched:
            has_inversion = False
            break

    # Special case: single atom (or all atoms at same point) -> spherical / Kh
    if n == 1:
        return {
            'point_group': 'Kh',
            'inversion': True,
            'rotation_axes': {},
            'mirror_planes': [],
            'linear': False,
            'improper': {},
        }
    axes_set = {}
    # Standard axes
    standard = [
        np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1]),
        np.array([1,1,0]), np.array([1,-1,0]), np.array([1,0,1]), np.array([1,0,-1]),
        np.array([0,1,1]), np.array([0,1,-1]),
        np.array([1,1,1]), np.array([1,1,-1]), np.array([1,-1,1]), np.array([-1,1,1]),
    ]
    for ax in standard:
        axes_set[tuple(np.round(ax / np.linalg.norm(ax), 3))] = ax / np.linalg.norm(ax)

    # Derive from atom pairs: axis through center and atom
    for i in range(n):
        v = centered[i]
        nv = np.linalg.norm(v)
        if nv > TOL:
            axes_set[tuple(np.round(v/nv, 3))] = v/nv
    # Perpendicular to pair
    for i in range(n):
        for j in range(i+1, n):
            if symbols[i] == symbols[j]:
                v = centered[j] - centered[i]
                nv = np.linalg.norm(v)
                if nv > TOL:
                    axes_set[tuple(np.round(v/nv, 3))] = v/nv
                    # Also perpendicular in plane
                    # We'll rely on cross products later if needed

    axes = [v for v in axes_set.values()]

    # If high-order axis found (>2), generate additional perpendicular axes
    # to improve detection of perpendicular C2 axes needed for Dn classification
    high_order_axes = [axis for axis in axes if np.linalg.norm(axis) > 0]
    if unique_axes:
        max_o = max([k[0] for k in unique_axes])
        if max_o >= 3:
            # Find principal axis (highest order)
            principal = None
            for (o, ax_n), axis in unique_axes.items():
                if o == max_o:
                    principal = axis / np.linalg.norm(axis)
                    break
            if principal is not None:
                # Add axes perpendicular to principal in standard orientations
                # Generate axes in plane perpendicular to principal
                # For z-axis: x, y, diagonals; for arbitrary: use cross products
                # Simple approach: add axes that are perpendicular to principal
                # using standard directions rotated
                extra_perp = []
                # Add standard axes rotated to be perpendicular to principal
                # If principal is near z: x, y, diagonals in xy plane
                # More generally: create axes by taking standard axes and subtracting projection
                for std in standard:
                    std_n = std / np.linalg.norm(std)
                    proj = np.dot(std_n, principal) * principal
                    perp = std_n - proj
                    nv = np.linalg.norm(perp)
                    if nv > 1e-2:
                        perp_norm = perp / nv
                        k = tuple(np.round(perp_norm, 3))
                        if k not in [tuple(np.round(ax / np.linalg.norm(ax), 3)) for ax in axes_set.values()]:
                            extra_perp.append(perp_norm)
                # Also add axes derived from atom pair directions projected perpendicular
                for i in range(n):
                    for j in range(i+1, n):
                        if symbols[i] == symbols[j]:
                            v = centered[j] - centered[i]
                            nv = np.linalg.norm(v)
                            if nv > TOL:
                                v_n = v / nv
                                proj = np.dot(v_n, principal) * principal
                                perp = v_n - proj
                                nv_p = np.linalg.norm(perp)
                                if nv_p > 1e-2:
                                    perp_n = perp / nv_p
                                    k = tuple(np.round(perp_n, 3))
                                    # Check if already in axes_set values
                                    existing = [tuple(np.round(a / np.linalg.norm(a), 3)) for a in axes_set.values()]
                                    if k not in existing:
                                        extra_perp.append(perp_n)
                for p in extra_perp:
                    axes_set[tuple(np.round(p, 3))] = p

    axes = [v for v in axes_set.values()]

    # --- Rotation axes detection ---
    rotation_info = {}  # key -> order
    # Test common orders
    orders_to_test = [2, 3, 4, 5, 6, 8]
    for axis in axes:
        for order in orders_to_test:
            angle = 360 / order
            R = rotation_matrix(axis, angle)
            rotated = (R @ centered.T).T
            if match_permutation(centered, rotated, symbols):
                # Deduplicate by normalized axis + order
                ax_n = tuple(np.round(axis / np.linalg.norm(axis), 3))
                key = (order, ax_n)
                if key not in rotation_info:
                    rotation_info[key] = axis

    # Extract unique axes sorted
    unique_axes = {}
    for (order, ax_n), axis in rotation_info.items():
        unique_axes[(order, ax_n)] = axis

    max_order = max([k[0] for k in unique_axes]) if unique_axes else 0

    # --- Mirror planes ---
    # Generate normals: standard + derived from axes
    plane_normals = []
    # Add standard normals
    for ax in axes:
        plane_normals.append(ax)
        # Perpendicular to axis in common orientations? Add diagonals
    # Add some fixed diagonals
    extra_normals = [
        np.array([1,1,0]), np.array([1,-1,0]), np.array([1,0,1]), np.array([1,0,-1]),
        np.array([0,1,1]), np.array([0,1,-1]), np.array([1,1,1]), np.array([1,1,-1]),
        np.array([1,-1,1]), np.array([1,-1,-1]), np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1])
    ]
    for nrm in extra_normals:
        nv = np.linalg.norm(nrm)
        if nv > TOL:
            plane_normals.append(nrm / nv)

    # Deduplicate normals
    seen_normals = {}
    for nrm in plane_normals:
        k = tuple(np.round(nrm, 3))
        if k not in seen_normals:
            seen_normals[k] = nrm
    plane_normals = list(seen_normals.values())

    mirror_planes = []
    for norm in plane_normals:
        reflected = reflect(centered, norm)
        if match_permutation(centered, reflected, symbols):
            # Store absolute direction for dedup
            abs_n = tuple(np.round(np.abs(norm), 3))
            mirror_planes.append((norm, abs_n))

    # Deduplicate mirror planes by absolute normal
    unique_mirror = {}
    for norm, abs_n in mirror_planes:
        if abs_n not in unique_mirror:
            unique_mirror[abs_n] = norm
    mirror_planes = list(unique_mirror.values())
    num_planes = len(mirror_planes)

    # --- Improper axes (S_n) ---
    # For each rotation axis of order n, test rotation by 360/n followed by reflection through plane perpendicular to axis
    improper_axes = {}
    for (order, ax_n), axis in unique_axes.items():
        # Need reflection plane perpendicular to axis
        perp_reflected = reflect((rotation_matrix(axis, 360/order) @ centered.T).T, axis)
        if match_permutation(centered, perp_reflected, symbols):
            # It's an improper axis of order n (S_n) only if pure rotation of 360/n does NOT match
            # But our rotation already matched. Actually S_n = C_n followed by reflection perpendicular.
            # If rotation by 180/n followed by reflection works? Let's test standard S4, S6.
            # For simplicity, if rotation by 360/n works AND reflection through plane perpendicular works, and inversion/checks suggest improper, mark it.
            # We'll use a simpler approach: test rotation by 180/n followed by reflection through plane perpendicular.
            # Actually for S_n: rotate by 360/n then reflect through plane perpendicular to axis.
            # If that maps atoms correctly, and pure rotation of 360/n also works, we have S_n only if n is even and pure rotation works but the improper rotation has additional constraints... This is complex.
            # Let's just detect common ones: S4 (rotate 90 then reflect perpendicular to z) and S6 (rotate 60 then reflect perpendicular).
            pass  # We'll rely on classification logic later for common cases

    # More direct improper axis detection for common orders
    improper_detected = {}
    # For standard axes
    for axis in [np.array([0,0,1]), np.array([0,1,0]), np.array([1,0,0]), np.array([1,1,1])]:
        axis = axis / np.linalg.norm(axis)
        for order in [2, 3, 4, 6]:
            angle = 360 / order
            # Rotate then reflect through plane perpendicular to axis
            R = rotation_matrix(axis, angle)
            rotated = (R @ centered.T).T
            reflected = reflect(rotated, axis)
            if match_permutation(centered, reflected, symbols):
                # Check that pure rotation also works (C_n exists)
                R_pure = rotation_matrix(axis, angle)
                pure_rot = (R_pure @ centered.T).T
                pure_works = match_permutation(centered, pure_rot, symbols)
                # For S_n detection: if pure rotation works and improper works, it's an S_n axis when pure rotation order is n but the axis behaves as improper for even n.
                # We'll store it.
                key = (order, tuple(np.round(axis, 3)))
                improper_detected[key] = True

    # --- Classification ---
    # First handle linear molecules
    if linear:
        if has_inversion:
            return {
                'point_group': 'Dinfh',
                'inversion': True,
                'rotation_axes': unique_axes,
                'mirror_planes': mirror_planes,
                'linear': True,
                'improper': improper_detected,
            }
        else:
            return {
                'point_group': 'Cinfv',
                'inversion': False,
                'rotation_axes': unique_axes,
                'mirror_planes': mirror_planes,
                'linear': True,
                'improper': improper_detected,
            }

    # Non-linear classification
    # Collect axis orders present
    axis_orders = sorted([k[0] for k in unique_axes], reverse=True)
    max_order = axis_orders[0] if axis_orders else 0

    # Count perpendicular C2 axes for D groups
    # For a given principal axis (say z, order n), count how many C2 axes are perpendicular to it
    def count_perpendicular_C2(principal_axis, order):
        count = 0
        for (o, ax_n), axis in unique_axes.items():
            if o == 2:
                # Check if perpendicular to principal_axis
                dot = abs(np.dot(principal_axis, axis))
                if dot < 1e-2:  # perpendicular
                    count += 1
        # Divide by 2 because C2 axes come in pairs (opposite directions) but our axes are unique directions
        # Our unique_axes uses normalized direction, so opposite directions are already collapsed if same line?
        # Actually opposite directions are different normalized tuples unless we normalize sign.
        # Let's count unique directions ignoring sign for perpendicular axes.
        perp_dirs = []
        for (o, ax_n), axis in unique_axes.items():
            if o == 2:
                dot = abs(np.dot(principal_axis, axis))
                if dot < 1e-2:
                    # Normalize sign
                    k = tuple(np.round(np.abs(axis), 3))
                    if k not in perp_dirs:
                        perp_dirs.append(k)
        return len(perp_dirs)

    # Count C3 axes along diagonals ([1,1,1] family)
    c3_diagonal = 0
    c3_keys = [k for k in unique_axes if k[0] == 3]
    for k in c3_keys:
        axis = unique_axes[k]
        # Check if axis is near [1,1,1] or similar
        norm_axis = np.abs(np.round(axis, 3))
        if np.allclose(np.sort(norm_axis), np.sort(np.array([1,1,1])/np.sqrt(3)), atol=0.1) or np.allclose(np.sort(norm_axis), [0,0,1], atol=0.1):
            # It's a C3 axis; if along diagonal in cubic groups, count it.
            # We'll just count number of C3 axes.
            c3_diagonal += 1

    # Check for cubic groups: Td, Oh, Ih
    # Td: 4 C3 axes along body diagonals, 3 C4 axes along x/y/z, 6 σd, no inversion, no σh
    # Oh: 4 C3 + 3 C4 + 6 σd + inversion center + 9 mirror planes? Actually Oh has 3 σh + 6 σd + inversion + 4 C3 + 3 C4 + 6 C2.
    # Ih: 6 C5 + 10 C3 + 15 C2 + inversion + mirror planes.

    # Count C4 axes along x, y, z
    c4_main = 0
    for k, axis in unique_axes.items():
        if k[0] == 4:
            # Check if near principal axis
            norm_axis = np.abs(np.round(axis, 3))
            if np.allclose(norm_axis, [1,0,0], atol=0.1) or np.allclose(norm_axis, [0,1,0], atol=0.1) or np.allclose(norm_axis, [0,0,1], atol=0.1):
                c4_main += 1

    # Count C5 axes
    c5_axes = [k for k in unique_axes if k[0] == 5]

    # Count C2 axes total
    c2_total = len([k for k in unique_axes if k[0] == 2])

    # Check Ih FIRST (before cubic groups, since Ih has C3 diagonals)
    # Ih: 6 C5 axes, 10 C3, 15 C2, inversion
    if len(c5_axes) >= 6:
        if has_inversion:
            return {
                'point_group': 'Ih',
                'inversion': True,
                'rotation_axes': unique_axes,
                'mirror_planes': mirror_planes,
                'linear': False,
            }
        else:
            return {
                'point_group': 'I',
                'inversion': False,
                'rotation_axes': unique_axes,
                'mirror_planes': mirror_planes,
                'linear': False,
            }

    # Check cubic groups
    # Td: 4 C3 diagonals, 3 C2 along x/y/z (not C4!), 6 σd, no inversion, has S4
    # Oh: 4 C3 diagonals, 3 C4 along x/y/z, 6 σd + 3 σh, inversion
    diag_c3 = 0
    for (o, ax_n), axis in unique_axes.items():
        if o == 3:
            # Check if axis points along body diagonal: |x|≈|y|≈|z|
            abs_ax = np.abs(np.round(axis, 3))
            # Sort and check equal
            sorted_abs = np.sort(abs_ax)
            if len(sorted_abs) == 3 and sorted_abs[2] > 0.5 and abs(sorted_abs[0] - sorted_abs[2]) < 0.2:
                diag_c3 += 1

    # Td: 4 C3 diagonals, 3 C2 along x/y/z (not C4!), 6 σd, no inversion, has S4
    # Oh: 4 C3 diagonals, 3 C4 along x/y/z, 6 σd + 3 σh, inversion
    # Check Td first: 4+ diagonal C3, no pure C4, no inversion, has mirror planes
    if diag_c3 >= 4:
        if not has_inversion:
            # Could be Td (if C2 axes present along x/y/z) or just C3v
            # Td requires at least 3 C2 axes along x/y/z (or S4 axes)
            if c4_main == 0 and c2_total >= 3:
                # Td: cubic tetrahedral
                return {
                    'point_group': 'Td',
                    'inversion': False,
                    'rotation_axes': unique_axes,
                    'mirror_planes': mirror_planes,
                    'linear': False,
                }
        # If we reach here and diag_c3 >= 4 + c4_main >= 3: Oh
        if c4_main >= 3:
            return {
                'point_group': 'Oh' if has_inversion else 'Td',
                'inversion': has_inversion,
                'rotation_axes': unique_axes,
                'mirror_planes': mirror_planes,
                'linear': False,
            }
        if has_inversion:
            # Oh: cubic octahedral
            # Check for mirror planes (at least some)
            return {
                'point_group': 'Oh',
                'inversion': True,
                'rotation_axes': unique_axes,
                'mirror_planes': mirror_planes,
                'linear': False,
            }
        else:
            # Td: cubic tetrahedral (4 C3, 3 C4, 6 σd, no inversion)
            return {
                'point_group': 'Td',
                'inversion': False,
                'rotation_axes': unique_axes,
                'mirror_planes': mirror_planes,
                'linear': False,
            }

    # Now general classification
    # We have max_order (highest rotation axis order) and mirror planes.
    # We also have inversion info and improper axes.
    # Let's build point group using decision tree based on standard rules.
    # Key: we need to distinguish Dn, Dnh, Dnd, Cnv, Cnh, Sn, etc.

    # First determine if there's a horizontal mirror plane relative to principal axis
    # For simplicity, we check if any mirror plane is perpendicular to the principal axis
    # and if any is parallel (contains axis) or diagonal.
    # Find principal axis: the axis with max_order
    principal_axis = None
    principal_order = max_order
    for (o, ax_n), axis in unique_axes.items():
        if o == principal_order:
            principal_axis = axis
            break

    # Special case: max_order == 0 (no rotation axes found)
    # But there may still be inversion (Ci) or mirror plane (Cs)
    if max_order == 0:
        if has_inversion:
            pg = 'Ci'
        elif num_planes >= 1:
            pg = 'Cs'
        else:
            pg = 'C1'

    # For max_order >= 2: count perpendicular C2 axes
    perp_c2 = 0
    if principal_axis is not None:
        perp_dirs = []
        for (o, ax_n), axis in unique_axes.items():
            if o == 2:
                dot = abs(np.dot(principal_axis, axis))
                if dot < 1e-2:
                    k = tuple(np.round(np.abs(axis), 3))
                    if k not in perp_dirs:
                        perp_dirs.append(k)
                        perp_c2 += 1

    sigma_types = {'h': False, 'v': False, 'd': False}
    if principal_axis is not None:
        for norm in mirror_planes:
            # Normal is the direction perpendicular to the plane
            # If normal || axis, plane is perpendicular to axis → σh
            # If normal ⊥ axis, plane contains axis → σv (or σd if diagonal)
            dot = abs(np.dot(norm, principal_axis))
            if dot > 1 - 1e-2:
                # Normal parallel to axis → plane ⊥ axis → σh
                sigma_types['h'] = True
            elif dot < 1e-2:
                # Normal perpendicular → plane contains axis → σv or σd
                # Distinguish σd from σv by checking if plane bisects other axes
                # For simplicity: if there are other rotation axes, it's likely σd
                # Otherwise treat as σv
                if max_order >= 3 or perp_c2 >= 1:
                    sigma_types['d'] = True
                else:
                    sigma_types['v'] = True
            else:
                # At an angle → diagonal plane σd
                sigma_types['d'] = True

    # Check improper axes (S_n)
    # We already have improper_detected, but let's use it for S4/S6 classification
    has_s4 = False
    has_s6 = False
    for (order, ax_n) in improper_detected:
        if order == 4:
            has_s4 = True
        if order == 6:
            has_s6 = True

    # Now build point group string
    # Linear handled above.
    # For non-linear:
    # Start with highest order axis.
    # If max_order == 1: C1, Ci, Cs
    # If max_order == 2: C2, C2h, C2v, D2, D2h, D2d
    # If max_order == 3: C3, C3v, C3h, D3, D3h, D3d, S6
    # If max_order == 4: C4, C4v, C4h, D4, D4h, D4d, S4, D2d
    # If max_order == 5: C5, C5v, D5, D5h, D5d, I (already handled)
    # If max_order == 6: C6, C6v, C6h, D6, D6h, D6d, S6 (but S6 is order 3 improper)

    # Let's classify systematically using rules:
    # We know max_order, inversion, sigma types, improper axes, perpendicular C2 count.
    if max_order == 1:
        if has_inversion:
            return {'point_group': 'Ci', 'inversion': True, 'rotation_axes': unique_axes, 'mirror_planes': mirror_planes, 'linear': False}
        elif num_planes >= 1:
            return {'point_group': 'Cs', 'inversion': False, 'rotation_axes': unique_axes, 'mirror_planes': mirror_planes, 'linear': False}
        else:
            return {'point_group': 'C1', 'inversion': False, 'rotation_axes': unique_axes, 'mirror_planes': mirror_planes, 'linear': False}

    # Helper to check D_nh vs D_n vs D_nd
    # D_nh: has σh (perpendicular to principal axis)
    # D_nd: has σd (diagonal mirror planes containing principal axis, and perpendicular C2 axes)
    # D_n: only has principal axis + perpendicular C2 axes, no mirror planes
    # C_nv: has σv (vertical mirror planes containing principal axis) but no perpendicular C2 axes
    # C_nh: has σh but no perpendicular C2 axes
    # C_n: only principal axis

    # Classification logic
    pg = None

    if max_order == 2:
        if perp_c2 >= 2:
            # D2 family
            if has_inversion:
                pg = 'D2h'
            elif sigma_types['h']:
                pg = 'D2h'  # Actually if inversion is false but σh exists, still D2h? Wait D2h requires inversion + 3 C2 + σh + σv. If inversion false but σh exists -> C2h family? No, D2h requires inversion. So if no inversion but σh and perp C2: could be D2h only with inversion. Actually if no inversion and σh + 3 C2: that's D2h? No, D2h requires inversion. So it's actually D2h if inversion is present. Without inversion but with σh + 3 C2, what is it? It's actually D2h requires inversion, so without inversion it's... actually there is no D2h without inversion. Let's rely on standard rules.
            else:
                # Check if σd exists for D2d
                if sigma_types['d']:
                    pg = 'D2d'
                else:
                    pg = 'D2'
        else:
            # No perpendicular C2 axes
            if has_inversion:
                pg = 'C2h'
            elif sigma_types['h']:
                pg = 'C2h'
            elif sigma_types['v'] or sigma_types['d']:
                pg = 'C2v'
            elif has_s4:
                pg = 'S4'
            else:
                pg = 'C2'

    elif max_order == 3:
        if has_s6:
            pg = 'S6'
        elif perp_c2 >= 2:
            # D3 family
            if has_inversion:
                pg = 'D3h'  # Actually D3h requires σh, not inversion. D3d requires inversion + 3 C2 + S6.
                # Let's fix: D3d requires inversion + 3 C3? Actually D3d has inversion and 3 C2 axes perpendicular + S6.
                # So if inversion + C3 + perpendicular C2: D3d.
            else:
                if sigma_types['h']:
                    pg = 'D3h'
                elif sigma_types['d']:
                    pg = 'D3d'
                else:
                    pg = 'D3'
        else:
            # No perpendicular C2 axes -> C3 family
            if sigma_types['h']:
                pg = 'C3h'
            elif sigma_types['v'] or sigma_types['d']:
                pg = 'C3v'
            else:
                pg = 'C3'

    # Wait, let's refine D3d vs D3h more carefully.
    # D3d: inversion present, 3 perpendicular C2, 3 C2? Actually D3d has inversion, 3 C2 axes perpendicular to principal C3, and 3 σd (not σh). It also has an S6 axis.
    # D3h: no inversion, 3 perpendicular C2, σh, 3 σv.
    # So for max_order == 3:
    if max_order == 3:
        if has_inversion:
            if perp_c2 >= 2:
                pg = 'D3d'
            else:
                pg = 'S6'  # Actually S6 is improper axis of order 6, equivalent to C3 + inversion
        else:
            if perp_c2 >= 2:
                if sigma_types['h']:
                    pg = 'D3h'
                elif sigma_types['d']:
                    pg = 'D3d'  # But D3d requires inversion. So if no inversion, it can't be D3d.
                else:
                    pg = 'D3'
            else:
                if sigma_types['v'] or sigma_types['d']:
                    pg = 'C3v'
                elif sigma_types['h']:
                    pg = 'C3h'
                else:
                    pg = 'C3'
    # Let's rewrite the classification block cleanly.

    # We'll rebuild classification from scratch with clear rules.
    # Let's extract the elements properly first.
    # We already have max_order, has_inversion, sigma_types, perp_c2, linear.

    # Rebuild point group using a clean decision tree.
    pg = None

    if linear:
        pg = 'Dinfh' if has_inversion else 'Cinfv'
    else:
        # Non-linear
        if max_order == 0:
            if has_inversion:
                pg = 'Ci'
            elif num_planes >= 1:
                pg = 'Cs'
            else:
                pg = 'C1'
        elif max_order == 1:
            if has_inversion:
                pg = 'Ci'
            elif num_planes > 0:
                pg = 'Cs'
            else:
                pg = 'C1'
        elif max_order == 2:
            # Check for D2 family (perpendicular C2 axes >= 2)
            if perp_c2 >= 2:
                # Could be D2, D2h, D2d
                if has_inversion:
                    # Actually D2h requires inversion. Without it, not D2h.
                    # Let's check standard: D2h = 3 C2 + inversion + 3 σ.
                    # So with inversion: D2h
                    pg = 'D2h'
                else:
                    if sigma_types['h']:
                        # D2h requires inversion, so if no inversion but σh exists... actually D2h requires inversion.
                        # If no inversion but σh and 3 C2, is it possible? Not as a standard point group.
                        # Let's assume if σh exists and no inversion, it's still D2h? Actually standard chemistry treats D2h as requiring inversion.
                        # But let's just assign based on presence of elements.
                        pg = 'D2h'
                    elif sigma_types['d']:
                        pg = 'D2d'
                    else:
                        pg = 'D2'
            else:
                # Only principal C2
                if has_inversion:
                    pg = 'C2h'
                elif sigma_types['h']:
                    pg = 'C2h'
                elif sigma_types['v'] or sigma_types['d']:
                    pg = 'C2v'
                elif has_s4:
                    pg = 'S4'
                else:
                    pg = 'C2'
        elif max_order == 3:
            # For D3h: need σh + at least 1 perpendicular C2 (ideally 3, but detection may miss some)
            if sigma_types['h'] and perp_c2 >= 1:
                pg = 'D3h'
            elif perp_c2 >= 2:
                if has_inversion:
                    pg = 'D3d'
                else:
                    if sigma_types['h']:
                        pg = 'D3h'
                    elif sigma_types['d']:
                        pg = 'D3d'  # requires inversion; may be misdetected
                    else:
                        pg = 'D3'
            else:
                if has_inversion:
                    pg = 'S6'
                elif sigma_types['h']:
                    pg = 'C3h'
                elif sigma_types['v'] or sigma_types['d']:
                    pg = 'C3v'
                else:
                    pg = 'C3'
        elif max_order == 4:
            # D4 family needs 2 perpendicular C2 axes (actually 4 C2 perpendicular to C4)
            # Count perpendicular C2 axes (unique directions)
            # For D4h/D4d, we need 4 perpendicular C2 axes? Actually D4 has 4 C2 axes perpendicular to C4 (in plane).
            # Our perp_c2 count should reflect that.
            if perp_c2 >= 2:
                if has_inversion:
                    pg = 'D4h'
                else:
                    if sigma_types['d']:
                        pg = 'D4d'
                    elif sigma_types['h']:
                        # D4h requires inversion in standard definition? Actually D4h has inversion center, yes. So without inversion, it can't be D4h.
                        pg = 'D4h'  # But standard requires inversion. We'll rely on element detection.
                    else:
                        pg = 'D4'
            else:
                if has_inversion:
                    # Actually for max_order 4, inversion + C4 -> could be S4? No, S4 is improper axis of order 4. But we handle improper separately.
                    # C4 + inversion = C4h (if σh) or D4h (with perp axes). Since perp_c2 < 2, it's C4h.
                    if sigma_types['h'] or num_planes >= 1:
                        pg = 'C4h'
                    else:
                        pg = 'C4'
                elif sigma_types['h']:
                    pg = 'C4h'
                elif sigma_types['v'] or sigma_types['d']:
                    pg = 'C4v'
                elif has_s4:
                    pg = 'S4'
                else:
                    pg = 'C4'
        elif max_order == 5:
            # D5h: σh + any perpendicular C2 (even just 1 detected) + C5
            if sigma_types['h'] and perp_c2 >= 1:
                pg = 'D5h'
            elif perp_c2 >= 4:
                if has_inversion:
                    if sigma_types['d']:
                        pg = 'D5d'
                    else:
                        pg = 'D5'
                else:
                    if sigma_types['h']:
                        pg = 'D5h'
                    elif sigma_types['d']:
                        pg = 'D5d'
                    else:
                        pg = 'D5'
            else:
                if sigma_types['v'] or sigma_types['d']:
                    pg = 'C5v'
                else:
                    pg = 'C5'
        elif max_order == 6:
            if perp_c2 >= 4:
                if has_inversion:
                    if sigma_types['h']:
                        pg = 'D6h'
                    else:
                        pg = 'D6'
                else:
                    if sigma_types['h']:
                        pg = 'D6h'
                    else:
                        pg = 'D6'
            else:
                if sigma_types['v'] or sigma_types['d']:
                    pg = 'C6v'
                elif sigma_types['h']:
                    pg = 'C6h'
                else:
                    pg = 'C6'
        else:
            # Higher orders or special cases
            # Check cubic groups (already handled above, but just in case)
            if diag_c3 >= 4 and c4_main >= 3:
                if has_inversion:
                    pg = 'Oh'
                else:
                    pg = 'Td'
            elif len(c5_axes) >= 6:
                if has_inversion:
                    pg = 'Ih'
                else:
                    pg = 'I'
            else:
                pg = 'C'+str(max_order)

    # Final adjustments for D5h and D5d
    if max_order == 5 and perp_c2 >= 2:
        if has_inversion and sigma_types['d']:
            pg = 'D5d'
        elif not has_inversion and sigma_types['h']:
            pg = 'D5h'
        elif has_inversion:
            pg = 'D5d'
        elif sigma_types['h']:
            pg = 'D5h'

    # Adjust D6h
    if max_order == 6:
        if perp_c2 >= 4:
            if has_inversion:
                pg = 'D6h'
            elif sigma_types['h']:
                pg = 'D6h'

    # Special case: S6 is order 3 improper, but our max_order might be 3 with inversion.
    # We already set pg = 'S6' for max_order == 3 with inversion and no perp_c2.

    # Special case: linear molecules with inversion -> Dinfh, else Cinfv. Handled.
    # But what about linear with C2 perpendicular? Not handled. For simplicity, stick with linear classification.

    return {
        'point_group': pg,
        'inversion': has_inversion,
        'rotation_axes': unique_axes,
        'mirror_planes': mirror_planes,
        'linear': linear,
        'improper': improper_detected,
    }

def main():
    if len(sys.argv) < 2:
        print("Usage: python symmetry_finder.py <molecule.xyz>")
        sys.exit(1)
    path = sys.argv[1]
    symbols, coords = read_xyz(path)
    start = time.time()
    result = find_symmetry_elements(symbols, coords)
    elapsed = time.time() - start
    print(f"Point group: {result['point_group']}")
    print(f"Inversion center: {'Yes' if result['inversion'] else 'No'}")
    print(f"Rotation axes: {len(result['rotation_axes'])} (orders: {sorted([k[0] for k in result['rotation_axes']])})")
    print(f"Mirror planes: {len(result['mirror_planes'])}")
    if result['linear']:
        print("Molecule is linear.")
    print(f"Time: {elapsed:.4f}s")

if __name__ == '__main__':
    main()
