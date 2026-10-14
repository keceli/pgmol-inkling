# Point Group Symmetry Benchmark

- **Run timestamp:** 2026-10-14T18:23:28.644583
- **Files tested:** 36
- **Library:** numpy only (`symmetry_finder.py`)
- **Environment:** Linux aarch64, python3-numpy 1.24.2

## Results Table

| File | Expected | Computed | Match | Time (s) |
|---|---|---|---|---|
| `Ar_Kh.xyz` | Kh | D6h | MISMATCH | 0.0666 |
| `B12H12_Ih.xyz` | Ih | C2v | MISMATCH | 0.0603 |
| `BF3_D3h.xyz` | D3h | C3v | MISMATCH | 0.0553 |
| `BrCClFH_C1.xyz` | C1 | C1 | MATCH | 0.0581 |
| `BrF5_C4v.xyz` | C4v | D4 | MISMATCH | 0.0655 |
| `C10FeH10_D5h.xyz` | D5h | C5v | MISMATCH | 0.0689 |
| `C10H8_D2h.xyz` | D2h | C2v | MISMATCH | 0.0573 |
| `C2Cl2F2H2_Ci.xyz` | Ci | C1 | MISMATCH | 0.0609 |
| `C2Cl2H2_C2h.xyz` | C2h | C2v | MISMATCH | 0.0703 |
| `C2H2_Dinfh.xyz` | Dinfh | Cinfv | MISMATCH | 0.0629 |
| `C2H4_D2h.xyz` | D2h | C2v | MISMATCH | 0.0612 |
| `C2H6_D3d.xyz` | D3d | C3v | MISMATCH | 0.0619 |
| `C2H6_D3h.xyz` | D3h | C3v | MISMATCH | 0.0652 |
| `C3H3P_C3.xyz` | C3 | C3 | MATCH | 0.0606 |
| `C3H4_D2d.xyz` | D2d | C2v | MISMATCH | 0.061 |
| `C4H6_C2h.xyz` | C2h | C2v | MISMATCH | 0.0586 |
| `C6H6_D6h.xyz` | D6h | C6v | MISMATCH | 0.0654 |
| `C8H8_Oh.xyz` | Oh | D4h | MISMATCH | 0.0665 |
| `CCl2H2_C2v.xyz` | C2v | C2v | MATCH | 0.0614 |
| `CCl3H_C3v.xyz` | C3v | Cs | MISMATCH | 0.062 |
| `CClFH2_Cs.xyz` | Cs | Cs | MATCH | 0.0628 |
| `CH4_Td.xyz` | Td | C2v | MISMATCH | 0.057 |
| `CHN_Cinfv.xyz` | Cinfv | Cinfv | MATCH | 0.0595 |
| `CO2_Dinfh.xyz` | Dinfh | Cinfv | MISMATCH | 0.0583 |
| `CO_Cinfv.xyz` | Cinfv | Cinfv | MATCH | 0.0566 |
| `Cl5P_D3h.xyz` | D3h | C3v | MISMATCH | 0.0599 |
| `ClHO_Cs.xyz` | Cs | Cs | MATCH | 0.0616 |
| `ClH_Cinfv.xyz` | Cinfv | Cinfv | MATCH | 0.0536 |
| `F4Xe_D4h.xyz` | D4h | D4h | MATCH | 0.0629 |
| `F6S_Oh.xyz` | Oh | D4h | MISMATCH | 0.0538 |
| `H2O2_C2.xyz` | C2 | C2 | MATCH | 0.0596 |
| `H2O_C2v.xyz` | C2v | C2v | MATCH | 0.054 |
| `H2_Dinfh.xyz` | Dinfh | Cinfv | MISMATCH | 0.0613 |
| `H3N_C3v.xyz` | C3v | C3v | MATCH | 0.053 |
| `O3_C2v.xyz` | C2v | C2v | MATCH | 0.0568 |
| `P4_Td.xyz` | Td | C2v | MISMATCH | 0.0548 |

## Summary

- **Matches:** 13 / 36 (36.1%)
- **Average time per file:** 0.0604 s
- **Total run time:** 2.1754 s

## Mismatches

- `Ar_Kh.xyz`: expected `Kh`, got `D6h`
- `B12H12_Ih.xyz`: expected `Ih`, got `C2v`
- `BF3_D3h.xyz`: expected `D3h`, got `C3v`
- `BrF5_C4v.xyz`: expected `C4v`, got `D4`
- `C10FeH10_D5h.xyz`: expected `D5h`, got `C5v`
- `C10H8_D2h.xyz`: expected `D2h`, got `C2v`
- `C2Cl2F2H2_Ci.xyz`: expected `Ci`, got `C1`
- `C2Cl2H2_C2h.xyz`: expected `C2h`, got `C2v`
- `C2H2_Dinfh.xyz`: expected `Dinfh`, got `Cinfv`
- `C2H4_D2h.xyz`: expected `D2h`, got `C2v`
- `C2H6_D3d.xyz`: expected `D3d`, got `C3v`
- `C2H6_D3h.xyz`: expected `D3h`, got `C3v`
- `C3H4_D2d.xyz`: expected `D2d`, got `C2v`
- `C4H6_C2h.xyz`: expected `C2h`, got `C2v`
- `C6H6_D6h.xyz`: expected `D6h`, got `C6v`
- `C8H8_Oh.xyz`: expected `Oh`, got `D4h`
- `CCl3H_C3v.xyz`: expected `C3v`, got `Cs`
- `CH4_Td.xyz`: expected `Td`, got `C2v`
- `CO2_Dinfh.xyz`: expected `Dinfh`, got `Cinfv`
- `Cl5P_D3h.xyz`: expected `D3h`, got `C3v`
- `F6S_Oh.xyz`: expected `Oh`, got `D4h`
- `H2_Dinfh.xyz`: expected `Dinfh`, got `Cinfv`
- `P4_Td.xyz`: expected `Td`, got `C2v`

---
*Generated automatically by benchmark script.*
