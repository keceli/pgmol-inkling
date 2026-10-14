# Point Group Symmetry Benchmark - Challenge XYZ

- **Run timestamp:** 2026-10-14T19:59:15.094180
- **Files tested:** 55
- **Library:** numpy only (`symmetry_finder.py`)

## Results Table

| File | Expected | Computed | Match | Time (s) |
|---|---|---|---|---|
| `B12H60_I.xyz` | I | I | MATCH | 0.5964 |
| `B2H6_D2h.xyz` | D2h | D2h | MATCH | 0.0643 |
| `BH3O3_C3h.xyz` | C3h | C3h | MATCH | 0.0622 |
| `C10FeH10_D5.xyz` | D5 | C5 | MISMATCH | 0.0808 |
| `C10FeH10_D5d.xyz` | D5d | D5d | MATCH | 0.0779 |
| `C10H16_Td.xyz` | Td | Td | MATCH | 0.095 |
| `C12CrH12_D6.xyz` | D6 | C6 | MISMATCH | 0.0952 |
| `C12CrH12_D6d.xyz` | D6d | C6v | MISMATCH | 0.0933 |
| `C12CrH12_D6h.xyz` | D6h | D6h | MATCH | 0.0935 |
| `C16H16U_D8d.xyz` | D8d | D4d | MISMATCH | 0.124 |
| `C16H16U_D8h.xyz` | D8h | D4h | MISMATCH | 0.1237 |
| `C20H20_Ih.xyz` | Ih | Oh | MISMATCH | 0.1554 |
| `C2Cl2H2_C2v.xyz` | C2v | C2v | MATCH | 0.0669 |
| `C2H2NO_C1.xyz` | C1 | C1 | MATCH | 0.0694 |
| `C2H2O2_D2h.xyz` | D2h | D2h | MATCH | 0.0697 |
| `C2H2_D2h.xyz` | D2h | D2h | MATCH | 0.0648 |
| `C2H4_D2.xyz` | D2 | D2 | MATCH | 0.0729 |
| `C2H6_D3.xyz` | D3 | C3 | MISMATCH | 0.0657 |
| `C3H6_D3h.xyz` | D3h | D3h | MATCH | 0.0633 |
| `C4H4P_C4.xyz` | C4 | C4 | MATCH | 0.0664 |
| `C4H4Si_S4.xyz` | S4 | S4 | MATCH | 0.085 |
| `C4H4_D4h.xyz` | D4h | D4h | MATCH | 0.0652 |
| `C4H8_D4h.xyz` | D4h | D4h | MATCH | 0.0733 |
| `C4N4Pt_C4h.xyz` | C4h | C4h | MATCH | 0.0656 |
| `C5H12_T.xyz` | T | Td | MISMATCH | 0.0766 |
| `C5H12_Td.xyz` | Td | Td | MATCH | 0.0775 |
| `C5H5NNiO_C5v.xyz` | C5v | C5v | MATCH | 0.0662 |
| `C5H5P_C5.xyz` | C5 | C5 | MATCH | 0.0645 |
| `C5H5_C5h.xyz` | C5h | C5 | MISMATCH | 0.0628 |
| `C5H5_D5h.xyz` | D5h | D5h | MATCH | 0.0643 |
| `C60_Ih.xyz` | Ih | Oh | MISMATCH | 1.2036 |
| `C6CrH6_C6v.xyz` | C6v | C6v | MATCH | 0.0684 |
| `C6H12_D3d.xyz` | D3d | D3d | MATCH | 0.0797 |
| `C6H6P_C6.xyz` | C6 | C6 | MATCH | 0.0836 |
| `C6H6_S6.xyz` | S6 | S6 | MATCH | 0.0681 |
| `C6O6_C6h.xyz` | C6h | C6h | MATCH | 0.0664 |
| `C8H8_D4.xyz` | D4 | D4 | MATCH | 0.0746 |
| `C8H8_S8.xyz` | S8 | C4 | MISMATCH | 0.072 |
| `C9CrH6O3_C3v.xyz` | C3v | C3v | MATCH | 0.0719 |
| `COS_Cinfv.xyz` | Cinfv | Cinfv | MATCH | 0.0602 |
| `Cl2OS_Cs.xyz` | Cs | Cs | MATCH | 0.0613 |
| `Cl3MoN3_C2v.xyz` | C2v | C2v | MATCH | 0.0634 |
| `Cl3MoN3_C3v.xyz` | C3v | C3v | MATCH | 0.0617 |
| `CoN6O12_Th.xyz` | Th | Oh | MISMATCH | 0.0891 |
| `F4S_C2v.xyz` | C2v | C2v | MATCH | 0.0624 |
| `F7I_D5h.xyz` | D5h | D5h | MATCH | 0.0662 |
| `F8Zr_D4d.xyz` | D4d | D4d | MATCH | 0.0682 |
| `Fe4S4_Td.xyz` | Td | Td | MATCH | 0.0647 |
| `H24N6Ni_O.xyz` | O | Td | MISMATCH | 0.1318 |
| `H3_Cinfv.xyz` | Cinfv | Cinfv | MATCH | 0.0616 |
| `H3_D3h.xyz` | D3h | D3h | MATCH | 0.0623 |
| `I3_Dinfh.xyz` | Dinfh | Dinfh | MATCH | 0.0645 |
| `N2_Dinfh.xyz` | Dinfh | Dinfh | MATCH | 0.0621 |
| `Na8_Oh.xyz` | Oh | Oh | MATCH | 0.0707 |
| `S8_D4d.xyz` | D4d | D4d | MATCH | 0.068 |

## Summary

- **Matches:** 42 / 55 (76.4%)
- **Average time per file:** 0.1056 s
- **Total run time:** 5.8083 s

## Mismatches

- `C10FeH10_D5.xyz`: expected `D5`, got `C5`
- `C12CrH12_D6.xyz`: expected `D6`, got `C6`
- `C12CrH12_D6d.xyz`: expected `D6d`, got `C6v`
- `C16H16U_D8d.xyz`: expected `D8d`, got `D4d`
- `C16H16U_D8h.xyz`: expected `D8h`, got `D4h`
- `C20H20_Ih.xyz`: expected `Ih`, got `Oh`
- `C2H6_D3.xyz`: expected `D3`, got `C3`
- `C5H12_T.xyz`: expected `T`, got `Td`
- `C5H5_C5h.xyz`: expected `C5h`, got `C5`
- `C60_Ih.xyz`: expected `Ih`, got `Oh`
- `C8H8_S8.xyz`: expected `S8`, got `C4`
- `CoN6O12_Th.xyz`: expected `Th`, got `Oh`
- `H24N6Ni_O.xyz`: expected `O`, got `Td`

---
*Benchmark generated automatically.*
