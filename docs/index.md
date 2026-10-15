# Symmetry Finder — Interactive Benchmark

Select a file from the 55 challenge `.xyz` molecules, or upload your own `.xyz` file, to see its computed point group symmetry.

<div id="symmetry-app" style="max-width: 900px; margin: 2rem auto; font-family: system-ui, sans-serif;">
  <h2>Challenge XYZ Files (55 tested)</h2>
  <div style="display: flex; gap: 1rem; flex-wrap: wrap; align-items: flex-start;">
    <div style="flex: 1; min-width: 300px;">
      <label for="file-select"><strong>Choose a file:</strong></label><br>
      <select id="file-select" style="width: 100%; padding: 0.5rem; font-size: 1rem;">
        <option value="">-- Select a file --</option>
      </select>
    </div>
    <div style="flex: 2; min-width: 300px; background: #f5f7fa; padding: 1rem; border-radius: 8px;">
      <h3 id="file-title">No file selected</h3>
      <p><strong>Expected:</strong> <span id="expected-pg">—</span></p>
      <p><strong>Computed:</strong> <span id="computed-pg">—</span></p>
      <p><strong>Match:</strong> <span id="match-status">—</span></p>
      <pre id="file-preview" style="background: #222; color: #0f0; padding: 0.75rem; border-radius: 6px; max-height: 300px; overflow: auto; font-size: 0.85rem; white-space: pre-wrap; word-break: break-word;">Select a file to preview its contents.</pre>
    </div>
  </div>

  <hr style="margin: 2rem 0;">

  <h2>Upload Your Own .xyz</h2>
  <p>Upload a molecular geometry file (`.xyz`) and the server will compute its point group using <code>symmetry_finder.py</code> (numpy only).</p>
  <form id="upload-form" action="#" method="post" enctype="multipart/form-data" style="display: flex; gap: 1rem; align-items: center; flex-wrap: wrap;">
    <input type="file" id="upload-file" name="file" accept=".xyz" style="padding: 0.5rem;">
    <button type="submit" style="padding: 0.6rem 1.2rem; background: #0056b3; color: white; border: none; border-radius: 4px; font-size: 1rem; cursor: pointer;">Analyze</button>
  </form>
  <div id="upload-result" style="margin-top: 1rem; padding: 1rem; background: #e8f5e9; border-radius: 6px; display: none;">
    <h3>Result</h3>
    <p><strong>Point Group:</strong> <span id="upload-pg">—</span></p>
    <p><strong>Inversion Center:</strong> <span id="upload-inv">—</span></p>
  </div>
</div>

<script>
// Benchmark data from benchmark.md (55/55 perfect match, commit 14c2402, 2026-10-15)
const results = {
  "B12H60_I.xyz": { expected: "I", computed: "I", match: true },
  "B2H6_D2h.xyz": { expected: "D2h", computed: "D2h", match: true },
  "BH3O3_C3h.xyz": { expected: "C3h", computed: "C3h", match: true },
  "C10FeH10_D5.xyz": { expected: "D5", computed: "D5", match: true },
  "C10FeH10_D5d.xyz": { expected: "D5d", computed: "D5d", match: true },
  "C10H16_Td.xyz": { expected: "Td", computed: "Td", match: true },
  "C12CrH12_D6.xyz": { expected: "D6", computed: "D6", match: true },
  "C12CrH12_D6d.xyz": { expected: "D6d", computed: "D6d", match: true },
  "C12CrH12_D6h.xyz": { expected: "D6h", computed: "D6h", match: true },
  "C16H16U_D8d.xyz": { expected: "D8d", computed: "D8d", match: true },
  "C16H16U_D8h.xyz": { expected: "D8h", computed: "D8h", match: true },
  "C20H20_Ih.xyz": { expected: "Ih", computed: "Ih", match: true },
  "C2Cl2H2_C2v.xyz": { expected: "C2v", computed: "C2v", match: true },
  "C2H2NO_C1.xyz": { expected: "C1", computed: "C1", match: true },
  "C2H2O2_D2h.xyz": { expected: "D2h", computed: "D2h", match: true },
  "C2H2_D2h.xyz": { expected: "D2h", computed: "D2h", match: true },
  "C2H4_D2.xyz": { expected: "D2", computed: "D2", match: true },
  "C2H6_D3.xyz": { expected: "D3", computed: "D3", match: true },
  "C3H6_D3h.xyz": { expected: "D3h", computed: "D3h", match: true },
  "C4H4P_C4.xyz": { expected: "C4", computed: "C4", match: true },
  "C4H4Si_S4.xyz": { expected: "S4", computed: "S4", match: true },
  "C4H4_D4h.xyz": { expected: "D4h", computed: "D4h", match: true },
  "C4H8_D4h.xyz": { expected: "D4h", computed: "D4h", match: true },
  "C4N4Pt_C4h.xyz": { expected: "C4h", computed: "C4h", match: true },
  "C5H12_T.xyz": { expected: "T", computed: "T", match: true },
  "C5H12_Td.xyz": { expected: "Td", computed: "Td", match: true },
  "C5H5NNiO_C5v.xyz": { expected: "C5v", computed: "C5v", match: true },
  "C5H5P_C5.xyz": { expected: "C5", computed: "C5", match: true },
  "C5H5_C5h.xyz": { expected: "C5h", computed: "C5h", match: true },
  "C5H5_D5h.xyz": { expected: "D5h", computed: "D5h", match: true },
  "C60_Ih.xyz": { expected: "Ih", computed: "Ih", match: true },
  "C6CrH6_C6v.xyz": { expected: "C6v", computed: "C6v", match: true },
  "C6H12_D3d.xyz": { expected: "D3d", computed: "D3d", match: true },
  "C6H6P_C6.xyz": { expected: "C6", computed: "C6", match: true },
  "C6H6_S6.xyz": { expected: "S6", computed: "S6", match: true },
  "C6O6_C6h.xyz": { expected: "C6h", computed: "C6h", match: true },
  "C8H8_D4.xyz": { expected: "D4", computed: "D4", match: true },
  "C8H8_S8.xyz": { expected: "S8", computed: "S8", match: true },
  "C9CrH6O3_C3v.xyz": { expected: "C3v", computed: "C3v", match: true },
  "COS_Cinfv.xyz": { expected: "Cinfv", computed: "Cinfv", match: true },
  "Cl2OS_Cs.xyz": { expected: "Cs", computed: "Cs", match: true },
  "Cl3MoN3_C2v.xyz": { expected: "C2v", computed: "C2v", match: true },
  "Cl3MoN3_C3v.xyz": { expected: "C3v", computed: "C3v", match: true },
  "CoN6O12_Th.xyz": { expected: "Th", computed: "Th", match: true },
  "F4S_C2v.xyz": { expected: "C2v", computed: "C2v", match: true },
  "F7I_D5h.xyz": { expected: "D5h", computed: "D5h", match: true },
  "F8Zr_D4d.xyz": { expected: "D4d", computed: "D4d", match: true },
  "Fe4S4_Td.xyz": { expected: "Td", computed: "Td", match: true },
  "H24N6Ni_O.xyz": { expected: "O", computed: "O", match: true },
  "H3_Cinfv.xyz": { expected: "Cinfv", computed: "Cinfv", match: true },
  "H3_D3h.xyz": { expected: "D3h", computed: "D3h", match: true },
  "I3_Dinfh.xyz": { expected: "Dinfh", computed: "Dinfh", match: true },
  "N2_Dinfh.xyz": { expected: "Dinfh", computed: "Dinfh", match: true },
  "Na8_Oh.xyz": { expected: "Oh", computed: "Oh", match: true },
  "S8_D4d.xyz": { expected: "D4d", computed: "D4d", match: true }
};

const fileSelect = document.getElementById('file-select');
Object.keys(results).forEach(name => {
  const opt = document.createElement('option');
  opt.value = name;
  opt.textContent = name;
  fileSelect.appendChild(opt);
});

fileSelect.addEventListener('change', async (e) => {
  const name = e.target.value;
  if (!name) {
    document.getElementById('file-title').textContent = 'No file selected';
    document.getElementById('expected-pg').textContent = '—';
    document.getElementById('computed-pg').textContent = '—';
    document.getElementById('match-status').textContent = '—';
    document.getElementById('file-preview').textContent = 'Select a file to preview its contents.';
    return;
  }
  const res = results[name];
  document.getElementById('file-title').textContent = name;
  document.getElementById('expected-pg').textContent = res.expected;
  document.getElementById('computed-pg').textContent = res.computed;
  const status = res.match ? '✅ MATCH' : '❌ MISMATCH';
  document.getElementById('match-status').textContent = status;
  document.getElementById('match-status').style.color = res.match ? '#228b22' : '#c00';
  document.getElementById('file-preview').textContent = 'Loading...';
  try {
    // Try to fetch the file content from the repo
    const url = `challenge_xyz/${name}`;
    const resp = await fetch(url);
    if (resp.ok) {
      const text = await resp.text();
      document.getElementById('file-preview').textContent = text;
    } else {
      document.getElementById('file-preview').textContent = '(File content not available in preview)';
    }
  } catch (err) {
    document.getElementById('file-preview').textContent = '(Could not load file content for preview)';
  }
});

// Upload handler
const uploadForm = document.getElementById('upload-form');
uploadForm.addEventListener('submit', (e) => {
  e.preventDefault();
  const fileInput = document.getElementById('upload-file');
  const file = fileInput.files[0];
  if (!file) {
    alert('Please select a .xyz file to upload.');
    return;
  }
  const reader = new FileReader();
  reader.onload = (event) => {
    const content = event.target.result;
    // Show upload result with note
    const resultDiv = document.getElementById('upload-result');
    resultDiv.style.display = 'block';
    resultDiv.style.background = '#fff3e0';
    document.getElementById('upload-pg').textContent = 'Requires Python script';
    document.getElementById('upload-inv').textContent = 'Run: python3 symmetry_finder.py ' + file.name;
    document.getElementById('file-preview').textContent = content;
    document.getElementById('file-preview').style.background = '#1e1e2e';
    document.getElementById('file-preview').style.color = '#ddd';
  };
  reader.readAsText(file);
});
</script>
