#!/usr/bin/env python3
"""
Generate sv_database.html — BIA reference structural variants database.
Run from the directory containing the xlsx files.
"""
import json, os
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
ANN  = os.path.join(BASE, "annotation_data")

# ── Chromosome maps ───────────────────────────────────────────────────────────
acc_to_chr = {}   # CM132547.1  -> Chr01
chr_to_acc = {}   # Chr01       -> CM132547.1
with open(os.path.join(ANN, "BIARef_vs_BIARef.Chrom.txt")) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) == 2:
            acc_to_chr[parts[1]] = parts[0]
            chr_to_acc[parts[0]] = parts[1]

cla_to_acc = {}   # chr1  -> CM132547.1 (for cross-reference display)
acc_to_cla = {}   # CM132547.1 -> chr1
with open(os.path.join(ANN, "CLA_vs_BIARef.mapids.txt")) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) == 2:
            cla_to_acc[parts[0]] = parts[1]
            acc_to_cla[parts[1]] = parts[0]

def chrom_sort_key(name):
    try:
        return int(name.replace("Chr", ""))
    except Exception:
        return 999

# ── BIA BLAST annotation ──────────────────────────────────────────────────────
blast = {}   # transcript_id -> best NR protein accession
with open(os.path.join(ANN, "BIARef.Proteomes.blast_nr.Max5.diamond.blastp")) as f:
    for line in f:
        parts = line.strip().split()
        if parts and parts[0] not in blast:
            blast[parts[0]] = parts[1]

# ── Load SVs ──────────────────────────────────────────────────────────────────
records = []

for category, fname in [
    ("Private", "BIA.PrivateSV.AllSV_WithCDSInfo.xlsx"),
    ("Shared",  "BIA.SharedSV.AllSV_WithCDSInfo.xlsx"),
]:
    print(f"  Reading {fname} ...")
    df = pd.read_excel(os.path.join(BASE, fname), sheet_name="Sheet1", dtype=str)
    # Fill NaN with empty string for string columns; keep booleans
    df = df.fillna("")
    print(f"    {len(df)} rows loaded")

    for _, row in df.iterrows():
        acc = str(row.get("Chrom", "")).strip()
        if not acc:
            continue

        chr_name = acc_to_chr.get(acc, acc)
        cla_name = acc_to_cla.get(acc, "")

        # SV length (read as str due to dtype=str above)
        raw_len = row.get("SVLength") or row.get("SV_Mean_Length") or "0"
        try:
            sv_len = round(float(raw_len))
        except Exception:
            sv_len = 0

        # Species
        sp_str  = str(row.get("Species", "")).strip()
        species = [s.strip() for s in sp_str.split(",") if s.strip()]

        raw_nsp = str(row.get("Number_species", "")).strip()
        try:
            n_sp = int(float(raw_nsp))
        except Exception:
            n_sp = len(species) if species else 1

        # CDS / gene overlap flags (stored as "True"/"False" strings after dtype=str)
        in_coding = str(row.get("in_coding", "")).strip().upper() == "TRUE"
        in_gene   = str(row.get("in_gene",   "")).strip().upper() == "TRUE"

        # Transcripts (NA when not in CDS)
        cds_str = str(row.get("Affected_CDS", "")).strip()
        if cds_str.upper() == "NA":
            cds_str = ""
        transcripts = [t.strip() for t in cds_str.split(",") if t.strip()]

        # BLAST protein hits (unique, best hit per transcript)
        proteins = []
        seen_p   = set()
        for tid in transcripts:
            p = blast.get(tid, "")
            if p and p not in seen_p:
                proteins.append(p)
                seen_p.add(p)

        records.append({
            "cat":  category,       # Private | Shared
            "chr":  chr_name,       # Chr01 … Chr24
            "acc":  acc,            # CM132547.1
            "cla":  cla_name,       # chr1 (may be empty)
            "pos":  int(float(row.get("Position") or 0)),
            "id":   str(row.get("ID", "")).strip(),
            "type": str(row.get("SVType", "")).strip(),
            "len":  sv_len,
            "nsp":  n_sp,
            "sp":   species,
            "cds":  in_coding,      # True = overlaps a CDS
            "gene": in_gene,        # True = overlaps a gene (exon or intron)
            "nt":   len(transcripts),
            "tx":   transcripts[:30],   # first 30 for detail panel
            "pr":   proteins[:8],       # top 8 unique protein hits
        })

print(f"  Total records: {len(records)}")

# ── Derived lists for filter dropdowns ───────────────────────────────────────
all_chroms  = sorted(set(r["chr"] for r in records), key=chrom_sort_key)
all_species = sorted(set(s for r in records for s in r["sp"]))

# ── Embed into HTML ───────────────────────────────────────────────────────────
data_js      = json.dumps(records,      separators=(",", ":"), ensure_ascii=False)
chroms_js    = json.dumps(all_chroms,   separators=(",", ":"))
species_js   = json.dumps(all_species,  separators=(",", ":"))

# The HTML template uses no Python f-strings; placeholders are replaced below.
HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Clownfish SV Database</title>
<link rel="stylesheet" href="https://cdn.datatables.net/1.13.8/css/dataTables.bootstrap5.min.css">
<link rel="stylesheet" href="https://cdn.datatables.net/buttons/2.4.2/css/buttons.bootstrap5.min.css">
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/css/bootstrap.min.css">
<style>
body { font-family: "Segoe UI", system-ui, sans-serif; background:#f8f9fa; }
.navbar { background: linear-gradient(135deg,#1a6e3c,#2a7a4a); }
.filter-card { background:#fff; border-radius:8px; box-shadow:0 1px 6px rgba(0,0,0,.08); padding:1rem 1.2rem; margin-bottom:1rem; }
.filter-card h6 { font-size:.72rem; text-transform:uppercase; letter-spacing:.07em; color:#666; margin-bottom:.4rem; }
.stat-bar { display:flex; gap:.6rem; flex-wrap:wrap; margin-bottom:.8rem; }
.stat-chip { background:#fff; border-radius:6px; padding:.3rem .8rem; font-size:.82rem; box-shadow:0 1px 4px rgba(0,0,0,.09); }
.stat-chip b { color:#333; }
/* SV type badges */
.sv-badge { display:inline-block; font-size:.72rem; font-weight:600; padding:.2em .55em; border-radius:4px; }
.sv-DEL   { background:#7B5EA7; color:#fff; }
.sv-INS   { background:#D87B7B; color:#fff; }
.sv-DUP   { background:#56BBF9; color:#1a4a6a; }
.sv-TRANS { background:#B7D86E; color:#3a4a1a; }
.sv-INV   { background:#F09837; color:#fff; }
/* SV type filter buttons */
.btn-sv-DEL   { border-color:#7B5EA7!important; color:#7B5EA7!important; }
.btn-sv-DEL:hover, .btn-sv-DEL.active   { background:#7B5EA7!important; color:#fff!important; }
.btn-sv-INS   { border-color:#D87B7B!important; color:#D87B7B!important; }
.btn-sv-INS:hover, .btn-sv-INS.active   { background:#D87B7B!important; color:#fff!important; }
.btn-sv-DUP   { border-color:#56BBF9!important; color:#1a4a6a!important; }
.btn-sv-DUP:hover, .btn-sv-DUP.active   { background:#56BBF9!important; color:#1a4a6a!important; }
.btn-sv-TRANS { border-color:#B7D86E!important; color:#3a4a1a!important; }
.btn-sv-TRANS:hover, .btn-sv-TRANS.active { background:#B7D86E!important; color:#3a4a1a!important; }
.btn-sv-INV   { border-color:#F09837!important; color:#F09837!important; }
.btn-sv-INV:hover, .btn-sv-INV.active   { background:#F09837!important; color:#fff!important; }
/* Category */
.cat-Private { background:#f0f0f0; color:#555; font-size:.7rem; padding:.15em .5em; border-radius:3px; font-weight:600; }
.cat-Shared  { background:#fff3e8; color:#b7580a; font-size:.7rem; padding:.15em .5em; border-radius:3px; font-weight:600; }
/* Misc */
.chrom-tag { font-family:monospace; font-size:.82rem; background:#eef; padding:.1em .4em; border-radius:3px; }
.sp-pill   { display:inline-block; background:#e8f5e9; color:#2e7d32; border-radius:8px; padding:.1em .5em; font-size:.7rem; margin:1px; }
.prot-tag  { font-family:monospace; font-size:.75rem; background:#f5f5f5; padding:.1em .4em; border-radius:3px; margin:1px; display:inline-block; }
.ov-cds    { background:#fde8e6; color:#c0392b; font-size:.7rem; font-weight:600; padding:.15em .5em; border-radius:3px; }
.ov-gene   { background:#fef3e8; color:#b7580a; font-size:.7rem; font-weight:600; padding:.15em .5em; border-radius:3px; }
.ov-inter  { background:#f0f0f0; color:#888;    font-size:.7rem; font-weight:600; padding:.15em .5em; border-radius:3px; }
/* Table */
.table-wrapper { background:#fff; border-radius:8px; box-shadow:0 1px 6px rgba(0,0,0,.08); padding:.8rem; }
#svTable tbody tr { cursor:pointer; }
#svTable tbody tr:hover { background:#f0f7ff !important; }
/* Modal detail sections */
.detail-section h6 { font-size:.72rem; text-transform:uppercase; letter-spacing:.07em; color:#888; margin-bottom:.3rem; }
.detail-section { margin-bottom:.9rem; }
.tx-box { font-family:monospace; font-size:.73rem; background:#fafafa; border:1px solid #e0e0e0; border-radius:4px; padding:.4rem; max-height:140px; overflow-y:auto; color:#444; line-height:1.6; }
</style>
</head>
<body>

<nav class="navbar navbar-dark mb-3">
  <div class="container-fluid">
    <span class="navbar-brand fw-bold">Clownfish SV Database
      <small class="fw-light ms-2" style="font-size:.75rem;opacity:.8;">
        <em>Premnas biaculeatus</em> reference &middot; GCA_053813585.1
      </small>
    </span>
  </div>
</nav>

<div class="container-fluid px-4">

  <!-- Stat bar -->
  <div class="stat-bar" id="statBar">
    <div class="stat-chip">Showing <b id="nShown">–</b> / <b id="nTotal">–</b></div>
    <div class="stat-chip">DEL <b id="n_DEL">–</b></div>
    <div class="stat-chip">INS <b id="n_INS">–</b></div>
    <div class="stat-chip">DUP <b id="n_DUP">–</b></div>
    <div class="stat-chip">TRANS <b id="n_TRANS">–</b></div>
    <div class="stat-chip">INV <b id="n_INV">–</b></div>
  </div>

  <!-- Filters -->
  <div class="filter-card">
    <div class="row g-3 align-items-end">

      <div class="col-auto">
        <h6>Category</h6>
        <div class="btn-group btn-group-sm" id="catGroup">
          <button class="btn btn-outline-secondary active" data-val="ALL"     >Both</button>
          <button class="btn btn-outline-secondary"        data-val="Private" >Private</button>
          <button class="btn btn-outline-warning"          data-val="Shared"  >Shared</button>
        </div>
      </div>

      <div class="col-auto">
        <h6>SV type</h6>
        <div class="btn-group btn-group-sm" id="typeGroup">
          <button class="btn btn-outline-secondary active" data-val="ALL"  >All</button>
          <button class="btn btn-outline-secondary btn-sv-DEL"   data-val="DEL"  >DEL</button>
          <button class="btn btn-outline-secondary btn-sv-INS"   data-val="INS"  >INS</button>
          <button class="btn btn-outline-secondary btn-sv-DUP"   data-val="DUP"  >DUP</button>
          <button class="btn btn-outline-secondary btn-sv-TRANS" data-val="TRANS">TRANS</button>
          <button class="btn btn-outline-secondary btn-sv-INV"   data-val="INV"  >INV</button>
        </div>
      </div>

      <div class="col-auto">
        <h6>Chromosome</h6>
        <select id="selChrom" class="form-select form-select-sm" style="min-width:120px">
          <option value="ALL">All</option>
        </select>
      </div>

      <div class="col-auto">
        <h6>Species</h6>
        <select id="selSpecies" class="form-select form-select-sm" style="min-width:100px">
          <option value="ALL">All</option>
        </select>
      </div>

      <div class="col-auto">
        <h6>Min length (bp)</h6>
        <input type="number" id="minLen" class="form-control form-control-sm" style="width:100px" placeholder="0" min="0">
      </div>
      <div class="col-auto">
        <h6>Max length (bp)</h6>
        <input type="number" id="maxLen" class="form-control form-control-sm" style="width:100px" placeholder="∞" min="0">
      </div>
      <div class="col-auto">
        <h6>Min # species</h6>
        <input type="number" id="minSp" class="form-control form-control-sm" style="width:80px" placeholder="1" min="1">
      </div>
      <div class="col-auto">
        <h6>Genomic overlap</h6>
        <div class="btn-group btn-group-sm" id="overlapGroup">
          <button class="btn btn-outline-secondary active" data-val="ALL"    >All</button>
          <button class="btn btn-outline-danger"           data-val="cds"    >In CDS</button>
          <button class="btn btn-outline-warning"          data-val="gene"   >In gene</button>
          <button class="btn btn-outline-secondary"        data-val="intergenic">Intergenic</button>
        </div>
      </div>

      <div class="col-auto ms-auto">
        <button class="btn btn-sm btn-outline-secondary" id="btnReset">Reset</button>
      </div>
    </div>
  </div>

  <!-- Table -->
  <div class="table-wrapper">
    <p class="text-muted mb-2" style="font-size:.8rem">Click any row to see full details.</p>
    <table id="svTable" class="table table-sm table-hover w-100">
      <thead class="table-light">
        <tr>
          <th>Chr</th>
          <th>Chr (CLA)</th>
          <th>Position</th>
          <th>SV ID</th>
          <th>Type</th>
          <th>Length</th>
          <th>Category</th>
          <th>Overlap</th>
          <th>Species</th>
          <th># CDS</th>
          <th>Top BLAST hits</th>
        </tr>
      </thead>
      <tbody></tbody>
    </table>
  </div>

</div><!-- /container -->

<!-- Detail modal -->
<div class="modal fade" id="detailModal" tabindex="-1">
  <div class="modal-dialog modal-lg modal-dialog-scrollable">
    <div class="modal-content">
      <div class="modal-header py-2">
        <h5 class="modal-title" id="modalTitle" style="font-size:1rem"></h5>
        <button type="button" class="btn-close" data-bs-dismiss="modal"></button>
      </div>
      <div class="modal-body" id="modalBody"></div>
    </div>
  </div>
</div>

<script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/js/bootstrap.bundle.min.js"></script>
<script src="https://cdn.datatables.net/1.13.8/js/jquery.dataTables.min.js"></script>
<script src="https://cdn.datatables.net/1.13.8/js/dataTables.bootstrap5.min.js"></script>
<script src="https://cdn.datatables.net/buttons/2.4.2/js/dataTables.buttons.min.js"></script>
<script src="https://cdn.datatables.net/buttons/2.4.2/js/buttons.bootstrap5.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/jszip/3.10.1/jszip.min.js"></script>
<script src="https://cdn.datatables.net/buttons/2.4.2/js/buttons.html5.min.js"></script>

<script>
// ── Data ─────────────────────────────────────────────────────────────────────
var RECORDS  = __RECORDS__;
var CHROMS   = __CHROMS__;
var SPECIES  = __SPECIES__;

// ── Filter state ──────────────────────────────────────────────────────────────
var fCat     = "ALL";
var fType    = "ALL";
var fChr     = "ALL";
var fSp      = "ALL";
var fMin     = null;
var fMax     = null;
var fMinSp   = null;
var fOverlap = "ALL";

// ── Build dropdowns ───────────────────────────────────────────────────────────
(function() {
  var sc = document.getElementById("selChrom");
  CHROMS.forEach(function(c) {
    var o = document.createElement("option"); o.value = c; o.textContent = c; sc.appendChild(o);
  });
  var ss = document.getElementById("selSpecies");
  SPECIES.forEach(function(s) {
    var o = document.createElement("option"); o.value = s; o.textContent = s; ss.appendChild(o);
  });
})();

// ── Button group wiring ───────────────────────────────────────────────────────
document.getElementById("catGroup").addEventListener("click", function(e) {
  var btn = e.target.closest("button");
  if (!btn) return;
  fCat = btn.dataset.val;
  this.querySelectorAll("button").forEach(function(b) { b.classList.toggle("active", b.dataset.val === fCat); });
  redraw();
});
document.getElementById("typeGroup").addEventListener("click", function(e) {
  var btn = e.target.closest("button");
  if (!btn) return;
  fType = btn.dataset.val;
  this.querySelectorAll("button").forEach(function(b) { b.classList.toggle("active", b.dataset.val === fType); });
  redraw();
});
document.getElementById("selChrom").addEventListener("change",   function() { fChr = this.value; redraw(); });
document.getElementById("selSpecies").addEventListener("change", function() { fSp  = this.value; redraw(); });
var dRedraw = debounce(redraw, 300);
document.getElementById("minLen").addEventListener("input",  function() { fMin   = this.value !== "" ? parseInt(this.value) : null; dRedraw(); });
document.getElementById("maxLen").addEventListener("input",  function() { fMax   = this.value !== "" ? parseInt(this.value) : null; dRedraw(); });
document.getElementById("minSp").addEventListener("input",   function() { fMinSp = this.value !== "" ? parseInt(this.value) : null; dRedraw(); });
document.getElementById("overlapGroup").addEventListener("click", function(e) {
  var btn = e.target.closest("button");
  if (!btn) return;
  fOverlap = btn.dataset.val;
  this.querySelectorAll("button").forEach(function(b) { b.classList.toggle("active", b.dataset.val === fOverlap); });
  redraw();
});
document.getElementById("btnReset").addEventListener("click", function() {
  fCat = fType = fChr = fSp = "ALL"; fMin = fMax = fMinSp = null; fOverlap = "ALL";
  document.getElementById("selChrom").value   = "ALL";
  document.getElementById("selSpecies").value = "ALL";
  document.getElementById("minLen").value = "";
  document.getElementById("maxLen").value = "";
  document.getElementById("minSp").value  = "";
  ["catGroup","typeGroup","overlapGroup"].forEach(function(id) {
    document.getElementById(id).querySelectorAll("button").forEach(function(b) { b.classList.toggle("active", b.dataset.val === "ALL"); });
  });
  redraw();
});

// ── Filter logic ──────────────────────────────────────────────────────────────
function matchRec(r) {
  if (fCat   !== "ALL" && r.cat  !== fCat)              return false;
  if (fType  !== "ALL" && r.type !== fType)             return false;
  if (fChr   !== "ALL" && r.chr  !== fChr)              return false;
  if (fSp    !== "ALL" && r.sp.indexOf(fSp) === -1)     return false;
  if (fMin   !== null  && r.len  <  fMin)               return false;
  if (fMax   !== null  && r.len  >  fMax)               return false;
  if (fMinSp !== null  && r.nsp  <  fMinSp)             return false;
  if (fOverlap === "cds"        && !r.cds)              return false;
  if (fOverlap === "gene"       && !r.gene)             return false;
  if (fOverlap === "intergenic" && r.gene)              return false;
  return true;
}

// ── Helpers ───────────────────────────────────────────────────────────────────
function fmtLen(n) {
  if (n >= 1000000) return (n/1000000).toFixed(2) + " Mb";
  if (n >= 1000)    return (n/1000).toFixed(1)    + " kb";
  return n + " bp";
}
function spPills(arr) {
  if (!arr || arr.length === 0) return "–";
  var show = arr.slice(0,4).map(function(s) { return '<span class="sp-pill">'+s+'</span>'; }).join("");
  if (arr.length > 4) show += '<span class="sp-pill">+' + (arr.length-4) + '</span>';
  return show;
}
function protTags(arr) {
  if (!arr || arr.length === 0) return '<span class="text-muted" style="font-size:.75rem">–</span>';
  return arr.slice(0,3).map(function(p) { return '<span class="prot-tag">'+p+'</span>'; }).join(" ") +
    (arr.length > 3 ? ' <span class="text-muted" style="font-size:.72rem">+' + (arr.length-3) + '</span>' : "");
}
function ovBadge(r) {
  if (r.cds)  return '<span class="ov-cds">CDS</span>';
  if (r.gene) return '<span class="ov-gene">Intronic</span>';
  return '<span class="ov-inter">Intergenic</span>';
}

// ── Debounce helper ───────────────────────────────────────────────────────────
function debounce(fn, ms) {
  var t;
  return function() { clearTimeout(t); t = setTimeout(fn, ms); };
}

// ── DataTable — initialised once, filtered via ext.search ────────────────────
var table = null;

// Custom search function: called by DataTables for each row, returns true = keep
$.fn.dataTable.ext.search.push(function(settings, data, dataIndex) {
  if (settings.nTable.id !== "svTable") return true;
  var r = RECORDS[dataIndex];
  return matchRec(r);
});

function initTable() {
  table = $("#svTable").DataTable({
    data: RECORDS,
    deferRender: true,          // only render rows actually visible
    columns: [
      { title:"Chr",         data:"chr",  orderable:true,
        render: function(d)    { return '<span class="chrom-tag">'+d+'</span>'; } },
      { title:"Chr (CLA)",   data:"cla",  orderable:true,
        render: function(d)    { return d ? '<span class="chrom-tag" style="background:#f0f0ff">'+d+'</span>' : '–'; } },
      { title:"Position",    data:"pos",  orderable:true,
        render: function(d)    { return d.toLocaleString(); } },
      { title:"SV ID",       data:"id",   orderable:true,
        render: function(d)    { return '<code style="font-size:.75rem">'+d+'</code>'; } },
      { title:"Type",        data:"type", orderable:true,
        render: function(d)    { return '<span class="sv-badge sv-'+d+'">'+d+'</span>'; } },
      { title:"Length",      data:"len",  orderable:true,
        render: function(d, type)  { return type !== "display" ? d : fmtLen(d); } },
      { title:"Category",    data:"cat",  orderable:true,
        render: function(d)    { return '<span class="cat-'+d+'">'+d+'</span>'; } },
      { title:"Overlap",     data:null,   orderable:true,
        render: function(d, t, r) { return ovBadge(r); } },
      { title:"Species",     data:"sp",   orderable:false,
        render: function(d, type)    {
          if (type !== "display") return d ? d.join(", ") : "–";
          return spPills(d);
        } },
      { title:"# CDS",       data:null,   orderable:true,
        render: function(d, t, r) { return r.cds ? '<strong>'+r.nt+'</strong>' : '–'; } },
      { title:"Top BLAST",   data:"pr",   orderable:false,
        render: function(d, type)    {
          if (type !== "display") return d ? d.join(", ") : "–";
          return protTags(d);
        } }
    ],
    createdRow: function(row, data, dataIndex) {
      $(row).attr("data-ridx", dataIndex);
    },
    pageLength: 25,
    lengthMenu: [[10,25,50,100,250],[10,25,50,100,250]],
    dom: '<"row mb-2"<"col-sm-4"l><"col-sm-4 text-center"B><"col-sm-4"f>>rtip',
    buttons: [{
      extend:      "csvHtml5",
      text:        "⬇ Export CSV",
      className:   "btn btn-sm btn-outline-secondary",
      filename:    "BIA_SV_database",
      exportOptions: { columns:[0,1,2,3,4,5,6,7,8,9,10], orthogonal:"export" }
    }],
    language: {
      emptyTable:  "No SVs match the current filters.",
      zeroRecords: "No SVs match the search."
    }
  });

  // Row click → modal (dataIndex is the index into RECORDS)
  $("#svTable tbody").on("click", "tr", function() {
    var idx = parseInt($(this).attr("data-ridx"));
    if (isNaN(idx)) return;
    openDetail(RECORDS[idx]);
  });
}

// ── updateStats: count across all currently-passing rows ─────────────────────
function updateStats() {
  var counts = {DEL:0, INS:0, DUP:0, TRANS:0, INV:0};
  var shown = 0;
  RECORDS.forEach(function(r) {
    if (!matchRec(r)) return;
    shown++;
    if (counts[r.type] !== undefined) counts[r.type]++;
  });
  document.getElementById("nShown").textContent = shown.toLocaleString();
  document.getElementById("nTotal").textContent = RECORDS.length.toLocaleString();
  ["DEL","INS","DUP","TRANS","INV"].forEach(function(t) {
    document.getElementById("n_"+t).textContent = counts[t].toLocaleString();
  });
}

// ── redraw: update stats then let DataTables re-filter in place ───────────────
function redraw() {
  updateStats();
  if (table) {
    table.draw();
  }
}

// ── Detail modal ──────────────────────────────────────────────────────────────
var bsModal = null;

function openDetail(rec) {
  if (!rec) return;

  document.getElementById("modalTitle").innerHTML =
    '<code>' + rec.id + '</code>' +
    '&ensp;<span class="sv-badge sv-' + rec.type + '">' + rec.type + '</span>' +
    '&ensp;<span class="chrom-tag">' + rec.chr + '</span>' +
    '&ensp;pos&thinsp;' + rec.pos.toLocaleString() + ' bp';

  var crossRef =
    '<div class="detail-section">' +
      '<h6>Chromosome</h6>' +
      '<span class="chrom-tag me-2">Chr: ' + rec.acc + '</span>' +
      (rec.cla ? '<span class="chrom-tag me-2" style="background:#f0f0ff">CLA: ' + rec.cla + '</span>' : '') +
      '<span class="chrom-tag" style="background:#eee">BIA: ' + rec.chr + '</span>' +
    '</div>' +
    '<div class="detail-section">' +
      '<h6>Genomic overlap</h6>' +
      ovBadge(rec) +
    '</div>';

  var spBlock =
    '<div class="detail-section">' +
      '<h6>Species (' + rec.nsp + ')</h6>' +
      rec.sp.map(function(s) { return '<span class="sp-pill">' + s + '</span>'; }).join(" ") +
    '</div>';

  var protBlock = "";
  if (rec.pr && rec.pr.length > 0) {
    protBlock =
      '<div class="detail-section">' +
        '<h6>BLAST hits – NR protein accessions (best hit per unique transcript, up to 10 shown)</h6>' +
        rec.pr.map(function(p) { return '<span class="prot-tag">' + p + '</span>'; }).join(" ") +
      '</div>';
  }

  var txBlock = "";
  if (rec.tx && rec.tx.length > 0) {
    var note = rec.nt > rec.tx.length ? " (showing first " + rec.tx.length + " of " + rec.nt + ")" : "";
    txBlock =
      '<div class="detail-section">' +
        '<h6>Affected transcripts' + note + '</h6>' +
        '<div class="tx-box">' + rec.tx.join(", ") + '</div>' +
      '</div>';
  }

  document.getElementById("modalBody").innerHTML =
    '<div class="row">' +
      '<div class="col-md-5">' + crossRef + spBlock + protBlock + '</div>' +
      '<div class="col-md-7">' + txBlock + '</div>' +
    '</div>';

  if (!bsModal) bsModal = new bootstrap.Modal(document.getElementById("detailModal"));
  bsModal.show();
}

// ── Start ─────────────────────────────────────────────────────────────────────
$(document).ready(function() {
  initTable();   // load all records once; ext.search handles filtering
  updateStats(); // populate stat bar for the unfiltered state
});
</script>
</body>
</html>
"""

# Inject data (replace placeholders — no f-string escaping needed)
html = HTML_TEMPLATE \
    .replace("__RECORDS__",  data_js) \
    .replace("__CHROMS__",   chroms_js) \
    .replace("__SPECIES__",  species_js)

out_path = os.path.join(BASE, "sv_database.html")
with open(out_path, "w", encoding="utf-8") as fh:
    fh.write(html)

size_mb = os.path.getsize(out_path) / 1e6
print(f"  Saved to: {out_path}")
print(f"  File size: {size_mb:.1f} MB")
