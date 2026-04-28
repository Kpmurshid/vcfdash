/* =========================================================
   vcfdash — viz.js
   D3.js v7 sparklines + table sort/filter logic
   All data is read from window.VCFDASH injected by report.py
   ========================================================= */

(function () {
  "use strict";

  /* ── Config from embedded data ──────────────────────── */
  const DATA        = window.VCFDASH || {};
  const VARIANTS    = DATA.variants    || [];
  const COVERAGE    = DATA.coverage    || [];
  const SPARK_DATA  = DATA.sparklines  || {};   // {chrom: {pos: depth}}
  const QC          = DATA.qc          || {};
  const THRESHOLDS  = DATA.thresholds  || { min_dp: 20, min_dp2: 10 };
  const MIN_DP      = THRESHOLDS.min_dp;
  const MIN_DP2     = THRESHOLDS.min_dp2;
  const HAS_SPARKS  = DATA.has_sparklines !== false;

  /* ── Colour helpers ─────────────────────────────────── */
  function depthColor(d) {
    if (d >= MIN_DP)  return "#22c55e";
    if (d >= MIN_DP2) return "#f59e0b";
    return "#ef4444";
  }

  function filterPillClass(f) {
    if (!f || f === "PASS") return "pill pill-pass";
    return "pill pill-warn";
  }

  function clinvarClass(cv) {
    if (!cv) return "";
    const lc = cv.toLowerCase();
    if (lc.includes("pathogenic") && !lc.includes("benign")) return "pill pill-path";
    if (lc.includes("benign"))  return "pill pill-benign";
    if (lc.includes("vus") || lc.includes("uncertain")) return "pill pill-vus";
    return "";
  }

  function formatAF(af) {
    if (af === null || af === undefined || af === 0 || af === "") return "—";
    const n = parseFloat(af);
    if (isNaN(n) || n === 0) return "—";
    if (n < 0.0001) return n.toExponential(2);
    return n.toFixed(4);
  }

  function formatVAF(vaf) {
    if (vaf === null || vaf === undefined || vaf === "") return "N/A";
    return (parseFloat(vaf) * 100).toFixed(1) + "%";
  }

  function dpClass(dp) {
    if (dp >= MIN_DP)  return "dp-good";
    if (dp >= MIN_DP2) return "dp-warn";
    return "dp-fail";
  }

  function varTypeBadge(vt) {
    const map = { SNV: "badge-snv", INDEL: "badge-indel", MNV: "badge-mnv" };
    return `<span class="badge ${map[vt] || 'badge-snv'}">${vt || "SNV"}</span>`;
  }

  /* ── Generic table sorter ───────────────────────────── */
  function makeTableSortable(tableEl) {
    const headers = tableEl.querySelectorAll("thead th[data-col]");
    let sortCol   = null;
    let sortDir   = 1;

    headers.forEach(th => {
      // inject sort icon span
      const icon = document.createElement("i");
      icon.className = "sort-icon";
      th.appendChild(icon);

      th.addEventListener("click", () => {
        const col = th.dataset.col;
        const type = th.dataset.type || "str";

        if (sortCol === col) {
          sortDir = -sortDir;
        } else {
          sortCol = col;
          sortDir = 1;
        }

        headers.forEach(h => h.classList.remove("sort-asc", "sort-desc"));
        th.classList.add(sortDir === 1 ? "sort-asc" : "sort-desc");

        const tbody = tableEl.querySelector("tbody");
        const rows  = Array.from(tbody.querySelectorAll("tr:not(.detail-row)"));

        rows.sort((a, b) => {
          const av = a.dataset[col] ?? a.cells[parseInt(th.dataset.idx || 0)]?.textContent ?? "";
          const bv = b.dataset[col] ?? b.cells[parseInt(th.dataset.idx || 0)]?.textContent ?? "";

          if (type === "num") {
            return (parseFloat(av) - parseFloat(bv)) * sortDir;
          }
          return av.localeCompare(bv, undefined, { sensitivity: "base" }) * sortDir;
        });

        // Reinsert rows (keep detail rows paired)
        rows.forEach(row => {
          tbody.appendChild(row);
          const detailId = row.dataset.detailFor;
          if (detailId) {
            const detail = document.getElementById(detailId);
            if (detail) tbody.appendChild(detail);
          }
          const nextDetail = row.nextElementSibling;
          if (nextDetail && nextDetail.classList.contains("detail-row")) {
            tbody.appendChild(nextDetail);
          }
        });
      });
    });
  }

  /* ── Gene coverage table ────────────────────────────── */
  function buildCoverageTable() {
    const tableEl = document.getElementById("coverage-table");
    if (!tableEl) return;

    const tbody  = tableEl.querySelector("tbody");
    const search = document.getElementById("gene-search");
    const filter = document.getElementById("status-filter");
    const count  = document.getElementById("cov-count");

    // Find max depth for bar scaling
    const maxDepth = Math.max(...COVERAGE.map(r => r.mean_depth || 0), 1);

    function renderRows() {
      const term     = (search?.value || "").toLowerCase();
      const statusF  = filter?.value || "all";

      let visible = 0;
      tbody.querySelectorAll("tr").forEach(tr => {
        const gene   = (tr.dataset.gene   || "").toLowerCase();
        const status = (tr.dataset.status || "").toUpperCase();
        const matchG = !term   || gene.includes(term);
        const matchS = statusF === "all" || status === statusF;
        const show   = matchG && matchS;
        tr.classList.toggle("hidden", !show);
        if (show) visible++;
      });
      if (count) {
        // Determine the right label: chromosomes vs genes vs regions
        const label = _covLabel();
        count.textContent = `${visible} / ${COVERAGE.length} ${label}`;
      }
    }

    // Build rows
    COVERAGE.forEach(r => {
      const status    = r.status || "FAIL";
      const rowClass  = { PASS: "row-pass", WARN: "row-warn", FAIL: "row-fail" }[status] || "";
      const statusLabel = {
        PASS: "✅ PASS", WARN: "⚠️ WARN", FAIL: "❌ FAIL"
      }[status] || status;
      const statusClass = {
        PASS: "status-pass", WARN: "status-warn", FAIL: "status-fail"
      }[status] || "";

      const barPct = Math.min(100, (r.mean_depth / maxDepth) * 100).toFixed(1);
      const pct20  = r.pct_20x != null ? r.pct_20x.toFixed(1) + "%" : "—";
      const pct10  = r.pct_10x != null ? r.pct_10x.toFixed(1) + "%" : "—";

      const tr = document.createElement("tr");
      tr.className   = rowClass;
      tr.dataset.gene   = r.gene || "";
      tr.dataset.status = status;
      tr.dataset.depth  = r.mean_depth || 0;
      tr.dataset.pct20  = r.pct_20x || 0;
      tr.dataset.pct10  = r.pct_10x || 0;

      tr.innerHTML = `
        <td><a href="#gene-${encodeURIComponent(r.gene)}" class="gene-link">${escHtml(r.gene)}</a></td>
        <td class="depth-bar-cell">
          <span class="${dpClass(r.mean_depth)}">${(r.mean_depth || 0).toFixed(1)}x</span>
          <div class="depth-bar-wrap" style="margin-top:4px">
            <div class="depth-bar-fill" style="width:${barPct}%"></div>
          </div>
        </td>
        <td>${pct20}</td>
        <td>${pct10}</td>
        <td><span class="${statusClass}">${statusLabel}</span></td>
      `;

      // clicking gene scrolls to its variants
      tr.querySelector(".gene-link")?.addEventListener("click", e => {
        e.preventDefault();
        const anchor = document.getElementById("gene-" + encodeURIComponent(r.gene));
        if (anchor) anchor.scrollIntoView({ behavior: "smooth", block: "start" });
        // highlight gene in variant filter
        const gf = document.getElementById("var-gene-filter");
        if (gf) {
          gf.value = r.gene;
          gf.dispatchEvent(new Event("input"));
        }
      });

      tbody.appendChild(tr);
    });

    makeTableSortable(tableEl);
    search?.addEventListener("input", renderRows);
    filter?.addEventListener("change", renderRows);
    renderRows();
  }

  /* ── Variant table ──────────────────────────────────── */
  const PAGE_SIZE = 200;   // render in batches for performance

  function buildVariantTable() {
    const tableEl   = document.getElementById("variant-table");
    if (!tableEl) return;

    const tbody     = tableEl.querySelector("tbody");
    const geneF     = document.getElementById("var-gene-filter");
    const csqF      = document.getElementById("var-csq-filter");
    const cvF       = document.getElementById("var-cv-filter");
    const textF     = document.getElementById("var-text-filter");
    const count     = document.getElementById("var-count");
    const loadMore  = document.getElementById("load-more-btn");

    let filtered    = [...VARIANTS];
    let rendered    = 0;

    // Build rows lazily
    const rowCache  = new Map();   // index -> {main, detail}

    function getRow(idx) {
      if (rowCache.has(idx)) return rowCache.get(idx);
      const v    = VARIANTS[idx];
      const pair = buildVariantRow(v, idx);
      rowCache.set(idx, pair);
      return pair;
    }

    function applyFilters() {
      const geneVal = (geneF?.value  || "").toLowerCase().trim();
      const csqVal  = (csqF?.value   || "").toLowerCase().trim();
      const cvVal   = (cvF?.value    || "").toLowerCase().trim();
      const textVal = (textF?.value  || "").toLowerCase().trim();

      filtered = VARIANTS.filter(v => {
        if (geneVal && !(v.gene || "").toLowerCase().includes(geneVal)) return false;
        if (csqVal  && !(v.consequence || "").toLowerCase().includes(csqVal)) return false;
        if (cvVal) {
          const cv = (v.clinvar || "").toLowerCase();
          if (cvVal === "pathogenic" && !(cv.includes("pathogenic") && !cv.includes("benign"))) return false;
          if (cvVal === "benign"     && !cv.includes("benign")) return false;
          if (cvVal === "vus"        && !(cv.includes("uncertain") || cv.includes("vus"))) return false;
          if (cvVal === "pass"       && !(v.filter === "PASS" || !v.filter)) return false;
        }
        if (textVal) {
          const hay = [v.gene, v.consequence, v.hgvs_c, v.hgvs_p,
                       v.chrom, String(v.pos), v.ref, v.alt].join(" ").toLowerCase();
          if (!hay.includes(textVal)) return false;
        }
        return true;
      });

      // Clear tbody
      while (tbody.firstChild) tbody.removeChild(tbody.firstChild);
      rendered = 0;
      renderBatch();

      if (count) count.textContent = `${filtered.length} / ${VARIANTS.length} variants`;
    }

    function renderBatch() {
      const end = Math.min(rendered + PAGE_SIZE, filtered.length);
      const frag = document.createDocumentFragment();

      for (let i = rendered; i < end; i++) {
        const origIdx = VARIANTS.indexOf(filtered[i]);
        const { main, detail } = getRow(origIdx >= 0 ? origIdx : i);
        frag.appendChild(main);
        frag.appendChild(detail);
      }
      tbody.appendChild(frag);
      rendered = end;

      if (loadMore) {
        loadMore.style.display = rendered < filtered.length ? "inline-block" : "none";
        loadMore.textContent = `Load more (${filtered.length - rendered} remaining)`;
      }
    }

    loadMore?.addEventListener("click", renderBatch);
    geneF?.addEventListener("input", applyFilters);
    csqF?.addEventListener("input", applyFilters);
    cvF?.addEventListener("change", applyFilters);
    textF?.addEventListener("input", applyFilters);

    makeTableSortable(tableEl);
    applyFilters();
  }

  /* ── Build a single variant row pair ───────────────── */
  function buildVariantRow(v, idx) {
    const detailId = `detail-${idx}`;
    const filt     = v.filter || "PASS";
    const cvClass  = clinvarClass(v.clinvar || "");

    // Main row
    const main = document.createElement("tr");
    main.id = `gene-${encodeURIComponent(v.gene || "")}`;  // anchor for gene links
    main.className = "gene-anchor";
    main.dataset.detailFor = detailId;
    main.dataset.gene      = v.gene || "";
    main.dataset.csq       = v.consequence || "";
    main.dataset.depth     = v.dp || 0;
    main.dataset.vaf       = v.vaf ?? "";

    const afDisplay = formatAF(v.gnomad_af);
    const afClass   = v.gnomad_af && v.gnomad_af > 0.01 ? "af-common" : "af-rare";

    main.innerHTML = `
      <td>
        <button class="expand-btn" data-target="${detailId}" title="Expand details">▶</button>
      </td>
      <td><strong>${escHtml(v.gene || "—")}</strong>${varTypeBadge(v.variant_type)}</td>
      <td class="mono" title="${escHtml(v.consequence || "")}">${truncate(escHtml(v.consequence || "—"), 30)}</td>
      <td class="mono" title="${escHtml(v.hgvs_c || "")}">${truncate(escHtml(v.hgvs_c || "—"), 25)}</td>
      <td class="mono" title="${escHtml(v.hgvs_p || "")}">${truncate(escHtml(v.hgvs_p || "—"), 20)}</td>
      <td class="${afClass}">${afDisplay}</td>
      <td>${v.clinvar ? `<span class="${cvClass} pill">${escHtml(v.clinvar)}</span>` : "—"}</td>
      <td>${v.cadd_phred ? parseFloat(v.cadd_phred).toFixed(1) : "—"}</td>
      <td class="mono">${escHtml(v.gt || ".")}</td>
      <td class="${dpClass(v.dp || 0)}">${v.dp ?? "—"}</td>
      <td class="${dpClass(v.gq || 0)}">${v.gq ?? "—"}</td>
      <td>${formatVAF(v.vaf)}</td>
      <td><span class="${filterPillClass(filt)}">${escHtml(filt)}</span></td>
    `;

    // Expand button handler
    main.querySelector(".expand-btn").addEventListener("click", function () {
      this.classList.toggle("open");
      const detailRow = document.getElementById(detailId);
      if (detailRow) {
        detailRow.classList.toggle("open");
        // Draw sparkline on first open
        if (detailRow.classList.contains("open") && !detailRow.dataset.drawn) {
          drawSparkline(detailRow, v);
          detailRow.dataset.drawn = "1";
        }
      }
    });

    // Detail row
    const detail = document.createElement("tr");
    detail.id = detailId;
    detail.className = "detail-row";

    const adStr = v.ad ? v.ad.join(",") : "N/A";
    const locStr = `${v.chrom}:${v.pos} ${v.ref}→${v.alt}`;

    detail.innerHTML = `
      <td colspan="13" class="detail-cell">
        <div class="detail-grid">
          <div>
            <div class="detail-section-title">Variant</div>
            <div class="detail-kv">
              <div class="kv-row"><span class="k">Location</span><span class="v">${escHtml(locStr)}</span></div>
              <div class="kv-row"><span class="k">Gene</span><span class="v">${escHtml(v.gene || "—")}</span></div>
              <div class="kv-row"><span class="k">Consequence</span><span class="v">${escHtml(v.consequence || "—")}</span></div>
              <div class="kv-row"><span class="k">HGVSc</span><span class="v">${escHtml(v.hgvs_c || "—")}</span></div>
              <div class="kv-row"><span class="k">HGVSp</span><span class="v">${escHtml(v.hgvs_p || "—")}</span></div>
              <div class="kv-row"><span class="k">Type</span><span class="v">${escHtml(v.variant_type || "—")}</span></div>
              <div class="kv-row"><span class="k">FILTER</span><span class="v">${escHtml(v.filter || "PASS")}</span></div>
            </div>
            <div class="detail-section-title" style="margin-top:.75rem">Annotation</div>
            <div class="detail-kv">
              <div class="kv-row"><span class="k">gnomAD AF</span><span class="v">${formatAF(v.gnomad_af)}</span></div>
              <div class="kv-row"><span class="k">ClinVar</span><span class="v">${escHtml(v.clinvar || "—")}</span></div>
              <div class="kv-row"><span class="k">CADD phred</span><span class="v">${v.cadd_phred ? parseFloat(v.cadd_phred).toFixed(2) : "—"}</span></div>
            </div>
          </div>
          <div>
            <div class="detail-section-title">Genotype</div>
            <div class="detail-kv">
              <div class="kv-row"><span class="k">GT</span><span class="v">${escHtml(v.gt || ".")}</span></div>
              <div class="kv-row"><span class="k">DP</span><span class="v">${v.dp ?? "N/A"}</span></div>
              <div class="kv-row"><span class="k">GQ</span><span class="v">${v.gq ?? "N/A"}</span></div>
              <div class="kv-row"><span class="k">AD</span><span class="v">${escHtml(adStr)}</span></div>
              <div class="kv-row"><span class="k">VAF</span><span class="v">${formatVAF(v.vaf)}</span></div>
            </div>
            <div class="sparkline-wrap" id="spark-wrap-${idx}">
              ${HAS_SPARKS && v.ad ? '<div class="sparkline-title">Coverage (200 bp window)</div><div class="sparkline-container"></div>' : '<div style="color:var(--c-muted);font-size:.8rem;margin-top:.5rem">Coverage sparkline not available (no per-base data or AD field)</div>'}
            </div>
          </div>
        </div>
      </td>
    `;

    return { main, detail };
  }

  /* ── D3 Sparkline ───────────────────────────────────── */
  function drawSparkline(detailRow, v) {
    if (!HAS_SPARKS || !v.ad) return;
    if (typeof d3 === "undefined") return;

    const container = detailRow.querySelector(".sparkline-container");
    if (!container) return;

    // Extract window data
    const chrom   = v.chrom;
    const pos     = v.pos;
    const WINDOW  = 100;  // ±100bp = 200bp total
    const chromData = SPARK_DATA[chrom] || {};

    const points = [];
    for (let p = Math.max(1, pos - WINDOW); p <= pos + WINDOW; p++) {
      points.push({ pos: p, depth: chromData[p] || 0 });
    }

    if (points.every(p => p.depth === 0)) {
      container.innerHTML = '<div style="color:var(--c-muted);font-size:.8rem">No per-base data for this region</div>';
      return;
    }

    const W    = 440;
    const H    = 90;
    const mL   = 38;
    const mR   = 10;
    const mT   = 8;
    const mB   = 22;
    const iW   = W - mL - mR;
    const iH   = H - mT - mB;

    const maxD = Math.max(d3.max(points, d => d.depth), MIN_DP * 1.1);

    const xScale = d3.scaleBand()
      .domain(points.map(d => d.pos))
      .range([0, iW])
      .padding(0.05);

    const yScale = d3.scaleLinear()
      .domain([0, maxD])
      .range([iH, 0])
      .nice();

    const svg = d3.select(container)
      .append("svg")
      .attr("width", W)
      .attr("height", H)
      .attr("class", "sparkline-svg")
      .attr("role", "img")
      .attr("aria-label", `Coverage at ${chrom}:${pos}`);

    const g = svg.append("g")
      .attr("transform", `translate(${mL},${mT})`);

    // Bars
    g.selectAll(".spark-bar")
      .data(points)
      .join("rect")
      .attr("class", "spark-bar")
      .attr("x", d => xScale(d.pos))
      .attr("y", d => yScale(d.depth))
      .attr("width", xScale.bandwidth())
      .attr("height", d => iH - yScale(d.depth))
      .attr("fill", d => depthColor(d.depth))
      .attr("opacity", 0.85);

    // Reference line at MIN_DP
    if (yScale(MIN_DP) >= 0) {
      g.append("line")
        .attr("x1", 0).attr("x2", iW)
        .attr("y1", yScale(MIN_DP)).attr("y2", yScale(MIN_DP))
        .attr("stroke", "#6c8ef5")
        .attr("stroke-dasharray", "3,2")
        .attr("stroke-width", 1)
        .attr("opacity", 0.7);

      g.append("text")
        .attr("x", -4)
        .attr("y", yScale(MIN_DP) + 3)
        .attr("text-anchor", "end")
        .attr("fill", "#6c8ef5")
        .attr("font-size", "9px")
        .text(MIN_DP + "x");
    }

    // Variant position marker
    const varX = xScale(pos);
    if (varX !== undefined) {
      g.append("line")
        .attr("x1", varX + xScale.bandwidth() / 2)
        .attr("x2", varX + xScale.bandwidth() / 2)
        .attr("y1", 0).attr("y2", iH)
        .attr("stroke", "#f59e0b")
        .attr("stroke-width", 1.5)
        .attr("opacity", 0.9);
    }

    // Y axis
    g.append("g")
      .call(d3.axisLeft(yScale).ticks(4).tickSize(2))
      .call(ax => ax.select(".domain").remove())
      .call(ax => ax.selectAll("text")
        .attr("fill", "#8892aa")
        .attr("font-size", "9px"))
      .call(ax => ax.selectAll("line").attr("stroke", "#2e3250"));

    // X axis labels (just start, variant, end)
    const xTickPositions = [points[0].pos, pos, points[points.length - 1].pos];
    g.append("g")
      .attr("transform", `translate(0,${iH})`)
      .selectAll("text")
      .data(xTickPositions)
      .join("text")
      .attr("x", d => (xScale(d) || 0) + xScale.bandwidth() / 2)
      .attr("y", 12)
      .attr("text-anchor", "middle")
      .attr("fill", "#8892aa")
      .attr("font-size", "9px")
      .text(d => d.toLocaleString());

    // Tooltip (shared)
    let tooltip = document.querySelector(".spark-tooltip");
    if (!tooltip) {
      tooltip = document.createElement("div");
      tooltip.className = "spark-tooltip";
      document.body.appendChild(tooltip);
    }

    g.selectAll(".spark-bar")
      .on("mousemove", function (event, d) {
        tooltip.style.display = "block";
        tooltip.style.left = (event.clientX + 12) + "px";
        tooltip.style.top  = (event.clientY - 28) + "px";
        tooltip.textContent = `pos ${d.pos.toLocaleString()}  depth ${d.depth}x`;
      })
      .on("mouseleave", () => { tooltip.style.display = "none"; });
  }

  /* ── Coverage label helper ──────────────────────────── */
  // Decide whether the coverage rows represent chromosomes, genes, or generic regions.
  function _covLabel() {
    if (!COVERAGE.length) return "regions";
    // If ≥80% of names look like standard chromosome names (chr1…chrY, 1…22, X, Y, MT)
    const chromRe = /^(chr)?(\d{1,2}|X|Y|MT?|Un)(_\w+)?$/i;
    const chromCount = COVERAGE.filter(r => chromRe.test(r.gene || "")).length;
    if (chromCount / COVERAGE.length >= 0.8) return "chromosomes";
    // If names look like coord-style "chr1:100-200"
    const coordCount = COVERAGE.filter(r => (r.gene || "").includes(":") && (r.gene || "").includes("-")).length;
    if (coordCount / COVERAGE.length >= 0.8) return "regions";
    return "genes";
  }

  /* ── Utility ────────────────────────────────────────── */
  function escHtml(s) {
    return String(s)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;");
  }

  function truncate(s, n) {
    if (!s) return s;
    return s.length > n ? s.slice(0, n) + "…" : s;
  }

  /* ── Initialise ─────────────────────────────────────── */
  document.addEventListener("DOMContentLoaded", () => {
    buildCoverageTable();
    buildVariantTable();
  });

})();
