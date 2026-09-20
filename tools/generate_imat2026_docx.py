from pathlib import Path
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from docx import Document
from docx.shared import Cm, Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT, WD_CELL_VERTICAL_ALIGNMENT
from docx.enum.section import WD_SECTION
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.enum.style import WD_STYLE_TYPE

OUTDIR = Path("conference/IMAT2026")
OUTDIR.mkdir(parents=True, exist_ok=True)
DOCX = OUTDIR / "IMAT2026_RMT_Extended_Abstract.docx"
FIG = OUTDIR / "rmt_mesh_dependence.png"

# ---------------------------------------------------------------------------
# Figure from the current RMT manuscript results
# ---------------------------------------------------------------------------
mesh_labels = ["108×36", "216×72", "432×144", "864×288"]
cells = [3888, 15552, 62208, 248832]
rmt = [0.293, 0.411, 0.393, 0.451]
v2 = [0.506, 0.509, 0.507, 0.504]
w2 = [0.509, 0.508, 0.504, 0.502]
v3 = [0.680, 0.613, 0.581, 0.561]

fig, ax = plt.subplots(figsize=(6.4, 3.0))
x = list(range(len(mesh_labels)))
ax.plot(x, rmt, marker="o", linewidth=1.5, label="RMT")
ax.plot(x, v2, marker="s", linewidth=1.2, label="V(2)-MG")
ax.plot(x, w2, marker="^", linewidth=1.2, label="W(2)-MG")
ax.plot(x, v3, marker="d", linewidth=1.2, label="V(3)-MG")
ax.set_xticks(x)
ax.set_xticklabels(mesh_labels)
ax.set_xlabel("Mesh")
ax.set_ylabel("Mean intergrid iterations per pressure solve")
ax.grid(True, alpha=0.25)
ax.legend(ncol=2, fontsize=8, frameon=False)
fig.tight_layout()
fig.savefig(FIG, dpi=320, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
doc = Document()
sec = doc.sections[0]
sec.page_width = Cm(21.0)
sec.page_height = Cm(29.7)
sec.top_margin = Cm(2.0)
sec.bottom_margin = Cm(2.0)
sec.left_margin = Cm(2.0)
sec.right_margin = Cm(2.0)

styles = doc.styles

normal = styles["Normal"]
normal.font.name = "Times New Roman"
normal._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
normal.font.size = Pt(10)
normal.paragraph_format.space_after = Pt(4)
normal.paragraph_format.line_spacing = 1.0

for sty_name in ["Title", "Subtitle"]:
    if sty_name in styles:
        s = styles[sty_name]
        s.font.name = "Times New Roman"
        s._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")

if "IMAT Section" not in styles:
    s = styles.add_style("IMAT Section", WD_STYLE_TYPE.PARAGRAPH)
else:
    s = styles["IMAT Section"]
s.font.name = "Times New Roman"
s._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
s.font.size = Pt(10)
s.font.bold = True
s.paragraph_format.space_before = Pt(7)
s.paragraph_format.space_after = Pt(3)
s.paragraph_format.keep_with_next = True

if "IMAT Caption" not in styles:
    s = styles.add_style("IMAT Caption", WD_STYLE_TYPE.PARAGRAPH)
else:
    s = styles["IMAT Caption"]
s.font.name = "Times New Roman"
s._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
s.font.size = Pt(9)
s.font.italic = True
s.paragraph_format.space_before = Pt(2)
s.paragraph_format.space_after = Pt(5)

def set_run_font(run, size=10, bold=None, italic=None):
    run.font.name = "Times New Roman"
    run._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
    run.font.size = Pt(size)
    if bold is not None:
        run.bold = bold
    if italic is not None:
        run.italic = italic

def add_p(text="", *, align=None, bold_lead=None, italic=False, after=4):
    p = doc.add_paragraph()
    if align is not None:
        p.alignment = align
    p.paragraph_format.space_after = Pt(after)
    if bold_lead and text.startswith(bold_lead):
        r1 = p.add_run(bold_lead)
        set_run_font(r1, bold=True)
        r2 = p.add_run(text[len(bold_lead):])
        set_run_font(r2)
    else:
        r = p.add_run(text)
        set_run_font(r, italic=italic)
    return p

def add_section(title):
    p = doc.add_paragraph(style="IMAT Section")
    r = p.add_run(title)
    set_run_font(r, bold=True)
    return p

def shade_cell(cell, fill="E7E6E6"):
    tcPr = cell._tc.get_or_add_tcPr()
    shd = tcPr.find(qn("w:shd"))
    if shd is None:
        shd = OxmlElement("w:shd")
        tcPr.append(shd)
    shd.set(qn("w:fill"), fill)

def set_cell_text(cell, text, bold=False, align=WD_ALIGN_PARAGRAPH.CENTER, size=9):
    cell.text = ""
    p = cell.paragraphs[0]
    p.alignment = align
    p.paragraph_format.space_after = Pt(0)
    r = p.add_run(str(text))
    set_run_font(r, size=size, bold=bold)
    cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER

def set_table_borders(table, color="808080", sz="4"):
    tbl = table._tbl
    tblPr = tbl.tblPr
    borders = tblPr.first_child_found_in("w:tblBorders")
    if borders is None:
        borders = OxmlElement("w:tblBorders")
        tblPr.append(borders)
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        tag = "w:" + edge
        el = borders.find(qn(tag))
        if el is None:
            el = OxmlElement(tag)
            borders.append(el)
        el.set(qn("w:val"), "single")
        el.set(qn("w:sz"), sz)
        el.set(qn("w:space"), "0")
        el.set(qn("w:color"), color)

def add_table(rows, headers, widths=None, font_size=9):
    table = doc.add_table(rows=1, cols=len(headers))
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = True
    hdr = table.rows[0].cells
    for i, h in enumerate(headers):
        set_cell_text(hdr[i], h, bold=True, size=font_size)
        shade_cell(hdr[i], "E7E6E6")
    for row in rows:
        cells_row = table.add_row().cells
        for i, val in enumerate(row):
            set_cell_text(cells_row[i], val, size=font_size)
    set_table_borders(table)
    table.rows[0]._tr.get_or_add_trPr().append(OxmlElement("w:tblHeader"))
    return table

# ---------------------------------------------------------------------------
# Front matter
# ---------------------------------------------------------------------------
p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
p.paragraph_format.space_after = Pt(4)
r = p.add_run("Robust Multigrid Pressure Correction for Parallel Opposed-Flow Reacting-Flow Simulation: Verification, Optimization, and Scaling")
set_run_font(r, size=14, bold=True)

p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
p.paragraph_format.space_after = Pt(2)
r = p.add_run("Ömer Uğur Zayıfoğlu")
set_run_font(r, size=10.5, bold=True)

p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
p.paragraph_format.space_after = Pt(1)
r = p.add_run("Department / Institute, University / Organization, Istanbul, Türkiye")
set_run_font(r, size=10)

p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
p.paragraph_format.space_after = Pt(7)
r = p.add_run("* Corresponding author: author@domain.edu")
set_run_font(r, size=10, italic=True)

add_section("Abstract")
abstract = (
"Opposed-flow reacting-flow simulations repeatedly solve an elliptic pressure-correction problem while advancing momentum, species, temperature, and chemistry, so pressure-solver robustness becomes increasingly important as the mesh and physical stiffness change. This study evaluates Martynenko’s Robust Multigrid Technique (RMT) as the pressure-correction strategy in a two-dimensional structured finite-volume reacting-flow solver. The implementation uses factor-three multiple shifted coarse grids, control-volume restriction of the fine-grid defect, a coarsest-to-finest sawtooth traversal, direct coarsest-grid solution, postsmoothing only, and full correction without conventional interpolation or a residual-dependent fallback. Verification combines four method-of-manufactured-solutions cases with a controlled reacting-flow comparison against single-grid red–black Gauss–Seidel, factor-two V- and W-cycle multigrid, and a factor-three V-cycle under identical physics, pressure tolerances, meshes, and OpenMP thread counts. At 16 threads in the low-flow sequence, the mean RMT cycles per pressure solve are 0.293, 0.411, 0.393, and 0.451 for 3,888, 15,552, 62,208, and 248,832 cells, respectively. Thus, a 64-fold increase in cell count produces only a limited change in the RMT intergrid iteration count. Implementation-level optimization of the same RMT cycle reduces pressure-solution time by 55.2% on a representative 432×144 case and by 58.5% on a representative 864×288 case, while preserving the cycle count, point-update count, and final pressure residual. The results indicate that RMT is primarily a robustness-oriented pressure strategy with weak mesh dependence and practically competitive total cost, rather than a method that universally minimizes wall-clock time at every operating point."
)
add_p(abstract, after=3)

p = doc.add_paragraph()
p.paragraph_format.space_after = Pt(2)
r = p.add_run("Keywords: ")
set_run_font(r, bold=True)
r = p.add_run("robust multigrid technique; reacting flow; pressure correction; computational fluid dynamics; OpenMP; method of manufactured solutions")
set_run_font(r)

p = doc.add_paragraph()
p.paragraph_format.space_after = Pt(6)
r = p.add_run("Track: ")
set_run_font(r, bold=True)
r = p.add_run("Track 1 — Fluid Mechanics and Computational Fluid Dynamics (High-Performance CFD)")
set_run_font(r)

# ---------------------------------------------------------------------------
# 1. Introduction
# ---------------------------------------------------------------------------
add_section("1. Introduction")
intro1 = (
"Reacting-flow computations couple momentum, pressure, scalar transport, thermochemistry, and chemical source terms over many physical time steps. In low-Mach formulations, the pressure-correction equation is solved repeatedly to enforce the velocity–pressure coupling, and its numerical cost can become significant as the spatial mesh is refined. Conventional multigrid can provide near-optimal complexity when smoothing, transfer, and cycling components are well matched to the operator, but these components may require tuning when coefficient variation, anisotropy, physical regime, or parallel granularity changes [1], [2]."
)
add_p(intro1)

intro2 = (
"Martynenko developed the Robust Multigrid Technique (RMT) for black-box scientific software, with the objective of reducing dependence on problem-specific numerical tuning while retaining efficient convergence over a family of boundary-value problems [1], [2]. The method differs structurally from a conventional V-cycle: in two dimensions the factor-three construction produces multiple shifted coarse grids; coarse-grid right-hand sides are formed from the fine-grid defect by control-volume averaging; the hierarchy is traversed from coarse to fine without presmoothing; and corrections remain associated with the fine-grid index space rather than being transferred by conventional interpolation. The multiple coarse-grid structure also exposes geometric parallelism, while red–black smoothing provides algebraic parallelism [3]."
)
add_p(intro2)

intro3 = (
"The present study asks whether this robustness-oriented construction can serve as the pressure solver of a parallel opposed-flow reacting-flow code without changing the underlying physics. The pressure algorithm is therefore varied while the momentum, species, chemistry, thermophysical treatment, mesh, time-step history, pressure tolerances, and thread counts are kept fixed. The objectives are to quantify mesh dependence of the RMT convergence rate, assess the computational cost relative to conventional in-house solvers, measure OpenMP scaling, and determine how much of the original RMT cost is implementation overhead rather than inherent numerical work."
)
add_p(intro3)

# ---------------------------------------------------------------------------
# 2. Materials and Methods
# ---------------------------------------------------------------------------
add_section("2. Materials and Methods")
m1 = (
"The target problem is a two-dimensional opposed-flow reacting-flow configuration in which fuel and oxidizer enter from opposing boundaries and form an interior reaction layer. A structured finite-volume discretization advances momentum, a pressure-correction equation, transported species, sensible energy, and a simplified Arrhenius chemical source. The production pressure convergence criteria are a relative tolerance of 10⁻⁴ and an absolute tolerance of 10⁻⁶. The controlled comparison uses four meshes, 108×36, 216×72, 432×144, and 864×288 cells, and 1, 2, 4, and 16 OpenMP threads. Six physical configurations span three inlet-flow intensities and three reaction-stiffness levels, giving 96 reacting-flow configurations for each pressure-solver formulation."
)
add_p(m1)

m2 = (
"The compared pressure algorithms are: (i) single-grid red–black Gauss–Seidel (SG-RBGS); (ii) factor-two V-cycle multigrid, V(2)-MG; (iii) factor-two W-cycle multigrid, W(2)-MG; (iv) factor-three V-cycle multigrid, V(3)-MG; and (v) RMT. OpenFOAM is retained as an application-level external reference under matched physical and spatial conditions, whereas the direct pressure-algorithm comparison is restricted to the in-house implementations so that software-architecture differences do not contaminate the interpretation."
)
add_p(m2)

m3 = (
"For RMT, the fine-grid pressure defect is rₕ = bₕ − Aₕpₕ. The factor-three hierarchy is formed from independent shifted grids. Coarse right-hand sides are obtained by control-volume averaging of the exact fine-grid defect, including boundary control volumes consistent with the coarse finite-volume operator. The cycle uses no presmoothing. The smallest shifted-grid systems are solved directly, followed by a coarsest-to-finest sawtooth traversal with parallel red–black Gauss–Seidel postsmoothing and a full correction update. No residual-monotonic line search, case-specific damping, mesh-specific smoothing count, or hidden single-grid fallback is used. A single globally fixed production postsmoothing count of 16 is retained."
)
add_p(m3)

m4 = (
"Verification precedes reacting-flow performance assessment. Four manufactured elliptic problems are used: a smooth Poisson-type problem, a high-contrast diffusion problem with coefficient contrast 10⁶, and two anisotropic problems with kᵧ = 0.03 and 0.003. In each case the exact field is prescribed first and the continuous forcing term is derived from the governing operator, following standard MMS practice [4], [5]. This separates discretization/solver verification from application-level comparison."
)
add_p(m4)

m5 = (
"After the numerically correct RMT cycle was frozen, implementation overhead was reduced without altering the mathematical iteration. The optimized implementation forms the fine-defect prefix integral once per RMT cycle and reuses it at all coarse levels, precomputes shifted-grid finite-volume coefficients, precomputes red/black point maps, and caches LU factorizations of invariant coarsest matrices while solving each current right-hand side directly. Regression checks require the optimized and reference implementations to produce equivalent corrections to roundoff, so any timing reduction is interpreted as lower implementation cost rather than retuning of the method."
)
add_p(m5)

# ---------------------------------------------------------------------------
# 3. Results
# ---------------------------------------------------------------------------
add_section("3. Results")
r1 = (
"The MMS study confirms convergence across smooth, discontinuous-coefficient, and strongly anisotropic operators. The most demanding representative case, kᵧ = 0.003 at n = 107, gives infinity-norm errors of 5.65×10⁻² for classical MG(3,3) and 5.84×10⁻² for the triple-coarsening RMT postsmoothing formulation. The corresponding serial solve times are 1.621 s and 0.332 s, respectively. These results illustrate why iteration count alone is not a sufficient performance metric: RMT may use more outer cycles while still reducing wall-clock cost through a different hierarchy and work distribution."
)
add_p(r1)

r2 = (
"In the reacting-flow application, the central robustness result is the weak dependence of the mean RMT cycle count on mesh resolution. Figure 1 compares RMT with the three conventional multigrid variants for the low-flow case at 16 threads. The RMT mean remains between 0.293 and 0.451 cycles per pressure solve while the number of cells increases by a factor of 64. The conventional multigrid methods also retain relatively stable cycle counts, as expected for effective multilevel methods, but the RMT values remain lower over this representative sequence."
)
add_p(r2)

p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
run = p.add_run()
run.add_picture(str(FIG), width=Cm(14.4))
p.paragraph_format.space_after = Pt(1)

p = doc.add_paragraph(style="IMAT Caption")
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = p.add_run("Figure 1. Mean intergrid iterations per pressure solve for the low-flow, 16-thread mesh sequence.")
set_run_font(r, size=9, italic=True)

table1_rows = [
    ["108×36", "3,888", "0.293", "23.66", "23.60", "1.00"],
    ["216×72", "15,552", "0.411", "50.13", "48.47", "1.03"],
    ["432×144", "62,208", "0.393", "348.44", "313.09", "1.11"],
    ["864×288", "248,832", "0.451", "4,826.37", "4,122.92", "1.17"],
]
p = doc.add_paragraph(style="IMAT Caption")
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = p.add_run("Table 1. RMT mesh dependence and total wall time relative to the mean of SG-RBGS, V(2)-MG, W(2)-MG, and V(3)-MG for the low-flow 16-thread sequence.")
set_run_font(r, size=9, italic=True)
add_table(
    table1_rows,
    ["Mesh", "Cells", "RMT cycles/solve", "RMT total [s]", "Reference mean [s]", "Ratio"],
    font_size=8.5,
)

r3 = (
"Table 1 shows that robustness does not imply minimum absolute time at every point. RMT is essentially equal to the mean in-house reference cost on the smallest mesh, 3% above it on 216×72, 11% above it on 432×144, and 17% above it on 864×288. Thus, the method remains in the same practical computational range while preserving a small, nearly mesh-independent intergrid iteration count. By contrast, SG-RBGS requires approximately 1.73, 6.06, 12.49, and 25.71 sweeps per pressure solve over the same mesh sequence, showing the expected growth of single-grid work with refinement."
)
add_p(r3)

table2_rows = [
    ["RMT strong scaling", "108×36", "84.81", "23.66", "3.58"],
    ["RMT strong scaling", "216×72", "393.67", "50.13", "7.85"],
    ["RMT strong scaling", "432×144", "3,865.14", "348.44", "11.09"],
    ["Pressure optimization", "432×144, representative low-flow", "148.15", "66.35", "55.2% reduction"],
    ["Pressure optimization", "864×288, representative high-flow", "2,711.94", "1,125.79", "58.5% reduction"],
]
p = doc.add_paragraph(style="IMAT Caption")
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = p.add_run("Table 2. Representative OpenMP scaling and implementation-level RMT optimization.")
set_run_font(r, size=9, italic=True)
add_table(
    table2_rows,
    ["Metric", "Case", "Before / T1 [s]", "After / T16 [s]", "Speedup / reduction"],
    font_size=8.2,
)

r4 = (
"Parallel efficiency improves with problem size. For the completed low-flow comparisons, the 1-to-16-thread speedup rises from 3.58 on 108×36 to 7.85 on 216×72 and 11.09 on 432×144. The same trend is consistent with the intended combination of geometric and algebraic parallelism in RMT [3]: small meshes expose synchronization and thread-management overhead, while larger meshes provide enough work to amortize those costs."
)
add_p(r4)

r5 = (
"The implementation optimization has a separate effect from numerical convergence. On a representative 432×144 low-flow case at 16 threads, pressure-solution time decreases from 148.15 s to 66.35 s, a 55.2% reduction. On a representative 864×288 high-flow case, it decreases from 2,711.94 s to 1,125.79 s, a 58.5% reduction. The corresponding complete-solver wall-clock reductions are approximately 19.4% and 23.5%. Because the optimized implementation preserves RMT cycle count, point-update count, and final pressure residual, these reductions isolate the cost of repeated indexing, coefficient construction, restriction support, and coarse factorization rather than a change in the RMT algorithm."
)
add_p(r5)

# ---------------------------------------------------------------------------
# 4. Discussion and Conclusion
# ---------------------------------------------------------------------------
add_section("4. Discussion and Conclusion")
d1 = (
"The results support a distinction between numerical robustness and absolute computational cost. The RMT cycle is more expensive than a conventional V-cycle because it operates on several shifted coarse grids, so a smaller cycle count does not automatically imply a smaller wall-clock time. Nevertheless, the key convergence metric, Ncycle(h), remains weakly dependent on mesh refinement in the reacting-flow pressure block. Increasing the problem size from 3,888 to 248,832 cells changes the low-flow 16-thread mean from 0.293 to only 0.451 RMT cycles per pressure solve. This behavior is consistent with the mesh-independent convergence objective of Martynenko’s framework [1], [2]."
)
add_p(d1)

d2 = (
"The practical limitation is therefore cost per RMT cycle rather than loss of convergence under refinement. The optimization study shows that more than half of the pressure-solution time in representative medium and large cases can be removed without altering the mathematical correction. This is important for thermofluid applications because pressure correction is invoked repeatedly inside the coupled transient calculation: even modest per-solve overhead accumulates over many time steps, whereas implementation-level savings propagate directly to the complete reacting-flow wall time."
)
add_p(d2)

d3 = (
"The present comparison is intentionally conservative. OpenFOAM is treated as an application-level external reference rather than a direct pressure-solver benchmark, and the in-house pressure methods share the same physics and convergence criteria. The study also does not claim that RMT is universally the fastest solver at every mesh or operating point. Its main advantage is the combination of a small and weakly mesh-dependent iteration count, one globally fixed postsmoothing parameter, and parallel behavior that improves with problem size. Future work will extend the same controlled methodology to broader reacting-flow conditions, quantify end-to-end comparisons with the external OpenFOAM workflow, and examine distributed-memory extensions in addition to the present shared-memory OpenMP implementation."
)
add_p(d3)

d4 = (
"In summary, a Martynenko-style RMT pressure correction has been verified and integrated into a parallel opposed-flow reacting-flow solver. Across the representative mesh ladder, the method maintains robust convergence while remaining within the same practical order of total computational cost as conventional in-house multigrid and single-grid references. Implementation optimization reduces representative pressure-solution times by approximately 55–59% without changing the numerical RMT cycle. These findings make RMT a relevant candidate for high-performance CFD applications in which consistent solver behavior across changes in mesh size and physical stiffness is as important as minimizing the run time of one individually tuned case."
)
add_p(d4)

# ---------------------------------------------------------------------------
# Acknowledgements and references
# ---------------------------------------------------------------------------
add_section("Acknowledgements")
add_p(
"The computational campaigns were designed for execution on the TRUBA high-performance computing infrastructure. The author gratefully acknowledges the scientific discussions on the Robust Multigrid Technique that informed the present implementation.",
after=5,
)

add_section("References")
refs = [
    '[1] S. I. Martynenko, “Robust multigrid technique for black box software,” Computational Methods in Applied Mathematics, vol. 6, no. 4, pp. 413–435, 2006, doi: 10.2478/cmam-2006-0026.',
    '[2] S. I. Martynenko, The Robust Multigrid Technique: For Black-Box Software. Berlin, Germany: De Gruyter, 2017, doi: 10.1515/9783110539264.',
    '[3] S. Martynenko, W. Zhou, I. Gökalp, V. Bakhtin, and P. Toktaliev, “Parallelization of Robust Multigrid Technique Using OpenMP Technology,” in Parallel Computing Technologies, PaCT 2021, Lecture Notes in Computer Science, vol. 12942, pp. 196–209, Springer, 2021, doi: 10.1007/978-3-030-86359-3_15.',
    '[4] P. J. Roache, “Code verification by the method of manufactured solutions,” Journal of Fluids Engineering, vol. 124, no. 1, pp. 4–10, 2002.',
    '[5] W. L. Oberkampf and C. J. Roy, Verification and Validation in Scientific Computing. Cambridge, U.K.: Cambridge University Press, 2010.',
    '[6] U. Trottenberg, C. W. Oosterlee, and A. Schüller, Multigrid. San Diego, CA, USA: Academic Press, 2001.',
]
for ref in refs:
    p = doc.add_paragraph()
    p.paragraph_format.left_indent = Cm(0.0)
    p.paragraph_format.first_line_indent = Cm(0.0)
    p.paragraph_format.space_after = Pt(2)
    r = p.add_run(ref)
    set_run_font(r, size=9)

# Footer
for section in doc.sections:
    fp = section.footer.paragraphs[0]
    fp.alignment = WD_ALIGN_PARAGRAPH.CENTER
    rr = fp.add_run("IMAT 2026 — Extended Abstract")
    set_run_font(rr, size=8, italic=True)

# Document properties
doc.core_properties.title = "Robust Multigrid Pressure Correction for Parallel Opposed-Flow Reacting-Flow Simulation"
doc.core_properties.subject = "IMAT 2026 Extended Abstract"
doc.core_properties.author = "Ömer Uğur Zayıfoğlu"
doc.core_properties.keywords = "RMT; CFD; reacting flow; multigrid; OpenMP; MMS"

doc.save(DOCX)
print(DOCX)
