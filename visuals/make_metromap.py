#!/usr/bin/env python3
"""
Generate the metaViraVerse metro map (nf-core workflow-schematic style).

Usage:
    python3 visuals/make_metromap.py [--font visuals/MavenPro.ttf] [--out visuals/metaviraverse_metromap.svg]

If --font is given, the font is embedded in the SVG (base64) so the image renders
identically everywhere (GitHub, browsers). Maven Pro is the font recommended by the
nf-core brand guidelines: https://fonts.google.com/specimen/Maven+Pro

Render the PNG with a headless browser, e.g.:
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" --headless \\
        --screenshot=visuals/metaviraverse_metromap.png --window-size=4820,1240 \\
        --hide-scrollbars visuals/metaviraverse_metromap.svg

Layout
------
Lines are the sequence types: viral sequences, prophages and plasmids run together
through pre-processing (the trunk); third-party input joins them at renaming after
its own preparation. At separation viral sequences and prophages become one line
(viral_sequences and prophages); after deduplication it takes the upper row and
plasmids the lower row; all lines meet again for the catalogue summary.
"""

import argparse
import base64
from pathlib import Path

H = 1240
SHIFT_Y = 40      # extra space at the top: everything moves down, grey boxes grow upwards
ICON_SCALE = 1.25 # size of the input/output file icons
STATION_ICON_SCALE = 1.0  # size of the file icons of additional inputs next to a station

# ---------------------------------------------------------------------------
# Lines, ordered top -> bottom inside the trunk
# ---------------------------------------------------------------------------
LINES = {
    "third":    {"color": "#0CA678", "label": "Third-party input (virus/plasmid)"},
    "viral":    {"color": "#1C7ED6", "label": "Viral sequences"},
    "prophage": {"color": "#FAB005", "label": "Prophages"},
    "viruses":  {"color": "#E03131", "label": "Viral sequences and prophages"},
    "plasmid":  {"color": "#AE3EC9", "label": "Plasmids"},
}
TRUNK = ["third", "viral", "prophage", "plasmid"]
MAIN = ["viral", "prophage", "plasmid"]   # before separation
VIRUSES = ["viruses"]                     # viral sequences and prophages, one line after separation
POOLED = ["viruses", "plasmid"]           # after separation

GAP = 14          # distance between parallel tracks
LW = 9            # line width
R = 13            # station radius

THIRD_Y = 250     # y of the third-party preparation row
TRUNK_TOP = 400   # y of the first track in the trunk
VIRAL_TOP = 262   # y of the first track of the viruses row
PLASMID_Y = 610   # y of the plasmids row

TY = {k: TRUNK_TOP + i * GAP for i, k in enumerate(TRUNK)}
TY["viruses"] = TY["prophage"]   # the pooled line continues on the prophage track
VY = {"viruses": VIRAL_TOP + GAP // 2}

X_START = 150
X_THIRD_JOIN0, X_THIRD_JOIN1 = 520, 600   # third-party row curves into the trunk
X_SEPARATE = 860                          # viral sequences and prophages become one line here
X_SPLIT0, X_SPLIT1 = 1420, 1500           # viruses up, plasmids down
X_MERGE0, X_MERGE1 = 4130, 4210           # back into the trunk
X_END = 4580

# viral clustering fork: vclust or vjdb-clustering (--viral_clustering)
FORK = (1560, 1620, 1820, 1880)   # out0, out1, in0, in1
FORK_X = 1720
FORK_DY = 52

W = X_END + 240

FONT = "'Maven Pro', 'Helvetica Neue', Arial, sans-serif"
TOOL_WEIGHT = 900  # tool names (Maven Pro is a variable font: 400-900)

# ---------------------------------------------------------------------------
# Stages (grey boxes)
# ---------------------------------------------------------------------------
STAGES = [
    (170, 560, ["1", "Input & third-party", "preparation"]),
    (575, 1430, ["2", "Pre-processing"]),
    (1445, 1910, ["3", "ANI-based clustering"]),
    (1925, 3000, ["4", "Taxonomy, host & mobility"]),
    (3015, 3640, ["5", "Protein annotation"]),
    (3655, 3960, ["6", "Protein clustering"]),
    (3975, 4568, ["7", "Catalogue & reporting"]),
]
BOX_TOP, BOX_BOTTOM = 60 - SHIFT_Y, 780

# ---------------------------------------------------------------------------
# Stations: (x, lines, tool name, description lines, optional, row)
#   tool name: a string, or a list of names stacked above the station (last one closest)
#   row: "trunk" | "third" | "viruses" | "plasmid" | "fork_top" | "fork_bottom"
# ---------------------------------------------------------------------------
AMR_TOOLS = ["AMRFinderPlus", "DeepARG", "RGI"]   # short name closest to the station

# Additional inputs of a station, drawn as a file icon above it with the caption on its right:
# (station tool name, x shift of the icon from the station, y of the icon top, icon label, caption lines, grey note)
STATION_INPUTS = [
    ("SpacePHARER", 0, BOX_TOP + 72, "FASTA", ["Custom CRISPR", "spacers & metadata"], "(--custom_spacers_*)"),
    ("metadata.py", -70, 278, "TSV", ["Metadata"], ["(--catalogues_metadata", "--analyses_metadata)"]),
]

STATIONS = [
    # 1. Third-party preparation
    (300, ["third"], "faLint", ["FASTA", "validation"], False, "third"),
    (450, ["third"], "Pyrodigal", ["Gene", "calling"], False, "third"),
    # 2. Pre-processing
    (680, TRUNK, "rename.py", ["Assign", "catalogue", "accessions"], False, "trunk"),
    (X_SEPARATE, MAIN, "separate.py", ["viral_sequences /", "prophages /", "plasmids"], False, "trunk"),
    (1040, VIRUSES, "CheckV", ["Quality", "evaluation"], True, "trunk"),
    (1200, VIRUSES, "barrnap", ["rRNA", "detection"], True, "trunk"),
    (1350, POOLED, "SHA-256", ["Deduplication", "& QC filter"], False, "trunk"),
    # 3. Clustering
    (FORK_X, VIRUSES, "vclust", ["ANI 95%, cov 85%"], False, "fork_top"),
    (FORK_X, VIRUSES, "vjdb-clustering", ["vclust + Leiden", "ANI 95%, cov 85%"], False, "fork_bottom"),
    (FORK_X, ["plasmid"], "vclust", ["Cluster into", "pOTUs", "ANI 70%, cov 50%"], False, "plasmid"),
    # 4. Taxonomy, host & mobility (cluster representatives)
    (2020, VIRUSES, "geNomad", ["Taxonomy"], True, "viruses"),
    (2210, VIRUSES, "VITAP", ["Taxonomy"], True, "viruses"),
    (2400, VIRUSES, "iPHoP", ["Host", "prediction"], True, "viruses"),
    (2600, VIRUSES, "SpacePHARER", ["Host prediction", "from database of", "CRISPR spacers"], True, "viruses"),
    (2830, VIRUSES, "BACPHLIP", ["Lifestyle"], False, "viruses"),
    (2160, ["plasmid"], "MOB-suite", ["Typing &", "mobility"], False, "plasmid"),
    (2600, ["plasmid"], "plaSquid", ["Replicons,", "Incompatibility", "& mobility"], False, "plasmid"),
    # 5. Protein annotation (cluster representatives)
    (3110, VIRUSES, ["HMMER", "ViPhOGs"], ["Protein-based", "taxonomy", "assignment"], True, "viruses"),
    (3310, VIRUSES, AMR_TOOLS, ["Antimicrobial", "resistance"], True, "viruses"),
    (3510, VIRUSES, "HMMER", ["Functional", "annotation"], True, "viruses"),
    (3310, ["plasmid"], AMR_TOOLS, ["Antimicrobial", "resistance"], True, "plasmid"),
    # 6. Protein clustering (cluster representatives)
    (3807, VIRUSES, "PhaMMseqs", ["Protein", "phamilies"], True, "viruses"),
    (3807, ["plasmid"], "MMseqs2", ["Protein", "clustering"], True, "plasmid"),
    # 7. Catalogue & reporting
    (4050, VIRUSES, ["Build", "final GFF"], ["Merge", "annotations"], False, "viruses"),
    (4050, ["plasmid"], ["Build", "final GFF"], ["Merge", "annotations"], False, "plasmid"),
    (4320, POOLED, "metadata.py", ["Statistics &", "metadata tables"], False, "trunk"),
    (4490, POOLED, "report.py", ["README &", "reports"], False, "trunk"),
]


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def path_through(points):
    """Polyline through waypoints; y changes are drawn as smooth S-curves."""
    d = f"M{points[0][0]},{points[0][1]}"
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if y0 == y1:
            d += f" L{x1},{y1}"
        else:
            xm = (x0 + x1) / 2
            d += f" C{xm},{y0} {xm},{y1} {x1},{y1}"
    return d


def line_paths():
    out = []
    # third-party: own row, then into the top trunk track; ends at renaming
    out.append(("third", path_through([(X_START, THIRD_Y), (X_THIRD_JOIN0, THIRD_Y),
                                       (X_THIRD_JOIN1, TY["third"]), (680, TY["third"])])))
    # plasmids: trunk -> lower row -> trunk
    k = "plasmid"
    out.append((k, path_through([(X_START, TY[k]), (X_SPLIT0, TY[k]), (X_SPLIT1, PLASMID_Y),
                                 (X_MERGE0, PLASMID_Y), (X_MERGE1, TY[k]), (X_END, TY[k])])))
    # viral sequences and prophages: separate lines up to separation, where they end
    for k in ("viral", "prophage"):
        out.append((k, path_through([(X_START, TY[k]), (X_SEPARATE, TY[k])])))
    # pooled viruses: from separation, trunk -> upper row with the clustering fork -> trunk
    f0, f1, f2, f3 = FORK
    for k in VIRUSES:
        y = VY[k]
        head = [(X_SEPARATE, TY[k]), (X_SPLIT0, TY[k]), (X_SPLIT1, y), (f0, y)]
        up = [(f0, y), (f1, y - FORK_DY), (f2, y - FORK_DY), (f3, y)]
        down = [(f0, y), (f1, y + FORK_DY), (f2, y + FORK_DY), (f3, y)]
        tail = [(f3, y), (X_MERGE0, y), (X_MERGE1, TY[k]), (X_END, TY[k])]
        out.append((k, path_through(head + up[1:] + tail[1:])))
        out.append((k, path_through(down)))
    return out


def text(x, y, s, size=22, weight=400, anchor="middle", fill="#000"):
    return (f'<text x="{x}" y="{y}" font-family="{FONT}" font-size="{size}" '
            f'font-weight="{weight}" text-anchor="{anchor}" fill="{fill}">{esc(s)}</text>')


def station_y(lines, row):
    if row == "trunk":
        ys = [TY[k] for k in lines]
    elif row == "viruses":
        ys = [VY[k] for k in lines]
    elif row == "fork_top":
        ys = [VY[k] - FORK_DY for k in lines]
    elif row == "fork_bottom":
        ys = [VY[k] + FORK_DY for k in lines]
    elif row == "third":
        ys = [THIRD_Y]
    else:
        ys = [PLASMID_Y]
    return min(ys), max(ys)


def label_anchors(lines, row):
    """y for the tool-name baseline (above) and first description line (below)."""
    y0, y1 = station_y(lines, row)
    if row == "trunk":
        # the name sits above the station's own top track (the trunk above it has ended: third-party
        # after renaming); the description clears the whole trunk, as lower tracks may pass under it
        return y0 - 30, TY[TRUNK[-1]] + 44
    return y0 - 30, y1 + 44


def station(x, lines, name, desc, optional, row):
    y0, y1 = station_y(lines, row)
    s = [f'<rect x="{x - R}" y="{y0 - R}" width="{2 * R}" height="{y1 - y0 + 2 * R}" '
         f'rx="{R}" fill="#fff" stroke="#000" stroke-width="4"/>']
    if not optional:   # mandatory steps carry a black dot
        r = 5.5
        s.append(f'<rect x="{x - r}" y="{y0 - r}" width="{2 * r}" height="{y1 - y0 + 2 * r}" '
                 f'rx="{r}" fill="#000"/>')
    top, bottom = label_anchors(lines, row)
    if isinstance(name, list):   # several tools stacked above the station
        for i, n in enumerate(name):
            s.append(text(x, top - 32 * (len(name) - 1 - i), n, 28, TOOL_WEIGHT))
        for i, d in enumerate(desc):
            s.append(text(x, bottom + i * 30, d, 26))
    elif row == "fork_top":
        # below the upper branch is the gap between the branches: description goes above, under the name
        s.append(text(x, top - 30 * len(desc), name, 28, TOOL_WEIGHT))
        for i, d in enumerate(desc):
            s.append(text(x, top - 30 * (len(desc) - 1 - i), d, 26))
    elif row == "fork_bottom":
        # the name goes below the lower branch, the description under it
        s.append(text(x, bottom, name, 28, TOOL_WEIGHT))
        for i, d in enumerate(desc):
            s.append(text(x, bottom + 30 * (i + 1), d, 26))
    else:
        s.append(text(x, top, name, 28, TOOL_WEIGHT))
        for i, d in enumerate(desc):
            s.append(text(x, bottom + i * 30, d, 26))
    return "\n".join(s)


def file_icon(x, y, label, caption, scale=1.0, caption_right=False, note=None):
    """nf-core style document icon; (x, y) = top-left corner.

    The caption goes below the icon, or to its right with caption_right (for icons placed
    next to a station); note is an extra grey line, or a list of lines, after the caption
    (e.g. the parameters)."""
    w, h, fold = 62 * scale, 80 * scale, 18 * scale
    p = f"M{x},{y} H{x + w - fold} L{x + w},{y + fold} V{y + h} H{x} Z"
    s = [f'<path d="{p}" fill="#fff" stroke="#000" stroke-width="5" stroke-linejoin="round"/>',
         f'<path d="M{x + w - fold},{y} V{y + fold} H{x + w}" fill="none" stroke="#000" '
         f'stroke-width="4" stroke-linejoin="round"/>',
         f'<rect x="{x - 8}" y="{y + h * 0.42}" width="{w + 16}" height="{26 * scale}" rx="4" fill="#000"/>',
         text(x + w / 2, y + h * 0.42 + 19 * scale, label, 18 * scale, 700, fill="#fff")]
    notes = [note] if isinstance(note, str) else (note or [])
    lines = [(c, 26, "#000") for c in caption] + [(n, 22, "#555") for n in notes]
    for i, (c, size, fill) in enumerate(lines):
        if caption_right:
            s.append(text(x + w + 18, y + h / 2 - 15 * (len(lines) - 1) + 8 + i * 30, c, size, 400, "start", fill))
        else:
            s.append(text(x + w / 2, y + h + 32 + i * 30, c, size, 400, fill=fill))
    return "\n".join(s)


def build(font_b64=None):
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">']
    if font_b64:
        o.append("<defs><style>@font-face{font-family:'Maven Pro';font-weight:400 900;"
                 f"src:url(data:font/ttf;base64,{font_b64}) format('truetype');}}</style></defs>")
    o.append(f'<rect width="{W}" height="{H}" fill="#fff"/>')
    o.append(f'<g transform="translate(0,{SHIFT_Y})">')

    # stage boxes
    for x0, x1, (num, *title) in STAGES:
        o.append(f'<rect x="{x0}" y="{BOX_TOP}" width="{x1 - x0}" height="{BOX_BOTTOM - BOX_TOP}" '
                 f'rx="24" fill="#E9E9E9"/>')
        title = [f"{num}. {title[0]}", *title[1:]]
        for i, t in enumerate(title):
            o.append(text(x0 + 24, BOX_TOP + 46 + i * 32, t, 28, 700, "start", "#8A8A8A"))

    # lines
    for k, d in line_paths():
        o.append(f'<path d="{d}" fill="none" stroke="{LINES[k]["color"]}" stroke-width="{LW}" '
                 f'stroke-linejoin="round"/>')
    # parameter choosing between the alternative viral clusterings, between the two branches
    o.append(text(FORK_X, (VY[VIRUSES[0]] + VY[VIRUSES[-1]]) / 2 + 9, "(--viral_clustering)", 22, 400, fill="#555"))

    # input termini (triangles) and output terminus (bar)
    ty0, ty1 = TY[MAIN[0]] - 12, TY[MAIN[-1]] + 12
    o.append(f'<path d="M{X_START},{ty0} L{X_START + 26},{(ty0 + ty1) / 2} L{X_START},{ty1} Z" fill="#000"/>')
    ty0, ty1 = TY[POOLED[0]] - 12, TY[POOLED[-1]] + 12
    o.append(f'<path d="M{X_START},{THIRD_Y - 12} L{X_START + 26},{THIRD_Y} L{X_START},{THIRD_Y + 12} Z" '
             f'fill="#000"/>')
    o.append(f'<rect x="{X_END}" y="{ty0}" width="11" height="{ty1 - ty0}" fill="#000"/>')

    # stations
    for st in STATIONS:
        o.append(station(*st))

    # additional inputs, above their station
    sw = 62 * STATION_ICON_SCALE
    for station_name, dx, y, label, caption, note in STATION_INPUTS:
        st = next(st for st in STATIONS if st[2] == station_name)
        o.append(file_icon(st[0] + dx - sw / 2, y, label, caption, STATION_ICON_SCALE, caption_right=True, note=note))

    # inputs / outputs
    iw = 62 * ICON_SCALE
    lx = 75.8 - iw / 2                    # widest caption ("samplesheet") starts at the left edge
    rx = W - 95 - iw / 2                  # widest caption ("catalogues") ends at the right edge
    o.append(file_icon(lx, 150, "CSV", ["Third-party", "samplesheet"], ICON_SCALE))
    o.append(file_icon(lx, 400, "CSV", ["MGnify", "samplesheet"], ICON_SCALE))
    ih = 80 * ICON_SCALE
    gap = 30
    y = 40
    for label, caption in [("FASTA", ["Viral & plasmid", "catalogues"]), ("GFF", ["Annotations"]),
                           ("TSV", ["Metadata"]), ("HTML", ["Reports"])]:
        o.append(file_icon(rx, y, label, caption, ICON_SCALE))
        y += ih + 40 + 30 * len(caption) + gap

    # ---------------- bottom: title + legend ----------------
    lt = 820
    o.append(text(STAGES[0][0], lt + 70, "ebi-metagenomics/", 46, 700, "start", "#555"))
    o.append(text(STAGES[0][0], lt + 140, "metaViraVerse", 86, 800, "start", "#000"))
    o.append(text(STAGES[0][0], lt + 190, "Viral and plasmid catalogues from MGnify and third-party predictions",
                  28, 400, "start", "#444"))

    # line legend
    bw = 1330
    bx = STAGES[-1][1] - bw   # right edge aligned with the last stage box
    o.append(f'<rect x="{bx}" y="{lt}" width="{bw}" height="350" rx="24" fill="#E9E9E9"/>')
    o.append(text(bx + 30, lt + 50, "Legend", 34, 700, "start"))
    for i, k in enumerate(LINES):
        y = lt + 100 + i * 42
        o.append(f'<line x1="{bx + 30}" y1="{y}" x2="{bx + 100}" y2="{y}" stroke="{LINES[k]["color"]}" '
                 f'stroke-width="{LW + 3}"/>')
        o.append(text(bx + 120, y + 8, LINES[k]["label"], 28, 400, "start"))

    # station legend
    sx = bx + 560
    items = [
        ("circle", "Mandatory"),
        ("optional", "Optional (--skip_* / databases given)"),
        ("in", "Input"),
        ("out", "Output"),
    ]
    for i, (kind, lab) in enumerate(items):
        y = lt + 100 + i * 42
        cx = sx + 20
        if kind == "circle":
            o.append(f'<circle cx="{cx}" cy="{y}" r="{R}" fill="#fff" stroke="#000" stroke-width="4"/>')
            o.append(f'<circle cx="{cx}" cy="{y}" r="5.5" fill="#000"/>')
        elif kind == "optional":
            o.append(f'<circle cx="{cx}" cy="{y}" r="{R}" fill="#fff" stroke="#000" stroke-width="4"/>')
        elif kind == "in":
            o.append(f'<path d="M{cx - 10},{y - 14} L{cx + 12},{y} L{cx - 10},{y + 14} Z" fill="#000"/>')
        else:
            o.append(f'<rect x="{cx - 5}" y="{y - 15}" width="10" height="30" fill="#000"/>')
        o.append(text(sx + 52, y + 9, lab, 28, 400, "start"))
    o.append(text(bx + 30, lt + 325, "Stages 4-7 run on cluster representatives only.", 24, 400, "start", "#555"))

    o.append("</g>")
    o.append("</svg>")
    return "\n".join(o)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--font", type=Path, help="TTF to embed (Maven Pro)")
    ap.add_argument("--out", type=Path, default=Path(__file__).parent / "metaviraverse_metromap.svg")
    a = ap.parse_args()
    font = base64.b64encode(a.font.read_bytes()).decode() if a.font else None
    a.out.write_text(build(font))
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
