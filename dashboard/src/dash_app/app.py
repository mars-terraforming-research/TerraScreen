# terrascreen_dash.py
import json
import math
import itertools
from functools import lru_cache

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import Dash, dcc, html, no_update
from dash.dependencies import Input, Output
from dash import html, dcc

# --- Splitter CSS and JS (inline) ---
SPLITTER_CSS = """
#grid{
  display:grid;
  grid-template-columns: var(--left,22rem) 10px minmax(0,1fr) 10px var(--right,26rem);
  grid-template-rows: 100vh;
  height:100vh; overflow:hidden;
}
.splitter{
  position:relative; z-index:5;
  background: linear-gradient(90deg,#e9ecef,#cfd8dc,#e9ecef);
  cursor: col-resize; user-select:none;
}
.splitter:hover{ background:#cfd8dc; }
/* Big hit area so dragging works anywhere near the divider */
.splitter::before{
  content:""; position:absolute; top:0; bottom:0; left:-8px; right:-8px;
}
#left,#center,#right{ overflow:auto; box-sizing:border-box; }
#left{ border-right:1px solid #eee; }
#right{ border-left:1px solid #eee; }
"""

SPLITTER_JS = """
document.addEventListener('DOMContentLoaded', () => {
  const grid = document.getElementById('grid');
  if (!grid) return;
  const s1 = document.getElementById('split-1');
  const s2 = document.getElementById('split-2');
  let dragging = 0;

  const rect = () => grid.getBoundingClientRect();

  function setLeft(px){ grid.style.setProperty('--left', Math.max(0,px) + 'px'); }
  function setRight(px){ grid.style.setProperty('--right', Math.max(0,px) + 'px'); }

  function onMove(clientX){
    if (!dragging) return;
    const r = rect();
    if (dragging === 1) setLeft(clientX - r.left);
    if (dragging === 2) setRight(r.right - clientX);
  }

  function onPointerMove(e){
    if (e.touches && e.touches[0]) onMove(e.touches[0].clientX);
    else onMove(e.clientX);
    if (dragging) e.preventDefault();
  }

  function onPointerUp(){
    dragging = 0;
    document.body.style.cursor = '';
  }

  function start(which){ return (e) => {
      dragging = which;
      document.body.style.cursor = 'col-resize';
      e.preventDefault();
    };
  }

  s1.addEventListener('mousedown', start(1));
  s2.addEventListener('mousedown', start(2));
  s1.addEventListener('touchstart', start(1), {passive:false});
  s2.addEventListener('touchstart', start(2), {passive:false});

  window.addEventListener('mousemove', onPointerMove, {passive:false});
  window.addEventListener('touchmove', onPointerMove, {passive:false});
  window.addEventListener('mouseup', onPointerUp);
  window.addEventListener('touchend', onPointerUp);
});
"""

# ------- HTTP (local Python or Pyodide) -------
try:
    from pyodide.http import open_url  # browser

    def fetch_text(url: str) -> str:
        return open_url(url).read()
except Exception:
    import requests  # local

    def fetch_text(url: str) -> str:
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        return r.text

@lru_cache(maxsize=256)
def fetch_text_cached(url: str) -> str:
    return fetch_text(url)

def fetch_json(url: str):
    return json.loads(fetch_text_cached(url))

# ------- Repo config -------
OWNER  = "mars-terraforming-research"
REPO   = "TerraScreen"
BRANCH = "main"
OUTPUT_DIR = "output"

RAW_BASE = f"https://raw.githubusercontent.com/{OWNER}/{REPO}/{BRANCH}/{OUTPUT_DIR}/"
API_URL  = f"https://api.github.com/repos/{OWNER}/{REPO}/contents/{OUTPUT_DIR}"

L_NSPECTI = 96
L_NSPECTV = 84
NLEV      = 24

IR_SCALE = 20  # <— fixed, no slider

# ------- Utilities -------
def _pretty(name: str) -> str:
    if name.endswith(".txt"):
        name = name[:-4]
    return name.replace("output_", "").replace("static_", "")

def _basename(name: str) -> str:
    return _pretty(name)

def list_repo_files(prefix: str):
    try:
        data = fetch_json(API_URL)
        names = [it["name"] for it in data if it.get("type") == "file"]
        return sorted([n for n in names if n.startswith(prefix) and n.endswith(".txt")])
    except Exception:
        return {
            "static_": ["static_dust1.5.txt", "static_Al_8um_60um.txt"],
            "output_": ["output_dust1.5.txt"],
        }[prefix]

def _floats_after_comma(line: str):
    if "," not in line:
        return []
    tail = line.split(",", 1)[1]
    vals = []
    for tok in tail.replace("  ", " ").split(","):
        tok = tok.strip()
        if not tok:
            continue
        try:
            vals.append(float(tok))
        except ValueError:
            pass
    return vals

def _find_line_idx_header(lines, *labels, start=0) -> int:
    labs = [l.lower() for l in labels]
    for i in range(start, len(lines)):
        s = lines[i].strip().lower()
        if any(s.startswith(lab) for lab in labs):
            return i
    return -1

def _find_line_idx_contains(lines, *substrs, start=0) -> int:
    subs = [s.lower() for s in substrs]
    for i in range(start, len(lines)):
        s = lines[i].lower()
        if any(sub in s for sub in subs):
            return i
    return -1

# ------- Parse: OUTPUT (Ts vs τ) -------
@lru_cache(maxsize=128)
def read_tau_ts(file_name: str) -> pd.DataFrame:
    text = fetch_text_cached(RAW_BASE + file_name)
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    hdr = _find_line_idx_header(lines, "it", "it   ,")
    tau, ts = [], []
    if hdr != -1:
        for ln in lines[hdr+1:]:
            if "," not in ln: break
            parts = [p.strip() for p in ln.split(",")]
            if len(parts) < 3: break
            try:
                tau.append(float(parts[1])); ts.append(float(parts[2]))
            except Exception:
                break
    df = pd.DataFrame({"tau": tau, "ts": ts})
    df["series"] = _pretty(file_name)
    if not df.empty:
        df = df.sort_values("tau")
    return df

# ------- Parse: STATIC (spectra at chosen τ) -------
class StaticCase:
    def __init__(self, meta, wl, Qext, Qscat, bwni, bwnv, solar_WN, OLR_WL, ASR_WL):
        self.meta = meta
        self.wl = wl
        self.Qext = Qext
        self.Qscat = Qscat
        self.bwni = bwni
        self.bwnv = bwnv
        self.solar_WN = solar_WN
        self.OLR_WL = OLR_WL
        self.ASR_WL = ASR_WL

@lru_cache(maxsize=128)
def read_static_case(file_name: str, tau_target_rounded: float) -> StaticCase:
    raw = fetch_text_cached(RAW_BASE + file_name)
    lines = [ln.rstrip() for ln in raw.splitlines()]

    # Core headers
    idx_bwni  = _find_line_idx_header(lines, "bwn ir")
    idx_bwnv  = _find_line_idx_header(lines, "bwn vis")
    idx_wl    = _find_line_idx_header(lines, "wavelenght", "wavelength")

    # IMPORTANT: Find spectral Qext/Qscat **after** the wavelength line
    idx_qext  = _find_line_idx_header(lines, "qext",  start=idx_wl + 1 if idx_wl != -1 else 0)
    idx_qscat = _find_line_idx_header(lines, "qscat", start=idx_qext + 1 if idx_qext != -1 else (idx_wl + 1))

    idx_pref  = _find_line_idx_header(lines, "p [mbar]")
    idx_tref  = _find_line_idx_header(lines, "temp initial")

    if min([idx_bwni, idx_bwnv, idx_wl, idx_qext, idx_qscat, idx_pref, idx_tref]) == -1:
        raise ValueError(f"Could not locate essential header blocks in {file_name}")

    # Optional per-band solar or scalar solar flux
    idx_sun_vis  = _find_line_idx_header(lines, "sun vis")
    idx_sun_flux = _find_line_idx_header(lines, "sun flux")
    solar_WN = None
    if idx_sun_vis != -1:
        arr = _floats_after_comma(lines[idx_sun_vis])
        solar_WN = np.array(arr) if len(arr) else None

    def floats_after(idx): return np.array(_floats_after_comma(lines[idx]))

    bwni   = floats_after(idx_bwni)
    bwnv   = floats_after(idx_bwnv)
    wl     = floats_after(idx_wl)
    Qext   = floats_after(idx_qext)
    Qscat  = floats_after(idx_qscat)
    _      = floats_after(idx_pref)
    _      = floats_after(idx_tref)

    # Iteration table (unchanged)
    idx_hdr = _find_line_idx_header(lines, "it", "it   ,")
    if idx_hdr == -1:
        idx_hdr = _find_line_idx_contains(lines, "it   ,")
    if idx_hdr == -1:
        raise ValueError("Could not find iteration table header")

    rows = []
    for ln in lines[idx_hdr+1:]:
        if "," not in ln: break
        parts = [p.strip() for p in ln.split(",")]
        try:
            vals = [float(x) for x in parts]
            rows.append(vals)
        except Exception:
            break
    if not rows:
        raise ValueError("No numeric rows found in iteration table")

    taus = np.array([r[1] for r in rows])
    itau = int(np.argmin(np.abs(taus - tau_target_rounded)))
    vals = rows[itau]

    base = 8
    p = base + NLEV
    OLR_WL = np.array(vals[p:p+L_NSPECTI]);  p += L_NSPECTI
    ASR_WL = np.array(vals[p:p+L_NSPECTV])

    meta = dict(tau=vals[1], ts=vals[2])
    return StaticCase(meta, wl=wl, Qext=Qext, Qscat=Qscat,
                      bwni=bwni, bwnv=bwnv, solar_WN=solar_WN,
                      OLR_WL=OLR_WL, ASR_WL=ASR_WL)

# ------- Band “tops” drawing helpers -------
def _segments_from_bounds(bounds_um: np.ndarray, values: np.ndarray):
    xs, ys = [], []
    for left, right, v in zip(bounds_um[:-1], bounds_um[1:], values):
        xs.extend([left, right, None])
        ys.extend([v,    v,    None])
    return xs, ys

def _vertical_sides(bounds_um: np.ndarray, values: np.ndarray):
    xs, ys = [], []
    for left, right, v in zip(bounds_um[:-1], bounds_um[1:], values):
        xs.extend([left, left, None,  right, right, None])
        ys.extend([0,    v,    None,  0,     v,     None])
    return xs, ys

def _add_bandline(fig, bounds_um, per_um, *, name, color=None, width=2.0, sides=False, side_width=0.5):
    x, y = _segments_from_bounds(bounds_um, per_um)
    fig.add_trace(go.Scatter(
        name=name, x=x, y=y, mode="lines",
        line=dict(width=width, color=color),
        hovertemplate="λ band<extra>"+name+"</extra>"
    ))
    if sides:
        xs, ys = _vertical_sides(bounds_um, per_um)
        fig.add_trace(go.Scatter(
            name=name+" (bands)", showlegend=False,
            x=xs, y=ys, mode="lines",
            line=dict(width=side_width, color=color),
            hoverinfo="skip"
        ))

# ------- Figures -------
def spectral_budget_figure(static_files, tau_target: float):
    if not static_files:
        return px.line(title="Select one or more cases that have spectral data")

    clear = read_static_case(static_files[0], round(0.0, 3))

    bnwlv = 1e4 / clear.bwnv[::-1]            # VIS edges
    bnwli = 1e4 / clear.bwni[::-1]            # IR edges
    all_bounds = np.append(bnwlv, bnwli[1:])
    res_wl     = np.diff(all_bounds)

    fig = go.Figure()

    # Optional solar VIS (per band)
    if (clear.solar_WN is not None) and (len(clear.solar_WN) >= L_NSPECTV):
        res_wl_vis = np.diff(bnwlv)
        solar_per_um = clear.solar_WN[::-1] / res_wl_vis
        _add_bandline(fig, bnwlv, solar_per_um, name="solar flux",
                      color="black", width=2.2, sides=False)

    # Clear (no aerosol), IR scaled by fixed constant
    net_clear = np.append(clear.ASR_WL[::-1], -IR_SCALE * clear.OLR_WL[::-1])
    clear_per_um = net_clear / res_wl
    _add_bandline(fig, all_bounds, clear_per_um, name="clear (no aerosol)",
                  color="orange", width=2.0, sides=True, side_width=0.7)

    # Selected static files at τ
    palette = itertools.cycle(["#17becf", "#d62728", "#9467bd", "#2ca02c"])
    for fname in static_files:
        case = read_static_case(fname, round(float(tau_target), 3))
        net_stack = np.append(case.ASR_WL[::-1], -IR_SCALE * case.OLR_WL[::-1])
        per_um    = net_stack / res_wl
        _add_bandline(fig, all_bounds, per_um, name=_pretty(fname),
                      color=next(palette), width=2.0, sides=False)

    fig.update_layout(
        title=f"Radiative budget at the top of the atmosphere (τ={tau_target:g}; IR scaled ×{IR_SCALE})",
        margin=dict(l=4, r=4, t=80, b=90),
        legend=dict(orientation="h", y=-0.25, yanchor="top", x=0, xanchor="left"),
        xaxis_title="Wavelength [μm]",
        yaxis_title=f"Irradiance [W/m²/μm]"
    )
    fig.update_xaxes(type="log", range=[math.log10(0.2), math.log10(60.0)], ticks="outside")
    fig.update_yaxes(zeroline=True)
    return fig

def optical_props_figure(cases, tau_target: float):
    """Bottom plot: Extinction only, log-scale Y with safe handling of zeros and length mismatches."""
    if not cases:
        return px.line(title="Select cases that have spectral data")

    fig = go.Figure()
    palette = itertools.cycle(["#d62728", "#17becf", "#9467bd", "#2ca02c"])
    ymin, ymax = np.inf, -np.inf

    for name in cases:
        case = read_static_case(name, round(float(tau_target), 3))

        x = np.asarray(case.wl,   dtype=float)
        y = np.asarray(case.Qext, dtype=float)

        # Align lengths defensively (some files can have off-by-one differences)
        n = int(min(len(x), len(y)))
        if n < 2:
            continue
        x, y = x[:n], y[:n]

        # Sort by wavelength
        order = np.argsort(x)
        x, y = x[order], y[order]

        # Log-safety: mask non-positive
        y_safe = y.copy()
        y_safe[y_safe <= 0] = np.nan
        if np.all(np.isnan(y_safe)):
            continue

        ymin = min(ymin, np.nanmin(y_safe))
        ymax = max(ymax, np.nanmax(y_safe))

        color = next(palette)
        fig.add_trace(go.Scatter(
            x=x, y=y_safe, mode="lines",
            name=f"{_pretty(name)} Extinction",
            line=dict(color=color, width=2)
        ))

    if not np.isfinite(ymin) or not np.isfinite(ymax):
        return px.line(title="No positive extinction values to plot on log scale")

    lo = max(ymin / 1.5, 1e-6)
    hi = ymax * 1.5

    fig.update_layout(
        title="Optical properties for the aerosols considered (Extinction only)",
        margin=dict(l=4, r=4, t=80, b=80),
        legend=dict(orientation="h", y=-0.25, yanchor="top", x=0, xanchor="left"),
        xaxis_title="Wavelength [μm]",
        yaxis_title="Extinction Efficiency",
    )
    fig.update_xaxes(type="log", range=[math.log10(0.2), math.log10(60.0)])
    fig.update_yaxes(type="log", range=[math.log10(lo), math.log10(hi)])

    return fig

def ts_vs_tau_figure(output_files):
    if not output_files:
        return px.line(title="Select cases that have temperature data (Ts vs τ)")
    frames = []
    for fname in output_files:
        df = read_tau_ts(fname)
        if not df.empty:
            frames.append(df)
    if not frames:
        return px.line(title="No temperature data")
    df = pd.concat(frames, ignore_index=True)
    fig = px.line(
        df, x="tau", y="ts", color="series",
        labels={"tau": "Visible opacity τ (670 nm)", "ts": "Surface temperature Ts [K]"},
        title="Ts vs τ"
    )
    fig.update_layout(margin=dict(l=0, r=0, t=60, b=0),
                      legend=dict(orientation="h", y=-0.25, yanchor="top", x=0, xanchor="left"))
    return fig

# ------- Build unified case list (no 'static'/'output' labels) -------
STATIC_FILES = list_repo_files("static_")
OUTPUT_FILES = list_repo_files("output_")

STATIC_MAP = { _basename(f): f for f in STATIC_FILES }
OUTPUT_MAP = { _basename(f): f for f in OUTPUT_FILES }
ALL_CASES  = sorted(set(STATIC_MAP) | set(OUTPUT_MAP))  # union

# ------- Layout (flexbox with resizable/collapsible sides) -------
app = Dash(__name__)

app.layout = html.Div(
    [
        # ------ Controls block (multi-column checklist + constrained slider) ------
        html.Div(
            [
                html.H3("TerraScreen dashboard", style={"margin": "0 0 0.5rem 0"}),

                html.H4("Cases", style={"margin": "0.5rem 0 0.25rem 0"}),

                # Checklist arranged into responsive "columns" by wrapping labels
                dcc.Checklist(
                    id="cases",
                    options=[{"label": name, "value": name} for name in ALL_CASES],
                    value=ALL_CASES[:2] if ALL_CASES else [],
                    # Make the label elements flow like tiles so they wrap into columns
                    style={
                        "display": "flex",
                        "flexWrap": "wrap",
                        "gap": "0.35rem 1.25rem",
                        "rowGap": "0.4rem",
                        # limit overall width so items form 3–5 columns depending on viewport
                        "maxWidth": "1200px",
                    },
                    # One style applied to every option label; fixed width makes columns
                    labelStyle={
                        "display": "inline-block",
                        "width": "230px",   # tweak: smaller -> more columns
                        "margin": "0",
                    },
                    inputStyle={"marginRight": "0.5rem"},
                ),

                html.Div(
                    [
                        html.Label("τ (for spectral plots)", style={"display": "block", "marginTop": "0.75rem"}),
                        dcc.Slider(
                            id="tau_target",
                            min=0.0,
                            max=2.0,
                            step=0.1,
                            value=0.5,
                            marks={0: "0", 0.5: "0.5", 1: "1", 2: "2"},
                        ),
                    ],
                    # Constrain slider width so it doesn't look silly on very wide pages
                    style={
                        "maxWidth": "480px",
                        "margin": "0.25rem 0 0.25rem 0",
                    },
                ),

                html.Div(id="status", style={"marginTop": "0.5rem", "fontSize": "0.9rem", "color": "#666"}),
                html.Hr(style={"margin": "1rem 0 0.75rem 0"}),
            ],
            # Center controls and give them breathing room
            style={
                "padding": "1rem",
                "margin": "0 auto",
                "maxWidth": "1400px",  # keeps the controls from stretching edge-to-edge
                "boxSizing": "border-box",
            },
        ),

        # ------ Plots (stacked vertically; each centered and width-limited) ------
        html.Div(
            [
                dcc.Graph(
                    id="spectral_budget",
                    style={
                        "height": "52vh",
                        "margin": "0 auto 1rem auto",
                        "maxWidth": "1400px",
                    },
                ),
                dcc.Graph(
                    id="optical_props",
                    style={
                        "height": "45vh",
                        "margin": "0 auto 1rem auto",
                        "maxWidth": "1400px",
                    },
                ),
                dcc.Graph(
                    id="ts_vs_tau",
                    style={
                        "height": "60vh",
                        "margin": "0 auto 1rem auto",
                        "maxWidth": "1400px",
                    },
                ),
            ],
            style={"padding": "0 1rem 1rem", "boxSizing": "border-box"},
        ),
    ],
    # Give the whole page a comfortable width and prevent horizontal scrollbars
    style={"maxWidth": "100%", "overflowX": "hidden"},
)

# ------- Callbacks -------
@app.callback(
    Output("spectral_budget","figure"),
    Output("optical_props","figure"),
    Output("ts_vs_tau","figure"),
    Output("status","children"),
    Input("cases","value"),
    Input("tau_target","value"),
)
def render(cases_selected, tau_target):
    bases = cases_selected or []
    static_files = [STATIC_MAP[b] for b in bases if b in STATIC_MAP]
    output_files = [OUTPUT_MAP[b] for b in bases if b in OUTPUT_MAP]

    fig_spec = spectral_budget_figure(static_files, float(tau_target))
    fig_opt  = optical_props_figure(static_files, float(tau_target))
    fig_ts   = ts_vs_tau_figure(output_files)
    status   = f"{len(bases)} case(s) • spectral={len(static_files)} • temperature={len(output_files)}"
    return fig_spec, fig_opt, fig_ts, status

if __name__ == "__main__":
    app.run(debug=True)
