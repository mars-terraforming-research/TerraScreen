# terrascreen_dash.py
import json
import math
import itertools
from functools import lru_cache

import numpy as np
import io
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

DARK_BG = "#1f222a"
DARK_FG = "#e6e6e6"

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

def _apply_dark_theme(fig):
    fig.update_layout(
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",   # transparent paper -> page background shows
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(color=DARK_FG, family='Helvetica'),
    )
    # Optional: slightly subtler gridlines
    fig.update_xaxes(gridcolor="rgba(255,255,255,0.08)")
    fig.update_yaxes(gridcolor="rgba(255,255,255,0.08)")
    return fig


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

def read_output_from_content(file_content, tau_target=None):
    '''
    Modified version of terrascreen_lib.read_output() that works with string content
    instead of file paths (for dashboard use with GitHub fetching)
    '''
    # Convert string content to file-like object
    f = io.StringIO(file_content)
    
    # Rest is identical to original read_output function
    Nlev = 24
    L_NSPECTI = 96
    L_NSPECTV = 84
    
    part_name = f.readline().split(',')[1]
    Ncase = int(f.readline().split(',')[1])
    dt = float(f.readline().split(',')[1])
    Qext670 = float(f.readline().split(',')[1])
    alb_sfc = float(f.readline().split(',')[1])
    conrath_nu = float(f.readline().split(',')[1])
    solar_flux = float(f.readline().split(',')[1])
    bwni = np.array([float(x) for x in f.readline().split(',')[1:]])
    bwnv = np.array([float(x) for x in f.readline().split(',')[1:]])
    solar_WN = np.array([float(x) for x in f.readline().split(',')[1:]])
    f.readline()  # skip ====
    wl = np.array([float(x) for x in f.readline().split(',')[1:]])
    Qext = np.array([float(x) for x in f.readline().split(',')[1:]])
    Qscat = np.array([float(x) for x in f.readline().split(',')[1:]])
    g = np.array([float(x) for x in f.readline().split(',')[1:]])
    f.readline()  # skip ====
    Pref = np.array([float(x) for x in f.readline().split(',')[1:]])
    Tref = np.array([float(x) for x in f.readline().split(',')[1:]])
    [f.readline() for _ in range(2)]  # skip 2 line
    
    # Read data arrays
    T = np.zeros((Ncase, Nlev))
    OLR_WL, SFC_dIR = [np.zeros((Ncase, L_NSPECTI)) for _ in range(2)]
    ASR_WL, SFC_dVIS = [np.zeros((Ncase, L_NSPECTV)) for _ in range(2)]
    tau, ts, net_top, net_bot, alb, OLR, ASR = [np.zeros((Ncase)) for _ in range(7)]
    
    for i in range(Ncase):
        vals = np.array([float(x) for x in f.readline().split(',')])
        tau[i] = vals[1]
        ts[i] = vals[2]
        alb[i] = vals[3]
        OLR[i] = vals[4]
        ASR[i] = vals[5]
        net_top[i] = vals[6]
        net_bot[i] = vals[7]
        T[i, :] = vals[8:8+Nlev]
        OLR_WL[i, :] = vals[9+Nlev:9+Nlev+L_NSPECTI]
        ASR_WL[i, :] = vals[9+Nlev+L_NSPECTI:9+Nlev+L_NSPECTI+L_NSPECTV]
        SFC_dIR[i, :] = vals[9+Nlev+L_NSPECTI+L_NSPECTV:9+Nlev+2*L_NSPECTI+L_NSPECTV]
        SFC_dVIS[i, :] = vals[9+Nlev+2*L_NSPECTI+L_NSPECTV:]
    f.close()

    # Compute center wavenumber
    res_wn_bwni = bwni[1:] - bwni[0:-1]
    center_wn_bwni = bwni[0:-1] + res_wn_bwni/2
    res_wn_bwnv = bwnv[1:] - bwnv[0:-1]
    center_wn_bwnv = bwnv[0:-1] + res_wn_bwnv/2
    
    # Initialize model
    class model(object):
        pass
    MOD = model()
    
    # Save all the static variables
    setattr(MOD, 'name', 'from_content')
    setattr(MOD, 'Qext670', Qext670)
    setattr(MOD, 'alb_sfc', alb_sfc)
    setattr(MOD, 'conrath_nu', conrath_nu)
    setattr(MOD, 'solar_flux', solar_flux)
    setattr(MOD, 'bwni', bwni)
    setattr(MOD, 'bwnv', bwnv)
    setattr(MOD, 'solar_WN', solar_WN)
    setattr(MOD, 'wni', center_wn_bwni)
    setattr(MOD, 'wnv', center_wn_bwnv)
    setattr(MOD, 'wli', 10**4/center_wn_bwni)
    setattr(MOD, 'wlv', 10**4/center_wn_bwnv)
    setattr(MOD, 'Qext', Qext)
    setattr(MOD, 'Qscat', Qscat)
    setattr(MOD, 'g', g)
    setattr(MOD, 'wl', wl)
    setattr(MOD, 'Tref', Tref)
    setattr(MOD, 'Pref', Pref)

    # Save variables
    var_list = ['tau', 'ts', 'alb', 'OLR', 'ASR', 'net_top', 'net_bot', 'OLR_WL', 'ASR_WL', 'SFC_dIR', 'SFC_dVIS', 'T']
    for ivar in var_list:
        if tau_target is None:
            setattr(MOD, ivar, eval(ivar))
        else:
            itau = np.argmin(np.abs(tau - tau_target))
            setattr(MOD, ivar, eval('%s[%i,...]' % (ivar, itau)))
    return MOD

def create_step_trace_original(bounds_wn, var, name, color, show_sides=True):
    """
    Exact copy of the working step trace function from original code.
    This creates the proper wavelength bins using center +/- half-width method.
    """
    x_vals, y_vals = [], []
    
    res_wn = bounds_wn[1:] - bounds_wn[:-1]
    center_wn = bounds_wn[:-1] + res_wn/2
    
    for i, val in enumerate(var):
        left = center_wn[i] - res_wn[i]/2
        right = center_wn[i] + res_wn[i]/2
        
        if show_sides:
            x_vals.extend([left, left, right, right])
            y_vals.extend([0, val, val, 0])
        else:
            x_vals.extend([left, right])
            y_vals.extend([val, val])
        
        if i < len(var) - 1:
            x_vals.append(None)  # Break line
            y_vals.append(None)
    
    return go.Scatter(x=x_vals, y=y_vals, name=name, 
                     line=dict(color=color, width=2), mode='lines',
                     hovertemplate="λ: %{x:.4g} μm<br>Value: %{y:.4g}<extra>"+name+"</extra>",
                     connectgaps=False)

def spectral_budget_figure(cases, tau_target: float):
    """
    CORRECTED VERSION - Uses read_output logic with dashboard fetching.
    This should eliminate all artifacts by using the exact same data parsing
    and plotting approach as the working original code.
    
    REPLACE your existing spectral_budget_figure function with this one.
    """
    fact_plot = 20  # IR scaling factor (matches original)
    fig = go.Figure()

    try:
        # Fetch and parse clear reference case using read_output logic
        clear_content = fetch_text_cached(RAW_BASE + "static_dust1.5.txt")
        clear = read_output_from_content(clear_content, 0.)
        
        # EXACT same wavelength calculation as working original
        bnwlv = 10**4 / clear.bwnv[::-1]
        bnwli = 10**4 / clear.bwni[::-1] 
        all_wl_bounds = np.append(bnwlv, bnwli[1:])
        res_wl = all_wl_bounds[1:] - all_wl_bounds[:-1]
        res_wl_vis = bnwlv[1:] - bnwlv[:-1]

        # Solar spectrum (exact same as original)
        if hasattr(clear, 'solar_WN') and clear.solar_WN is not None and len(clear.solar_WN):
            solar_trace = create_step_trace_original(
                bnwlv, clear.solar_WN[::-1]/res_wl_vis, 
                'solar flux', '#ffffff', show_sides=False
            )
            fig.add_trace(solar_trace)

        # Clear reference (exact same as original)
        net_clear = np.append(clear.ASR_WL[::-1], -fact_plot * clear.OLR_WL[::-1])
        clear_trace = create_step_trace_original(
            all_wl_bounds, net_clear/res_wl,
            'clear (no aerosol)', '#f6a21a', show_sides=True
        )
        fig.add_trace(clear_trace)

        # Selected cases - use same read_output logic
        for fname in (cases or []):
            try:
                case_content = fetch_text_cached(RAW_BASE + fname)
                case = read_output_from_content(case_content, tau_target)
                label = _basename(fname).replace("_", " ")
                color = color_for(_basename(fname))

                # Exact same stacking as original
                net_case = np.append(case.ASR_WL[::-1], -fact_plot * case.OLR_WL[::-1])
                case_trace = create_step_trace_original(
                    all_wl_bounds, net_case/res_wl,
                    label, color, show_sides=False
                )
                fig.add_trace(case_trace)
            except Exception as e:
                print(f"Error loading {fname}: {e}")
                continue

    except Exception as e:
        print(f"Error loading reference: {e}")
        return px.line(title="Could not load data")

    # Same layout as before
    ticks = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5,2,3,4,5,6,7,8,9,10,12,14,16,20,25,30,40,50,60]
    fig.update_xaxes(
        type="log", range=[np.log10(0.2), np.log10(60.0)],
        tickmode="array", tickvals=ticks, ticktext=[str(t) for t in ticks],
        title="Wavelength [μm]", showgrid=True
    )
    fig.update_yaxes(
        title=f"Irradiance [W/m²/μm] (×{fact_plot} scaled in IR)",
        showgrid=True
    )
    fig.update_layout(
        title=f"Radiative budget at the top of the atmosphere (τ = {float(tau_target):.2f})",
        margin=dict(l=8, r=8, t=80, b=60),
        legend=dict(orientation="h", x=0, y=-0.25, xanchor="left", yanchor="top"),
    )

    return _apply_dark_theme(fig)

def optical_props_figure(cases, tau_target: float, show_scattering: bool):
    """
    Bottom plot: Extinction (solid) and optional scattering (dashed),
    log-scale Y with safe handling of zeros and length mismatches.
    """
    if not cases:
        return px.line(title="Select aerosols that have spectral data")

    fig = go.Figure()
    ymin, ymax = np.inf, -np.inf

    for name in cases:
        case = read_static_case(name, round(float(tau_target), 3))

        base = _basename(name)  # key for color and legend group
        color = color_for(base)

        # Ensure arrays are consistent and sorted by wavelength
        x = np.asarray(case.wl,   dtype=float)
        y_ext = np.asarray(case.Qext, dtype=float)
        y_scat = np.asarray(case.Qscat, dtype=float)

        n = int(min(len(x), len(y_ext), len(y_scat)))
        if n < 2:
            continue
        x, y_ext, y_scat = x[:n], y_ext[:n], y_scat[:n]

        order = np.argsort(x)
        x, y_ext, y_scat = x[order], y_ext[order], y_scat[order]

        # Safe values for log scale
        y_ext_safe = y_ext.copy()
        y_ext_safe[y_ext_safe <= 0] = np.nan

        y_scat_safe = y_scat.copy()
        y_scat_safe[y_scat_safe <= 0] = np.nan

        # Skip if no positive values
        if np.all(np.isnan(y_ext_safe)) and (not show_scattering or np.all(np.isnan(y_scat_safe))):
            continue

        # Update y-range based on what we will plot
        if np.any(~np.isnan(y_ext_safe)):
            ymin = min(ymin, np.nanmin(y_ext_safe))
            ymax = max(ymax, np.nanmax(y_ext_safe))
        if show_scattering and np.any(~np.isnan(y_scat_safe)):
            ymin = min(ymin, np.nanmin(y_scat_safe))
            ymax = max(ymax, np.nanmax(y_scat_safe))

        # Extinction (solid)
        fig.add_trace(go.Scatter(
            x=x, y=y_ext_safe, mode="lines",
            name=base.replace("_", " ") + " — extinction",
            line=dict(color=color, width=2, dash="solid"),
            legendgroup=base,
        ))

        # Scattering (dashed) - optional
        if show_scattering and np.any(~np.isnan(y_scat_safe)):
            fig.add_trace(go.Scatter(
                x=x, y=y_scat_safe, mode="lines",
                name=base.replace("_", " ") + " — scattering",
                line=dict(color=color, width=2, dash="dash"),
                legendgroup=base,
            ))

    if not np.isfinite(ymin) or not np.isfinite(ymax):
        return px.line(title="No positive values to plot on log scale")

    lo = max(ymin / 1.5, 1e-6)
    hi = ymax * 1.5

    fig.update_layout(
        title="Optical properties for the aerosols considered (Extinction & Scattering)",
        margin=dict(l=4, r=4, t=80, b=80),
        legend=dict(orientation="h", y=-0.25, yanchor="top", x=0, xanchor="left"),
        xaxis_title="Wavelength [μm]",
        yaxis_title="Efficiency",
    )
    ticks = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5,2,3,4,5,6,7,8,9,10,12,14,16,20,25,30,40,50,60]
    fig.update_xaxes(
        type="log", range=[np.log10(0.2), np.log10(60.0)],
        tickmode="array", tickvals=ticks, ticktext=[str(t) for t in ticks],
        title="Wavelength [μm]", showgrid=True
    )
    # fig.update_xaxes(type="log", range=[math.log10(0.2), math.log10(60.0)])
    fig.update_yaxes(type="log", range=[math.log10(lo), math.log10(hi)])

    return _apply_dark_theme(fig)


def ts_vs_tau_figure(output_files):
    if not output_files:
        return px.line(title="Select aerosols that have temperature data (Ts vs τ)")

    frames = []
    for fname in output_files:
        df = read_tau_ts(fname)
        if not df.empty:
            frames.append(df)

    if not frames:
        return px.line(title="No temperature data")

    df = pd.concat(frames, ignore_index=True)

    # color_discrete_map must be keyed by the 'series' values (which are base names)
    series_keys = df["series"].dropna().unique().tolist()
    color_map = {k: color_for(k) for k in series_keys}

    fig = px.line(
        df, x="tau", y="ts", color="series",
        color_discrete_map=color_map,
        labels={"tau": "Visible opacity τ (670 nm)", "ts": "Surface temperature Ts [K]"},
        title="Ts vs τ"
    )
    fig.update_layout(margin=dict(l=0, r=0, t=60, b=0),
                      legend=dict(orientation="h", y=-0.25, yanchor="top", x=0, xanchor="left"))
    return _apply_dark_theme(fig)


# ------- Build unified case list (no 'static'/'output' labels) -------
STATIC_FILES = list_repo_files("static_")
OUTPUT_FILES = list_repo_files("output_")

STATIC_MAP = { _basename(f): f for f in STATIC_FILES }
OUTPUT_MAP = { _basename(f): f for f in OUTPUT_FILES }
ALL_CASES  = sorted(set(STATIC_MAP) | set(OUTPUT_MAP))  # union
BASES = sorted(ALL_CASES)

# A long, stable palette; repeat if needed
COLOR_POOL = (
    px.colors.qualitative.Plotly
    + px.colors.qualitative.D3
    + px.colors.qualitative.Safe
    + px.colors.qualitative.Set3
)

COLOR_BY_BASE = {base: COLOR_POOL[i % len(COLOR_POOL)] for i, base in enumerate(BASES)}

def color_for(base_name: str) -> str:
    """Return the stable color for a given base case (e.g., 'Al_8um_60um')."""
    return COLOR_BY_BASE.get(base_name, "#7f7f7f")

# ------- Layout (flexbox with resizable/collapsible sides) -------
app = Dash(__name__)

# Wrap everything in a dark container
app.layout = html.Div(
    [
        html.Div(
            [
                html.H2("TerraScreen Dashboard", style={"margin": "0 0 0.5rem 0"}),
                html.H4("Candidate Particles", style={"margin": "0.5rem 0 0.25rem 0"}),
                dcc.Checklist(
                    id="cases",
                    options=[{"label": name, "value": name} for name in ALL_CASES],
                    value=ALL_CASES[:2] if ALL_CASES else [],
                    style={
                        "display": "flex",
                        "flexWrap": "wrap",
                        "gap": "0.35rem 1.25rem",
                        "rowGap": "0.4rem",
                        "maxWidth": "1200px",
                    },
                    labelStyle={"display": "inline-block", "width": "230px", "margin": "0"},
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
                    style={"maxWidth": "480px", "margin": "0.25rem 0 0.25rem 0"},
                ),
                html.Div(
                    [
                        dcc.Checklist(
                            id="scattering_toggle",
                            options=[{"label": "Show scattering (dashed)", "value": "scat"}],
                            value=["scat"],  # default ON; set [] if you prefer OFF by default
                            style={"marginTop": "0.5rem"},
                        ),
                    ],
                    style={"maxWidth": "480px", "margin": "0.25rem 0 0.25rem 0"},
                ),
                html.Div(id="status", style={"marginTop": "0.5rem", "fontSize": "0.9rem", "color": "#B8BDC9"}),
                html.Hr(style={"margin": "1rem 0 0.75rem 0", "borderColor": "#2a2f3a"}),
            ],
            style={"padding": "1rem", "margin": "0 auto", "maxWidth": "1400px", "boxSizing": "border-box"},
        ),

        html.Div(
            [
                dcc.Graph(
                    id="spectral_budget",
                    style={"height": "52vh", "margin": "0 auto 1rem auto", "maxWidth": "1400px"},
                ),
                dcc.Graph(
                    id="optical_props",
                    style={"height": "45vh", "margin": "0 auto 1rem auto", "maxWidth": "1400px"},
                ),
                dcc.Graph(
                    id="ts_vs_tau",
                    style={"height": "60vh", "margin": "0 auto 1rem auto", "maxWidth": "1400px"},
                ),
            ],
            style={"padding": "0 1rem 1rem", "boxSizing": "border-box"},
        ),

        # ----- Top-right image overlay (non-clickable, tooltip attribution) -----
        html.Img(
            src="assets/terraformed_mars.png", 
            alt="Terraformed Mars (illustration)",
            title="Image credit: Science.org — 'Terraforming Mars could be easier than scientists thought'",
            style={
                "position": "fixed",
                "top": "24px",
                "right": "36px",
                "width": "180px",    
                "height": "auto",
                "opacity": "0.85",
                "borderRadius": "8px",
                # "boxShadow": "0 2px 12px rgba(0,0,0,0.6)",
                "zIndex": 1000,
                "pointerEvents": "auto",   # allow hover tooltip
            },
        ),
    ],
    # Dark page background + default text color
    style={"backgroundColor": DARK_BG, "color": DARK_FG, "minHeight": "100vh", "maxWidth": "100%", "overflowX": "hidden", "fontFamily": "Helvetica"},
)

# ------- Callbacks -------
@app.callback(
    Output("spectral_budget","figure"),
    Output("optical_props","figure"),
    Output("ts_vs_tau","figure"),
    Output("status","children"),
    Input("cases","value"),
    Input("tau_target","value"),
    Input("scattering_toggle","value"),
)
def render(cases_selected, tau_target, scattering_toggle):
    bases = cases_selected or []
    static_files = [STATIC_MAP[b] for b in bases if b in STATIC_MAP]
    output_files = [OUTPUT_MAP[b] for b in bases if b in OUTPUT_MAP]

    show_scattering = "scat" in (scattering_toggle or [])

    fig_spec = spectral_budget_figure(static_files, float(tau_target))
    fig_opt  = optical_props_figure(static_files, float(tau_target), show_scattering)
    fig_ts   = ts_vs_tau_figure(output_files)
    status   = f"{len(bases)} case(s) • spectral={len(static_files)} • temperature={len(output_files)}"
    return fig_spec, fig_opt, fig_ts, status

if __name__ == "__main__":
    app.run(debug=True)
