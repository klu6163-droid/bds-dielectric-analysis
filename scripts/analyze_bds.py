"""Read-only BDS inputs -> separately saved, explicitly provisional analysis.

Python 3.10+; numpy, scipy, matplotlib. See references/input-and-run.md.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import re
import shutil
import sys
import zipfile

import numpy as np
from scipy.signal import find_peaks, peak_widths
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import bds_enhancements as enhancements
import result_lock

VERSION = "2.4"
EPS0 = 8.8541878128e-12
WINDOWS = (7, 9, 13)
DEFAULTS = dict(min_span_decades=1.0, min_points=5, max_abs_slope=0.15,
                max_ratio=1.5, max_log_residual=0.05, min_relative_prominence=0.15,
                min_width_decades=0.15, min_height_fraction=0.20,
                peak_match_decades=0.16, edge_warning_decades=0.5,
                max_frequency_gap_decades=0.15)
COLMAP = {"Cp": "Cp_F", "D": "tan_delta", "εr": "epsilon_r_export",
          "εr'": "epsilon_prime_export", "εr''": "epsilon_doubleprime_export"}
PEAK_COLUMNS = ["sample", "temperature_C", "candidate_id", "frequency_Hz",
                "height", "frequency_min_Hz", "frequency_max_Hz",
                "tau_peak_s_NOT_assigned_tau_s", "low_edge_warning", "high_edge_warning"]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def resolve(base, path):
    p = Path(path)
    return (p if p.is_absolute() else base / p).resolve()


def write_json(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False),
                    encoding="utf-8")


def write_csv(path, rows, fallback=()):
    fields = list(dict.fromkeys(k for r in rows for k in r)) or list(fallback)
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: "" if isinstance(v, (float, np.floating)) and not np.isfinite(v)
                             else v for k, v in row.items()})


def geometry_value(text, unit, scale):
    if text is None:
        return None
    match = re.fullmatch(r"\s*([0-9.eE+-]+)\s*" + unit + r"\s*", text)
    require(match is not None, f"Unknown geometry unit: {text}")
    value = float(match[1]) * scale
    require(np.isfinite(value) and value > 0, f"Invalid geometry: {text}")
    return value


def read_fsm(path):
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    meta, blocks, temperatures = {}, {}, None
    for i, line in enumerate(lines):
        if not line.lstrip().startswith("---"):
            if not blocks:
                parts = line.split()
                if len(parts) == 2:
                    meta[parts[0]] = parts[1]
            continue
        name = line.strip().strip("-").strip().replace("′", "'").replace("″", "''")
        require(name not in blocks, f"Duplicate section: {name}")
        header = lines[i + 1]
        require("频率(Hz)" in header, f"Frequency unit not verified in {path}")
        temps = [float(t) for t in re.findall(r"温度_([+-]?\d+(?:\.\d+)?)℃", header)]
        require(temps and len(set(temps)) == len(temps), "Missing or duplicate temperature columns")
        if temperatures is None:
            temperatures = temps
        require(temps == temperatures, "Temperature columns differ between sections")
        rows = []
        for following in lines[i + 2:]:
            if not following.strip() or following.lstrip().startswith("---"):
                break
            parts = following.split()
            require(len(parts) == len(temps) + 2 and parts[0].isdigit(),
                    f"Malformed row in {path}, section {name}: {following[:100]}")
            rows.append([float(v) for v in parts])
        require(rows, f"Empty section {name}")
        blocks[name] = np.asarray(rows)
    require(blocks, "No FSM data sections")
    base = blocks.get("εr'", next(iter(blocks.values())))
    for name, values in blocks.items():
        require(values.shape == base.shape and np.array_equal(values[:, :2], base[:, :2]),
                f"Index/frequency mismatch in section {name}")
        require(np.isfinite(values[:, :2]).all(), f"Nonfinite index/frequency in section {name}")
    records = []
    for j, temperature in enumerate(temperatures):
        for i in range(len(base)):
            record = dict(temperature_C=temperature, point_index=int(base[i, 0]),
                          frequency_Hz=float(base[i, 1]))
            for key, arr in blocks.items():
                record[COLMAP.get(key, f"source_section_{key}")] = float(arr[i, j + 2])
            records.append(record)
    return meta, records


def read_long_csv(path, sample):
    numeric = set(COLMAP.values()) | {"temperature_C", "frequency_Hz"}
    with path.open(encoding="utf-8-sig", newline="") as stream:
        source = list(csv.DictReader(stream))
    require(source, "Empty CSV")
    require({"temperature_C", "frequency_Hz"} <= set(source[0]), "Missing frequency/temperature columns")
    records = []
    for i, row in enumerate(source):
        if "sample" in row:
            require(row["sample"] == sample, f"CSV sample differs from mapping at row {i + 2}")
        record = {k: float(v) for k, v in row.items() if k in numeric and v is not None and v.strip()}
        require({"temperature_C", "frequency_Hz"} <= set(record), f"Missing frequency/temperature at row {i + 2}")
        record["point_index"] = i
        record.update({f"source_{k}": v for k, v in row.items() if k not in numeric})
        records.append(record)
    return {"format": "normalized long CSV"}, records


def load_datasets(config, base):
    loaded, samples = [], set()
    require(config.get("datasets"), "No datasets")
    for spec in config["datasets"]:
        sample = str(spec["sample"])
        require(sample and sample not in samples, "Sample labels must be unique per configured input")
        samples.add(sample)
        path = resolve(base, spec["path"])
        fmt = spec.get("format", "fsm_txt")
        if fmt == "fsm_txt":
            meta, rows = read_fsm(path)
        elif fmt == "long_csv":
            meta, rows = read_long_csv(path, sample)
        else:
            raise ValueError(f"Unsupported format {fmt}; inspect and adapt explicitly")
        if "expected_source_sample" in spec:
            require(meta.get("样品名称") == spec["expected_source_sample"], "Source sample mapping mismatch")
        area = geometry_value(meta.get("样品面积"), r"mm(?:2|²)", 1e-6)
        thickness = geometry_value(meta.get("样品厚度"), "mm", 1e-3)
        for key, old in (("area_m2", area), ("thickness_m", thickness)):
            if key in spec and old is not None and not np.isclose(float(spec[key]), old, rtol=1e-8, atol=0):
                require(spec.get("geometry_note"), "Geometry override requires geometry_note")
        area = float(spec["area_m2"]) if "area_m2" in spec else area
        thickness = float(spec["thickness_m"]) if "thickness_m" in spec else thickness
        for v in (area, thickness):
            require(v is None or (np.isfinite(v) and v > 0), "Invalid geometry")
        temps = sorted(set(r["temperature_C"] for r in rows))
        require(all(np.isfinite(temps)) and min(temps) > -273.15, "Invalid Celsius temperature")
        groups = []
        for temp in temps:
            part = [r for r in rows if r["temperature_C"] == temp]
            f = np.asarray([r["frequency_Hz"] for r in part])
            require(len(f) >= 3 and np.isfinite(f).all() and (f > 0).all(), "Need >=3 finite positive frequencies")
            require((np.diff(f) > 0).all(), "Frequency must be strictly increasing; no silent sorting/deduplication")
            groups.append(dict(sample=sample, temperature_C=temp, f=f, source_rows=part))
        loaded.append(dict(sample=sample, path=path, metadata=meta, area_m2=area,
                           thickness_m=thickness, groups=groups, specification=spec))
    return loaded


def convert_group(group, dataset, mode="exported_epsilon"):
    require(mode == "exported_epsilon", "Only original epsilon columns are supported")
    rows, f = group["source_rows"], group["f"]
    require(all("epsilon_prime_export" in r and "epsilon_doubleprime_export" in r for r in rows),
            "Direct epsilon analysis needs both original epsilon columns; no reconstruction fallback")
    ep = np.asarray([r["epsilon_prime_export"] for r in rows])
    el = np.asarray([r["epsilon_doubleprime_export"] for r in rows])
    for name, values in (("epsilon_prime_export", ep), ("epsilon_doubleprime_export", el)):
        bad = np.flatnonzero(~np.isfinite(values) | (values <= 0))
        require(not len(bad), f"{group['sample']} {group['temperature_C']} C: invalid {name}, indices {bad[:8]}; no silent filtering")
    group.update(ep=ep, el=el, sigma=2*np.pi*f*EPS0*el/100,
        area_m2=dataset["area_m2"], thickness_m=dataset["thickness_m"],
        conversion_audit=dict(input_mapping="epsilon_r_prime -> epsilon_prime; epsilon_r_doubleprime -> epsilon_doubleprime",
                              direct_values_unchanged=True))
    return group


def local_derivative(x, values, window):
    out = np.full(len(x), np.nan)
    half = window // 2
    for i in range(half, len(x) - half):
        sl = slice(i-half, i+half+1)
        u = x[sl] - x[i]
        scale = max(abs(u))
        coeff = np.polynomial.polynomial.polyfit(u / scale, values[sl], 2)
        out[i] = -np.pi / 2 * coeff[1] / scale
    return out


def numeric_peaks(f, values, window, p):
    ids = np.flatnonzero(np.isfinite(values))
    if len(ids) < 3:
        return []
    # Span windows can be unavailable inside a large sampling gap. Never bridge it.
    splits = np.split(ids, np.flatnonzero(np.diff(ids) > 1)+1)
    if len(splits) > 1:
        combined = []
        for segment in splits:
            if len(segment) < 3:
                continue
            a, b = int(segment[0]), int(segment[-1])+1
            part = numeric_peaks(f[a:b], values[a:b], window, p)
            for record in part:
                record["index"] += a
            combined.extend(part)
        scale = max([v["height"] for v in combined if v["height"] > 0], default=np.inf)
        for record in combined:
            if record["height"] < p["min_height_fraction"]*scale:
                record["shape_screen_passed"] = False
                if "small" not in record["rejection_reasons"]:
                    record["rejection_reasons"] += ";small"
        return combined
    v = values[ids]
    peaks, props = find_peaks(v, prominence=0)
    if not len(peaks):
        return []
    _, _, left, right = peak_widths(v, peaks, rel_height=0.5)
    scale = max([v[k] for k in peaks if v[k] > 0], default=np.inf)
    records = []
    for j, k in enumerate(peaks):
        width = float(np.interp(right[j], np.arange(len(ids)), np.log10(f[ids])) -
                      np.interp(left[j], np.arange(len(ids)), np.log10(f[ids])))
        prominence = float(props["prominences"][j])
        rel = prominence / v[k] if v[k] > 0 else 0
        reasons = []
        if v[k] <= 0: reasons.append("nonpositive")
        if rel < p["min_relative_prominence"]: reasons.append("low_prominence")
        if width < p["min_width_decades"]: reasons.append("narrow")
        if v[k] < p["min_height_fraction"] * scale: reasons.append("small")
        records.append(dict(window_id=str(window), window_points=window if isinstance(window,int) else None,
                            index=int(ids[k]), frequency_Hz=float(f[ids[k]]),
                            height=float(v[k]), prominence=prominence, width_decades=width,
                            shape_screen_passed=not reasons, rejection_reasons=";".join(reasons)))
    return records


def stable_peaks(f, derivatives, p, primary=9):
    all_peaks = {w: numeric_peaks(f, derivatives[w], w, p) for w in derivatives}
    selected = []
    for peak in all_peaks[primary]:
        if not peak["shape_screen_passed"]:
            continue
        matches = {primary: peak}
        for w in (w for w in derivatives if w != primary):
            possible = [v for v in all_peaks[w] if v["shape_screen_passed"] and
                        abs(np.log10(v["frequency_Hz"]/peak["frequency_Hz"])) <= p["peak_match_decades"]]
            if possible:
                matches[w] = min(possible, key=lambda v: abs(np.log10(v["frequency_Hz"]/peak["frequency_Hz"])))
        if len(matches) != len(derivatives):
            continue
        frequencies = [matches[w]["frequency_Hz"] for w in derivatives]
        selected.append(dict(candidate_id=f"P{len(selected)+1}", index=peak["index"],
            frequency_Hz=peak["frequency_Hz"], height=peak["height"],
            primary_window=str(primary),
            frequencies_by_window=json.dumps({str(w):matches[w]["frequency_Hz"] for w in derivatives}),
            frequency_min_Hz=min(frequencies), frequency_max_Hz=max(frequencies),
            tau_peak_s_NOT_assigned_tau_s=float(1/(2*np.pi*peak["frequency_Hz"])),
            low_edge_warning=bool(np.log10(peak["frequency_Hz"]/f[0]) < p["edge_warning_decades"]),
            high_edge_warning=bool(np.log10(f[-1]/peak["frequency_Hz"]) < p["edge_warning_decades"]),
            physical_assignment="Unassigned; not confirmed segmental alpha relaxation"))
    return selected, [v for values in all_peaks.values() for v in values]


def interval_stats(f, sig, el, i, j):
    x = np.log10(f[i:j+1]); xc = x-x.mean()
    y = np.log10(sig[i:j+1]); e = np.log10(el[i:j+1])
    slope = float(np.sum(xc*(y-y.mean()))/np.sum(xc*xc))
    slope_e = float(np.sum(xc*(e-e.mean()))/np.sum(xc*xc))
    require(abs(slope_e-(slope-1)) < 1e-9, "Conductivity/loss slope identity failed")
    return dict(i_start=i, i_end=j, n_points=j-i+1, span_decades=float(x[-1]-x[0]),
                f_start_Hz=float(f[i]), f_end_Hz=float(f[j]), beta_sigma=slope,
                beta_epsilon_loss=slope_e, sigma_max_min_ratio=float(10**(y.max()-y.min())),
                sigma_geometric_mean_S_cm=float(10**y.mean()),
                largest_frequency_gap_decades=float(np.diff(x).max()),
                max_log10_residual=float(np.max(abs(y-(y.mean()+slope*xc)))))


def passes(interval, p, slope=None, ratio=None):
    return (abs(interval["beta_sigma"]) <= (p["max_abs_slope"] if slope is None else slope)
            and interval["sigma_max_min_ratio"] <= (p["max_ratio"] if ratio is None else ratio)
            and interval["max_log10_residual"] <= p["max_log_residual"])


def choose(intervals, p, slope=None, ratio=None):
    candidates = [v for v in intervals if passes(v, p, slope, ratio)]
    return max(candidates, key=lambda v: (v["span_decades"], -abs(v["beta_sigma"]))) if candidates else None


def conduction_scan(f, sig, el, p):
    intervals, shortest = [], []
    x = np.log10(f)
    gaps = np.diff(x)
    gap_limit = p["max_frequency_gap_decades"]
    gap_indices = np.flatnonzero(gaps > gap_limit)
    gap_audit = [dict(left_index=int(k), right_index=int(k+1),
                      left_frequency_Hz=float(f[k]), right_frequency_Hz=float(f[k+1]),
                      gap_decades=float(gaps[k]), limit_decades=float(gap_limit),
                      action="split_scan_no_cross_gap_window") for k in gap_indices]
    segments = np.r_[0, np.cumsum(gaps > gap_limit)]
    for i in range(len(f)-1):
        first = True
        for j in range(i + p["min_points"] - 1, len(f)):
            if segments[j] != segments[i]:
                break
            if x[j]-x[i] < p["min_span_decades"]:
                continue
            a = interval_stats(f, sig, el, i, j)
            intervals.append(a)
            if first:
                shortest.append(a)
                first = False
    picked = choose(intervals, p)
    accepted = picked is not None
    if picked is None and shortest:
        picked = min(shortest, key=lambda a: (abs(a["beta_sigma"]), a["sigma_max_min_ratio"]))
    if picked is None:
        status = "Insufficient_contiguous_frequency_span_or_points" if len(gap_indices) else "Insufficient_frequency_span_or_points"
        result = dict(status=status, numerical_candidate=False)
    else:
        result = dict(picked, numerical_candidate=accepted,
            status="Candidate_IC_not_confirmed" if accepted else "Rejected_comparison_not_IC",
            touches_low_frequency_limit=picked["i_start"] == 0,
            touches_high_frequency_limit=picked["i_end"] == len(f)-1)
        i,j = picked["i_start"],picked["i_end"]
        result["touches_internal_frequency_gap"] = bool((i > 0 and gaps[i-1] > gap_limit) or
                                                        (j < len(f)-1 and gaps[j] > gap_limit))
    result.update(frequency_gap_count=len(gap_audit), max_frequency_gap_decades=gap_limit,
                  gap_audit=gap_audit)
    sensitivity = []
    for slope in (.10, .15, .20):
        for ratio in (1.25, 1.5, 2.0):
            q = choose(intervals, p, slope, ratio)
            sensitivity.append(dict(slope_limit=slope, ratio_limit=ratio, candidate_found=q is not None,
                f_start_Hz=q["f_start_Hz"] if q else None, f_end_Hz=q["f_end_Hz"] if q else None))
    return result, shortest, sensitivity


def analyze_group(group, p, derivative_config=None):
    x, ep = np.log(group["f"]), group["ep"]
    group["raw_der"] = -np.pi/2*np.gradient(ep, x, edge_order=2)
    config = enhancements.derivative_settings(derivative_config or {"mode":"legacy_points"})
    ds, primary, diagnostics, label = enhancements.build_derivatives(group["f"], ep, config, local_derivative)
    group.update(derivatives=ds, primary_window=primary, derivative_diagnostics=diagnostics,
                 derivative_config=config, derivative_label=label)
    group["peaks"], group["all_peaks"] = stable_peaks(group["f"], ds, p, primary)
    group["IC"], group["scan"], group["sensitivity"] = conduction_scan(group["f"], group["sigma"], group["el"], p)
    group["window_spans_decades"] = {str(w): float(np.median([r["actual_span_decades"] for r in records if r["status"]=="ok"]))
        if any(r["status"]=="ok" for r in records) else None for w,records in diagnostics.items()}
    return group


def freq(value):
    return f"{value/1000:.3g} kHz" if value >= 1000 else f"{value:.3g} Hz"


def common_axes(ax, g):
    ax.set_xscale("log")
    ax.set_xlim(g["f"][0]/1.15, g["f"][-1]*1.15)
    ax.set_xlabel("Frequency (Hz)")
    ax.grid(which="major", color="#DEE3E8", lw=.6)
    ax.tick_params(which="both", direction="in", top=True, right=True)


def readable_signed_ticks(ax, threshold):
    lo, hi = ax.get_ylim()
    first = int(np.ceil(np.log10(threshold*2)))
    last = int(np.ceil(np.log10(max(abs(lo), abs(hi), threshold*2))))
    positive = [10.0**k for k in range(first, last+1)]
    ticks = sorted([-v for v in positive if lo <= -v] + [0] + [v for v in positive if v <= hi])
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{v:g}" for v in ticks])


def derivative_panel(ax, g, color):
    f = g["f"]
    ax.plot(f, g["raw_der"], color="#AFB7C1", lw=.7, alpha=.7)
    primary = g["primary_window"]
    secondary = [w for w in g["derivatives"] if w != primary]
    for k,w in enumerate(secondary):
        ax.plot(f,g["derivatives"][w],color=color,lw=.8,ls="--" if k%2==0 else ":",alpha=.5)
    ax.plot(f,g["derivatives"][primary],color=color,lw=1.5)
    amplitude = max(float(np.max(abs(g["raw_der"]))), 1e-6)
    threshold = max(min(amplitude*.005, .5), 1e-8)
    ax.set_yscale("symlog", linthresh=threshold, linscale=.5)
    # Reserve positive space above the curves for candidate labels.
    ax.set_ylim(min(float(g["raw_der"].min())*1.3, -threshold),
                max(float(g["raw_der"].max())*5, threshold*10))
    readable_signed_ticks(ax, threshold)
    ax.axhline(0, lw=.6, ls="--", color="#7C8996")
    ax.set_ylabel(r"$\varepsilon''_{\mathrm{der}}$")
    common_axes(ax, g)
    if not g["peaks"]:
        ax.text(.98, .96, "No stable candidate peak", transform=ax.transAxes, ha="right", va="top", fontsize=8)
    for i, peak in enumerate(g["peaks"]):
        fp, yp = peak["frequency_Hz"], peak["height"]
        ax.plot(fp, yp, "D", ms=5, mfc="white", mec="#26313D", zorder=8)
        ax.plot([peak["frequency_min_Hz"], peak["frequency_max_Hz"]], [yp, yp], color="#26313D", lw=1)
        label = f"{peak['candidate_id']}: {freq(fp)}\nWindow range: {freq(peak['frequency_min_Hz'])}–{freq(peak['frequency_max_Hz'])}"
        if peak["low_edge_warning"] or peak["high_edge_warning"]:
            label += "\nMeasurement-edge warning"
        # Full labels for the first three by frequency; all other peaks retain IDs.
        if i < 3:
            ax.annotate(label, (fp, yp), xytext=(.42, .97-i*.28), textcoords="axes fraction",
                fontsize=7.4, ha="left", va="top",
                bbox=dict(boxstyle="round,pad=.2", fc="white", ec="#D1D7DD", alpha=.95),
                arrowprops=dict(arrowstyle="->", color="#657585", lw=.7))
        else:
            ax.annotate(peak["candidate_id"], (fp, yp), xytext=(4, 6),
                        textcoords="offset points", fontsize=7)
    ax.set_title("Derivative spectrum", loc="left", fontsize=10, weight="bold")


def conduction_panel(ax, g, color, key):
    y, f, region = g[key], g["f"], g["IC"]
    ax.plot(f, y, color=color, lw=1.4, marker="o", ms=1.7, mew=0)
    ax.set_yscale("log")
    common_axes(ax, g)
    ax.set_ylabel(r"$\varepsilon''$" if key == "el" else r"$\sigma'$ (S cm$^{-1}$)")
    ax.set_title("Dielectric loss" if key == "el" else "AC conductivity",
                 loc="left", fontsize=10, weight="bold")
    # Annotation below data; no hidden low-frequency endpoints.
    ax.set_ylim(float(y.min())/8, float(y.max())*2)
    if "i_start" not in region:
        label = "No IC window: insufficient span / points"
    else:
        i, j = region["i_start"], region["i_end"]
        ff, yy = f[i:j+1], y[i:j+1]
        colorbox = "#CC7A00" if region["numerical_candidate"] else "#747E87"
        y0, y1 = float(yy.min()/1.15), float(yy.max()*1.15)
        ax.add_patch(Rectangle((ff[0], y0), ff[-1]-ff[0], y1-y0,
                               facecolor=colorbox, edgecolor="none", alpha=.10))
        ax.add_patch(Rectangle((ff[0], y0), ff[-1]-ff[0], y1-y0,
                               facecolor="none", edgecolor=colorbox, lw=1.3, ls="--"))
        value = region["sigma_geometric_mean_S_cm"]
        guide = value*100/(2*np.pi*EPS0*ff) if key == "el" else np.full(len(ff), value)
        ax.plot(ff, guide, ":", color="#27323C", lw=1)
        label = "Candidate Ionic Conduction" if region["numerical_candidate"] else "No clear IC; gray = rejected comparison"
        label += f"\n{freq(ff[0])}–{freq(ff[-1])}; βσ={region['beta_sigma']:+.3f}"
        label += f"; βε″={region['beta_epsilon_loss']:+.3f}"
        if region["touches_low_frequency_limit"] or region["touches_high_frequency_limit"]:
            label += "\nTouches measurement limit"
    ax.text(.02, .025, label, transform=ax.transAxes, fontsize=7.2, va="bottom",
            bbox=dict(boxstyle="round,pad=.2", fc="white", ec="#D1D7DD", alpha=.96))


def safe_filename(value):
    name = re.sub(r'[<>:"/\\|?*\x00-\x1f]', "_", str(value)).strip(". ")
    return name or "sample"


def save_figure(fig, out, name, formats):
    for ext in formats:
        fig.savefig(out/f"{name}.{ext}", dpi=180)
    plt.close(fig)


def plot_results(groups, out, formats, mode):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.linewidth": 1,
                         "svg.fonttype": "none", "pdf.fonttype": 42})
    temps = sorted(set(g["temperature_C"] for g in groups))
    cmap = plt.get_cmap("viridis")
    colors = {t: cmap(.12+.75*i/max(1, len(temps)-1)) for i, t in enumerate(temps)}
    samples = list(dict.fromkeys(g["sample"] for g in groups))
    note = "Original relative permittivity real and loss values; used directly"
    def rows_figure(rows, name, title):
        fig, axes = plt.subplots(len(rows), 3, figsize=(15.5, 3.55*len(rows)+1.1), squeeze=False)
        for i, g in enumerate(rows):
            derivative_panel(axes[i, 0], g, colors[g["temperature_C"]])
            conduction_panel(axes[i, 1], g, colors[g["temperature_C"]], "el")
            conduction_panel(axes[i, 2], g, colors[g["temperature_C"]], "sigma")
            for ax in axes[i]:
                ax.set_title(f"{g['sample']} | {g['temperature_C']:g} °C — "+ax.get_title(loc="left"), loc="left", fontsize=9.7)
        fig.suptitle(title, fontsize=15, weight="bold", y=.99)
        fig.text(.5, .99-.28/fig.get_figheight(), note, ha="center", fontsize=8.5, color="#52606D")
        fig.text(.5, .017, "Diamonds: unassigned peaks. Gray: direct derivative; colored: "+rows[0]["derivative_label"]+". Orange: IC candidate; gray box: rejected.",
                 ha="center", fontsize=8)
        fig.subplots_adjust(left=.06, right=.985, bottom=.65/fig.get_figheight(),
                            top=1-.9/fig.get_figheight(), hspace=.5, wspace=.3)
        save_figure(fig, out, name, formats)
    for index, sample in enumerate(samples, 1):
        rows = [g for g in groups if g["sample"] == sample]
        stem = f"{index:02d}_{safe_filename(sample)}"
        rows_figure(rows, stem+"_annotated", f"{sample} | Candidate spectra at all temperatures")
        fig, axes = plt.subplots(2, 2, figsize=(11, 8.4))
        for key, ax, title, ylabel in [
            ("ep", axes[0, 0], "Real permittivity", r"$\varepsilon'$"),
            ("el", axes[0, 1], "Dielectric loss", r"$\varepsilon''$"),
            ("sigma", axes[1, 0], "AC conductivity", r"$\sigma'$ (S cm$^{-1}$)"),
            ("raw_der", axes[1, 1], "Derivative spectrum", r"$\varepsilon''_{\mathrm{der}}$")]:
            for g in rows:
                y = g["derivatives"][g["primary_window"]] if key == "raw_der" else g[key]
                ax.plot(g["f"], y, color=colors[g["temperature_C"]], label=f"{g['temperature_C']:g} °C", lw=1.3)
                if key == "raw_der":
                    ax.plot(g["f"], g["raw_der"], color=colors[g["temperature_C"]], lw=.6, alpha=.25)
            ax.set_xscale("log")
            if key == "raw_der":
                ax.set_yscale("symlog", linthresh=.5, linscale=.5)
                ax.axhline(0, ls="--", color="#89939C", lw=.6)
                readable_signed_ticks(ax, .5)
            else:
                ax.set_yscale("log")
            ax.set(title=title, xlabel="Frequency (Hz)", ylabel=ylabel)
            ax.grid(which="major", color="#DEE3E8", lw=.6)
            ax.legend(frameon=False, fontsize=8)
        fig.suptitle(f"{sample} | BDS spectra", fontsize=16, weight="bold")
        fig.text(.5, .94, note, ha="center", fontsize=8.5)
        fig.text(.5, .02, "Derivative: direct differences + "+enhancements.window_label(rows[0]["primary_window"])+" local quadratic; signed-log axis. Original values retained.", ha="center", fontsize=8.5)
        fig.subplots_adjust(left=.09, right=.98, bottom=.1, top=.87, wspace=.28, hspace=.35)
        save_figure(fig, out, stem+"_four_spectra", formats)
    for t in temps:
        rows_figure([g for g in groups if g["temperature_C"] == t],
                    f"comparison_{safe_filename(f'{t:g}')}C", f"BDS comparison | {t:g} °C")


def deliver_tables(groups, out):
    peaks, regions, all_peaks, scans, sensitivity, points, audit, window_records = [], [], [], [], [], [], [], []
    gap_records = []
    for g in groups:
        identity = dict(sample=g["sample"], temperature_C=g["temperature_C"],conversion_branch=g.get("conversion_branch","unspecified"))
        peaks += [dict(identity, **v) for v in g["peaks"]]
        regions.append(dict(identity, **{k:v for k,v in g["IC"].items() if k != "gap_audit"}))
        gap_records += [dict(identity, **v) for v in g["IC"]["gap_audit"]]
        all_peaks += [dict(identity, **v) for v in g["all_peaks"]]
        scans += [dict(identity, **v) for v in g["scan"]]
        sensitivity += [dict(identity, **v) for v in g["sensitivity"]]
        audit.append(dict(identity, conversion_checks=g["conversion_audit"],
                          derivative_window_span_decades=g["window_spans_decades"]))
        for w,records in g["derivative_diagnostics"].items():
            window_records += [dict(identity,window_id=str(w),**r) for r in records]
        peak_ids = {v["index"]: v["candidate_id"] for v in g["peaks"]}
        ic = g["IC"]
        for i, source in enumerate(g["source_rows"]):
            marked = "i_start" in ic and ic["i_start"] <= i <= ic["i_end"]
            row = dict(source, **{"sample":g["sample"], "area_m2":g["area_m2"], "thickness_m":g["thickness_m"]})
            row.update(epsilon_prime_used=float(g["ep"][i]), epsilon_doubleprime_used=float(g["el"][i]),
                       sigma_prime_S_m=float(g["sigma"][i]*100), sigma_prime_S_cm=float(g["sigma"][i]),
                       epsilon_derivative_raw=float(g["raw_der"][i]),
                       epsilon_derivative_primary=float(g["derivatives"][g["primary_window"]][i]),
                       conversion_branch=g.get("conversion_branch","unspecified"),
                       primary_derivative_window=str(g["primary_window"]), candidate_peak_id=peak_ids.get(i, ""),
                       candidate_track_id=next((p.get("track_id","") for p in g["peaks"] if p["index"]==i),""),
                       primary_window_status=ic["status"], inside_marked_window=marked,
                       inside_numerical_IC_candidate=bool(marked and ic["numerical_candidate"]))
            row.update({f"epsilon_derivative_local{w}": float(g["derivatives"][w][i]) for w in g["derivatives"]})
            points.append(row)
    write_csv(out/"all_data.csv", points)
    write_csv(out/"candidate_peaks.csv", peaks, PEAK_COLUMNS)
    write_csv(out/"candidate_IC_regions.csv", regions)
    write_csv(out/"IC_frequency_gap_audit.csv", gap_records,
              ["sample","temperature_C","left_frequency_Hz","right_frequency_Hz","gap_decades","limit_decades","action"])
    write_csv(out/"all_numerical_peaks.csv", all_peaks, ["sample", "temperature_C", "frequency_Hz", "rejection_reasons"])
    write_csv(out/"shortest_windows_full_frequency_scan.csv", scans, ["sample", "temperature_C", "f_start_Hz", "f_end_Hz"])
    write_csv(out/"threshold_sensitivity.csv", sensitivity)
    write_json(out/"conversion_audit.json", audit)
    write_csv(out/"derivative_window_support.csv",window_records)
    return peaks, regions, len(points)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--inspect", action="store_true")
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = json.loads(config_path.read_text(encoding="utf-8-sig"))
    base = config_path.parent
    mode = config.get("conversion", "exported_epsilon")
    require(mode == "exported_epsilon", "Only original epsilon columns are supported; update the legacy configuration")
    retired = {"compare_conversions", "cpd_definition_confirmed", "cpd_definition_evidence", "user_selection"} & set(config)
    require(not retired, f"Retired two-batch configuration keys: {sorted(retired)}; use the direct-epsilon configuration")
    loaded = load_datasets(config, base)
    summary = [dict(sample=d["sample"], source=str(d["path"]), metadata=d["metadata"],
                    area_m2=d["area_m2"], thickness_m=d["thickness_m"],
                    temperatures=[dict(temperature_C=g["temperature_C"], points=len(g["f"]),
                                       frequency_min_Hz=float(g["f"][0]), frequency_max_Hz=float(g["f"][-1]))
                                  for g in d["groups"]]) for d in loaded]
    if args.inspect:
        print(json.dumps(summary, ensure_ascii=False, indent=2))
        return
    require(args.out is not None, "--out required for analysis")
    out = args.out.resolve()
    require(not out.exists(), "Output path already exists; choose a new directory")
    derivative_config = enhancements.derivative_settings(config.get("derivative",{}))
    tracking_config = enhancements.track_settings(config.get("tracking",{}))
    supplied = config.get("screening", {})
    require(not (set(supplied)-set(DEFAULTS)), f"Unknown screening keys: {set(supplied)-set(DEFAULTS)}")
    params = dict(DEFAULTS, **supplied)
    require(all(isinstance(v, (float, int)) and np.isfinite(v) and v > 0 for v in params.values()), "Invalid threshold")
    require(isinstance(params["min_points"], int) and params["min_points"] >= 3, "min_points must be integer >=3")
    require(params["max_ratio"] >= 1, "max_ratio must be >=1")
    formats = config.get("formats", ["png", "pdf", "svg"])
    require(formats and set(formats) <= {"png", "pdf", "svg"}, "Allowed formats: png/pdf/svg")
    protected = {d["path"] for d in loaded} | {config_path, Path(__file__).resolve(),Path(enhancements.__file__).resolve(),Path(result_lock.__file__).resolve()}
    for raw in config.get("protect_paths", []):
        item = resolve(base, raw)
        require(item.exists(), f"Protected path missing: {item}")
        require(not out.is_relative_to(item), "Output must be outside protected source directories")
        protected.update(p for p in item.rglob("*") if p.is_file()) if item.is_dir() else protected.add(item)
    before = {str(p): sha(p) for p in sorted(protected)}
    groups = []
    for dataset in loaded:
        for group in dataset["groups"]:
            group.update(conversion_mode=mode,conversion_branch=enhancements.branch_name(mode))
            groups.append(analyze_group(convert_group(group,dataset,mode),params,derivative_config))
    observations,track_summary,track_links = enhancements.temperature_tracks(groups,tracking_config)
    out.mkdir(parents=True, exist_ok=False)
    # Persist exact inputs/configuration and executable, without copying raw data.
    recorded = dict(config, datasets=[dict(d["specification"], path=str(d["path"])) for d in loaded],
                    protect_paths=[str(resolve(base, p)) for p in config.get("protect_paths", [])])
    recorded.update(conversion=mode,skill_version=VERSION,derivative=derivative_config,
                    tracking=tracking_config,screening=params)
    write_json(out/"run_config.json", recorded)
    write_json(out/"input_summary.json", summary)
    shutil.copy2(__file__, out/"analyze_bds.py")
    shutil.copy2(enhancements.__file__,out/"bds_enhancements.py")
    shutil.copy2(result_lock.__file__,out/"result_lock.py")
    peaks, regions, count = deliver_tables(groups, out)
    write_csv(out/"candidate_peak_tracks.csv",observations,["sample","conversion_branch","track_id","temperature_C","frequency_Hz"])
    write_csv(out/"candidate_track_summary.csv",track_summary,["sample","conversion_branch","track_id","n_temperatures"])
    write_csv(out/"tracking_link_audit.csv",track_links,["sample","conversion_branch","status"])
    plot_results(groups, out, formats, mode)
    enhancements.plot_tracks(groups,observations,out,formats,save_figure,safe_filename)
    note = f"""# BDS 候选谱分析

本批直接使用原始文件εr′作为ε′、εr″作为ε″，不对这两列再次归一化。
电导率σ′=2πfε0ε″（S/m），除以100为S/cm；导数谱由原始ε′求得。
共 {len(loaded)} 个样品、{len(groups)} 条温度谱、{count} 个数据点。原始点不删除、不重采样。
导数公式：ε″der=−(π/2)dε′/dlnω；直接差分与{groups[0]['derivative_label']}局部二次求导均保存。
菱形为未归属候选峰；窗口间峰位范围不是置信区间。τ_peak不等于已确定的τs。
橙框为数值近平台候选；灰框为未通过的最近对照区。两图的斜率关系是同一换算的结果，不是独立证据。
框接触测量边界时不能声称过程已完全分离。无候选不等于无离子导电。
电导率窗口几何均值仅为代表交流值，未指定最终σdc；未计算D、离子浓度、解离率或σTτs。

## 本次筛查参数（探索性规则，不是文献标准）

{json.dumps(params, ensure_ascii=False, indent=2)}

## 候选峰

|样品|温度/℃|候选|频率/Hz|窗口范围/Hz|边界警示|
|---|---:|---|---:|---|---|
"""
    for q in peaks:
        note += f"|{q['sample']}|{q['temperature_C']:g}|{q['candidate_id']}|{q['frequency_Hz']:.6g}|{q['frequency_min_Hz']:.6g}–{q['frequency_max_Hz']:.6g}|{q['low_edge_warning'] or q['high_edge_warning']}|\n"
    note += "\n## 导电区筛查\n\n|样品|温度/℃|状态|频率/Hz|βσ|\n|---|---:|---|---|---:|\n"
    for q in regions:
        span = f"{q['f_start_Hz']:.6g}–{q['f_end_Hz']:.6g}" if "f_start_Hz" in q else "不足"
        slope = f"{q['beta_sigma']:.4f}" if "beta_sigma" in q else "—"
        note += f"|{q['sample']}|{q['temperature_C']:g}|{q['status']}|{span}|{slope}|\n"
    note += "\n详见 conversion_audit.json、candidate_peaks.csv、candidate_IC_regions.csv 和 threshold_sensitivity.csv。PNG供检查，SVG保留文字和路径，PDF为矢量。\n"
    note += "\n## 跨温度候选峰跟踪\n\n"
    note += f"版本 {VERSION}。所有谱和候选均基于原始文件中的εr′、εr″。\n\n"
    note += "跟踪仅在同一样品内进行，跨样品轨迹编号不表示同一物理过程。相近匹配存在歧义时断开，缺峰不插值、不跨越中间空温度。\n\n"
    note += "|样品|换算|轨迹|温度数|起止温度/℃|边缘警示数|\n|---|---|---|---:|---|---:|\n"
    for tr in track_summary:
        note += f"|{tr['sample']}|{tr['conversion_branch']}|{tr['track_id']}|{tr['n_temperatures']}|{tr['first_temperature_C']:g}–{tr['last_temperature_C']:g}|{tr['boundary_warning_count']}|\n"
    note += "\n轨迹连线为候选关联，不构成共同弛豫机制或α峰归属的证据；峰位误差条为窗口敏感性，不是置信区间。全部配置、每点窗口支持情况及歧义链接均可查表。\n"
    note += "\n导电窗口不跨越超过所设上限的频率缺口，详见IC_frequency_gap_audit.csv。后续取值和重绘须核验result_manifest.json，不混用其他结果版本。\n"
    (out/"分析说明.md").write_text(note, encoding="utf-8")
    require(all(sha(Path(p)) == digest for p, digest in before.items()), "Protected input changed during analysis")
    with (out/"all_data.csv").open(encoding="utf-8-sig", newline="") as stream:
        require(sum(1 for _ in csv.DictReader(stream)) == count, "Saved row count mismatch")
    write_json(out/"verification.json", dict(protected_files_unchanged=True, protected_sha256=before,
        source_points_removed=0, output_rows=count, candidate_peaks=len(peaks),
        numerical_IC_candidates=sum(r["numerical_candidate"] for r in regions),
        epsilon0_F_m=EPS0, conversion=mode, screening=params,skill_version=VERSION,
        derivative=derivative_config,tracking=tracking_config,
        tracked_sequences=len(track_summary),tracked_candidate_observations=len(observations),
        versions=dict(python=sys.version.split()[0], numpy=np.__version__,
                      scipy=__import__("scipy").__version__, matplotlib=matplotlib.__version__)))
    science = dict(purpose="fresh_analysis", skill_version=VERSION, conversion=mode,
                   column_mapping={"epsilon_prime":"epsilon_prime_export", "epsilon_doubleprime":"epsilon_doubleprime_export"},
                   derivative=derivative_config, screening=params, tracking=tracking_config,
                   samples=[dict(sample=d["sample"], source=str(d["path"]),
                                 temperatures_C=[g["temperature_C"] for g in d["groups"]]) for d in loaded],
                   units={"frequency":"Hz", "temperature_input":"degC", "product_temperature":"K",
                          "epsilon":"relative_dimensionless", "sigma":"S/cm", "tau":"s"},
                   tau_convention="tau_peak=1/(2*pi*f_peak); physical assignment unconfirmed")
    payload = [p.name for p in out.iterdir() if p.is_file() and p.suffix in {".csv",".json",".py",".md"}]
    analysis_sources = {str(p):before[str(p)] for p in ({d["path"] for d in loaded} | {config_path})}
    manifest = result_lock.seal_run(out, science, analysis_sources, payload)
    result_lock.verify_run(out)
    archive = out/"BDS_results.zip"
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for file in sorted(out.rglob("*")):
            if file.is_file() and file != archive:
                z.write(file, file.relative_to(out))
    with zipfile.ZipFile(archive) as z:
        require(z.testzip() is None, "ZIP CRC failure")
    print(json.dumps(dict(output=str(out), samples=len(loaded), temperature_spectra=len(groups),
                         rows=count, candidate_peaks=len(peaks),
                         numerical_IC_candidates=sum(r["numerical_candidate"] for r in regions),
                         candidate_sequences=len(track_summary),
                         analysis_id=manifest["analysis_id"],
                         protected_files_unchanged=True), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    try:
        main()
    except (ValueError, KeyError, OSError) as exc:
        print(f"BDS analysis stopped: {exc}", file=sys.stderr)
        raise SystemExit(2)
