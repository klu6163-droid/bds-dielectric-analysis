"""BDS: frequency-span derivatives and provisional temperature tracks."""
import copy
import numpy as np
import matplotlib.pyplot as plt

DERIVATIVE_DEFAULTS = dict(mode="frequency_span", spans_decades=[0.25, 0.35, 0.50],
    primary_span_decades=0.35, windows_points=[7, 9, 13], primary_window_points=9,
    min_points=5, min_points_each_side=2, max_gap_decades=0.15)
TRACK_DEFAULTS = dict(enabled=True, max_shift_decades_per_C=0.12,
    max_temperature_gap_C=10.0, ambiguity_margin_decades=0.10, use_prediction=True)


def checked_config(user, defaults, name):
    unknown = set(user)-set(defaults)
    if unknown:
        raise ValueError(f"Unknown {name} keys: {unknown}")
    return dict(copy.deepcopy(defaults), **user)


def derivative_settings(user):
    c = checked_config(user, DERIVATIVE_DEFAULTS, "derivative")
    if c["mode"] not in ("frequency_span", "legacy_points"):
        raise ValueError("derivative.mode must be frequency_span or legacy_points")
    windows = c["windows_points"]
    if not isinstance(windows, list) or len(windows) < 2:
        raise ValueError("Need >=2 distinct point windows for sensitivity screening")
    if not all(type(w) is int and w >= 5 and w % 2 == 1 for w in windows):
        raise ValueError("Point windows must be odd integers >=5")
    if len(set(windows)) != len(windows):
        raise ValueError("Point windows must be distinct")
    if type(c["primary_window_points"]) is not int or c["primary_window_points"] not in windows:
        raise ValueError("primary_window_points must occur in windows_points")
    spans = c["spans_decades"]
    if not isinstance(spans, list) or len(spans) < 2 or len(set(spans)) != len(spans):
        raise ValueError("Need >=2 distinct derivative sensitivity spans")
    if not all(isinstance(s, (float, int)) and np.isfinite(s) and s > 0 for s in spans):
        raise ValueError("Derivative spans must be finite and positive")
    if c["primary_span_decades"] not in spans:
        raise ValueError("Primary span must occur in spans_decades")
    for k, minimum in (("min_points", 5), ("min_points_each_side", 1)):
        if not isinstance(c[k], int) or c[k] < minimum:
            raise ValueError(f"Invalid derivative {k}")
    if not np.isfinite(c["max_gap_decades"]) or c["max_gap_decades"] <= 0:
        raise ValueError("Invalid derivative max_gap_decades")
    return c


def span_id(span):
    return f"span_{span:.8g}dec"


def span_derivative(f, ep, span, config):
    """Fit inside a centered FULL log10-frequency width, without resampling."""
    x = np.log10(f)
    values = np.full(len(f), np.nan)
    diagnostics = []
    radius = span/2
    for i in range(len(f)):
        ids = np.flatnonzero((x >= x[i]-radius-1e-12) & (x <= x[i]+radius+1e-12))
        left = int(np.sum(ids < i)); right = int(np.sum(ids > i))
        actual = float(x[ids[-1]]-x[ids[0]])
        status = "ok"
        if x[i]-radius < x[0]-1e-12 or x[i]+radius > x[-1]+1e-12:
            status = "measurement_edge"
        elif len(ids) < config["min_points"] or min(left, right) < config["min_points_each_side"]:
            status = "insufficient_support"
        elif np.max(np.diff(x[ids])) > config["max_gap_decades"]:
            status = "frequency_gap"
        if status == "ok":
            # Convert actual log10 frequencies to ln coordinates before differentiation.
            u = np.log(f[ids])-np.log(f[i])
            scale = max(abs(u))
            coef = np.polynomial.polynomial.polyfit(u/scale, ep[ids], 2)
            values[i] = -np.pi/2*coef[1]/scale
        diagnostics.append(dict(index=i, frequency_Hz=float(f[i]), target_span_decades=float(span),
            points_used=len(ids), actual_span_decades=actual,
            support_start_Hz=float(f[ids[0]]), support_end_Hz=float(f[ids[-1]]), status=status))
    return values, diagnostics


def build_derivatives(f, ep, config, legacy_function):
    if config["mode"] == "legacy_points":
        ds = {w: legacy_function(np.log(f), ep, w) for w in sorted(config["windows_points"])}
        diagnostics = {}
        for w in ds:
            records = []
            for i in range(len(f)):
                left, right = max(0, i-w//2), min(len(f)-1, i+w//2)
                status = "ok" if np.isfinite(ds[w][i]) else "measurement_edge"
                if status == "ok" and np.max(np.diff(np.log10(f[left:right+1]))) > config["max_gap_decades"]:
                    ds[w][i] = np.nan
                    status = "frequency_gap"
                records.append(dict(index=i, frequency_Hz=float(f[i]), target_span_decades=None,
                    points_used=right-left+1, actual_span_decades=float(np.log10(f[right]/f[left])),
                    support_start_Hz=float(f[left]), support_end_Hz=float(f[right]),
                    status=status))
            diagnostics[w] = records
        primary = config["primary_window_points"]
        return ds, primary, diagnostics, "/".join(str(w) for w in sorted(ds))+f"-point (primary {primary})"
    ds, diagnostics = {}, {}
    for span in sorted(config["spans_decades"]):
        key = span_id(span)
        ds[key], diagnostics[key] = span_derivative(f, ep, span, config)
    return ds, span_id(config["primary_span_decades"]), diagnostics, "/".join(f"{v:g}" for v in sorted(config["spans_decades"]))+" decade"


def window_label(key):
    return f"{key} pt" if isinstance(key, int) else key.replace("span_", "").replace("dec", " dec")


def track_settings(user):
    c = checked_config(user, TRACK_DEFAULTS, "tracking")
    if not isinstance(c["enabled"], bool) or not isinstance(c["use_prediction"], bool):
        raise ValueError("Tracking enabled/use_prediction must be boolean")
    for k in ("max_shift_decades_per_C", "max_temperature_gap_C", "ambiguity_margin_decades"):
        if not isinstance(c[k], (int, float)) or not np.isfinite(c[k]) or c[k] <= 0:
            raise ValueError(f"Invalid tracking {k}")
    return c


def branch_name(mode):
    if mode != "exported_epsilon":
        raise ValueError("Only original epsilon columns are supported")
    return "exported"


def temperature_tracks(groups, config):
    """Conservative mutual-best links. Ambiguous links are recorded, never forced."""
    observations, summaries, links = [], [], []
    if not config["enabled"]:
        return observations, summaries, links
    keys = sorted(set((g["sample"], g["conversion_branch"]) for g in groups))
    for sample, branch in keys:
        sequence = sorted([g for g in groups if (g["sample"], g["conversion_branch"]) == (sample, branch)],
                          key=lambda g: g["temperature_C"])
        tracks, previous = [], []
        for g in sequence:
            temp = g["temperature_C"]
            current = g["peaks"]
            possible = []
            # previous is reset at every measured temperature, including spectra with no peaks.
            for tr in previous:
                last = tr[-1]
                dt = temp-last["temperature_C"]
                if dt <= 0 or dt > config["max_temperature_gap_C"]:
                    continue
                lastlog = np.log10(last["frequency_Hz"])
                pred = lastlog
                if config["use_prediction"] and len(tr) >= 2:
                    older = tr[-2]
                    slope = (lastlog-np.log10(older["frequency_Hz"]))/(last["temperature_C"]-older["temperature_C"])
                    pred += slope*dt
                gate = config["max_shift_decades_per_C"]*dt
                for j, peak in enumerate(current):
                    now = np.log10(peak["frequency_Hz"])
                    eligible = abs(now-lastlog) <= gate and abs(now-pred) <= gate
                    record = dict(sample=sample, conversion_branch=branch, from_track_id=last["track_id"],
                        from_temperature_C=last["temperature_C"], to_temperature_C=temp,
                        from_frequency_Hz=last["frequency_Hz"], to_candidate_id=peak["candidate_id"],
                        to_frequency_Hz=peak["frequency_Hz"], predicted_log10_frequency=float(pred),
                        cost_decades=float(abs(now-pred)), eligible=bool(eligible),
                        status="eligible_not_selected" if eligible else "outside_shift_gate")
                    links.append(record)
                    if eligible:
                        possible.append((tr, j, record))
            ambiguous_tracks, ambiguous_peaks = set(), set()
            for tr in previous:
                opts = sorted([v for v in possible if v[0] is tr], key=lambda v: v[2]["cost_decades"])
                if len(opts) > 1 and opts[1][2]["cost_decades"]-opts[0][2]["cost_decades"] <= config["ambiguity_margin_decades"]:
                    ambiguous_tracks.add(tr[-1]["track_id"])
                    ambiguous_peaks.update(v[1] for v in opts if v[2]["cost_decades"] <= opts[0][2]["cost_decades"]+config["ambiguity_margin_decades"])
            for j in range(len(current)):
                opts = sorted([v for v in possible if v[1] == j], key=lambda v: v[2]["cost_decades"])
                if len(opts) > 1 and opts[1][2]["cost_decades"]-opts[0][2]["cost_decades"] <= config["ambiguity_margin_decades"]:
                    ambiguous_peaks.add(j)
                    ambiguous_tracks.update(v[0][-1]["track_id"] for v in opts if v[2]["cost_decades"] <= opts[0][2]["cost_decades"]+config["ambiguity_margin_decades"])
            matches = {}
            for tr, j, record in possible:
                tid = tr[-1]["track_id"]
                if tid in ambiguous_tracks or j in ambiguous_peaks:
                    record["status"] = "ambiguous_not_linked"
                    continue
                rowbest = min([v for v in possible if v[0] is tr], key=lambda v: v[2]["cost_decades"])
                colbest = min([v for v in possible if v[1] == j], key=lambda v: v[2]["cost_decades"])
                if rowbest[1] == j and colbest[0] is tr:
                    matches[j] = tr
                    record["status"] = "provisional_link"
            active = []
            for j, peak in enumerate(current):
                if j in matches:
                    tr = matches[j]
                    tid, status = tr[-1]["track_id"], "provisional_continuation"
                else:
                    tid = f"{sample}-{branch}-T{len(tracks)+1:02d}"
                    tr = []
                    tracks.append(tr)
                    status = "ambiguous_new_track" if j in ambiguous_peaks else "new_track"
                peak["track_id"], peak["track_link_status"] = tid, status
                ob = dict(sample=sample, conversion_branch=branch, track_id=tid,
                    temperature_C=temp, candidate_id=peak["candidate_id"], frequency_Hz=peak["frequency_Hz"],
                    frequency_min_Hz=peak["frequency_min_Hz"], frequency_max_Hz=peak["frequency_max_Hz"],
                    tau_peak_s_NOT_assigned_tau_s=peak["tau_peak_s_NOT_assigned_tau_s"],
                    low_edge_warning=peak["low_edge_warning"], high_edge_warning=peak["high_edge_warning"],
                    link_status=status, physical_assignment="Unassigned candidate sequence")
                tr.append(ob); observations.append(ob); active.append(tr)
            previous = active
        for tr in tracks:
            summaries.append(dict(sample=sample, conversion_branch=branch, track_id=tr[0]["track_id"],
                n_temperatures=len(tr), first_temperature_C=tr[0]["temperature_C"], last_temperature_C=tr[-1]["temperature_C"],
                first_frequency_Hz=tr[0]["frequency_Hz"], last_frequency_Hz=tr[-1]["frequency_Hz"],
                covers_all_measured_temperatures=len(tr)==len(sequence),
                boundary_warning_count=sum(v["low_edge_warning"] or v["high_edge_warning"] for v in tr),
                sequence_status="continuous_candidate_sequence" if len(tr)>=3 else "short_candidate_sequence",
                physical_assignment="Unassigned; tracking is not proof of a common physical mechanism"))
    return observations, summaries, links


def plot_tracks(groups, observations, out, formats, save, safe):
    if not observations:
        return
    samples=list(dict.fromkeys(g["sample"] for g in groups))
    branches=list(dict.fromkeys(g["conversion_branch"] for g in groups))
    fig,axes=plt.subplots(len(samples),len(branches),figsize=(6.3*len(branches),4.1*len(samples)+.6),squeeze=False,sharey=True,sharex=True)
    all_temps=sorted(set(g["temperature_C"] for g in groups))
    padding=max((max(all_temps)-min(all_temps))*.08,1)
    minfreq=min(v["frequency_min_Hz"] for v in observations)
    maxfreq=max(v["frequency_max_Hz"] for v in observations)
    for r,sample in enumerate(samples):
        for c,branch in enumerate(branches):
            ax=axes[r,c]
            rows=[v for v in observations if v["sample"]==sample and v["conversion_branch"]==branch]
            for k,tid in enumerate(dict.fromkeys(v["track_id"] for v in rows)):
                seq=sorted([v for v in rows if v["track_id"]==tid],key=lambda v:v["temperature_C"])
                color=plt.get_cmap("tab10")(k%10)
                x=np.array([v["temperature_C"] for v in seq]); y=np.array([v["frequency_Hz"] for v in seq])
                low=y-np.array([v["frequency_min_Hz"] for v in seq]);high=np.array([v["frequency_max_Hz"] for v in seq])-y
                ax.errorbar(x,y,yerr=[low,high],marker="o",ls="--",lw=1.2,capsize=3,color=color,label=tid)
                for v in seq:
                    if v["low_edge_warning"] or v["high_edge_warning"]:
                        ax.plot(v["temperature_C"],v["frequency_Hz"],"x",color="#C63736",ms=10,mew=1.5)
            if not rows:
                ax.text(.5,.5,"No stable candidate peaks",transform=ax.transAxes,ha="center")
            ax.set_yscale("log");ax.set_xlabel("Temperature (°C)");ax.set_ylabel("Candidate peak frequency (Hz)")
            ax.set_xlim(min(all_temps)-padding,max(all_temps)+padding)
            ax.set_ylim(minfreq*.6,maxfreq*1.6)
            if len(all_temps)<=10:ax.set_xticks(all_temps)
            measured=[g["temperature_C"] for g in groups if g["sample"]==sample and g["conversion_branch"]==branch]
            missing=[t for t in measured if not any(v["temperature_C"]==t for v in rows)]
            if missing:
                ax.text(.02,.03,"No stable candidate: "+", ".join(f"{t:g}" for t in missing)+" °C",
                        transform=ax.transAxes,fontsize=8,color="#725139")
            ax.set_title(f"{sample} | {branch}",loc="left",weight="bold")
            ax.grid(color="#DFE5E9",lw=.6)
            if rows:ax.legend(frameon=False,fontsize=8)
    fig.suptitle("Provisional peak sequences across temperature",fontsize=13,weight="bold",y=.99)
    fig.text(.5,.945,"Each sample is tracked separately using original epsilon values.\nMissing peaks remain blank.",ha="center",fontsize=8.5)
    fig.text(.5,.025,"Dashed links: numerical continuity, not physical assignment.\nBars: derivative-window sensitivity, not confidence intervals.\nRed x: measurement-edge warning.",ha="center",fontsize=8.5)
    fig.subplots_adjust(left=.16 if len(branches)==1 else .08,right=.98,bottom=.13,top=.89,hspace=.4,wspace=.25)
    save(fig,out,"candidate_peak_temperature_tracks",formats)
