"""Plot sigma*T*tau versus tau from one locked, explicitly selected BDS table.

This is a replot/calculation entry point, not a conductivity-window or peak selector.
Dependencies: NumPy and Matplotlib. Origin handoff files are data/mappings, not OPJU.
"""
from pathlib import Path
import argparse
import copy
import csv
import json
import math
import shutil
import textwrap

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter

import result_lock as lock

VERSION = "2.5"
PRODUCT = "product_S_cm_inv_K_s"
SIGMA_STATES = {"confirmed_dc", "candidate_plateau", "ac_reference", "missing"}
TAU_STATES = {"assigned_segmental", "unassigned", "missing"}
COLORS = ["#2E80C5", "#E77C23", "#29836B", "#8D5AB5", "#B24B63", "#687A8A"]
MARKERS = ["s", "D", "o", "^", "v", "P"]


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows, fields=None):
    if fields is None:
        fields = list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open("x", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: "" if isinstance(v, float) and not math.isfinite(v) else v
                             for k, v in row.items()})


def write_json(path, value):
    with Path(path).open("x", encoding="utf-8") as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)


def number(row, key):
    value = row.get(key, "")
    if str(value).strip().lower() in {"", "nan", "na", "none"}:
        return float("nan")
    return float(value)


def close(a, b):
    return math.isclose(a, b, rel_tol=1e-8, abs_tol=0)


def warning(row):
    return any("warning" in k and str(v).strip().lower() in {"true", "1", "yes"}
               for k, v in row.items())


def prepare(rows, expected_id, reference_temperature=25.0):
    required = {"sample", "temperature_C", "sigma_S_cm", "f_peak_Hz", "tau_peak_s",
                "f_start_Hz", "f_end_Hz", "sigma_status", "tau_status", "selection_reason"}
    if not rows or not required <= set(rows[0]):
        raise ValueError(f"Selected table must contain {sorted(required)}")
    seen, groups, out = set(), {}, []
    for line, source in enumerate(rows, 2):
        r = dict(source)
        sample = r["sample"].strip()
        t = number(r, "temperature_C")
        if not sample or not math.isfinite(t) or t + 273.15 <= 0:
            raise ValueError(f"Invalid sample/temperature at CSV line {line}")
        r["sample"] = sample
        # One selected process per sample and temperature in a comparison table.
        if (sample, t) in seen:
            raise ValueError(f"Duplicate sample/temperature: {sample}, {t}")
        seen.add((sample, t))
        present_ids = [r.get(k, "") for k in ("analysis_id", "sigma_analysis_id", "tau_analysis_id")]
        if any(present_ids) and (not all(present_ids) or any(v != expected_id for v in present_ids)):
            raise ValueError(f"Mixed or unbound sigma/tau analysis IDs at CSV line {line}")
        if not any(present_ids):
            for k in ("analysis_id", "sigma_analysis_id", "tau_analysis_id"):
                r[k] = expected_id
        if r["sigma_status"] not in SIGMA_STATES or r["tau_status"] not in TAU_STATES:
            raise ValueError(f"Unknown evidence status at CSV line {line}")
        if not r["selection_reason"].strip():
            raise ValueError(f"Missing selection reason at CSV line {line}")
        if r["sigma_status"] == "confirmed_dc" and not r.get("sigma_evidence", "").strip():
            raise ValueError("confirmed_dc requires documented sigma_evidence")
        if r["tau_status"] == "assigned_segmental" and not r.get("tau_evidence", "").strip():
            raise ValueError("assigned_segmental requires documented tau_evidence")
        tk = t + 273.15
        if math.isfinite(number(r, "temperature_K")) and not close(number(r, "temperature_K"), tk):
            raise ValueError(f"Kelvin conversion mismatch at CSV line {line}")
        r["temperature_K"] = tk
        sigma, tau, peak = (number(r, k) for k in ("sigma_S_cm", "tau_peak_s", "f_peak_Hz"))
        lo, hi = number(r, "f_start_Hz"), number(r, "f_end_Hz")
        if math.isfinite(sigma) and sigma > 0 and r["sigma_status"] != "missing":
            if not (math.isfinite(lo) and math.isfinite(hi) and 0 < lo <= hi):
                raise ValueError(f"Invalid sigma support interval at CSV line {line}")
        if math.isfinite(tau) and tau > 0 and r["tau_status"] != "missing":
            if not (math.isfinite(peak) and peak > 0 and close(tau, 1/(2*math.pi*peak))):
                raise ValueError(f"tau != 1/(2*pi*f_peak) at CSV line {line}")
        reasons = []
        if not (math.isfinite(sigma) and sigma > 0) or r["sigma_status"] == "missing":
            reasons.append("missing/nonpositive sigma")
        if not (math.isfinite(tau) and tau > 0) or r["tau_status"] == "missing":
            reasons.append("missing/nonpositive tau")
        valid = not reasons
        p = sigma * tk * tau if valid else float("nan")
        if valid and not (math.isfinite(p) and p > 0):
            raise ValueError(f"Nonfinite/underflow product at CSV line {line}")
        old = number(r, PRODUCT)
        if valid and math.isfinite(old) and not close(old, p):
            raise ValueError(f"Saved product disagrees with selected sigma/T/tau at CSV line {line}")
        # Preserve an already verified saved product exactly; never move a point for appearance.
        r[PRODUCT] = old if valid and math.isfinite(old) else p
        r["plot_valid"] = valid
        r["plot_exclusion_reason"] = "; ".join(reasons)
        r["source_csv_line"] = line
        track = r.get("comparison_group", "").strip() or r.get("track_id", "").strip() or "selected"
        key = (sample, r["sigma_status"], r["tau_status"], track)
        if key not in groups:
            groups[key] = []
        groups[key].append(r)
        out.append(r)
    stats = []
    for index, (key, group) in enumerate(groups.items(), 1):
        gid = f"G{index:02d}"
        good = [r for r in group if r["plot_valid"]]
        values = np.array([r[PRODUCT] for r in good])
        times = np.array([number(r, "tau_peak_s") for r in good])
        baseline = next((r[PRODUCT] for r in good if number(r, "temperature_C") == reference_temperature), None)
        for r in group:
            r["plot_group"] = gid
            r["reference_temperature_C"] = reference_temperature
            r["reference_product"] = baseline if baseline is not None else float("nan")
            r["relative_to_reference"] = r[PRODUCT]/baseline if r["plot_valid"] and baseline else float("nan")
            r["normalization_status"] = "available" if baseline else "reference_temperature_missing_in_group"
        n = len(good)
        mean = float(values.mean()) if n else float("nan")
        sd = float(values.std(ddof=1)) if n > 1 else float("nan")
        stats.append(dict(plot_group=gid, sample=key[0], sigma_status=key[1], tau_status=key[2],
            comparison_group=key[3], n_measured=len(group), n_valid=n,
            excluded_temperatures_C=";".join(str(r["temperature_C"]) for r in group if not r["plot_valid"]),
            temperature_min_C=min(number(r,"temperature_C") for r in group),
            temperature_max_C=max(number(r,"temperature_C") for r in group),
            tau_span_decades=float(np.log10(times.max()/times.min())) if n else float("nan"),
            geometric_mean=float(np.exp(np.log(values).mean())) if n else float("nan"),
            arithmetic_mean=mean, sample_std_ddof1=sd, CV_percent=100*sd/mean if n > 1 else float("nan"),
            max_min_ratio=float(values.max()/values.min()) if n else float("nan"),
            edge_warning_count=sum(warning(r) for r in group),
            interpretation="descriptive; cannot assess constancy" if n < 2 else "descriptive; not proof of coupling"))
    return out, stats


def style_groups(rows, stats, sample_order=None):
    samples = list(dict.fromkeys(r["sample"] for r in rows))
    if sample_order:
        if len(sample_order) != len(set(sample_order)) or set(sample_order) != set(samples):
            raise ValueError("--sample-order must name each measured sample exactly once")
        samples = sample_order
    styles = {}
    for s in stats:
        i = samples.index(s["sample"])
        multiple = sum(g["sample"] == s["sample"] for g in stats) > 1
        styles[s["plot_group"]] = dict(order=i,color=COLORS[i % len(COLORS)],
            marker=MARKERS[(i if not multiple else list(g["plot_group"] for g in stats).index(s["plot_group"])) % len(MARKERS)],
            label=s["sample"] + (f' / {s["sigma_status"]} / {s["comparison_group"]}' if multiple else ""))
    return styles


def labels(rows):
    ss = {r["sigma_status"] for r in rows if r["plot_valid"]}
    ts = {r["tau_status"] for r in rows if r["plot_valid"]}
    sub = "dc" if ss == {"confirmed_dc"} else "window" if ss == {"ac_reference"} else "candidate" if ss == {"candidate_plateau"} else "selected"
    tau = "s" if ts == {"assigned_segmental"} else "peak"
    return sub, tau


def draw(rows, stats, styles, destination, two_panel, reference_temperature, dpi):
    plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":11, "axes.linewidth":1.3,
        "svg.fonttype":"none", "pdf.fonttype":42, "xtick.direction":"in", "ytick.direction":"in"})
    fig, axes = plt.subplots(1, 2 if two_panel else 1, figsize=(14,6.8) if two_panel else (8.4,6.8))
    axes = np.atleast_1d(axes);ax=axes[0]
    sub,tau_name=labels(rows)
    for s in sorted(stats,key=lambda item:styles[item["plot_group"]]["order"]):
        group=sorted([r for r in rows if r["plot_group"]==s["plot_group"]],key=lambda r:number(r,"temperature_C"))
        valid=[r for r in group if r["plot_valid"]]
        st=styles[s["plot_group"]]
        x=[number(r,"tau_peak_s") for r in valid];y=[r[PRODUCT] for r in valid]
        if valid:
            ax.plot(x,y,linestyle="none",marker=st["marker"],mfc="white",mec=st["color"],mew=1.7,ms=8,label=st["label"],zorder=3)
        if len(valid)>1:
            ax.hlines(s["geometric_mean"],min(x)/1.15,max(x)*1.15,color=st["color"],linestyle=":",lw=1.4,zorder=1)
        for j,r in enumerate(valid):
            offset=(8 if j%2==0 else -19)
            ax.annotate(f'{number(r,"temperature_C"):g} °C'+('*' if warning(r) else ''),
                (number(r,"tau_peak_s"),r[PRODUCT]),xytext=(0,offset),textcoords="offset points",
                ha="center",color=st["color"],fontsize=9)
        if two_panel:
            axes[1].plot([number(r,"temperature_C") for r in group],[r["relative_to_reference"] for r in group],
                color=st["color"],marker=st["marker"],mfc="white",mew=1.5,ms=7,lw=1.3,label=st["label"])
            for r in group:
                relative=r["relative_to_reference"]
                if math.isfinite(relative) and number(r,"temperature_C") != reference_temperature:
                    axes[1].annotate(f'{relative:.3f}',(number(r,"temperature_C"),relative),
                        xytext=(0,9 if st["order"]%2==0 else -17),textcoords="offset points",
                        ha="center",fontsize=8.5,color=st["color"])
    ax.set(xscale="log",yscale="log",xlabel=rf'$\tau_{{\mathrm{{{tau_name}}}}}$ (s)',
           ylabel=rf'$\sigma_{{\mathrm{{{sub}}}}}T\tau_{{\mathrm{{{tau_name}}}}}$ (S cm$^{{-1}}$ K s)')
    valid=[r for r in rows if r["plot_valid"]]
    if valid:
        xs=[number(r,"tau_peak_s") for r in valid];ys=[r[PRODUCT] for r in valid]
        ax.set_xlim(min(xs)/1.7,max(xs)*1.7);ax.set_ylim(min(ys)/1.6,max(ys)*1.8)
        ax.legend(frameon=False,fontsize=10,loc="best")
    else:
        ax.set_xlim(1e-6,1);ax.set_ylim(1e-12,1e-6)
        ax.text(.5,.5,"No valid positive products",transform=ax.transAxes,ha="center")
    ax.yaxis.set_major_locator(LogLocator(base=10,subs=(1.0,)))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_title("Product versus selected peak time",loc="left",fontsize=13)
    if two_panel:
        axes[1].axhline(1,color="#8995A1",lw=1,ls=":")
        axes[1].set(xlabel="Temperature (°C)",ylabel=rf'$(\sigma T\tau)/(\sigma T\tau)_{{{reference_temperature:g}\,^\circ C}}$')
        axes[1].set_title(f"Relative to {reference_temperature:g} °C",loc="left",fontsize=13)
        axes[1].legend(frameon=False,fontsize=10,loc="best")
        temperatures=sorted({number(r,"temperature_C") for r in rows})
        if len(temperatures)<=12:axes[1].set_xticks(temperatures)
        ratios=[r["relative_to_reference"] for r in rows if math.isfinite(r["relative_to_reference"])]
        if ratios:
            pad=max((max(ratios)-min(ratios))*.15,.03)
            axes[1].set_ylim(min(min(ratios),1)-pad,max(max(ratios),1)+pad)
    for a in axes:
        a.tick_params(which="both",top=True,right=True)
        a.grid(axis="y",which="major",color="#DDE3E8",lw=.65)
    statuses="; ".join(dict.fromkeys(f'{s["sample"]}: {s["sigma_status"]}, {s["tau_status"]}' for s in stats))
    note=statuses+". Dotted lines: geometric means; *: retained source warning."
    note+=" Descriptive comparison; flatness alone does not establish segmental coupling."
    if two_panel and any(r["normalization_status"]!="available" for r in rows):
        note+=" Groups without the reference temperature have no normalized series."
    fig.suptitle("Conductivity-relaxation product",fontsize=16,y=.97)
    fig.subplots_adjust(left=.10 if two_panel else .15,right=.97,bottom=.23,top=.86,wspace=.32)
    fig.text(.10 if two_panel else .15,.10,textwrap.fill(note,150 if two_panel else 88),fontsize=8.5,color="#536170",va="center")
    stem="sigmaTtau_two_panel" if two_panel else "sigmaTtau_main"
    for ext in ("png","svg","pdf"):
        fig.savefig(destination/(stem+"."+ext),dpi=dpi)
    plt.close(fig)


def origin_handoff(rows, stats, styles, destination):
    wide, references, mappings = [], [], []
    by_group={s["plot_group"]:[r for r in rows if r["plot_group"]==s["plot_group"]] for s in stats}
    size=max(map(len,by_group.values()))
    for i in range(size):
        record={}
        for gid,group in by_group.items():
            r=group[i] if i<len(group) else None
            record[gid+"_tau_s"]=number(r,"tau_peak_s") if r and r["plot_valid"] else float("nan")
            record[gid+"_product"]=r[PRODUCT] if r and r["plot_valid"] else float("nan")
        wide.append(record)
    for i,s in enumerate(stats):
        gid=s["plot_group"];valid=[r for r in by_group[gid] if r["plot_valid"]]
        mappings.append(dict(group=gid,**styles[gid],x_column_0based=2*i,y_column_0based=2*i+1,
            x_field=gid+"_tau_s",y_field=gid+"_product",plot="scatter",open_symbol=True,connect=False))
        if len(valid)>1:
            for x in [min(number(r,"tau_peak_s") for r in valid)/1.15,max(number(r,"tau_peak_s") for r in valid)*1.15]:
                references.append(dict(plot_group=gid,tau_s=x,product=s["geometric_mean"],role="descriptive_geometric_mean"))
    write_csv(destination/"origin_points_wide.csv",wide)
    write_csv(destination/"origin_reference_lines.csv",references,["plot_group","tau_s","product","role"])
    sub,tau=labels(rows)
    write_json(destination/"origin_handoff.json",dict(schema="bds-sigma-tau-origin-handoff/1",
        is_render_plan=False,creates_opju=False,points="origin_points_wide.csv",references="origin_reference_lines.csv",
        labels_from="sigmaTtau_values.csv",axes={"x_scale":"log10","y_scale":"log10",
        "x_label":f"τ\\-({tau}) (s)","y_label":f"σ\\-({sub})Tτ\\-({tau}) (S cm\\+(-1) K s)"},series=mappings,
        instructions="Use originpro-cli-plot scatter workflow or a separately implemented native template. Preserve each XY pairing and separate dotted reference plots. Run its normal source-bound delivery and inspect every exported graph."))


def run(args):
    if not args.selection_note.strip():raise ValueError("A selection/grouping note is required")
    if not math.isfinite(args.reference_temperature) or args.dpi < 72:
        raise ValueError("Reference temperature must be finite; DPI must be >=72")
    checked=lock.locked_inputs(args.run,[args.table],check_sources=not args.artifacts_only)
    units=checked["science"].get("units",{})
    if units.get("sigma") != "S/cm" or units.get("tau") != "s":
        raise ValueError("Parent manifest must identify sigma in S/cm and tau in s; convert explicitly in a new version first")
    source=Path(checked["paths"][0]);manifest_path=Path(args.run).resolve()/"result_manifest.json"
    expected=checked["science"].get("selected_input_analysis_id",checked["analysis_id"])
    rows,stats=prepare(read_csv(source),expected,args.reference_temperature)
    styles=style_groups(rows,stats,args.sample_order)
    output=Path(args.out).resolve()
    # Do not put new results inside a locked input version.
    if output.is_relative_to(Path(args.run).resolve()):raise ValueError("Output must be outside the locked input result directory")
    if output.exists():raise ValueError("Output already exists; choose a new directory")
    sources={str(p.resolve()):lock.file_sha(p) for p in [source,manifest_path,Path(__file__),Path(lock.__file__)]}
    output.mkdir(parents=True,exist_ok=False)
    shutil.copy2(source,output/"selected_input.csv")
    shutil.copy2(manifest_path,output/"parent_manifest.json")
    shutil.copy2(__file__,output/"plot_sigma_tau.py")
    shutil.copy2(lock.__file__,output/"result_lock.py")
    write_csv(output/"sigmaTtau_values.csv",rows)
    write_csv(output/"descriptive_statistics.csv",stats)
    config=dict(source_run=str(Path(args.run).resolve()),table=args.table,parent_analysis_id=checked["analysis_id"],
        selected_input_analysis_id=expected,selection_note=args.selection_note,layout=args.layout,
        reference_temperature_C=args.reference_temperature,source_hashes_checked=checked["source_hashes_checked"],
        dpi=args.dpi,styles=styles,calculation="sigma_S_cm * (temperature_C + 273.15) * tau_peak_s",
        uncertainty="No error bars invented; window sensitivity is not a confidence interval")
    write_json(output/"plot_config.json",config)
    if args.layout in {"main","both"}:draw(rows,stats,styles,output,False,args.reference_temperature,args.dpi)
    if args.layout in {"two-panel","both"}:draw(rows,stats,styles,output,True,args.reference_temperature,args.dpi)
    origin_handoff(rows,stats,styles,output)
    # Verify exact sources again, including the manifest-bound acquisition files by default.
    lock.verify_run(args.run,check_sources=not args.artifacts_only)
    science=copy.deepcopy(checked["science"])
    science.update(purpose="historical_replot",skill_version=VERSION,
        historical_source_note="Replot the explicitly selected manifest-bound table: "+args.selection_note,
        operation="sigmaTtau plot only; no sigma window changes or new peak selection",
        selected_input_analysis_id=expected,parent_result_id=checked["analysis_id"],plot=config)
    result=lock.seal_run(output,science,sources,[p.name for p in output.iterdir() if p.is_file()])
    lock.verify_run(output)
    return dict(output=str(output),analysis_id=result["analysis_id"],measured_rows=len(rows),
                valid_rows=sum(r["plot_valid"] for r in rows),source_hashes_checked=checked["source_hashes_checked"],
                visual_review_required=True,origin_project_created=False)


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run",required=True,type=Path,help="Existing result_lock result directory")
    p.add_argument("--table",required=True,help="Manifest-bound selected CSV path relative to --run")
    p.add_argument("--out",required=True,type=Path,help="New directory outside --run")
    p.add_argument("--selection-note",required=True,help="Actual selection history and within-sample grouping rationale")
    p.add_argument("--layout",choices=["main","two-panel","both"],default="main")
    p.add_argument("--reference-temperature",type=float,default=25.0)
    p.add_argument("--sample-order",nargs="+")
    p.add_argument("--dpi",type=int,default=300)
    p.add_argument("--artifacts-only",action="store_true",help="Explicit offline use only; report source hashes not rechecked")
    return p


if __name__=="__main__":
    try:print(json.dumps(run(parser().parse_args()),ensure_ascii=False,indent=2))
    except (ValueError,KeyError,OSError) as exc:raise SystemExit(f"sigmaTtau plot failed: {exc}")
