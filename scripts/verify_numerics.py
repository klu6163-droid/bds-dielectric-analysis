"""Numerical invariants and CLI safety checks, using synthetic spectra only."""
from pathlib import Path
import csv
import json
import subprocess
import sys
import tempfile
import numpy as np
sys.dont_write_bytecode = True
from analyze_bds import EPS0, DEFAULTS, WINDOWS, local_derivative, stable_peaks, conduction_scan, convert_group, sha
import bds_enhancements as ext
import result_lock


def expect_value_error(action):
    try:
        action()
    except ValueError:
        return
    raise AssertionError("Expected a rejected invalid request")


def verify_v24():
    # A user-selected 13-point primary actually reaches peak extraction.
    f=np.logspace(-2,5,211)
    ep=4+30/(1+(2*np.pi*f*.02)**2)
    cfg=ext.derivative_settings(dict(mode="legacy_points",primary_window_points=13))
    ds,primary,diagnostics,label=ext.build_derivatives(f,ep,cfg,local_derivative)
    assert primary==13 and "primary 13" in label
    assert np.array_equal(ds[13],local_derivative(np.log(f),ep,13),equal_nan=True)
    peaks,_=stable_peaks(f,ds,DEFAULTS,primary)
    assert peaks and all(p["primary_window"]=="13" for p in peaks)
    assert ext.build_derivatives(f,ep,ext.derivative_settings({"mode":"legacy_points"}),local_derivative)[1]==9
    for invalid in [{"windows_points":[7,9,12]}, {"windows_points":[7,9,9]},
                    {"windows_points":[13]}, {"windows_points":[3,7,9]},
                    {"primary_window_points":11}, {"primary_window_points":True}]:
        expect_value_error(lambda: ext.derivative_settings(dict(mode="legacy_points",**invalid)))
    gap_f=10**np.array([0,.05,.10,.15,.20,2,2.05,2.10,2.15,2.20])
    sigma=np.full(len(gap_f),2e-9)
    region,scan,sens=conduction_scan(gap_f,sigma,sigma*100/(2*np.pi*gap_f*EPS0),DEFAULTS)
    assert not region["numerical_candidate"] and region["status"]=="Insufficient_contiguous_frequency_span_or_points"
    assert region["frequency_gap_count"]==1 and abs(region["gap_audit"][0]["gap_decades"]-1.8)<1e-12
    assert not scan and not any(s["candidate_found"] for s in sens)
    # Two independently adequate segments may contribute a window, but never one spanning the gap.
    split_f=10**np.r_[np.linspace(0,1.2,25),np.linspace(3,4.2,25)]
    sigma=np.full(len(split_f),2e-9)
    region,scan,sens=conduction_scan(split_f,sigma,sigma*100/(2*np.pi*split_f*EPS0),DEFAULTS)
    assert region["numerical_candidate"] and region["span_decades"]<=1.2+1e-12
    assert region["touches_internal_frequency_gap"]
    assert all(s["largest_frequency_gap_decades"]<=DEFAULTS["max_frequency_gap_decades"] for s in scan)
    for s in sens:
        if s["candidate_found"]:assert not (s["f_start_Hz"]<100 and s["f_end_Hz"]>1000)
    _,_,diag,_=ext.build_derivatives(split_f,30-np.log10(split_f),cfg,local_derivative)
    assert any(r["status"]=="frequency_gap" for rows in diag.values() for r in rows)
    with tempfile.TemporaryDirectory(prefix="bds_lock_verify_") as td:
        root=Path(td);source=root/"raw.csv";source.write_text("original\n",encoding="utf-8")
        out=root/"result";out.mkdir();table=out/"selected.csv";table.write_text("sample,temperature_C\nTest,25\n",encoding="utf-8")
        science=dict(purpose="fresh_analysis",skill_version="2.4",conversion="exported_epsilon",
            column_mapping={"epsilon_prime":"epsilon_prime_export"},derivative=cfg,screening=DEFAULTS,
            samples=["Test"],units={"sigma":"S/cm"},tau_convention="1/(2*pi*f_peak)")
        manifest=result_lock.seal_run(out,science,{str(source):sha(source)},["selected.csv"])
        result=result_lock.locked_inputs(out,["selected.csv"])
        assert result["analysis_id"]==manifest["analysis_id"] and result["source_hashes_checked"]
        ident=result["analysis_id"]
        result_lock.require_same_analysis(ident,ident,ident)
        expect_value_error(lambda:result_lock.require_same_analysis(ident,"other_run",ident))
        expect_value_error(lambda:result_lock.locked_inputs(out,["unbound.csv"]))
        expect_value_error(lambda:result_lock.artifact_path(out,"../raw.csv"))
        expect_value_error(lambda:result_lock.seal_run(out,science,{str(source):sha(source)},["selected.csv"]))
        original=table.read_bytes();table.write_text("changed\n",encoding="utf-8")
        expect_value_error(lambda:result_lock.verify_run(out));table.write_bytes(original)
        source.write_text("changed source\n",encoding="utf-8")
        expect_value_error(lambda:result_lock.verify_run(out))
        assert not result_lock.verify_run(out,check_sources=False)["source_hashes_checked"]
        source.write_text("original\n",encoding="utf-8")
        mf=out/"result_manifest.json";saved=mf.read_bytes()
        altered=dict(manifest,analysis_id="wrong");mf.write_text(json.dumps(altered),encoding="utf-8")
        expect_value_error(lambda:result_lock.verify_run(out));mf.write_bytes(saved)
        other=root/"history";other.mkdir();(other/"selected.csv").write_bytes(original)
        expect_value_error(lambda:result_lock.seal_run(other,dict(science,purpose="historical_replot"),{str(source):sha(source)},["selected.csv"]))
        assert result_lock.verify_run(out)["pass_"]
    print("PASS v2.4: selectable 13-point primary, legacy default 9, point-window validation, IC/derivative gaps, result identity/tamper rejection, mixed-version rejection.")


def verify_v2():
    cfg = ext.derivative_settings({})
    # Same quadratic on different/nonuniform grids: same physical derivative.
    support_counts = []
    for n in (141, 281, 421):
        x = np.linspace(-2,5,n)+.0005*np.sin(np.arange(n))
        f = 10**x
        ln = np.log(f)
        ep = 8+2*ln+.3*ln*ln
        ds, primary, diagnostics, _ = ext.build_derivatives(f,ep,cfg,local_derivative)
        for w, values in ds.items():
            good = np.isfinite(values)
            assert np.allclose(values[good],(-np.pi/2*(2+.6*ln))[good],rtol=1e-9,atol=1e-9)
            assert all(r["actual_span_decades"]<=r["target_span_decades"]+1e-10 for r in diagnostics[w])
        support_counts.append(diagnostics[primary][n//2]["points_used"])
        tau=.02
        ep=4+30/(1+(2*np.pi*f*tau)**2)
        ds,primary,_,_=ext.build_derivatives(f,ep,cfg,local_derivative)
        peaks,_=stable_peaks(f,ds,DEFAULTS,primary)
        assert len(peaks)==1 and abs(np.log10(peaks[0]["frequency_Hz"]/(1/(2*np.pi*tau))))<.06
    assert support_counts[0]<support_counts[1]<support_counts[2]
    # Do not smooth across a large unmeasured frequency gap.
    f=np.r_[np.logspace(0,1,30),np.logspace(1.4,2.4,30)]
    _,_,diagnostics,_=ext.build_derivatives(f,10-np.log10(f),dict(cfg,spans_decades=[.5,1.0],primary_span_decades=.5),local_derivative)
    assert any(r["status"]=="frequency_gap" for records in diagnostics.values() for r in records)
    # Partial-window edge values remain absent.
    values,_=ext.span_derivative(np.logspace(0,3,101),np.linspace(20,5,101),.5,cfg)
    assert np.isnan(values[0]) and np.isnan(values[-1])
    def group(temp,freqs,sample="Test"):
        return dict(sample=sample,conversion_branch="exported",temperature_C=temp,
            peaks=[dict(candidate_id=f"P{i+1}",frequency_Hz=v,height=10,frequency_min_Hz=v*.95,
                frequency_max_Hz=v*1.05,tau_peak_s_NOT_assigned_tau_s=1/(2*np.pi*v),
                low_edge_warning=False,high_edge_warning=False) for i,v in enumerate(freqs)])
    tracking=ext.track_settings({})
    rows=[group(t,[100*1.5**i,10000*.9**i]) for i,t in enumerate((25,30,35,40,45))]
    obs,tracks,links=ext.temperature_tracks(rows,tracking)
    assert len(tracks)==2 and all(t["n_temperatures"]==5 for t in tracks)
    assert len([l for l in links if l["status"]=="provisional_link"])==8
    # Ambiguous split-like candidates are not forcibly joined.
    _,tracks,links=ext.temperature_tracks([group(25,[100]),group(30,[90,110])],tracking)
    assert len(tracks)==3 and any(l["status"]=="ambiguous_not_linked" for l in links)
    assert not any(l["status"]=="provisional_link" for l in links)
    # No bridge across an explicitly measured temperature without a peak.
    _,tracks,_=ext.temperature_tracks([group(25,[100]),group(30,[]),group(35,[120])],tracking)
    assert len(tracks)==2
    _,tracks,_=ext.temperature_tracks([group(25,[100]),group(45,[120])],tracking)
    assert len(tracks)==2
    _,tracks,_=ext.temperature_tracks([group(25,[100]),group(30,[120]),
        group(25,[100],"Other"),group(30,[120],"Other")],tracking)
    assert len(tracks)==2 and len(set(t["sample"] for t in tracks))==2
    print("PASS v2: frequency-span sampling-density invariance, gaps/edges, continuous tracks, ambiguous links, missing temperatures, sample isolation.")


def main():
    # Irregular ln(f) grid: quadratic derivative has a closed-form result.
    x = np.linspace(0, 11, 121) + .005*np.sin(np.arange(121))
    y = 3 + 2*x + .3*x*x
    for w in WINDOWS:
        actual = local_derivative(x, y, w)
        h = w//2
        expected = -np.pi/2*(2+.6*x)
        assert np.allclose(actual[h:-h], expected[h:-h], rtol=1e-10, atol=1e-10)
        assert np.isnan(actual[:h]).all() and np.isnan(actual[-h:]).all()
    # Debye epsilon' derivative peaks at omega*tau=1.
    f = np.logspace(-2, 5, 211)
    tau = .02
    ep = 4 + 30/(1+(2*np.pi*f*tau)**2)
    ds = {w: local_derivative(np.log(f), ep, w) for w in WINDOWS}
    peaks, _ = stable_peaks(f, ds, DEFAULTS)
    assert len(peaks) == 1
    assert abs(np.log10(peaks[0]["frequency_Hz"]/(1/(2*np.pi*tau)))) < .06
    # Constant sigma must be a candidate; a power-law loss background is not.
    sigma = np.full(len(f), 2e-9)
    el = sigma*100/(2*np.pi*f*EPS0)
    region, _, _ = conduction_scan(f, sigma, el, DEFAULTS)
    assert region["numerical_candidate"] and abs(region["beta_sigma"]) < 1e-10
    sigma2 = 2e-9*f**.36
    region, _, _ = conduction_scan(f, sigma2, sigma2*100/(2*np.pi*f*EPS0), DEFAULTS)
    assert not region["numerical_candidate"] and abs(region["beta_sigma"]-.36) < 1e-10
    region, _, _ = conduction_scan(f[:5], sigma[:5], el[:5], DEFAULTS)
    assert region["status"] == "Insufficient_frequency_span_or_points"
    empty, _ = stable_peaks(f[:5], {w: local_derivative(np.log(f[:5]), ep[:5], w) for w in WINDOWS}, DEFAULTS)
    assert empty == []
    # Direct epsilon is independent of Cp and geometry and preserves source numbers.
    direct_rows = [dict(epsilon_prime_export=7., epsilon_doubleprime_export=3., Cp_F=-1.) for _ in f]
    direct = convert_group(dict(f=f, source_rows=direct_rows, sample="Direct", temperature_C=25),
                           dict(area_m2=None, thickness_m=None), "exported_epsilon")
    assert np.array_equal(direct["ep"], np.full(len(f), 7.))
    assert np.array_equal(direct["el"], np.full(len(f), 3.))
    assert np.allclose(direct["sigma"], 2*np.pi*f*EPS0*3/100, rtol=1e-12, atol=0)
    # Read-only inspect, empty candidates, short span, existing-output refusal.
    with tempfile.TemporaryDirectory(prefix="bds_skill_verify_") as td:
        root = Path(td)
        source = root/"short.csv"
        with source.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(["temperature_C", "frequency_Hz", "Cp_F", "tan_delta", "epsilon_prime_export", "epsilon_doubleprime_export"])
            for freq in (10, 11, 13, 16, 20):
                writer.writerow([31.5, freq, -1, -1, 10, 1])
        cfg = root/"config.json"
        config = dict(formats=["png"],
            datasets=[dict(sample="Synthetic", path="short.csv", format="long_csv")])
        cfg.write_text(json.dumps(config), encoding="utf-8")
        script = Path(__file__).with_name("run_bds.py")
        command = [sys.executable, "-X", "utf8", "-B", str(script), "--config", str(cfg)]
        before = sha(source)
        result = subprocess.run(command+["--inspect"], capture_output=True, text=True, encoding="utf-8")
        assert result.returncode == 0, result.stderr
        assert set(p.name for p in root.iterdir()) == {"short.csv", "config.json"}
        out = root/"results"
        result = subprocess.run(command+["--out", str(out)], capture_output=True, text=True, encoding="utf-8")
        assert result.returncode == 0, result.stderr
        assert sha(source) == before
        with (out/"all_data.csv").open(encoding="utf-8-sig") as stream:
            saved = list(csv.DictReader(stream))
        assert all(np.isclose(float(r["epsilon_prime_used"]), 10) for r in saved)
        assert all(np.isclose(float(r["epsilon_doubleprime_used"]), 1) for r in saved)
        assert all(float(r["epsilon_prime_export"]) == 10 for r in saved)
        assert "epsilon_prime_CpD" not in saved[0] and "epsilon_doubleprime_CpD" not in saved[0]
        assert not (out/"user_selection.json").exists()
        assert not (out/"alternate_conversion").exists()
        assert not list(out.glob("conversion_comparison*"))
        assert not (out/"comparison_status.json").exists()
        # Invalid modes and old two-route entry points must fail before output.
        for override in ({"conversion":"parallel_cp_d"},{"conversion":"conditional_cp_d"},{"compare_conversions":True},{"cpd_definition_confirmed":True},{"user_selection":None}):
            badcfg = root/"retired.json"
            badcfg.write_text(json.dumps(dict(config,**override)),encoding="utf-8")
            rejected = root/"must_not_exist"
            result = subprocess.run([sys.executable,"-B",str(script),"--config",str(badcfg),
                                     "--out",str(rejected)],capture_output=True,text=True,encoding="utf-8")
            assert result.returncode != 0 and not rejected.exists()
        # Missing epsilon must stop; Cp and D cannot silently substitute.
        missing = root/"missing.csv"
        missing.write_text("temperature_C,frequency_Hz,Cp_F,tan_delta\n25,10,1e-11,1\n25,20,1e-11,1\n25,30,1e-11,1\n",encoding="utf-8")
        badcfg = root/"missing_config.json"
        badcfg.write_text(json.dumps(dict(formats=["png"],datasets=[dict(sample="Missing",path="missing.csv",format="long_csv")])),encoding="utf-8")
        missing_out=root/"must_not_reconstruct"
        result=subprocess.run([sys.executable,"-X","utf8","-B",str(script),"--config",str(badcfg),"--out",str(missing_out)],capture_output=True,text=True,encoding="utf-8")
        assert result.returncode != 0 and not missing_out.exists()
        with (out/"candidate_peaks.csv").open(encoding="utf-8-sig") as stream:
            assert list(csv.DictReader(stream)) == []
        ver = json.loads((out/"verification.json").read_text(encoding="utf-8"))
        assert ver["output_rows"] == 5 and ver["candidate_peaks"] == 0
        locked=result_lock.locked_inputs(out,["all_data.csv","candidate_peaks.csv","candidate_IC_regions.csv"])
        assert locked["pass_"] and locked["science"]["skill_version"]=="2.4"
        saved_config=json.loads((out/"run_config.json").read_text(encoding="utf-8"))
        assert saved_config["screening"]["max_frequency_gap_decades"]==.15
        assert (out/"IC_frequency_gap_audit.csv").is_file()
        output_before = sha(out/"verification.json")
        rerun = subprocess.run(command+["--out", str(out)], capture_output=True, text=True, encoding="utf-8")
        assert rerun.returncode != 0 and sha(out/"verification.json") == output_before
    verify_v2()
    verify_v24()
    print("PASS: irregular-grid derivative, Debye peak, plateau/nonplateau, units, empty/short spectra, read-only inspect, no overwrite.")


if __name__ == "__main__":
    main()
