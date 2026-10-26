"""Behavior checks for the selected-value sigmaTtau plot entry point; synthetic only."""
from pathlib import Path
from tempfile import TemporaryDirectory
from argparse import Namespace
import copy
import json
import math
import unittest

import plot_sigma_tau as plot
import result_lock as lock


def row(t=25, sigma=2e-8, peak=100):
    return dict(sample="SyntheticA",temperature_C=str(t),sigma_S_cm=str(sigma),
        f_peak_Hz=str(peak),tau_peak_s=str(1/(2*math.pi*peak)),f_start_Hz="100",f_end_Hz="1000",
        sigma_status="ac_reference",tau_status="unassigned",selection_reason="Synthetic fixed support and peak",
        comparison_group="P1",edge_warning="False")


class SelectedTableTests(unittest.TestCase):
    def test_formula_and_descriptive_statistics(self):
        rows=[row(),row(30,3e-8,200),row(35,4e-8,300)]
        actual,stats=plot.prepare(rows,"version")
        p=2e-8*298.15/(2*math.pi*100)
        self.assertAlmostEqual(actual[0][plot.PRODUCT]/p,1)
        self.assertEqual(actual[0]["relative_to_reference"],1)
        vals=[float(r["sigma_S_cm"])*(float(r["temperature_C"])+273.15)*float(r["tau_peak_s"]) for r in rows]
        self.assertAlmostEqual(stats[0]["geometric_mean"],math.prod(vals)**(1/3))
        self.assertEqual(stats[0]["n_valid"],3)
        self.assertEqual([r["sigma_S_cm"] for r in actual],[r["sigma_S_cm"] for r in rows])
        self.assertEqual([r["tau_peak_s"] for r in actual],[r["tau_peak_s"] for r in rows])

    def test_mixed_identity_rejected(self):
        r=row();r.update(analysis_id="A",sigma_analysis_id="A",tau_analysis_id="B")
        with self.assertRaisesRegex(ValueError,"Mixed"):plot.prepare([r],"A")

    def test_formula_mismatch_rejected(self):
        for key,value in [("temperature_K","25"),("tau_peak_s","0.2"),(plot.PRODUCT,"1.0")]:
            r=row();r[key]=value
            with self.assertRaises(ValueError):plot.prepare([r],"version")

    def test_missing_and_nonpositive_remain(self):
        rows=[row(),row(30,0),row(35)]
        rows[2]["tau_peak_s"]="";rows[2]["tau_status"]="missing"
        actual,stats=plot.prepare(rows,"version")
        self.assertEqual(len(actual),3)
        self.assertFalse(actual[1]["plot_valid"])
        self.assertTrue(math.isnan(actual[1][plot.PRODUCT]))
        self.assertEqual(actual[1]["sigma_S_cm"],"0")
        self.assertEqual(sum(s["n_valid"] for s in stats),1)

    def test_missing_reference_not_substituted(self):
        actual,_=plot.prepare([row(30),row(35)],"version")
        self.assertTrue(all(math.isnan(r["relative_to_reference"]) for r in actual))

    def test_evidence_states_separate_groups(self):
        rows=[row(),row(30)]
        rows[1]["sigma_status"]="candidate_plateau"
        actual,stats=plot.prepare(rows,"version")
        self.assertEqual(len(stats),2)
        self.assertTrue(math.isnan(actual[1]["relative_to_reference"]))
        self.assertTrue(all(s["n_valid"]==1 for s in stats))

    def test_claim_requires_evidence(self):
        for key,value in [("sigma_status","confirmed_dc"),("tau_status","assigned_segmental")]:
            r=row();r[key]=value
            with self.assertRaises(ValueError):plot.prepare([r],"version")

    def test_duplicate_rejected(self):
        with self.assertRaisesRegex(ValueError,"Duplicate"):plot.prepare([row(),row()],"version")

    def test_requested_sample_order(self):
        first=row();first["sample"]="B"
        second=row();second["sample"]="A"
        actual,stats=plot.prepare([first,second],"version")
        styles=plot.style_groups(actual,stats,["A","B"])
        self.assertEqual(styles[stats[0]["plot_group"]]["order"],1)
        self.assertEqual(styles[stats[1]["plot_group"]]["order"],0)

    def test_full_export_protection_and_round_trip(self):
        with TemporaryDirectory(prefix="bds_sigma_tau_") as temp:
            root=Path(temp);source=root/"measurement.txt";source.write_text("Synthetic source, no real experiment")
            original=lock.file_sha(source)
            parent=root/"parent";parent.mkdir()
            rows=[row(),row(30,3e-8,200),row(35,4e-8,300)]
            rows[0]["edge_warning"]="True"
            plot.write_csv(parent/"chosen.csv",rows)
            science=dict(purpose="historical_replot",skill_version="test",conversion="exported_epsilon",
                column_mapping={"sigma":"sigma_S_cm"},derivative={"mode":"synthetic"},screening={},
                samples=["SyntheticA"],units={"sigma":"S/cm","tau":"s"},
                tau_convention="1/(2*pi*f_peak)",historical_source_note="Synthetic fixture")
            lock.seal_run(parent,science,{str(source):original},["chosen.csv"])
            args=Namespace(run=parent,table="chosen.csv",out=root/"result",selection_note="Synthetic selected sequence",
                layout="both",reference_temperature=25.0,sample_order=None,dpi=100,artifacts_only=False)
            result=plot.run(args)
            self.assertEqual(result["valid_rows"],3)
            self.assertFalse(result["origin_project_created"])
            self.assertEqual(lock.file_sha(source),original)
            self.assertTrue(lock.verify_run(args.out)["pass_"])
            for stem in ("sigmaTtau_main","sigmaTtau_two_panel"):
                for ext in ("png","svg","pdf"):
                    self.assertGreater((args.out/(stem+"."+ext)).stat().st_size,1000)
            saved=plot.read_csv(args.out/"sigmaTtau_values.csv")
            self.assertEqual([r["tau_peak_s"] for r in saved],[r["tau_peak_s"] for r in rows])
            handoff=json.loads((args.out/"origin_handoff.json").read_text(encoding="utf-8"))
            self.assertFalse(handoff["is_render_plan"])
            self.assertEqual(handoff["series"][0]["x_field"],"G01_tau_s")
            self.assertEqual(handoff["series"][0]["y_field"],"G01_product")
            with self.assertRaisesRegex(ValueError,"already exists"):plot.run(args)
            bad=copy.copy(args);bad.out=parent/"nested"
            with self.assertRaisesRegex(ValueError,"outside"):plot.run(bad)
            source.write_text("tampered")
            with self.assertRaisesRegex(ValueError,"Changed source"):plot.run(args)


if __name__=="__main__":unittest.main(verbosity=2)
