"""Bind a BDS result to source columns, processing settings and exact saved tables.

Integrity/version checking, not authentication or proof of a physical assignment.
No original files are changed; the manifest is written once into the result folder.
"""
from pathlib import Path
import argparse
import hashlib
import json

SCHEMA = "bds-result-lock/1"


def file_sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest(value):
    raw = json.dumps(value, sort_keys=True, ensure_ascii=False, allow_nan=False,
                     separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def artifact_path(directory, name):
    directory = Path(directory).resolve()
    rel = Path(name)
    path = (directory / rel).resolve()
    if rel.is_absolute() or not path.is_relative_to(directory) or path == directory:
        raise ValueError(f"Result artifact must stay inside result directory: {name}")
    return path


def seal_run(directory, science, source_sha256, artifacts):
    """Also usable for an explicitly identified historical snapshot; never infer its route."""
    directory = Path(directory).resolve()
    target = directory / "result_manifest.json"
    if target.exists():
        raise ValueError("Existing result manifest is immutable; create a separately named result version")
    required = {"purpose","skill_version","conversion","column_mapping","derivative","screening","samples","units","tau_convention"}
    if not required <= set(science):
        raise ValueError(f"Missing scientific provenance: {sorted(required-set(science))}")
    if science["purpose"] not in {"fresh_analysis","historical_replot"}:
        raise ValueError("Purpose must distinguish new analysis from historical replot")
    if science["purpose"] == "historical_replot" and not science.get("historical_source_note"):
        raise ValueError("Historical snapshot requires its actual route/version evidence")
    if not source_sha256 or not artifacts:
        raise ValueError("Both source hashes and result artifacts are required")
    for path, expected in source_sha256.items():
        if file_sha(path) != expected:
            raise ValueError(f"Source changed before locking: {path}")
    records = {}
    for name in artifacts:
        path = artifact_path(directory, name)
        if path == target:
            raise ValueError("A manifest cannot include its own hash")
        records[path.relative_to(directory).as_posix()] = file_sha(path)
    body = dict(schema=SCHEMA, science=science, source_sha256=source_sha256,
                artifacts_sha256=records)
    manifest = dict(body, analysis_id=digest(body))
    with target.open("x", encoding="utf-8") as stream:
        json.dump(manifest, stream, ensure_ascii=False, indent=2, allow_nan=False)
    return manifest


def verify_run(directory, *, check_sources=True):
    directory = Path(directory).resolve()
    manifest = json.loads((directory / "result_manifest.json").read_text(encoding="utf-8"))
    body = {k:v for k,v in manifest.items() if k != "analysis_id"}
    if body.get("schema") != SCHEMA or digest(body) != manifest.get("analysis_id"):
        raise ValueError("Result manifest identity does not match its contents")
    if not body.get("artifacts_sha256") or not body.get("source_sha256"):
        raise ValueError("Incomplete result manifest")
    for name, expected in body["artifacts_sha256"].items():
        if file_sha(artifact_path(directory, name)) != expected:
            raise ValueError(f"Changed result artifact: {name}")
    if check_sources:
        for path, expected in body["source_sha256"].items():
            if file_sha(path) != expected:
                raise ValueError(f"Changed source: {path}")
    return dict(pass_=True, analysis_id=manifest["analysis_id"],
                source_hashes_checked=check_sources,
                artifacts_checked=len(body["artifacts_sha256"]), science=body["science"])


def locked_inputs(directory, filenames, *, check_sources=True):
    """Resolve sigma/peak tables from one verified result version, never by latest-file glob."""
    result = verify_run(directory, check_sources=check_sources)
    manifest = json.loads((Path(directory)/"result_manifest.json").read_text(encoding="utf-8"))
    paths = []
    for name in filenames:
        path = artifact_path(directory, name)
        key = path.relative_to(Path(directory).resolve()).as_posix()
        if key not in manifest["artifacts_sha256"]:
            raise ValueError(f"Table is not bound to this result version: {name}")
        paths.append(str(path))
    return dict(result, paths=paths)


def require_same_analysis(sigma_analysis_id, tau_analysis_id, expected_analysis_id):
    if not expected_analysis_id or sigma_analysis_id != expected_analysis_id or tau_analysis_id != expected_analysis_id:
        raise ValueError("Sigma and tau must belong to the same explicitly selected analysis version")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--files", nargs="*", default=[])
    parser.add_argument("--artifacts-only", action="store_true",
                        help="Explicit offline check; reports that original source hashes were not rechecked")
    args = parser.parse_args()
    result = locked_inputs(args.run, args.files, check_sources=not args.artifacts_only)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, OSError) as exc:
        raise SystemExit(f"BDS result verification failed: {exc}")
