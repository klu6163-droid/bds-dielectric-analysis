"""Entry point for direct original-epsilon BDS analysis."""
import sys
sys.dont_write_bytecode = True
from analyze_bds import main

if __name__ == "__main__":
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    try:
        main()
    except (ValueError, KeyError, OSError) as exc:
        print(f"BDS analysis stopped: {exc}", file=sys.stderr)
        raise SystemExit(2)
