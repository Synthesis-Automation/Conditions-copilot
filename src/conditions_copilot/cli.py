from __future__ import annotations
import argparse, json, os
from .controller import build_discovery_payload, run_once

def main():
    ap = argparse.ArgumentParser(description="Conditions Copilot CLI")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p1 = sub.add_parser("discover", help="Print universal discovery payload")
    p1.add_argument("--rxn", required=True, help="Reaction SMILES")

    p2 = sub.add_parser("run", help="Run one discovery→propose→validate loop (interactive LLM by default)")
    p2.add_argument("--rxn", required=True, help="Reaction SMILES")
    p2.add_argument("--dataset", default=os.path.join(os.path.dirname(__file__), "..", "..", "data", "demo_neighbors.csv"))
    p2.add_argument("--dictdir", default=os.path.join(os.path.dirname(__file__), "..", "..", "data", "dicts"))
    p2.add_argument("--system", default=os.path.join(os.path.dirname(__file__), "llm_prompting", "system.txt"))

    args = ap.parse_args()
    if args.cmd == "discover":
        print(json.dumps(build_discovery_payload(args.rxn), indent=2))
    elif args.cmd == "run":
        out = run_once(args.rxn, args.dataset, args.dictdir, args.system)
        print(json.dumps(out, indent=2))

if __name__ == "__main__":
    main()
