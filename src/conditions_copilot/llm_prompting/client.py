from __future__ import annotations
import os, json, subprocess, tempfile, sys
from typing import Dict, Any

def call_llm(payload: Dict[str, Any], system_path: str) -> Dict[str, Any]:
    """Call external LLM via command specified in $LLM_CMD.
    If not set, fall back to interactive mode: print payload and ask user to paste JSON.
    """
    cmd = os.environ.get("LLM_CMD")
    system_txt = open(system_path, "r", encoding="utf-8").read()
    if not cmd:
        print("\n=== SYSTEM ===\n" + system_txt)
        print("\n=== USER (payload) ===\n" + json.dumps(payload, indent=2))
        print("\nPaste the LLM JSON reply then press Enter (Ctrl-D to finish):\n")
        data = sys.stdin.read()
        return json.loads(data)

    # Non-interactive: pipe 'SYSTEM:\n...\n\nUSER:\n<payload>' to the command
    prompt = f"SYSTEM:\n{system_txt}\n\nUSER:\n{json.dumps(payload)}"
    proc = subprocess.run(cmd, input=prompt.encode("utf-8"), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        raise RuntimeError(f"LLM command failed: {proc.stderr.decode('utf-8', 'ignore')}")
    out = proc.stdout.decode("utf-8", "ignore").strip()
    # Try to extract JSON if extra text present
    first = out.find("{"); last = out.rfind("}")
    if first >= 0 and last >= 0:
        out = out[first:last+1]
    return json.loads(out)
