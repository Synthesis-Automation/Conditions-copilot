from __future__ import annotations
import os, json, subprocess, sys
from typing import Dict, Any, Optional

# Built-in test defaults (requested hardcoded fallback)
DEFAULT_TEST_API_KEY = "sk-f4207c24fda4439f8d8fd711da49cf8c"
DEFAULT_TEST_MODEL = "deepseek-v3.1"
DEFAULT_TEST_BASE_URL = "https://dashscope.aliyuncs.com/compatible-mode/v1"


def _json_extract(text: str) -> str:
    first = text.find("{")
    last = text.rfind("}")
    return text[first:last + 1] if first >= 0 and last >= 0 else text


def call_openai(payload: Dict[str, Any], system_txt: str, *,
                model: str, base_url: Optional[str] = None, json_mode: bool = True,
                temperature: float = 0.2) -> Dict[str, Any]:
    from openai import OpenAI
    # Support both OpenAI and Aliyun DashScope (OpenAI-compatible) envs, with hardcoded test fallback
    api_key = (
        os.environ.get("OPENAI_API_KEY")
        or os.environ.get("DASHSCOPE_API_KEY")
        or DEFAULT_TEST_API_KEY
    )
    # Resolve base URL precedence: explicit arg -> OPENAI_BASE_URL -> DASHSCOPE_BASE_URL -> test default
    env_base = os.environ.get("OPENAI_BASE_URL") or os.environ.get("DASHSCOPE_BASE_URL")
    resolved_base = base_url or env_base or DEFAULT_TEST_BASE_URL
    client = OpenAI(api_key=api_key, base_url=resolved_base)
    messages = [
        {"role": "system", "content": system_txt + "\nReturn strictly JSON only."},
        {"role": "user", "content": json.dumps(payload)},
    ]
    kwargs: Dict[str, Any] = {
        "model": model,
        "messages": messages,
        "temperature": temperature,
    }
    if json_mode:
        kwargs["response_format"] = {"type": "json_object"}
    resp = client.chat.completions.create(**kwargs)
    content = (resp.choices[0].message.content or "").strip()
    return json.loads(_json_extract(content))


def call_llm(payload: Dict[str, Any], system_path: str) -> Dict[str, Any]:
    """LLM adapter priority:
    1) $LLM_CMD → shell command (non-interactive)
    2) $OPENAI_API_KEY → OpenAI Chat Completions (asks for confirmation)
    3) Interactive paste fallback (prints payload and waits for JSON on stdin)
    """
    cmd = os.environ.get("LLM_CMD")
    system_txt = open(system_path, "r", encoding="utf-8").read()

    # 1) External command
    if cmd:
        prompt = f"SYSTEM:\n{system_txt}\n\nUSER:\n{json.dumps(payload)}"
        proc = subprocess.run(cmd, input=prompt.encode("utf-8"), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if proc.returncode != 0:
            raise RuntimeError(f"LLM command failed: {proc.stderr.decode('utf-8', 'ignore')}")
        out = proc.stdout.decode("utf-8", "ignore").strip()
        return json.loads(_json_extract(out))

    # 2) OpenAI-compatible API with confirmation (OpenAI or DashScope or built-in test key)
    model = os.environ.get("OPENAI_MODEL") or os.environ.get("DASHSCOPE_MODEL") or DEFAULT_TEST_MODEL
    base_url = (
        os.environ.get("OPENAI_BASE_URL")
        or os.environ.get("DASHSCOPE_BASE_URL")
        or DEFAULT_TEST_BASE_URL
    )
    # show preview and ask for confirmation
    preview = json.dumps(payload)[:1000]
    print("\n=== SYSTEM ===\n" + system_txt)
    print("\n=== USER (payload preview) ===\n" + preview + ("..." if len(json.dumps(payload)) > 1000 else ""))
    dest = f"model '{model}'" + (f" at '{base_url}'" if base_url else "")
    ans = input(f"\nSend to API {dest}? [y/N]: ").strip().lower()
    if ans not in ("y", "yes"):
        print("Aborted by user. Falling back to paste mode.")
    else:
        try:
            return call_openai(payload, system_txt, model=model, base_url=base_url, json_mode=True)
        except Exception as e:
            raise RuntimeError(f"OpenAI call failed: {e}")

    # 3) Interactive paste fallback
    print("\n=== SYSTEM ===\n" + system_txt)
    print("\n=== USER (payload) ===\n" + json.dumps(payload, indent=2))
    print("\nPaste the LLM JSON reply then press Enter (Ctrl-D to finish):\n")
    data = sys.stdin.read()
    return json.loads(data)
