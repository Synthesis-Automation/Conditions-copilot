from __future__ import annotations
import os, sys, json, platform

def mask(s: str, show: int = 4) -> str:
    if not s:
        return "<missing>"
    return s[:show] + "…" + str(len(s))

def main():
    print("Python:", sys.version.split(" ")[0], platform.platform())
    key = os.environ.get("OPENAI_API_KEY")
    model = os.environ.get("OPENAI_MODEL", "gpt-4o-mini")
    base = os.environ.get("OPENAI_BASE_URL")
    cmd = os.environ.get("LLM_CMD")
    print("LLM_CMD:", cmd or "<empty>")
    print("OPENAI_API_KEY:", mask(key))
    print("OPENAI_MODEL:", model)
    print("OPENAI_BASE_URL:", base or "<default>")

    if not key:
        print("\nNo OPENAI_API_KEY in current process environment.")
        print("Restart your shell/app or set it for this session.")
        sys.exit(2)

    try:
        from openai import OpenAI
    except Exception as e:
        print("\nFailed to import openai package:", e)
        print("Install with: pip install -e .  (or pip install openai)")
        sys.exit(3)

    try:
        client = OpenAI(api_key=key, base_url=base) if base else OpenAI(api_key=key)
        resp = client.chat.completions.create(
            model=model,
            messages=[
                {"role": "system", "content": "Return strictly json only. Reply with a json object."},
                {"role": "user", "content": json.dumps({"ping": True})},
            ],
            response_format={"type": "json_object"},
            temperature=0.0,
        )
        content = (resp.choices[0].message.content or "").strip()
        print("\nAPI call succeeded. First 200 chars of JSON:")
        print(content[:200])
        sys.exit(0)
    except Exception as e:
        print("\nOpenAI API call failed:")
        print(repr(e))
        sys.exit(4)

if __name__ == "__main__":
    main()
