# ts_embed_from_py.py
import sys
from pathlib import Path

def strip_main_block(source: str) -> str:
    lines = source.splitlines()
    out = []
    cut = False
    for ln in lines:
        s = ln.strip()
        if s.startswith('if __name__ == "__main__"') or s.startswith("if __name__ == '__main__'"):
            cut = True
        if not cut:
            out.append(ln)
    return "\n".join(out).rstrip() + "\n"

def to_ts_embed(py_source: str, var_name: str = "window.dashApp") -> str:
    py_escaped = py_source.replace("`", "\\`")  # escape backticks
    return (
        "// Auto-generated from Python by ts_embed_from_py.py\n"
        "declare global { interface Window { dashApp: string } }\n\n"
        f"{var_name} = `\n"
        f"{py_escaped}"
        "`;\n"
    )

def main():
    if len(sys.argv) < 2:
        print("Usage: python ts_embed_from_py.py <path-to-python-file> [--keep-main]", file=sys.stderr)
        sys.exit(1)
    path = Path(sys.argv[1])
    keep_main = "--keep-main" in sys.argv[2:]
    src = path.read_text()
    if not keep_main:
        src = strip_main_block(src)
    ts = to_ts_embed(src)
    sys.stdout.write(ts)

if __name__ == "__main__":
    main()