#!/usr/bin/env python3
"""
semantic_commit.py — generate Conventional Commits messages via local Ollama.

Usage:
    python tools/semantic_commit.py [--model MODEL] [--commit]

Requirements:
    - Ollama running locally (https://ollama.com)
    - A pulled model, e.g.: ollama pull llama3.2
"""

import argparse
import json
import subprocess
import sys
import textwrap
import urllib.error
import urllib.request

OLLAMA_URL = "http://localhost:11434/api/chat"
DEFAULT_MODEL = "qwen3-coder:30b-a3b-q4_K_M"
MAX_DIFF_CHARS = 12_000

SYSTEM_PROMPT = """\
You are an expert at writing Conventional Commits v1.0.0 messages.

## Format

<type>[(<scope>)][!]: <description>

[body]

[footers]

## Type selection

| Type     | Use when                                   |
|----------|--------------------------------------------|
| feat     | A new feature is introduced                |
| fix      | A bug is corrected                         |
| perf     | Performance improvement, no API change     |
| refactor | Code restructure, no feature or fix        |
| docs     | Documentation changes only                 |
| style    | Formatting, whitespace (no logic change)   |
| test     | Test additions or updates only             |
| chore    | Build, tooling, dependency bumps           |
| ci       | CI/CD pipeline configuration changes       |
| build    | Build system or compiler option changes    |
| revert   | Reverts a previous commit                  |

## Rules

- Subject line: max 72 chars · imperative mood · lowercase · no trailing period
- Scope: optional, identifies the module/component most affected; omit if cross-cutting
- Breaking change: append ! after type/scope AND add a BREAKING CHANGE: footer
- Body: explain the *why*, not the *what*; wrap at 90 chars; blank line before body
- Footers: each on its own line after a blank line (BREAKING CHANGE:, Closes #N, Refs #N)
- Never add co-authors

Output ONLY the raw commit message — no explanations, no markdown fences.\
"""


def run_git(*args: str) -> str:
    result = subprocess.run(["git", *args], capture_output=True, text=True)
    if result.returncode != 0:
        print(f"git error: {result.stderr.strip()}", file=sys.stderr)
        sys.exit(1)
    return result.stdout.strip()


def staged_stat() -> str:
    return run_git("diff", "--cached", "--stat")


def staged_diff() -> str:
    diff = run_git("diff", "--cached", "--unified=3")
    if len(diff) > MAX_DIFF_CHARS:
        diff = diff[:MAX_DIFF_CHARS] + "\n\n[diff truncated — showing first 12 000 chars]"
    return diff


def recent_commits(n: int = 10) -> str:
    return run_git("log", "--oneline", f"-{n}")


def ask_ollama(model: str, user_prompt: str) -> str:
    payload = json.dumps({
        "model": model,
        "messages": [
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user",   "content": user_prompt},
        ],
        "stream": True,
    }).encode()

    req = urllib.request.Request(
        OLLAMA_URL,
        data=payload,
        headers={"Content-Type": "application/json"},
    )

    collected = []
    try:
        with urllib.request.urlopen(req, timeout=120) as resp:
            for raw_line in resp:
                line = raw_line.decode().strip()
                if not line:
                    continue
                chunk = json.loads(line)
                token = chunk.get("message", {}).get("content", "")
                if token:
                    print(token, end="", flush=True)
                    collected.append(token)
                if chunk.get("done"):
                    break
    except urllib.error.URLError as exc:
        print(f"\nError: cannot reach Ollama at {OLLAMA_URL}\n{exc}", file=sys.stderr)
        sys.exit(1)

    print()  # final newline
    return "".join(collected).strip()


def wrap_message(message: str, width: int = 90) -> str:
    """Wrap each line of a commit message to at most `width` columns."""
    lines = message.splitlines()
    wrapped = []
    for line in lines:
        if len(line) <= width:
            wrapped.append(line)
            continue
        indent = len(line) - len(line.lstrip())
        prefix = line[:indent]
        wrapped.append(textwrap.fill(
            line,
            width=width,
            initial_indent=prefix,
            subsequent_indent=prefix,
            break_long_words=False,
            break_on_hyphens=False,
        ))
    return "\n".join(wrapped)


def build_prompt(stat: str, diff: str, commits: str) -> str:
    return (
        "## Staged changes (stat)\n\n"
        f"{stat}\n\n"
        "## Staged diff\n\n"
        f"{diff}\n\n"
        "## Recent commits (style reference)\n\n"
        f"{commits}\n\n"
        "Generate the semantic commit message for the staged changes above."
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate a Conventional Commits message via local Ollama."
    )
    parser.add_argument(
        "--model", default=DEFAULT_MODEL,
        help=f"Ollama model to use (default: {DEFAULT_MODEL})"
    )
    parser.add_argument(
        "--commit", action="store_true",
        help="Run `git commit` with the generated message after review"
    )
    args = parser.parse_args()

    stat = staged_stat()
    if not stat:
        print("Nothing staged. Run `git add <files>` first.", file=sys.stderr)
        sys.exit(1)

    diff    = staged_diff()
    commits = recent_commits()
    prompt  = build_prompt(stat, diff, commits)

    print(f"[{args.model}] Generating commit message...\n", file=sys.stderr)
    message = wrap_message(ask_ollama(args.model, prompt))

    if args.commit:
        print("\n--- Confirm commit? [y/N] ", end="", flush=True)
        if input().strip().lower() == "y":
            subprocess.run(["git", "commit", "-m", message])
        else:
            print("Aborted.", file=sys.stderr)


if __name__ == "__main__":
    main()
