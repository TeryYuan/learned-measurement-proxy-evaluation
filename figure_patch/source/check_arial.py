#!/usr/bin/env python3
"""Find actual installed Arial faces without downloading or copying font files.

Writes a private path-bearing configuration and a path-free review summary.
Exit 0: required faces found; 2: faces missing; 3: dependency/setup failure.
This is a preflight only: final SVG/PDF fonts and glyphs still need inspection.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import tempfile


FACE_ORDER = ("regular", "bold", "italic", "bold_italic")
EXPECTED_NAMES = {
    "regular": "arial.ttf", "bold": "arialbd.ttf",
    "italic": "ariali.ttf", "bold_italic": "arialbi.ttf",
}
CHECK_GLYPHS = {
    "minus": "−", "beta": "β", "lambda": "λ", "gamma": "γ",
    "delta": "Δ", "mu": "μ", "greater_equal": "≥", "plus_minus": "±",
    "right_arrow": "→", "combining_hat": "\u0302",
}


def save_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def font_roots(extra: list[Path]) -> list[Path]:
    roots = list(extra)
    if os.name == "nt":
        windir = os.environ.get("WINDIR") or os.environ.get("SystemRoot")
        if windir:
            roots.append(Path(windir) / "Fonts")
        local_app_data = os.environ.get("LOCALAPPDATA")
        if local_app_data:
            roots.append(Path(local_app_data) / "Microsoft" / "Windows" / "Fonts")
    else:
        roots.extend([
            Path("/usr/share/fonts"), Path("/usr/local/share/fonts"),
            Path.home() / ".fonts", Path.home() / ".local/share/fonts",
            Path("/Library/Fonts"), Path("/System/Library/Fonts"),
            Path.home() / "Library/Fonts",
            Path("/mnt/c/Windows/Fonts"),
        ])
    return [p for p in roots if p.is_dir()]


def classify_face(style_name: str) -> str | None:
    style = style_name.casefold().replace("-", " ")
    if any(t in style for t in ("light", "black", "narrow", "condensed")):
        return None
    bold = "bold" in style
    italic = "italic" in style or "oblique" in style
    if bold and italic:
        return "bold_italic"
    if bold:
        return "bold"
    if italic:
        return "italic"
    if style.strip() in ("regular", "normal", "roman", "book", ""):
        return "regular"
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--font-dir", action="append", type=Path, default=[],
                        help="Additional existing local font folder; may repeat.")
    parser.add_argument("--private-config", type=Path,
                        default=Path("private/arial_fonts.local.json"))
    parser.add_argument("--public-summary", type=Path,
                        default=Path("validation/arial_preflight_summary.json"))
    args = parser.parse_args()
    if args.private_config.resolve() == args.public_summary.resolve():
        parser.error("Private config and public summary must be different files.")

    # Avoid creating a user-profile font cache or exposing that path in warnings.
    font_cache = tempfile.TemporaryDirectory(prefix="patterns_font_check_")
    os.environ["MPLCONFIGDIR"] = font_cache.name
    try:
        from matplotlib import font_manager, ft2font
        import matplotlib
    except ImportError:
        summary = {
            "status": "DEPENDENCY_UNAVAILABLE", "dependency": "matplotlib",
            "action": "Run in the project's existing plotting environment.",
            "private_paths_included": False, "font_files_copied": False,
        }
        save_json(args.public_summary, summary)
        print(json.dumps(summary))
        return 3

    candidates: set[Path] = set()
    try:
        candidates.update(Path(p) for p in font_manager.findSystemFonts())
    except (OSError, RuntimeError):
        pass
    for root in font_roots(args.font_dir):
        try:
            candidates.update(p for p in root.rglob("*")
                              if p.is_file() and p.suffix.casefold() in (".ttf", ".otf", ".ttc"))
        except OSError:
            continue

    found: dict[str, list[tuple[Path, dict]]] = {k: [] for k in FACE_ORDER}
    unreadable = 0
    for path in sorted(candidates):
        try:
            font = ft2font.FT2Font(str(path))
            family = str(font.family_name)
            if family.casefold().strip() != "arial":
                continue
            face = classify_face(str(font.style_name))
            if face is None:
                continue
            cmap = font.get_charmap()
            info = {
                "family": family, "style": str(font.style_name),
                "postscript_name": str(font.postscript_name),
                "missing_requested_glyphs": [name for name, glyph in CHECK_GLYPHS.items()
                                             if ord(glyph) not in cmap],
            }
            found[face].append((path.resolve(), info))
        except (OSError, RuntimeError, ValueError):
            unreadable += 1

    chosen = {}
    for face, choices in found.items():
        if choices:
            choices.sort(key=lambda v: (v[0].name.casefold() != EXPECTED_NAMES[face], str(v[0])))
            chosen[face] = choices[0]
    required = ("regular", "bold", "italic")
    missing = [face for face in required if face not in chosen]
    status = "ARIAL_FACES_FOUND" if not missing else "ARIAL_REQUIRED_FACES_MISSING"
    private_config = {
        "status": status,
        "private_do_not_publish": True,
        "faces": {face: {"path": str(path), **info} for face, (path, info) in chosen.items()},
        "missing_required_faces": missing,
        "integration": "Register each path with font_manager.fontManager.addfont; use family Arial and custom Arial mathtext. Final PDF font inspection remains required.",
        "no_redistribution": "Font files remain in their existing installation; do not copy this private configuration or font files to publication_patch.",
    }
    summary = {
        "status": status, "matplotlib_version": matplotlib.__version__,
        "actual_family_required": "Arial", "missing_required_faces": missing,
        "faces": {face: info for face, (_, info) in chosen.items()},
        "font_files_copied": False, "font_files_downloaded": False,
        "private_paths_included": False, "unreadable_font_count": unreadable,
        "scope": "Installed-font discovery and requested-glyph coverage only; not a final export validation or redistribution license grant.",
        "next_step": "Use selected installed faces and inspect exported PDF/SVG fonts and glyphs." if not missing else
                     "Continue text/layout preparation. Obtain access to an existing lawful Arial installation before claiming Arial-compliant final exports; do not silently substitute another family.",
    }
    save_json(args.private_config, private_config)
    save_json(args.public_summary, summary)
    print(json.dumps({"status": status, "found_faces": list(chosen),
                      "missing_required_faces": missing,
                      "private_config_written": True,
                      "public_summary_written": True}, ensure_ascii=False))
    return 0 if not missing else 2


if __name__ == "__main__":
    sys.exit(main())
