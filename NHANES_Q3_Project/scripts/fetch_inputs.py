"""Recover missing public NCHS inputs listed in the supplied source manifest.

Existing inputs with the recorded size and expected format marker are reused.
No hashes, credentials, alternate data sources, or outcome-dependent choices.
"""
from pathlib import Path
import argparse
import json
import urllib.request

ROOT = Path(__file__).resolve().parents[1]


def looks_complete(path, entry):
    if not path.is_file() or path.stat().st_size != entry["bytes"]:
        return False
    with path.open("rb") as handle:
        head = handle.read(160)
    if path.suffix.lower() == ".xpt":
        return head.startswith(b"HEADER RECORD")
    return bool(head.strip()) and b"<html" not in head.lower()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true", help="List sources without network requests")
    args = parser.parse_args()
    manifest = json.loads((ROOT / "sources/source_manifest.json").read_text())
    for entry in manifest["data_files"]:
        dest = ROOT / entry["local_path"]
        if args.list:
            print(entry["file"], entry["data_url"])
            continue
        if looks_complete(dest, entry):
            print("REUSE", entry["file"])
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        request = urllib.request.Request(entry["data_url"], headers={"User-Agent": "NHANES-public-research-input-fetch/1.0"})
        partial = dest.with_suffix(dest.suffix + ".part")
        try:
            with urllib.request.urlopen(request, timeout=60) as response, partial.open("wb") as handle:
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    handle.write(chunk)
            if not looks_complete(partial, dict(entry, file=entry["file"])):
                # .part changes the suffix; check the XPT marker explicitly too.
                raise RuntimeError("Downloaded file differs from the recorded size or has an HTML response. Check official source before replacing inputs.")
            if dest.suffix.lower() == ".xpt":
                with partial.open("rb") as handle:
                    if not handle.read(80).startswith(b"HEADER RECORD"):
                        raise RuntimeError("Expected a SAS XPORT file, received another format")
            partial.replace(dest)
            print("DOWNLOADED", entry["file"])
        except Exception:
            partial.unlink(missing_ok=True)
            raise


if __name__ == "__main__":
    main()
