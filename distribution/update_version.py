"""
Set or get the mfusg version in meson.build, pixi.toml, and src/mfusg.f.

The version banner in src/mfusg.f also gets the current UTC date.
"""

import argparse
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# file, pattern, and the group that holds the version
TARGETS = [
    (ROOT / "meson.build", re.compile(r"^(  version: ')([^']*)(',)", re.M)),
    (ROOT / "pixi.toml", re.compile(r'^(version = ")([^"]*)(")', re.M)),
    (ROOT / "src" / "mfusg.f", re.compile(r"(PARAMETER \(VERSION=')([^']*)('\))")),
]
VERSION = re.compile(r"^\d+\.\d+\.\d+(\.dev\d+)?$")


def read(fpth):
    with open(fpth, newline="") as f:
        return f.read()


def write(fpth, text):
    with open(fpth, "w", newline="") as f:
        f.write(text)


def get_version():
    fpth, pattern = TARGETS[0]
    match = pattern.search(read(fpth))
    if match is None:
        sys.exit(f"no version found in {fpth.name}")
    return match.group(2)


def next_version(version):
    major, minor = version.split(".")[:2]
    return f"{major}.{int(minor) + 1}.0.dev0"


def set_version(version):
    if VERSION.match(version) is None:
        sys.exit(f"{version} is not major.minor.patch or major.minor.patch.devN")
    date = datetime.now(timezone.utc).strftime("%m/%d/%Y")
    for fpth, pattern in TARGETS:
        text = read(fpth)
        value = f"{version} {date}" if fpth.suffix == ".f" else version
        text, n = pattern.subn(lambda m: m.group(1) + value + m.group(3), text)
        if n != 1:
            sys.exit(f"expected one version in {fpth.name}, found {n}")
        write(fpth, text)
        print(f"{fpth.relative_to(ROOT)}: {value}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-v", "--version", help="version to set")
    group.add_argument("-g", "--get", action="store_true", help="print the version")
    group.add_argument(
        "-n", "--next", action="store_true", help="print the next minor dev version"
    )
    args = parser.parse_args()
    if args.get:
        print(get_version())
    elif args.next:
        print(next_version(get_version()))
    else:
        set_version(args.version)
