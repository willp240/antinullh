#!/usr/bin/env python3
"""Audit and optionally resubmit incomplete fixed-oscillation Condor jobs."""

import argparse
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


FIT_VALID_RE = re.compile(rb"(?m)^\s*FitValid:\s*[01]\s*$")
EXEC_RE = re.compile(r"^\s*(?!#).*bin/fixedosc_fit(?:\s|$)")
LOG_EVENT_RE = re.compile(r"(?m)^(\d{3})\s+\(")


@dataclass
class Job:
    name: str
    submit: Path
    expected: int
    completed: int
    active: bool


def count_expected(script: Path) -> int:
    with script.open(errors="replace") as handle:
        return sum(bool(EXEC_RE.search(line)) for line in handle)


def count_completed(output: Path) -> int:
    if not output.is_file():
        return 0
    return len(FIT_VALID_RE.findall(output.read_bytes()))


def appears_active(log: Path) -> bool:
    """Treat a job as active unless its latest event is terminated or removed."""
    if not log.is_file():
        return False
    events = LOG_EVENT_RE.findall(log.read_text(errors="replace"))
    # 005: terminated; 009: removed (for example, via condor_rm).
    return bool(events) and events[-1] not in {"005", "009"}


def audit(fit_dir: Path) -> list[Job]:
    submit_dir = fit_dir / "submit"
    if not submit_dir.is_dir():
        raise ValueError(f"submit directory not found: {submit_dir}")

    jobs = []
    for submit in sorted(submit_dir.glob("*.submit")):
        name = submit.stem
        script = fit_dir / "sh" / f"{name}.sh"
        if not script.is_file():
            raise ValueError(f"job script not found: {script}")
        expected = count_expected(script)
        if expected == 0:
            raise ValueError(f"no fixedosc_fit commands found in: {script}")
        jobs.append(
            Job(
                name=name,
                submit=submit,
                expected=expected,
                completed=count_completed(fit_dir / "output" / f"{name}.output"),
                active=appears_active(fit_dir / "log" / f"{name}.log"),
            )
        )
    if not jobs:
        raise ValueError(f"no submit files found in: {submit_dir}")
    return jobs


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Find fixed-oscillation jobs whose output has fewer FitValid records "
            "than fixedosc_fit commands, and optionally resubmit their saved Condor files."
        )
    )
    parser.add_argument("fit_dir", type=Path, help="top-level fixed-oscillation fit directory")
    parser.add_argument(
        "--submit",
        action="store_true",
        help="actually call condor_submit (default: audit/dry-run only)",
    )
    parser.add_argument(
        "--include-active",
        action="store_true",
        help="also select jobs whose latest Condor log event is not termination",
    )
    parser.add_argument(
        "--batch-name",
        help="Condor batch name (default: <fit-directory-name>_recovery)",
    )
    args = parser.parse_args()

    fit_dir = args.fit_dir.expanduser().resolve()
    try:
        jobs = audit(fit_dir)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))

    incomplete = [job for job in jobs if job.completed != job.expected]
    active = [job for job in incomplete if job.active]
    selected = incomplete if args.include_active else [job for job in incomplete if not job.active]

    print(f"Fit directory: {fit_dir}")
    print(f"Jobs: {len(jobs)} total, {len(jobs) - len(incomplete)} complete, "
          f"{len(incomplete)} incomplete, {len(active)} apparently active")
    for job in incomplete:
        state = "active; skipped" if job.active and not args.include_active else "selected"
        print(f"  {job.name}: {job.completed}/{job.expected} ({state})")

    if not args.submit:
        print(f"Dry run: {len(selected)} job(s) would be resubmitted. Add --submit to proceed.")
        return 0

    if not selected:
        print("No jobs selected for resubmission.")
        return 0

    batch_name = args.batch_name or f"{fit_dir.name}_recovery"
    failures = 0
    for job in selected:
        command = ["condor_submit", "-batch-name", batch_name, str(job.submit)]
        print("Submitting:", " ".join(command), flush=True)
        result = subprocess.run(command, check=False)
        if result.returncode:
            failures += 1

    if failures:
        print(f"ERROR: condor_submit failed for {failures} of {len(selected)} jobs.", file=sys.stderr)
        return 1
    print(f"Submitted {len(selected)} job(s) in batch {batch_name!r}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
