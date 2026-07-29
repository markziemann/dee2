#!/usr/bin/env python3
"""
Parse an NCBI SRA efetch XML (EXPERIMENT_PACKAGE_SET) into a flat CSV.

Usage:
    python3 parse_sra_xml.py input.xml output.csv

Only uses the Python standard library (xml.etree.ElementTree, csv) —
no pysradb / rentrez / biopython.
"""

import sys
import os
import glob
import csv
import xml.etree.ElementTree as ET


def expand_input_paths(raw_paths):
    """
    Turn the input args into a concrete, sorted list of XML file paths.
    Each raw path can be:
      - a direct path to an .xml file
      - a directory (all *.xml files inside it are used)
      - a glob pattern (only if the shell didn't already expand it, e.g.
        because it was quoted) -- handled here as a fallback
    """
    resolved = []
    missing = []

    for raw in raw_paths:
        if os.path.isdir(raw):
            found = sorted(glob.glob(os.path.join(raw, "*.xml")))
            if not found:
                missing.append(f"{raw} (directory, but no *.xml files inside)")
            resolved.extend(found)
        elif os.path.isfile(raw):
            resolved.append(raw)
        else:
            # Might be an unexpanded glob pattern
            matches = sorted(glob.glob(raw))
            if matches:
                resolved.extend(matches)
            else:
                missing.append(raw)

    if missing:
        print("ERROR: the following input path(s) could not be found:", file=sys.stderr)
        for m in missing:
            print(f"  - {m}", file=sys.stderr)
        print(
            "Check that you're pointing at the right directory (e.g. the "
            "sra_batches_<date>/ folder from the fetch script) and that the "
            "path is correct relative to your current working directory.",
            file=sys.stderr,
        )
        sys.exit(1)

    if not resolved:
        print("ERROR: no input XML files resolved from the given arguments.", file=sys.stderr)
        sys.exit(1)

    return resolved


def get_text(elem, path, default=""):
    node = elem.find(path)
    return node.text.strip() if node is not None and node.text else default


def get_external_ids(elem, namespace_wanted):
    """
    Pull EXTERNAL_ID values matching a given namespace (e.g. 'GEO')
    from an IDENTIFIERS block. A single record can have zero, one,
    or (rarely) multiple matching IDs; return them semicolon-joined.
    """
    ids = []
    for ext_id in elem.findall("./IDENTIFIERS/EXTERNAL_ID"):
        ns = ext_id.get("namespace", "")
        if ns.upper() == namespace_wanted.upper():
            if ext_id.text:
                ids.append(ext_id.text.strip())
    return ";".join(ids)


def parse_package(pkg):
    """Extract one row per RUN within an EXPERIMENT_PACKAGE."""
    rows = []

    study = pkg.find("./STUDY")
    experiment = pkg.find("./EXPERIMENT")
    sample = pkg.find("./SAMPLE")
    submission = pkg.find("./SUBMISSION")
    runset = pkg.find("./RUN_SET")

    # --- STUDY level ---
    study_acc = study.get("accession", "") if study is not None else ""
    study_title = get_text(study, "./DESCRIPTOR/STUDY_TITLE")
    geo_series = get_external_ids(study, "GEO") if study is not None else ""

    # --- SUBMISSION level ---
    submission_acc = submission.get("accession", "") if submission is not None else ""

    # --- SAMPLE level ---
    sample_acc = sample.get("accession", "") if sample is not None else ""
    geo_sample = get_external_ids(sample, "GEO") if sample is not None else ""
    organism = get_text(sample, "./SAMPLE_NAME/SCIENTIFIC_NAME")
    sample_title = get_text(sample, "./TITLE")

    # --- EXPERIMENT level ---
    experiment_acc = experiment.get("accession", "") if experiment is not None else ""
    experiment_title = get_text(experiment, "./TITLE")
    library_strategy = get_text(experiment, "./DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY")
    library_source = get_text(experiment, "./DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SOURCE")
    library_selection = get_text(experiment, "./DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SELECTION")
    layout_elem = experiment.find("./DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT") if experiment is not None else None
    if layout_elem is not None and len(list(layout_elem)) > 0:
        library_layout = layout_elem[0].tag  # SINGLE or PAIRED
    else:
        library_layout = ""
    platform = ""
    if experiment is not None:
        platform_elem = experiment.find("./PLATFORM")
        if platform_elem is not None and len(list(platform_elem)) > 0:
            instrument_model_elem = platform_elem[0].find("INSTRUMENT_MODEL")
            platform = instrument_model_elem.text.strip() if instrument_model_elem is not None and instrument_model_elem.text else platform_elem[0].tag

    # --- RUN level (one row per run) ---
    runs = runset.findall("./RUN") if runset is not None else []
    if not runs:
        # Still emit one row even if no RUN element is present
        runs = [None]

    for run in runs:
        run_acc = run.get("accession", "") if run is not None else ""
        run_date = run.get("run_date", "") if run is not None else ""
        spots = run.get("total_spots", "") if run is not None else ""
        bases = run.get("total_bases", "") if run is not None else ""
        run_size = run.get("size", "") if run is not None else ""

        rows.append({
            "geo_series": geo_series,
            "geo_sample": geo_sample,
            "study_accession": study_acc,
            "study_title": study_title,
            "submission_accession": submission_acc,
            "sample_accession": sample_acc,
            "sample_title": sample_title,
            "organism": organism,
            "experiment_accession": experiment_acc,
            "experiment_title": experiment_title,
            "library_strategy": library_strategy,
            "library_source": library_source,
            "library_selection": library_selection,
            "library_layout": library_layout,
            "platform": platform,
            "run_accession": run_acc,
            "run_date": run_date,
            "total_spots": spots,
            "total_bases": bases,
            "run_size_bytes": run_size,
        })

    return rows


def main():
    if len(sys.argv) < 3:
        print(
            f"Usage: {sys.argv[0]} input1.xml [input2.xml ...] output.csv\n"
            f"   or: {sys.argv[0]} batch_directory/ output.csv",
            file=sys.stderr,
        )
        sys.exit(1)

    # Last arg is the output CSV; everything before it is input
    # (a mix of files, directories, or glob patterns is fine)
    out_path = sys.argv[-1]
    raw_in = sys.argv[1:-1]
    in_paths = expand_input_paths(raw_in)
    print(f"Resolved {len(in_paths)} input XML file(s):", file=sys.stderr)
    for p in in_paths:
        print(f"  - {p}", file=sys.stderr)

    fieldnames = [
        "geo_series", "geo_sample",
        "study_accession", "study_title",
        "submission_accession",
        "sample_accession", "sample_title", "organism",
        "experiment_accession", "experiment_title",
        "library_strategy", "library_source", "library_selection", "library_layout",
        "platform",
        "run_accession", "run_date", "total_spots", "total_bases", "run_size_bytes",
    ]

    all_rows = []
    total_packages = 0

    for in_path in in_paths:
        try:
            tree = ET.parse(in_path)
        except ET.ParseError as e:
            print(f"WARNING: {in_path} is not valid XML ({e}) -- skipping", file=sys.stderr)
            continue
        except OSError as e:
            print(f"WARNING: could not open {in_path} ({e}) -- skipping", file=sys.stderr)
            continue

        root = tree.getroot()
        packages = root.findall("./EXPERIMENT_PACKAGE")
        if not packages:
            print(f"WARNING: no EXPERIMENT_PACKAGE elements found in {in_path}", file=sys.stderr)
            continue

        total_packages += len(packages)
        for pkg in packages:
            all_rows.extend(parse_package(pkg))

    if not all_rows:
        print("No rows extracted from any input file.", file=sys.stderr)
        sys.exit(1)

    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_rows)

    print(
        f"Wrote {len(all_rows)} rows ({total_packages} experiment packages "
        f"from {len(in_paths)} input file(s)) to {out_path}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
