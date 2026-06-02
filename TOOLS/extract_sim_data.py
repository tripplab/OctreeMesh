#!/usr/bin/env python3
"""Extract summary statistics from GiD ``.post.res`` simulation output.

The script focuses on the OctreeMesh solver results most commonly needed for
post-processing: nodal Displacement vectors and Von Mises scalar values on
Gauss points. It intentionally uses only the Python standard library so it can
run in a freshly checked-out repository without additional package installs.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import shlex
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO

DEFAULT_RESULTS = ("Displacement", "Von Mises")
DEFAULT_PERCENTILES = (5.0, 25.0, 50.0, 75.0, 95.0)


@dataclass(frozen=True)
class ResultHeader:
    name: str
    analysis: str
    step: str
    result_type: str
    location: str
    gauss_point_set: str | None = None


@dataclass(frozen=True)
class GaussPointSet:
    name: str
    element_type: str | None
    count: int


class ParseError(RuntimeError):
    """Raised when the input result file is malformed in a fatal way."""


class WarningCollector:
    def __init__(self, strict: bool) -> None:
        self.strict = strict
        self.count = 0

    def warn(self, line_number: int, message: str) -> None:
        if self.strict:
            raise ParseError(f"line {line_number}: {message}")
        self.count += 1
        print(f"warning: line {line_number}: {message}", file=sys.stderr)


def parse_result_header(line: str) -> ResultHeader | None:
    """Parse a GiD Result header, preserving quoted fields."""
    stripped = line.strip()
    if not stripped.startswith("Result "):
        return None

    try:
        tokens = shlex.split(stripped)
    except ValueError as exc:
        raise ParseError(f"invalid Result header: {exc}") from exc

    if len(tokens) < 6 or tokens[0] != "Result":
        raise ParseError(f"invalid Result header: {line.rstrip()}")

    return ResultHeader(
        name=tokens[1],
        analysis=tokens[2],
        step=tokens[3],
        result_type=tokens[4],
        location=tokens[5],
        gauss_point_set=tokens[6] if len(tokens) > 6 else None,
    )


def parse_gauss_points_header(line: str) -> tuple[str, str | None] | None:
    stripped = line.strip()
    if not stripped.startswith("GaussPoints "):
        return None

    try:
        tokens = shlex.split(stripped)
    except ValueError as exc:
        raise ParseError(f"invalid GaussPoints header: {exc}") from exc

    if len(tokens) < 2 or tokens[0] != "GaussPoints":
        raise ParseError(f"invalid GaussPoints header: {line.rstrip()}")

    element_type = None
    for index, token in enumerate(tokens):
        if token == "ElemType" and index + 1 < len(tokens):
            element_type = tokens[index + 1]
            break

    return tokens[1], element_type


def is_int_token(token: str) -> bool:
    try:
        int(token)
    except ValueError:
        return False
    return True


def parse_float_tokens(
    tokens: Iterable[str], line_number: int, warnings: WarningCollector
) -> list[float]:
    values: list[float] = []
    for token in tokens:
        try:
            values.append(float(token))
        except ValueError:
            warnings.warn(line_number, f"could not parse float value {token!r}")
    return values


def percentile(sorted_values: list[float], percent: float) -> float:
    if not sorted_values:
        raise ValueError("cannot compute percentile of an empty data set")
    if len(sorted_values) == 1:
        return sorted_values[0]

    position = (len(sorted_values) - 1) * (percent / 100.0)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return sorted_values[int(position)]
    fraction = position - lower
    return sorted_values[lower] * (1.0 - fraction) + sorted_values[upper] * fraction


def summarize_values(values: list[float], percentiles: tuple[float, ...]) -> dict[str, float | int]:
    if not values:
        return {"count": 0}

    count = len(values)
    mean = math.fsum(values) / count
    variance = math.fsum((value - mean) ** 2 for value in values) / count
    sorted_values = sorted(values)

    summary: dict[str, float | int] = {
        "count": count,
        "min": sorted_values[0],
        "max": sorted_values[-1],
        "mean": mean,
        "stddev_population": math.sqrt(variance),
        "median": percentile(sorted_values, 50.0),
    }

    for percent in percentiles:
        label = f"p{int(percent):02d}" if percent.is_integer() else f"p{percent:g}"
        summary[label] = percentile(sorted_values, percent)

    return summary


def parse_percentiles(raw: str) -> tuple[float, ...]:
    percentiles: list[float] = []
    for item in raw.split(","):
        item = item.strip()
        if not item:
            continue
        try:
            value = float(item)
        except ValueError as exc:
            raise argparse.ArgumentTypeError(f"invalid percentile {item!r}") from exc
        if value < 0.0 or value > 100.0:
            raise argparse.ArgumentTypeError("percentiles must be between 0 and 100")
        percentiles.append(value)
    if not percentiles:
        raise argparse.ArgumentTypeError("at least one percentile is required")
    return tuple(percentiles)


def parse_displacement_block(
    lines: TextIO,
    start_line_number: int,
    warnings: WarningCollector,
    percentiles: tuple[float, ...],
) -> tuple[dict[str, object], int]:
    ux_values: list[float] = []
    uy_values: list[float] = []
    uz_values: list[float] = []
    magnitudes: list[float] = []
    min_location: dict[str, float | int] | None = None
    max_location: dict[str, float | int] | None = None
    min_magnitude = math.inf
    max_magnitude = -math.inf
    line_number = start_line_number

    for line_number, line in enumerate(lines, start=start_line_number):
        stripped = line.strip()
        if not stripped:
            continue
        if stripped == "End Values":
            break

        tokens = stripped.split()
        if len(tokens) < 4:
            warnings.warn(line_number, "Displacement row has fewer than 4 fields")
            continue
        if not is_int_token(tokens[0]):
            warnings.warn(line_number, f"invalid Displacement node id {tokens[0]!r}")
            continue

        node_id = int(tokens[0])
        components = parse_float_tokens(tokens[1:4], line_number, warnings)
        if len(components) != 3:
            warnings.warn(line_number, "Displacement row does not contain 3 valid components")
            continue
        if len(tokens) > 4:
            warnings.warn(line_number, "ignoring extra fields in Displacement row")

        ux, uy, uz = components
        magnitude = math.sqrt(ux * ux + uy * uy + uz * uz)
        ux_values.append(ux)
        uy_values.append(uy)
        uz_values.append(uz)
        magnitudes.append(magnitude)

        location = {
            "node_id": node_id,
            "ux": ux,
            "uy": uy,
            "uz": uz,
            "magnitude": magnitude,
        }
        if magnitude < min_magnitude:
            min_magnitude = magnitude
            min_location = location
        if magnitude > max_magnitude:
            max_magnitude = magnitude
            max_location = location
    else:
        raise ParseError("Displacement Values block ended before 'End Values'")

    return (
        {
            "type": "Vector",
            "location": "OnNodes",
            "count": len(magnitudes),
            "components": {
                "x": summarize_values(ux_values, percentiles),
                "y": summarize_values(uy_values, percentiles),
                "z": summarize_values(uz_values, percentiles),
            },
            "magnitude": summarize_values(magnitudes, percentiles),
            "min_location": min_location,
            "max_location": max_location,
        },
        line_number,
    )


def parse_von_mises_block(
    lines: TextIO,
    start_line_number: int,
    warnings: WarningCollector,
    percentiles: tuple[float, ...],
    gauss_point_count: int,
) -> tuple[dict[str, object], int]:
    values: list[float] = []
    element_count = 0
    current_element: int | None = None
    current_gauss_index = 0
    min_location: dict[str, float | int] | None = None
    max_location: dict[str, float | int] | None = None
    min_value = math.inf
    max_value = -math.inf
    line_number = start_line_number

    def add_value(value: float) -> None:
        nonlocal current_gauss_index, min_value, max_value, min_location, max_location
        if current_element is None:
            raise ParseError(f"line {line_number}: scalar value found before an element id")
        current_gauss_index += 1
        if current_gauss_index > gauss_point_count:
            raise ParseError(
                f"line {line_number}: too many Von Mises values for element {current_element}"
            )
        values.append(value)
        location = {
            "element_id": current_element,
            "gauss_point": current_gauss_index,
            "value": value,
        }
        if value < min_value:
            min_value = value
            min_location = location
        if value > max_value:
            max_value = value
            max_location = location

    for line_number, line in enumerate(lines, start=start_line_number):
        stripped = line.strip()
        if not stripped:
            continue
        if stripped == "End Values":
            if current_element is not None and current_gauss_index != gauss_point_count:
                warnings.warn(
                    line_number,
                    f"element {current_element} has {current_gauss_index} "
                    f"Von Mises values; expected {gauss_point_count}",
                )
            break

        tokens = stripped.split()
        if current_element is None or current_gauss_index == gauss_point_count:
            if len(tokens) < 2:
                warnings.warn(line_number, "Von Mises element row has fewer than 2 fields")
                current_element = None
                current_gauss_index = 0
                continue
            if not is_int_token(tokens[0]):
                warnings.warn(line_number, f"invalid Von Mises element id {tokens[0]!r}")
                current_element = None
                current_gauss_index = 0
                continue
            current_element = int(tokens[0])
            current_gauss_index = 0
            element_count += 1
            scalar_tokens = tokens[1:]
        else:
            scalar_tokens = tokens

        scalar_values = parse_float_tokens(scalar_tokens, line_number, warnings)
        for scalar_value in scalar_values:
            add_value(scalar_value)
    else:
        raise ParseError("Von Mises Values block ended before 'End Values'")

    return (
        {
            "type": "Scalar",
            "location": "OnGaussPoints",
            "gauss_points_per_element": gauss_point_count,
            "element_count": element_count,
            "count": len(values),
            "values": summarize_values(values, percentiles),
            "min_location": min_location,
            "max_location": max_location,
        },
        line_number,
    )


def skip_values_block(lines: TextIO, start_line_number: int) -> int:
    line_number = start_line_number
    for line_number, line in enumerate(lines, start=start_line_number):
        if line.strip() == "End Values":
            return line_number
    raise ParseError("Values block ended before 'End Values'")


def parse_result_file(
    path: Path,
    requested_results: set[str],
    percentiles: tuple[float, ...],
    strict: bool,
) -> dict[str, object]:
    warnings = WarningCollector(strict)
    gauss_points: dict[str, GaussPointSet] = {}
    results: dict[str, object] = {}

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        iterator = iter(handle)
        line_number = 0
        for line_number, line in enumerate(iterator, start=1):
            gp_header = parse_gauss_points_header(line)
            if gp_header is not None:
                name, element_type = gp_header
                count: int | None = None
                for line_number, gp_line in enumerate(iterator, start=line_number + 1):
                    stripped = gp_line.strip()
                    if stripped.startswith("Number Of Gauss Points:"):
                        raw_count = stripped.split(":", 1)[1].strip()
                        try:
                            count = int(raw_count)
                        except ValueError as exc:
                            raise ParseError(
                                f"line {line_number}: invalid Gauss point count {raw_count!r}"
                            ) from exc
                    if stripped == "End GaussPoints":
                        break
                else:
                    raise ParseError(f"GaussPoints {name!r} ended before 'End GaussPoints'")
                if count is None:
                    warnings.warn(line_number, f"GaussPoints {name!r} has no count")
                else:
                    gauss_points[name] = GaussPointSet(name, element_type, count)
                continue

            header = parse_result_header(line)
            if header is None:
                continue

            for line_number, values_line in enumerate(iterator, start=line_number + 1):
                if values_line.strip() == "Values":
                    break
            else:
                raise ParseError(f"Result {header.name!r} ended before 'Values'")

            if header.name not in requested_results:
                line_number = skip_values_block(iterator, line_number + 1)
                continue

            if header.name == "Displacement":
                summary, line_number = parse_displacement_block(
                    iterator, line_number + 1, warnings, percentiles
                )
                results[header.name] = {
                    **summary,
                    "analysis": header.analysis,
                    "step": header.step,
                }
            elif header.name == "Von Mises":
                if header.gauss_point_set is None:
                    raise ParseError("Von Mises result does not reference a Gauss point set")
                gp_set = gauss_points.get(header.gauss_point_set)
                if gp_set is None:
                    raise ParseError(
                        f"Von Mises references unknown Gauss point set {header.gauss_point_set!r}"
                    )
                summary, line_number = parse_von_mises_block(
                    iterator, line_number + 1, warnings, percentiles, gp_set.count
                )
                results[header.name] = {
                    **summary,
                    "analysis": header.analysis,
                    "step": header.step,
                    "gauss_point_set": header.gauss_point_set,
                }
            else:
                line_number = skip_values_block(iterator, line_number + 1)

    missing = sorted(requested_results - set(results))
    if missing:
        raise ParseError(f"requested result(s) not found: {', '.join(missing)}")

    return {
        "file": str(path),
        "gauss_points": {
            name: {
                "element_type": gp_set.element_type,
                "number_of_gauss_points": gp_set.count,
            }
            for name, gp_set in gauss_points.items()
        },
        "results": results,
        "warnings": warnings.count,
    }


def format_float(value: object) -> str:
    if isinstance(value, float):
        return f"{value:.8e}"
    return str(value)


def format_summary_block(summary: dict[str, object], indent: str = "    ") -> list[str]:
    order = ["count", "min", "max", "mean", "stddev_population", "median"]
    ordered_keys = [key for key in order if key in summary]
    ordered_keys.extend(key for key in summary if key not in ordered_keys)
    return [f"{indent}{key}: {format_float(summary[key])}" for key in ordered_keys]


def format_text(report: dict[str, object]) -> str:
    lines = [f"File: {report['file']}", ""]
    results = report["results"]
    if not isinstance(results, dict):
        return ""

    displacement = results.get("Displacement")
    if isinstance(displacement, dict):
        lines.extend(
            [
                "Result: Displacement",
                f"  Analysis: {displacement.get('analysis')}",
                f"  Step: {displacement.get('step')}",
                f"  Location: {displacement.get('location')}",
                f"  Count: {displacement.get('count')}",
                "  Components:",
            ]
        )
        components = displacement.get("components", {})
        if isinstance(components, dict):
            for axis in ("x", "y", "z"):
                axis_summary = components.get(axis)
                if isinstance(axis_summary, dict):
                    lines.append(f"    {axis}:")
                    lines.extend(format_summary_block(axis_summary, indent="      "))
        magnitude = displacement.get("magnitude")
        if isinstance(magnitude, dict):
            lines.append("  Magnitude:")
            lines.extend(format_summary_block(magnitude, indent="    "))
        for label in ("min_location", "max_location"):
            location = displacement.get(label)
            if isinstance(location, dict):
                lines.append(f"  {label}:")
                for key, value in location.items():
                    lines.append(f"    {key}: {format_float(value)}")
        lines.append("")

    von_mises = results.get("Von Mises")
    if isinstance(von_mises, dict):
        lines.extend(
            [
                "Result: Von Mises",
                f"  Analysis: {von_mises.get('analysis')}",
                f"  Step: {von_mises.get('step')}",
                f"  Location: {von_mises.get('location')}",
                f"  Gauss point set: {von_mises.get('gauss_point_set')}",
                f"  Gauss points per element: {von_mises.get('gauss_points_per_element')}",
                f"  Elements: {von_mises.get('element_count')}",
                f"  Count: {von_mises.get('count')}",
                "  Values:",
            ]
        )
        value_summary = von_mises.get("values")
        if isinstance(value_summary, dict):
            lines.extend(format_summary_block(value_summary, indent="    "))
        for label in ("min_location", "max_location"):
            location = von_mises.get(label)
            if isinstance(location, dict):
                lines.append(f"  {label}:")
                for key, value in location.items():
                    lines.append(f"    {key}: {format_float(value)}")
        lines.append("")

    warnings = report.get("warnings", 0)
    if warnings:
        lines.append(f"Warnings: {warnings}")
    return "\n".join(lines).rstrip() + "\n"


def format_csv(report: dict[str, object]) -> str:
    output = io.StringIO()
    fieldnames = [
        "file",
        "result",
        "series",
        "metric",
        "value",
        "analysis",
        "step",
        "location",
        "gauss_point_set",
        "gauss_points_per_element",
        "element_count",
        "extreme_node_id",
        "extreme_element_id",
        "extreme_gauss_point",
    ]
    writer = csv.DictWriter(output, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()

    def write_summary_rows(
        result_name: str,
        result_data: dict[str, object],
        series: str,
        summary: dict[str, object],
    ) -> None:
        base_row = {
            "file": report.get("file"),
            "result": result_name,
            "series": series,
            "analysis": result_data.get("analysis"),
            "step": result_data.get("step"),
            "location": result_data.get("location"),
            "gauss_point_set": result_data.get("gauss_point_set"),
            "gauss_points_per_element": result_data.get("gauss_points_per_element"),
            "element_count": result_data.get("element_count"),
        }
        for metric, value in summary.items():
            writer.writerow({**base_row, "metric": metric, "value": value})

    def write_location_rows(
        result_name: str,
        result_data: dict[str, object],
        series: str,
        location_data: dict[str, object],
    ) -> None:
        base_row = {
            "file": report.get("file"),
            "result": result_name,
            "series": series,
            "analysis": result_data.get("analysis"),
            "step": result_data.get("step"),
            "location": result_data.get("location"),
            "gauss_point_set": result_data.get("gauss_point_set"),
            "gauss_points_per_element": result_data.get("gauss_points_per_element"),
            "element_count": result_data.get("element_count"),
            "extreme_node_id": location_data.get("node_id"),
            "extreme_element_id": location_data.get("element_id"),
            "extreme_gauss_point": location_data.get("gauss_point"),
        }
        for metric, value in location_data.items():
            if metric in {"node_id", "element_id", "gauss_point"}:
                continue
            writer.writerow({**base_row, "metric": metric, "value": value})

    results = report.get("results")
    if not isinstance(results, dict):
        return output.getvalue()

    displacement = results.get("Displacement")
    if isinstance(displacement, dict):
        components = displacement.get("components", {})
        if isinstance(components, dict):
            for axis in ("x", "y", "z"):
                axis_summary = components.get(axis)
                if isinstance(axis_summary, dict):
                    write_summary_rows("Displacement", displacement, f"component_{axis}", axis_summary)
        magnitude = displacement.get("magnitude")
        if isinstance(magnitude, dict):
            write_summary_rows("Displacement", displacement, "magnitude", magnitude)
        for series in ("min_location", "max_location"):
            location = displacement.get(series)
            if isinstance(location, dict):
                write_location_rows("Displacement", displacement, series, location)

    von_mises = results.get("Von Mises")
    if isinstance(von_mises, dict):
        value_summary = von_mises.get("values")
        if isinstance(value_summary, dict):
            write_summary_rows("Von Mises", von_mises, "values", value_summary)
        for series in ("min_location", "max_location"):
            location = von_mises.get(series)
            if isinstance(location, dict):
                write_location_rows("Von Mises", von_mises, series, location)

    return output.getvalue()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        add_help=False,
        description="Extract Displacement and Von Mises statistics from a GiD .post.res file.",
        epilog=(
            "examples:\n"
            "  python TOOLS/extract_sim_data.py octreemesh.post.res\n"
            "  python TOOLS/extract_sim_data.py octreemesh.post.res --format json --output sim_stats.json\n"
            "  python TOOLS/extract_sim_data.py octreemesh.post.res --format csv --output sim_stats.csv\n"
            "  python TOOLS/extract_sim_data.py octreemesh.post.res --percentiles 1,5,50,95,99\n\n"
            "notes:\n"
            "  Displacement rows are treated as node_id ux uy uz. The script reports component\n"
            "  statistics plus vector magnitude statistics. Von Mises values are grouped by\n"
            "  the referenced GaussPoints definition, for example 8 values per hexahedral element."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-h", "--help", action="help", help="show this usage guide and exit")
    parser.add_argument("input_file", type=Path, help="Path to octreemesh.post.res")
    parser.add_argument(
        "--results",
        nargs="+",
        default=list(DEFAULT_RESULTS),
        help="Result names to extract. Default: Displacement 'Von Mises'.",
    )
    parser.add_argument(
        "--format",
        choices=("text", "json", "csv"),
        default="text",
        help="Output format. Default: text.",
    )
    parser.add_argument("--output", type=Path, help="Write output to this file instead of stdout.")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Treat malformed rows as fatal parse errors instead of warnings.",
    )
    parser.add_argument(
        "--percentiles",
        type=parse_percentiles,
        default=DEFAULT_PERCENTILES,
        help="Comma-separated percentile list. Default: 5,25,50,75,95.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.input_file.exists():
        parser.error(f"input file does not exist: {args.input_file}")
    if not args.input_file.is_file():
        parser.error(f"input path is not a file: {args.input_file}")

    try:
        report = parse_result_file(
            args.input_file,
            requested_results=set(args.results),
            percentiles=args.percentiles,
            strict=args.strict,
        )
    except ParseError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    if args.format == "json":
        output = json.dumps(report, indent=2, sort_keys=True) + "\n"
    elif args.format == "csv":
        output = format_csv(report)
    else:
        output = format_text(report)

    if args.output:
        args.output.write_text(output, encoding="utf-8")
    else:
        print(output, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
