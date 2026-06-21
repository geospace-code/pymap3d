#!/usr/bin/env python3
"""
pymap3d CLI

Run with: python -m pymap3d <command> [options]

Example:
    python -m pymap3d geodetic2ecef 40.7128 -74.0060 10
"""

import argparse
from datetime import datetime
import inspect
import json
from typing import Any

import pymap3d as pm


def _is_datetime_param(param: inspect.Parameter) -> bool:
    """Detect datetime-like parameters from annotation or conventional parameter names."""
    ann = param.annotation
    name = param.name.lower()
    return ann is datetime or name in {"t", "time"}


def _normalize_dest(name: str) -> str:
    return name.replace("-", "_")


def _coerce_value(raw: Any, param: inspect.Parameter) -> Any:
    """Convert string CLI values into practical Python values for pymap3d functions."""
    if raw is None:
        return None

    if isinstance(raw, bool):
        return raw

    if _is_datetime_param(param):
        return pm.str2dt(raw)

    if isinstance(raw, str):
        text = raw.strip()

        if "," in text:
            try:
                return [float(x.strip()) for x in text.split(",")]
            except ValueError:
                return text

        try:
            return float(text)
        except ValueError:
            return text

    return raw


def _add_parameter_arg(parser: argparse.ArgumentParser, param: inspect.Parameter) -> None:
    """Create argparse options for one function parameter."""
    name = param.name
    dest = _normalize_dest(name)

    if name == "ell":
        return

    default = param.default
    has_default = default is not inspect.Parameter.empty

    if has_default and isinstance(default, bool):
        flag = name.replace("_", "-")
        if default:
            parser.add_argument(
                f"--{flag}",
                dest=dest,
                action="store_true",
                default=True,
                help=f"Enable {name} (default: enabled)",
            )
            parser.add_argument(
                f"--no-{flag}",
                dest=dest,
                action="store_false",
                help=f"Disable {name}",
            )
        else:
            parser.add_argument(
                f"--{flag}",
                dest=dest,
                action="store_true",
                default=False,
                help=f"Enable {name}",
            )
            parser.add_argument(
                f"--no-{flag}",
                dest=dest,
                action="store_false",
                help=f"Disable {name} (default: disabled)",
            )
        return

    if has_default:
        parser.add_argument(
            f"--{name.replace('_', '-')}",
            dest=dest,
            default=str(default),
            help=f"{name} (default: {default})",
        )
    else:
        parser.add_argument(name, help=name)


def _register_commands(subparsers: argparse._SubParsersAction) -> dict[str, Any]:
    """Build command parsers from practical callable exports in pymap3d.__all__."""
    excluded = {"Ellipsoid"}
    commands: dict[str, Any] = {}

    for name in sorted(pm.__all__):
        if name in excluded:
            continue

        obj = getattr(pm, name, None)
        if obj is None or not callable(obj) or inspect.isclass(obj):
            continue

        sig = inspect.signature(obj)

        # Skip functions that need variadic arguments; these are not practical via CLI.
        if any(
            p.kind in (inspect.Parameter.VAR_POSITIONAL, inspect.Parameter.VAR_KEYWORD)
            for p in sig.parameters.values()
        ):
            continue

        cmd_parser = subparsers.add_parser(name, help=(obj.__doc__ or "").strip().splitlines()[0])
        cmd_parser.set_defaults(_func_name=name)

        for param in sig.parameters.values():
            _add_parameter_arg(cmd_parser, param)

        commands[name] = obj

    return commands


def _call_selected_function(
    args: argparse.Namespace, commands: dict[str, Any], ell: pm.Ellipsoid
) -> Any:
    """Build call arguments from argparse namespace and invoke the selected pymap3d function."""
    func_name = args.command
    func = commands[func_name]
    sig = inspect.signature(func)

    kwargs: dict[str, Any] = {}

    for param in sig.parameters.values():
        name = param.name
        if name == "ell":
            kwargs[name] = ell
            continue

        dest = _normalize_dest(name)
        raw = getattr(args, dest)
        kwargs[name] = _coerce_value(raw, param)

    return func(**kwargs)


def _format_scalar(value: Any, precision: int | None) -> str:
    if precision is not None and isinstance(value, float):
        return f"{value:.{precision}f}"

    return str(value)


def _to_jsonable(value: Any) -> Any:
    """Convert results into JSON-serializable structures."""
    match value:
        case tuple():
            return [_to_jsonable(v) for v in value]
        case list():
            return [_to_jsonable(v) for v in value]
        case datetime():
            return value.isoformat()

    return value


def _emit_result(result: Any, output_mode: str, precision: int | None) -> None:
    """Render results according to output mode."""
    if output_mode == "json":
        print(json.dumps(_to_jsonable(result)))
        return

    match result:
        case tuple():
            print(" ".join(_format_scalar(x, precision) for x in result))
            return
        case list():
            print(" ".join(_format_scalar(x, precision) for x in result))
            return

    print(_format_scalar(result, precision))


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="pymap3d",
        description="pymap3d CLI - Geographic coordinate conversions",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--ellipsoid",
        default="wgs84",
        choices=[name.lower() for name in pm.Ellipsoid.models],
        help="Ellipsoid model to use",
    )
    parser.add_argument(
        "--output",
        default="text",
        choices=["text", "json"],
        help="Output format",
    )
    parser.add_argument(
        "--precision",
        type=int,
        default=None,
        help="Decimal places for floating-point text output",
    )

    subparsers = parser.add_subparsers(dest="command", required=True, help="Conversion command")

    commands = _register_commands(subparsers)

    args = parser.parse_args()

    ell = pm.Ellipsoid.from_name(args.ellipsoid)

    if args.command not in commands:
        parser.error(f"Command {args.command} unknown or not yet implemented")

    result = _call_selected_function(args, commands, ell)
    _emit_result(result, args.output, args.precision)


if __name__ == "__main__":
    main()
