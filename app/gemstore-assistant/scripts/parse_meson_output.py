import argparse
import json
import math
from pathlib import Path


def load_output(path):
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)

    if not isinstance(data, dict):
        raise ValueError("Top-level output must be a JSON object")

    states = data.get("states")
    if not isinstance(states, list):
        raise ValueError("Output JSON must contain a 'states' array")

    parsed_states = []
    for item in states:
        if not isinstance(item, dict):
            raise ValueError("Each state entry must be a JSON object")

        index = item.get("index")
        mass = item.get("mass")
        rms_radius = item.get("rms_radius")
        eigenvector = item.get("eigenvector")

        if not isinstance(index, int):
            raise ValueError("State 'index' must be an integer")
        if not isinstance(mass, (int, float)):
            raise ValueError("State 'mass' must be numeric")
        if not isinstance(rms_radius, (int, float)):
            raise ValueError("State 'rms_radius' must be numeric")
        if not isinstance(eigenvector, list):
            raise ValueError("State 'eigenvector' must be an array")

        parsed_states.append(
            {
                "index": index,
                "mass": float(mass),
                "rms_radius": float(rms_radius),
                "eigenvector": [float(value) for value in eigenvector],
            }
        )

    return {
        "generated": data.get("generated"),
        "project": data.get("project"),
        "task": data.get("task"),
        "model": data.get("model"),
        "system": data.get("system"),
        "basis": data.get("basis"),
        "states": parsed_states,
    }


def summarize_output(data, limit=None):
    states = data["states"]
    if limit is not None:
        states = states[:limit]

    lines = []
    lines.append(f"project: {data.get('project')}")
    lines.append(f"generated: {data.get('generated')}")
    lines.append(f"states: {len(data['states'])}")

    for state in states:
        eigenvector = state["eigenvector"]
        max_coeff = max((abs(value) for value in eigenvector), default=0.0)
        norm_sq = sum(value * value for value in eigenvector)
        norm = math.sqrt(norm_sq)
        lines.append(
            "state {index}: mass={mass:.6f} GeV, rms={rms:.6f} fm, "
            "max|c|={maxc:.6f}, ||c||={norm:.6f}".format(
                index=state["index"],
                mass=state["mass"],
                rms=state["rms_radius"],
                maxc=max_coeff,
                norm=norm,
            )
        )

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Parse gemstore meson JSON output")
    parser.add_argument("output_file", help="Path to gemstore JSON output file")
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Only summarize the first N states",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Print normalized parsed JSON instead of a text summary",
    )
    args = parser.parse_args()

    data = load_output(Path(args.output_file))

    if args.json:
        print(json.dumps(data, indent=2))
    else:
        print(summarize_output(data, limit=args.limit))


if __name__ == "__main__":
    main()
