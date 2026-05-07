import json


template = {
    "project": "name",
    "task": "SPECTRA",
    "system": {
        "type": "MESON",
        "f1": 3,
        "f2": 3,
        "S": 1,
        "L": 0,
        "J": 1,
    },
    "model": {
        "type": "GISCREEN",
        "param": "GISCREEN_MESON",
    },
    "basis": {
        "type": "GEM",
        "nmax": 16,
        "rmax": 30.0,
        "rmin": 0.1,
    },
    "print": {
        "pot": "false",
        "wfn": "false"
    }
}

states = [
    ("1S0", 0, 0, 0),
    ("3S1", 1, 0, 1),
    ("1P1", 0, 1, 1),
    ("3P0", 1, 1, 0),
    ("3P1", 1, 1, 1),
    ("3P2", 1, 1, 2),
    ("1D2", 0, 2, 2),
    ("3D1", 1, 2, 1),
    ("3D2", 1, 2, 2),
    ("3D3", 1, 2, 3),
]


def write_input_file(data, filename):
    with open(filename, "w", encoding="ascii") as f:
        json.dump(data, f, indent=2)
        f.write("\n")


def set_custom_param_file(data, model_type, file_path):
    data["model"]["type"] = model_type
    data["model"]["param"] = f"{model_type}_CUSTOM"
    data["model"]["file"] = file_path


def generate_family_inputs(prefix, f1, f2, param_name):
    for name, S, L, J in states:
        content = json.loads(json.dumps(template))
        content["project"] = f"{prefix}_{name}"
        content["system"]["f1"] = f1
        content["system"]["f2"] = f2
        content["system"]["S"] = S
        content["system"]["L"] = L
        content["system"]["J"] = J
        content["model"]["param"] = param_name

        fname = f"{prefix}_{name}.json"
        write_input_file(content, fname)
        print(f"Created {fname}")


if __name__ == "__main__":
    generate_family_inputs("charmonium", 3, 3, "GISCREEN_CCBAR")
    # generate_family_inputs("bottomonium", 4, 4, "GISCREEN_BBBAR")
    # generate_family_inputs("Bc", 3, 4, "GISCREEN_MESON")
    # generate_family_inputs("Bs", 2, 4, "GISCREEN_MESON")
    # generate_family_inputs("Ds", 2, 3, "GISCREEN_MESON")
    # generate_family_inputs("B", 1, 4, "GISCREEN_MESON")
    # generate_family_inputs("D", 1, 3, "GISCREEN_MESON")
    print("All input files generated successfully.")
