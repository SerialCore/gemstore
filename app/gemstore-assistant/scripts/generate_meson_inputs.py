import os

template = """&GLOBAL
  project = name
  task = SPECTRA
&END
&SYSTEM
  model = GI_SCREEN
  system = MESON
&END
&PARAMS
  params = GIScreen_meson
&END
&QUANTUM
  f1 = 3
  f2 = 3
  S = 1
  L = 0
  J = 1
&END
&GAUSS
  nmax = 20
  rmax = 20.0
  rmin = 0.1
&END
"""

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


def generate_ccbar_inputs():
    # For charmonium f1=f2=3
    content = template.replace("f1 = 3", f"f1 = 3")
    content = content.replace("f2 = 3", f"f2 = 3")
    content = content.replace("params = GIScreen_meson", f"params = GIScreen_ccbar")
    for name, S, L, J in states:
        content = content.replace("project = name", f"project = charmonium_{name}")
        content = content.replace("S = 1", f"S = {S}")
        content = content.replace("L = 0", f"L = {L}")
        content = content.replace("J = 1", f"J = {J}")
        fname = f"charmonium_{name}.inp"
        with open(fname, "w") as f:
            f.write(content)
        print(f"Created {fname}")

def generate_bbbar_inputs():
    # For bottomonium
    content = template.replace("f1 = 3", f"f1 = 4")
    content = content.replace("f2 = 3", f"f2 = 4")
    content = content.replace("params = GIScreen_meson", f"params = GIScreen_bbbar")
    for name, S, L, J in states:
        content = content.replace("project = name", f"project = bottomonium_{name}")
        content = content.replace("S = 1", f"S = {S}")
        content = content.replace("L = 0", f"L = {L}")
        content = content.replace("J = 1", f"J = {J}")
        fname = f"bottomonium_{name}.inp"
        with open(fname, "w") as f:
            f.write(content)
        print(f"Created {fname}")

def generate_cbbar_inputs():
    # For Bc
    content = template.replace("f1 = 3", f"f1 = 3")
    content = content.replace("f2 = 3", f"f2 = 4")
    for name, S, L, J in states:
        content = content.replace("project = name", f"project = Bc_{name}")
        content = content.replace("S = 1", f"S = {S}")
        content = content.replace("L = 0", f"L = {L}")
        content = content.replace("J = 1", f"J = {J}")
        fname = f"Bc_{name}.inp"
        with open(fname, "w") as f:
            f.write(content)
        print(f"Created {fname}")

def generate_sbbar_inputs():
    # For Bs
    content = template.replace("f1 = 3", f"f1 = 2")
    content = content.replace("f2 = 3", f"f2 = 4")
    for name, S, L, J in states:
        content = content.replace("project = name", f"project = Bs_{name}")
        content = content.replace("S = 1", f"S = {S}")
        content = content.replace("L = 0", f"L = {L}")
        content = content.replace("J = 1", f"J = {J}")
        fname = f"Bs_{name}.inp"
        with open(fname, "w") as f:
            f.write(content)
        print(f"Created {fname}")

def generate_scbar_inputs():
    # For Ds
    content = template.replace("f1 = 3", f"f1 = 2")
    content = content.replace("f2 = 3", f"f2 = 3")
    for name, S, L, J in states:
        content = content.replace("project = name", f"project = Ds_{name}")
        content = content.replace("S = 1", f"S = {S}")
        content = content.replace("L = 0", f"L = {L}")
        content = content.replace("J = 1", f"J = {J}")
        fname = f"Ds_{name}.inp"
        with open(fname, "w") as f:
            f.write(content)
        print(f"Created {fname}")

def generate_nbbar_inputs():
    # For B
    content = template.replace("f1 = 3", f"f1 = 1")
    content = content.replace("f2 = 3", f"f2 = 4")
    for name, S, L, J in states:
        content = content.replace("project = name", f"project = B_{name}")
        content = content.replace("S = 1", f"S = {S}")
        content = content.replace("L = 0", f"L = {L}")
        content = content.replace("J = 1", f"J = {J}")
        fname = f"B_{name}.inp"
        with open(fname, "w") as f:
            f.write(content)
        print(f"Created {fname}")

def generate_ncbar_inputs():
    # For D
    content = template.replace("f1 = 3", f"f1 = 1")
    content = content.replace("f2 = 3", f"f2 = 3")
    for name, S, L, J in states:
        content = content.replace("project = name", f"project = D_{name}")
        content = content.replace("S = 1", f"S = {S}")
        content = content.replace("L = 0", f"L = {L}")
        content = content.replace("J = 1", f"J = {J}")
        fname = f"D_{name}.inp"
        with open(fname, "w") as f:
            f.write(content)
        print(f"Created {fname}")

if __name__ == "__main__":
    generate_ccbar_inputs()
    #generate_bbbar_inputs()
    #generate_cbbar_inputs()
    #generate_sbbar_inputs()
    #generate_scbar_inputs()
    #generate_nbbar_inputs()
    #generate_ncbar_inputs()
    print("All input files generated successfully.")
