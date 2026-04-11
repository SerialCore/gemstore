# Gemstore Meson Spectra Input Template

This template defines the input format for the gemstore meson spectra computation.

---

## Charmonium Spectra

```
&GLOBAL
  project = {file name}
  task = SPECTRA
  # task = RADIUS
  # task = DECAY3P0
  # task = COUPLCHN
  # task = SCATTER
&END
&SYSTEM
  # model = GI_STRING
  model = GI_SCREEN
  # model = GI_QUADRA
  system = MESON
  # system = BARYON
  # system = MOLECULE
&END
&PARAMS
  # params = GIString_meson
  # params = GIScreen_meson
  # params = GIQuadra_meson
  params = GIScreen_ccbar
  # params = GIScreen_bbbar
  # params = GIQuadra_ccbar
  # params = GIQuadra_bbbar
  # mn = 0.4713455847642
  # ms = 0.6283121820133
  # mc = 1.810505119204
  # mb = 5.156014766761
  # b1 = 0.2575467075473
  # mu = 0.1453562021339
  # c = -0.658943240626
  # sigma_0 = 1.884145499156
  # s = 1.113514380624
  # epsilon_cont = -0.32452949845
  # epsilon_sov = -0.5404734834836
  # epsilon_sos = 0.9999999508829
  # epsilon_tens = -0.4999502878773
&END
&QUANTUM
  f1 = 3
  f2 = 3
  # f3 = 1
  # f4 = 1
  S = 1
  L = 0
  # jl = 0.5
  J = 1
&END
&GAUSS
  nmax = 20
  rmax = 20.0
  rmin = 0.1
&END
```

## Parameter Definition

--f1, f2, quark flavor, 1 for n, 2 for s, 3 for c, 4 for b.
--S for spin, L for orbital angular momentum, J for total angular momentum.
--jl=s1+L for Jj coupling
--nmax, rmax, rmin will be usually fixed
