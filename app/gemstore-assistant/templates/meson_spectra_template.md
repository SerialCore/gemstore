# Gemstore Meson Spectra Input Template

This template defines the input format for the gemstore meson spectra computation.

---

## Charmonium Spectra

```
&GLOBAL
  project = charmonium
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
  params = GIScreen_meson
  # params = GIQuadra_meson
  # params = GIScreen_bbbar
  # params = GIQuadra_bbbar
  # params = GIScreen_ccbar
  # params = GIQuadra_ccbar
  # mn = 0.4433275191676
  # ms = 0.606437654085
  # mc = 1.797694082834
  # mb = 5.142643566529
  # b1 = 0.2530091124005
  # b2 = 0.02
  # mu = 0.1401631404922
  # c = -0.6300631479611
  # sigma_0 = 1.776663048762
  # s = 1.180712032887
  # epsilon_cont = -0.3068965865533
  # epsilon_sov = -0.3715004921572
  # epsilon_sos = 0.9276232568028
  # epsilon_tens = -0.5056764482414
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
  rmin = 0.01
&END
```

## Parameter Definition

--f1, f2, quark flavor, 1 for n, 2 for s, 3 for c, 4 for b.
--S for spin, L for orbital angular momentum, J for total angular momentum.
--jl=s1+L for Jj coupling
--nmax, rmax, rmin will be usually fixed
--GIScreen_meson is the default parameter set
