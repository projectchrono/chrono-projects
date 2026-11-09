DEM validation tests
====================

Programs that replicate laboratory tests on granular materials with Chrono::Multicore, to validate Chrono's discrete
element models against measurements, and Python scripts that make figures, tables and clips from their results.

- [Building](#building)
- [Element tests](#element-tests): direct shear in a square box and in a Jenike cell, triaxial test, steady shear in an
  annular cell
- [Runs, states and outputs](#runs-states-and-outputs)
- [Rendering](#rendering)
- [Post-processing scripts](#post-processing-scripts)
- [Other programs](#other-programs)
- [References](#references)

## Building

The programs are built with chrono-projects when `ENABLE_PROJECTS` and `ENABLE_PROJECTS_DEM_VALIDATION` are on (the
second is on by default) and Chrono was built with the Multicore module. Chrono::VSG is optional: it provides run-time
visualization (`--vsg`) and the rendering of saved states (`--render`).

## Element tests

Three programs replicate element tests of soil mechanics on spheres, and a fourth measures the μ(I) rheology in steady
shear. They use soft-sphere (SMC) contact with Hertz normal forces. The tangential displacement model (`multistep`, the
default, `onestep` or `none`) is set with `--history`; the direct shear and triaxial programs can set a different one
for the test (`--history_shear`, `--history_test`). The programs take their parameters from the command line, in SI
units; `--help` lists every option with its default.

### Direct shear in a square box: `test_VDEM_directShear_squareBox`

The lower box moves at constant speed under the fixed upper box, and a lid that can only move vertically applies the
normal stress with its weight. Two samples are predefined (`--sample`), and every parameter of a sample can be changed
on the command line:

- `hartl` (default): a model of the direct shear tests of Härtl and Ooi (2008), 5000 glass beads of 6 mm in a 12 cm
  square box;
- `history`: the sample of the contact-history study of Fleischmann et al. (2016), 1800 spheres of 5 mm in a 6 cm
  square box, compacted with a low friction (`--dense`).

The spheres fall into the box from a jittered grid, the lid compacts them, and the lower box then shears the sample.
The time series gives the shear stress from the force on the lower box (`tau_bin`) and from the force on the upper box
and the lid (`tau_upper`).

```
test_VDEM_directShear_squareBox --sigma 3100 --out hartl_3.1kPa

# The history study: one sample, sheared with each tangential displacement model
# (a tiny travel ends the first run right after the sample is saved)
test_VDEM_directShear_squareBox --sample history --save_state history.bin --travel 1e-6 --out history_prep
test_VDEM_directShear_squareBox --sample history --load_state history.bin --history_shear multistep --out history_multistep
test_VDEM_directShear_squareBox --sample history --load_state history.bin --history_shear onestep --out history_onestep
test_VDEM_directShear_squareBox --sample history --load_state history.bin --history_shear none --out history_none
```

### Direct shear in a Jenike cell: `test_VDEM_directShear_Jenike`

A replication of the direct shear tests of Härtl and Ooi (2008) on 6 mm glass beads, in a Jenike cell of 143 mm inner
diameter: a bottom ring closed by a base, a top ring resting on it, and a lid on the sample inside the top ring. The
beads rain into the cell in layers; those above the rim are then retired, and the lid is placed on the sample, so that
its weight and the weight of the beads above the shear plane give the normal stress (`--sigma`). The top ring and the
lid then move at constant speed. As in Härtl and Ooi's simulations, the top ring starts 3 mm off center, so the shear
displacement starts at -3 mm. The shear stress (`tau`) is the force of the beads on the top ring and the lid over the
cell area.

Other options set the grooved lid and base of the experiment (`--grooves 1`), rolling friction (`--mu_roll`), a lid
free to tilt (`--lid_tilt`), the filling (`--rain_*`, and the bead friction while filling, `--mu_fill`) and gravity
(`--gravity`). One filling serves several normal stresses:

```
test_VDEM_directShear_Jenike --seed 1 --fill_only --save_fill fill_seed1.bin --out fill_seed1
test_VDEM_directShear_Jenike --seed 1 --load_fill fill_seed1.bin --sigma 3100 --out sigma3100_seed1
test_VDEM_directShear_Jenike --seed 1 --load_fill fill_seed1.bin --sigma 24200 --out sigma24200_seed1
```

### Triaxial test: `test_VDEM_triaxial_membrane`

A replication of the drained triaxial tests of Cui, O'Sullivan and O'Neill (2007) on chrome steel spheres: a
cylindrical specimen 101 mm across and about 203 mm high, between rigid platens, under a cell pressure of 80 kPa
(`--sigma3`). The specimen has 15,400 spheres of 5 mm by default; `--radii 2,2.5,3 --num_spheres 14350` gives the
experiment's mixture of 4, 5 and 6 mm spheres in equal numbers. The density, Poisson ratio and friction coefficients
are the experiment's; Young's modulus is 0.2 GPa by default, a thousandth of steel's.

The lateral boundary is a stress-controlled membrane: forces toward the axis on the spheres at the surface of the
specimen, which carry the cell pressure. With Cui et al.'s rule (`--membrane voronoi`, the default), each sphere near
the surface that touches no sphere farther out carries the pressure on its Voronoi cell of the unrolled membrane; with
`--membrane bins`, the outermost sphere of each bin of angle and height carries the pressure on the bin. The platens
are wider than the specimen, so that no sphere can leave past a platen edge as the specimen bulges.

The spheres settle in a rigid mold. The top platen is then placed, the membrane and a stress servo on the top platen
consolidate the specimen, and the mold is lifted clear. The test is either

- triaxial (`--test triaxial`, the default): the top platen moves down at a constant axial strain rate
  (`--strain_rate`, 5 %/s by default) to `--strain_end`, with the membrane holding the cell pressure; or
- isotropic (`--test isotropic`): the platen servo and the membrane follow a list of pressures (`--path`), with a hold
  at each, for a compression with unload-reload loops.

Stresses and strains are positive in compression. The axial stress is computed on the current cross-section (the
volume enclosed by the membrane over the height) and, as in the laboratory, on the cross-section at the start of the
test (`sigma1_top_a0`).

```
# Spheres deposited without friction; the friction is --mu from the placement of the top platen on
test_VDEM_triaxial_membrane --mu_fill 0 --save_fill fill_mono.bin --out mono_s80
test_VDEM_triaxial_membrane --load_fill fill_mono.bin --mu_fill 0 --test isotropic --ramp 4e5 --out iso_mono
test_VDEM_triaxial_membrane --radii 2,2.5,3 --num_spheres 14350 --mu_fill 0 --out poly_s80
```

### Steady shear in an annular cell: `test_VDEM_annularShear`

μ(I), the friction coefficient of steady dense flow against the inertial number I, on the glass beads of the direct
shear programs. Chrono::Multicore has no periodic boundaries, so an annulus closes the flow direction on itself: 14,400
beads fill the ring between a frictionless inner cylinder and a frictionless outer ring (radii 0.09 and 0.15 m) and are
sheared between two rough plates, flat disks covered with a random layer of glued spheres of the bead size. The plates
turn about the axis in opposite directions, each at half the relative angular speed. A stress servo on the height of
the top plate holds the normal pressure (`--pressure`): the plate moves at `--servo_gain` times the pressure error, at
most `--servo_vmax`.

The beads settle under gravity. At `--t_place`, gravity is switched off and the top plate is placed above the beads; it
moves down until it carries a tenth of the pressure, and its servo then holds the pressure. At `--t_shear`, the plates
start to turn. The relative angular speed follows from a nominal inertial number at the mid radius (`--I_nominal`):
Ω = I √(p/ρ) H / (d r_m), with H the distance between the centers of the two glued layers at the start of the shear,
d the bead diameter and r_m the mid radius. It ramps up over a nominal strain of `--strain_ramp`, and the run ends at
the nominal strain `--strain_end`. The time series gives the pressure and the shear stress on both plates and their
ratio (`mu_top`, `mu_bottom`); the local shear rate, stresses and solid fraction come from the particle and contact
snapshots, by coarse-graining.

One filling, and one consolidated sample at each pressure, serve every inertial number:

```
test_VDEM_annularShear --seed 1 --fill_only --save_fill fill_seed1.bin --out fill_seed1
test_VDEM_annularShear --load_fill fill_seed1.bin --pressure 3100 --prepare_only --save_state cons_p3100.bin --out cons_p3100
test_VDEM_annularShear --load_state cons_p3100.bin --pressure 3100 --I_nominal 0.01 --strain_ramp 0.5 --servo_gain 1e-5 --out p3100_I0.01
```

## Runs, states and outputs

Each run of an element test writes to its output directory (`--out`):

| File | Contents |
|---|---|
| `meta.json` | all parameters, run information and a summary |
| `timeseries.csv` | one row per output interval (`--out_dt`; in the annular cell, every `--out_strain` of nominal strain during the shear), written as the run goes; its comment lines define the columns |
| `particles/*.csv`, `contacts/*.csv` | particle and contact snapshots at chosen shear displacements or axial strains (`--snaps`), and frames at a fixed interval (`--frame_disp`, `--frame_strain`); in the annular cell, every `--snap_every` of nominal strain from `--snap_from` |
| `states/*.bin` | the full state at each snapshot (`--snap_states`) or frame (`--frame_states`) |

A run can start from a saved state: a filled sample (`--save_fill` and `--load_fill`, in the Jenike cell, the triaxial
test and the annular cell), the sample at the start of the test (`--save_state` and `--load_state`), or a snapshot or
frame state (`--load_state`). State files are read and written at the paths given, relative to the working directory
and not to `--out`. A run that loads a state must use the sample options of the run that saved it.

`scripts/run_case.py` launches a run in its own folder. It copies the executable there (so that later rebuilds do not
affect the run), records in `meta.json` the commit and a fingerprint of any uncommitted changes of the Chrono and
chrono-projects source trees (through `--provenance`), and runs the program with the run folder as working directory,
its console output in `stdout.txt` and its exit code in `exit.txt`. On Windows, the run has below-normal priority.
`--path` puts a folder first on the PATH of the run (on Windows, the folder of the Chrono DLLs).

```
python run_case.py --exe <build>/bin/Release/test_VDEM_directShear_Jenike.exe --run sigma3100_seed1 --chrono <Chrono source> --projects <chrono-projects source> -- --load_fill ../fill_seed1.bin --sigma 3100 --out .
```

## Rendering

With Chrono::VSG, the square-box and triaxial programs render saved states to images and exit, without simulating
(`--render`, with a state file, or a directory whose state files are all rendered). The sample options must match those
of the run that saved the states. The camera is fixed, and the sample is cut at its vertical mid-plane. Each state
gives one image per layer:

- `grains/`: the spheres behind the mid-plane, light gray or in gray stripes (`--render_stripes`) that show how the
  sample deforms;
- `network/`: the contacts in a slab behind the mid-plane (`--render_slab`) whose normal force is above a fraction
  (`--render_fmin`) of a reference force (`--render_fref`, by default the mean normal force of the first state), drawn
  as cylinders whose radius and color (light to dark orange) are proportional to the normal force, up to 3 times the
  reference;
- `ratio/`: the same network colored by the Coulomb ratio |f_t| / (mu f_n) (with `--render_ratio`);
- `walls/`: the walls or platens alone, to separate the network from them when compositing.

`scripts/ds_render.py` composites the network over the grains, frame by frame, into a loopable MP4 clip with a poster
frame:

```
# frames/: the frame states of a triaxial run (states/frame_*.bin)
test_VDEM_triaxial_membrane --render frames --render_out images --render_stripes 12
python ds_render.py --render images --out forcenet
```

## Post-processing scripts

The scripts in `scripts/` read run folders and write figures (PNG at 300 dpi and SVG), Markdown tables and clips. They
need Python 3 with NumPy and Matplotlib; `ds_networks.py` and `ds_render.py` also need Pillow, and `ds_render.py` needs
FFmpeg on the PATH. The figures use the Calibri font where it is installed. `ds_curves.py`, `ds_plateau.py` and
`tri_curves.py` draw runs that differ only in their seed as their mean, with a band over their range. The usage comment
at the top of each script lists its options.

| Script | Output |
|---|---|
| `ds_curves.py` | shear stress against shear displacement at several normal stresses, on the axes of Härtl and Ooi (2008), Fig. 2; table of peak and residual stress ratios |
| `ds_plateau.py` | plateau stress ratio against normal stress for groups of runs, next to Härtl and Ooi's Table 2 |
| `ds_porosity.py` | plateau stress ratio against the porosity before shearing, for samples filled in different ways, next to Härtl and Ooi's three filling methods |
| `ds_gravity.py` | stress ratio against shear displacement under different gravities; table of plateau ratios, porosity, lid rise and the split of the normal load between the lid and the beads |
| `ds_history.py` | stress ratio against shear displacement for one sample sheared with each tangential displacement model; table of ratios, dilation, sliding fraction and coordination |
| `ds_networks.py` | force-network images of several runs side by side |
| `ds_render.py` | MP4 clip, poster and stills from the image layers of `--render` |
| `tri_curves.py` | stress ratio, dilation, q/p or void ratio against axial strain, next to the laboratory curves; table of void ratios, peak and end values |
| `tri_csl.py` | stress paths and end states in the (p, q) plane, the isotropic compression in the (ln p, v) plane, and the critical state line in both; table of the fitted parameters |
| `tri_love_weber.py` | table of the Love-Weber stress of each snapshot against the platen and membrane stresses |
| `dsplot.py` | shared helpers: reading run folders, figure style |
| `run_case.py` | launches a run in its own folder (see above) |

`scripts/data/cui2007_fig4_lab.csv` holds the laboratory curves of Cui et al. (2007), Fig. 4: the stress ratio against
axial strain at 80 kPa, without area correction, read from the vector paths of the published figure.
`tri_curves.py --lab` draws them.

In these examples, `<scripts>` is the `scripts/` folder, and the commands run in the folder that holds the run folders:

```
python <scripts>/ds_curves.py --out ds_curves --tau tau sigma3100_seed1 sigma6400_seed1 sigma12500_seed1 sigma24200_seed1
python <scripts>/tri_curves.py --out tri_curves --by_specimen --area0 --lab <scripts>/data/cui2007_fig4_lab.csv --xmax 12.5 mono_s80 poly_s80
python <scripts>/tri_csl.py --pq tri_csl_pq --lnp tri_ncl_lnp --iso iso_mono mono_s20 mono_s40 mono_s80 mono_s160 mono_s320
```

## Other programs

These programs set their parameters in constants at the top of the source. `USE_SMC` selects SMC contact (defined) or
NSC contact (commented out), and `problem` selects the stage to run. The first stage (`SETTLING`) settles the granular
material and saves a checkpoint; each later stage starts from the checkpoint of the stage before it. In the last three
programs of the table, `TESTING` runs the mechanism with a single test ball. Results go to `../<NAME>_SMC` or
`../<NAME>_NSC`, relative to the working directory, with POV-Ray data in `POVRAY/`.

| Program | Test | Units |
|---|---|---|
| `test_VDEM_demo_crater` | a 25.4 mm sphere dropped onto a bed of 1 mm beads (low-velocity cratering) | SI |
| `test_VDEM_demo_massflow` | 0.5 mm beads flowing out of a container through a gap that an angled insert opens, and the mass collected below | SI |
| `test_VDEM_demo_penetrometer` | a sphere or a cone dropped into a bed of 3 mm beads | SI |
| `test_VDEM_directShear_actuator` | direct shear of 6 mm beads, the shear box driven at constant speed by a linear actuator | cm, g, s |
| `test_VDEM_pressureSinkage` | a 16 by 7 cm plate pushed at constant speed into a bin of 8 mm beads | cm, g, s |
| `test_VDEM_singleWheel` | a wheel driven through a bin of granular material at a set slip | cm, g, s |

## References

- L. Cui, C. O'Sullivan, S. O'Neill, "An analysis of the triaxial apparatus using a mixed boundary three-dimensional
  discrete element model", *Géotechnique* 57(10), 831-844, 2007. https://doi.org/10.1680/geot.2007.57.10.831
- J. Fleischmann, R. Serban, D. Negrut, P. Jayakumar, "On the importance of displacement history in soft-body contact
  models", *ASME Journal of Computational and Nonlinear Dynamics* 11(4), 044502, 2016.
- J. Härtl, J. Y. Ooi, "Experiments and simulations of direct shear tests: porosity, contact friction and bulk
  friction", *Granular Matter* 10(4), 263-271, 2008. https://doi.org/10.1007/s10035-008-0085-3
