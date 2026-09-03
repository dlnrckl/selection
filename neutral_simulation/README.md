# Neutral Simulation

This directory contains the SLiM model and helper files used to generate neutral 50-kb genomic simulations under the demographic model implemented in `neutral.slim`.

## Files

### `neutral.slim`

Main SLiM simulation script.

The model simulates a **50-kb genomic segment** with:

- mutation rate: `1.2 × 10^-8`
- recombination rate: `1 × 10^-8`
- neutral mutation types with selection coefficient `s = 0`
- a focal neutral mutation (`m2`) introduced at the center of the simulated region

The demographic history implemented in the script includes West Asian, Asian, Caucasus hunter-gatherer, Anatolian Farmer, Bronze Age Anatolian, and present-day Anatolian populations.

### `simulation_freqs.csv`

Input table used by `slimModifier.py`.

The file contains two columns:

```text
frequency,id
```

The first column provides the values used by the modifier script, while the second column is used to name the generated SLiM files.

### `slimModifier.py`

Python helper script intended to generate multiple SLiM scripts from a template using the values in `simulation_freqs.csv`.

It creates a directory named:

```text
slim_scripts/
```

and writes one `.slim` file for each row of the CSV file.

## Demographic model

The population history implemented in `neutral.slim` is summarized below.

| Generation | Event |
|---|---|
| 1 | `p1` (West Asian population) is initialized with a population size of 7,500. |
| 75,000 | `p1` expands to 15,000. A focal neutral mutation is introduced at the center of the 50-kb region. |
| 75,001 | `p2` (Asian population) splits from `p1` with a population size of 10,000. |
| 75,700 | `p3` (Caucasus hunter-gatherer) and `p4` (Anatolian Farmer) split from `p1`, with population sizes of 3,000 and 5,000, respectively. |
| 75,700 | The focal mutation frequency is evaluated in `p4`. Simulations outside the specified `0.45–0.55` range are restarted from the saved population state. |
| 76,330 | `p5` (Bronze Age Anatolians) is created with a population size of 20,000 and ancestry contributions of 25% from `p3` and 75% from `p4`. |
| 76,330 | Migration is stopped and the source populations `p3` and `p4` are removed. |
| 76,500 | `p6` (present-day Anatolians) is created with a population size of 30,000 and ancestry contributions of 10% from `p2` and 90% from `p5`. |
| 76,500 | Migration is stopped and the source populations `p2` and `p5` are removed. |
| 76,501 | Simulation terminates. |

## MS-format outputs

The simulation produces MS-format samples at three time points:

```text
p_1_1000.ms   # Neolithic / Anatolian Farmer population (p4)
p_1_2000.ms   # Bronze Age Anatolian population (p5)
p_1_3000.ms   # Present-day Anatolian population (p6)
```

Each output is generated using:

```text
outputMSSample(200, replace=F)
```

## Running the simulation

With SLiM installed, the main model can be run directly as:

```bash
slim neutral.slim
```

## Generating modified SLiM scripts

`slimModifier.py` is called as:

```bash
python3 slimModifier.py neutral.slim simulation_freqs.csv
```

Python dependency:

```bash
pip install pandas
```

The generated scripts are written to:

```text
slim_scripts/
```

## Important note about `slimModifier.py`

The current `slimModifier.py` searches the input SLiM template for the markers:

```text
FIRST TARGET
SECOND TARGET
```

The current version of `neutral.slim` in this directory does **not** contain these markers. Therefore, `slimModifier.py` cannot be used with the currently uploaded `neutral.slim` without the corresponding template markers.

The main `neutral.slim` simulation itself is independent of this helper script and can be run directly with SLiM.
