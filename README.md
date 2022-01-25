This document outlines the methods required for SARS-CoV2 antibody ELISA
validation testing for SeroNet (see the [SARS-CoV-2 Antibody ELISA
Validation Testing
Plan](https://abcs-amp.cancer.gov/uploads/singleton/27302/648c626297fe66fabf658854bb150d68e1f5bb91)
for more details).

## Tests

### Measures

The validation testing plan requires the following measures (defined
below) for each experiment. Due to the log-normal nature of these data,
the geometric mean and standard deviation are used.

<table>
<thead>
<tr class="header">
<th style="text-align: left;">Experiment</th>
<th style="text-align: center;"><span class="math inline"><em>x̄</em><sub><em>g</em></sub></span></th>
<th style="text-align: center;">RSD</th>
<th style="text-align: center;"><span class="math inline"><em>δ</em></span></th>
<th style="text-align: left;">For each</th>
<th style="text-align: left;">Acceptance</th>
<th style="text-align: center;">Value</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">LLOQ</td>
<td style="text-align: center;">x</td>
<td style="text-align: center;">x</td>
<td style="text-align: center;">x</td>
<td style="text-align: left;">dilution</td>
<td style="text-align: left;"><span class="math inline"><em>δ</em></span> ≤ 50%; RSD ≤ 30%</td>
<td style="text-align: center;"></td>
</tr>
<tr class="even">
<td style="text-align: left;">ULOQ</td>
<td style="text-align: center;">x</td>
<td style="text-align: center;">x</td>
<td style="text-align: center;">x</td>
<td style="text-align: left;">dilution</td>
<td style="text-align: left;"><span class="math inline"><em>δ</em></span> ≤ 50%; RSD ≤ 30%</td>
<td style="text-align: center;"></td>
</tr>
<tr class="odd">
<td style="text-align: left;">Linearity</td>
<td style="text-align: center;">x</td>
<td style="text-align: center;">x</td>
<td style="text-align: center;">x</td>
<td style="text-align: left;">dilution</td>
<td style="text-align: left;"><span class="math inline"><em>δ</em></span> ≤ 50%; RSD ≤ 30%</td>
<td style="text-align: center;"></td>
</tr>
<tr class="even">
<td style="text-align: left;">Cut point</td>
<td style="text-align: center;">x</td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
<td style="text-align: left;"></td>
<td style="text-align: left;"></td>
<td style="text-align: center;">95th Pctl</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Sensitivity</td>
<td style="text-align: center;">x</td>
<td style="text-align: center;">x</td>
<td style="text-align: center;">x</td>
<td style="text-align: left;">concentration level</td>
<td style="text-align: left;"><span class="math inline"><em>δ</em></span> ≤ 50%; RSD ≤ 30%</td>
<td style="text-align: center;"></td>
</tr>
<tr class="even">
<td style="text-align: left;">Accuracy</td>
<td style="text-align: center;">x</td>
<td style="text-align: center;"></td>
<td style="text-align: center;">x</td>
<td style="text-align: left;">inter/intra-day/analyst</td>
<td style="text-align: left;"><span class="math inline"><em>δ</em></span> ≤ 25%</td>
<td style="text-align: center;"></td>
</tr>
<tr class="odd">
<td style="text-align: left;">Precision</td>
<td style="text-align: center;">x</td>
<td style="text-align: center;">x</td>
<td style="text-align: center;"></td>
<td style="text-align: left;">inter/intra-day/analyst</td>
<td style="text-align: left;">RSD ≤ 25%</td>
<td style="text-align: center;"></td>
</tr>
<tr class="even">
<td style="text-align: left;">Carry-over</td>
<td style="text-align: center;">x</td>
<td style="text-align: center;">x</td>
<td style="text-align: center;"></td>
<td style="text-align: left;"></td>
<td style="text-align: left;">RSD ≤ 25% or Neg ≤ LLOQ</td>
<td style="text-align: center;"></td>
</tr>
<tr class="odd">
<td style="text-align: left;">Stability</td>
<td style="text-align: center;">x</td>
<td style="text-align: center;"></td>
<td style="text-align: center;">x</td>
<td style="text-align: left;"></td>
<td style="text-align: left;"><span class="math inline"><em>δ</em></span> ≤ 25%</td>
<td style="text-align: center;"></td>
</tr>
<tr class="even">
<td style="text-align: left;">Specificity</td>
<td style="text-align: center;">x</td>
<td style="text-align: center;"></td>
<td style="text-align: center;">x</td>
<td style="text-align: left;"></td>
<td style="text-align: left;"></td>
<td style="text-align: center;"></td>
</tr>
</tbody>
</table>

#### Geometric Mean (*x̄*<sub>*g*</sub>)

Given a vector of values, *x*<sub>*i*</sub>, of length *n*,

$$ \\log(\\bar{x}\_g) = \\frac{1}{n} \\sum\_{i=1}^n \\log(x\_i). $$

#### Geometric Standard Deviation (*s*<sub>*g*</sub>)

Given a vector of values *x*<sub>*i*</sub>, of length *n*,

$$ \\log(s\_g) = \\frac{1}{\\sqrt(n)} \\sqrt{\\sum\_{i=1}^n \\left\[\\log(x\_i) - \\log(\\bar{x}\_g)\\right\]^2}. $$

#### Relative Standard Deviation (RSD)

Given the sample standard deviation, *s*<sub>*g*</sub>, and the sample
mean, *x̄*<sub>*g*</sub>,

$$ RSD = \\frac{100 \* s\_g}{\\bar{x}\_g}. $$

#### Coefficient of Variation (CV)

Given the sample standard deviation, *s*<sub>*g*</sub>, and the sample
mean, *x̄*<sub>*g*</sub>,

$$ CV = \\frac{s\_g}{\\bar{x}\_g}. $$

##### Inter-assay variability

Inter-assay variability is calculated as the Coefficient of Variation
(CV) within the group.

##### Intra-assay variability

Intra-assay variability is calculated as the mean CV across groups.
[Estimating intra- and inter-assay variability in salivary
cortisol](https://pubmed.ncbi.nlm.nih.gov/21498487/) has some discussion
on this topic.

#### Percent Error (δ)

Given a sample mean, *x̄*<sub>*g*</sub>, and an expected experimental
value, *x*<sub>*e**x**p*</sub>,

$$ δ = \\left| \\frac{\\bar{x}\_g - x\_{exp}}{x\_{exp}} \\right| \* 100% $$

To ascertain the expected experimental value, the dilution for the first
result in each series (e.g. `ACC_01`, `LLOQ_C1`) are used to calculate
the theoretical concentrations for the remaining results. For accuracy,
each following dilution is 3-fold, and for sensitivity analysis, each
following dilution is 2-fold. Thus, if the average dilution for `ACC_01`
measures is 300, the expected experimental value would be 100. Likewise,
if the value for `LLOQ_C1` were 3300, the expected experimental value
would be 1650.

Expected values for sensitivity are treated the same, and the dilution
factor is included in the data set.

Expected values for stability and specificity are assumed to be the same
as the control group.

Expected values for accuracy are taken as the median of all values for
the group (assay and sample ID). This allows for non-normality if some
replicates fail badly.

For lot-to-lot comparisons, the expected concentrations for the *new*
group are assumed to be the same as the observed concentrations in the
*old* group.

For ULOQ and linearity, many or most of the first results in each
dilution series are out of range, so we use which ever value is closest
to 10 as the base value and calculate the expected value from that.

#### Percentile (Pctl)

The 95th percentile is the value, *x*, at which 95% of the distribution
is less than or equal to *x*.

## Data format

Each tab of the excel spreadsheet provided is formatted similarly with
the following variables:

-   `Assay` - Assay being tested
-   `Sample_ID` - Differentiates each sample for a specific day. Sample
    IDs are reused each new day.
    -   This includes the replicate number along with the sample ID for
        the `LLOQ`, `ULOQ`, and `LINEARITY` tabs. Remove the dash and
        number at the end of the sample ID for these tabs.
    -   Range of `Sample_ID` in the `STABILITY_HEAT` tab is
        `{Heat MED, Heat HIGH, LOW, MED, HIGH}`. The low/med/high
        designation refers to the concentration level, and heat/no heat
        is the treatment. This applies to other stability tabs as well.
-   `Treatment` - Treatment applied to the sample. This appears to be
    constant within each tab and can be ignored.
    -   Only included in Stability tabs.
-   `Lot Status: Old or New` - Differentiates between old and new
    antigen lots.
    -   Only included in `Lot-to-Lot Antigen` and `Lot-to-Lot Conjugate`
        tabs.
    -   Not case sensitive.
-   `Analyst` - Analyst ID (initials?) to differentiate between samples
    run by different analysts.
-   `Dil_Factor` - Dilution factor
    -   Only included in `LLOQ`, `ULOQ` and `LINEARITY` tabs.
-   `OD` - Optical Density. This is a raw instrument value and can be
    ignored.
    -   Only included in `LLOQ`, `ULOQ` and `LINEARITY` tabs.
-   `Antigen Lot Number` - Lot number. These are constant for each assay
    in each tab, with the exception of the `Lot-to-Lot Antigen` and
    `Lot-to-Lot Conjugate` tabs.
-   `Conjugate Lot Number` - Additional lot number that will help
    differentiate between lots.
-   `Day` - Day the plate was run. Starts at 1 for a series of tests.
-   `Plate_per_Day` - Plate number for the day. Starts over at 1 for
    each new day.
-   `Sample_per_plate` - Appears to be associated with a group of
    replicates for a specific sample on a given plate. For example the
    first replicate of a group of dilution factors for a given sample on
    a given plate is labeled 1. The next group of dilution factors for
    either a new replicate of the same sample or a new sample is
    labeled 2. This starts over at 1 for each plate.
-   `Replicate_per_sample` - Replicate number for each sample/treatment
    combination (e.g. LLOQ1 at dilution factor of 150 appears 4 times on
    plate 1 and has replicate numbers 1:4).
-   `Calculated Result from Titer (AU/mL)` - Calculated concentration.
    This is the primary end point.
    -   “Range?” indicates a sample is beyond the range of detection.
        Most are below the range of detection, but at least some are
        above. *Will need to consider how to handle these values.*
    -   Some cells are blank. These will be treated as missing data.
    -   Called `Result (Calculated Con.)` in the `LLOQ` and `ULOQ` tabs.
    -   Called `Result (AU/mL)` in the `LINEARITY` tab.
-   `ACon (AU/mL)` - Arbitrary units per mL. *Not sure what to do with
    this. Ignore.*
    -   Only included in `LLOQ`, `ULOQ` and `LINEARITY` tabs.
    -   Called `ACon` in the `ULOQ` tab.
    -   Called `Calculated Results (AU/mL)` in the `LINEARITY` tab.
-   `Result File` - Ignore for analysis purposes
-   `Plate Map File` - Ignore for analysis purposes
