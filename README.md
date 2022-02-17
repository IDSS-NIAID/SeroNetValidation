This document outlines the methods required for SARS-CoV2 antibody ELISA
validation testing for SeroNet (see the [Coronavirus Luinex-based
15-PLEX Immunoasssay Validation Testing
Plan](https://abcs-amp.cancer.gov/uploads/external/300/544f6f97076053733ff7f6c3762619ea23b943f2)
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
<td style="text-align: center;">Upper 95% CI Bound</td>
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

For estimating the geometric mean of groups of replicates performed by
different analysts (e.g. LLOQ and ULOQ summary results), we use the
[`meta::metamean()`](https://www.rdocumentation.org/packages/meta/versions/4.9-6/topics/metamean)
function, with the log means (log *x̄*<sub>*g*</sub>) for each analyst,
standard deviation of the log mean for each analyst
(log *s*<sub>*g*</sub>), and the number of replicates for each analyst.

#### Geometric Standard Deviation (*s*<sub>*g*</sub>)

Given a vector of values *x*<sub>*i*</sub>, of length *n*,

$$ \\log(s\_g) = \\frac{1}{\\sqrt(n)} \\sqrt{\\sum\_{i=1}^n \\left\[\\log(x\_i) - \\log(\\bar{x}\_g)\\right\]^2}. $$

#### Relative Standard Deviation (RSD) / Coefficient of Variation (CV%)

Given the sample standard deviation, *s*, and the sample mean, *x̄*,

$$ RSD = \\frac{100 \* s}{\\bar{x}}. $$

For log-normally distributed data, RSD is calculated as a function of
the standard deviation of the log mean, *s*<sub>*g*</sub> (for more
discussion see [Koopmans, Owen, and
Rosenblatt](https://doi.org/10.1093%2Fbiomet%2F51.1-2.25)).

$$ RSD = \\sqrt{e^{s\_g^2} - 1} $$

##### Inter-assay variability

Inter-assay variability is calculated as the Coefficient of Variation
(CV) within the group (e.g. the group of all replicates of a specific
sample by a specific analyst).

##### Intra-assay variability

Intra-assay variability is calculated as the mean CV across groups,
weighting by sample size. [Estimating intra- and inter-assay variability
in salivary cortisol](https://pubmed.ncbi.nlm.nih.gov/21498487/) has
some discussion on this topic.

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

## Data format

Each tab of the excel spreadsheet provided is formatted similarly with
the following variables:

-   `Assay` - Assay being tested
-   `Sample_ID` - Differentiates each sample for a specific day. Sample
    IDs are reused each new day.
    -   This sometimes includes the replicate number along with the
        sample ID. Remove the dash and number at the end of the sample
        ID for these cases.
    -   Range of `Sample_ID` for heat tests in the `STABILITY` tab is
        `{High (1X), High (5X), High (10X), High (BIL), High (HB), High (Ctrl), ..., Low (1X), ...}`.
        The low/high designation refers to the concentration/heat level.
        Use the `1X` values as control for heat stability and `Ctrl` as
        control for matrix components.
-   `Treatment` - Only included in STABILITY tab. Treatment applied to
    the sample.
-   `Lot Status: Old or New` - Differentiates between old and new
    antigen lots.
    -   Only included in `Lot-to-Lot Antigen` and `Lot-to-Lot Conjugate`
        tabs.
    -   Not case sensitive.
-   `Analyst` - Analyst ID (initials?) to differentiate between samples
    run by different analysts.
-   `Dil_Factor` - Dilution factor
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
        above. Throw these out.
    -   Some cells are blank. These will be treated as missing data.
    -   Often called `Result (Calculated Con.)` or `Result (AU/mL)`.
-   `ACon (AU/mL)` - Arbitrary units per mL. *Not sure what to do with
    this. Ignore.*
    -   Sometimes called `ACon`, `ACon (AU/mL)`, or
        `Calculated Results (AU/mL)`.
-   `Result File` - Ignore for analysis purposes
-   `Plate Map File` - Ignore for analysis purposes
-   `Interpretation` - Contains “Negative” for some individuals in the
    CoV2 N assay. These are individuals who have been vaccinated but
    shouldn’t have any other antibodies. Drop these rows from precision
    and lot-to-lot calculations.
