# TALOS
### Temporal Analysis of Locomotion & Orthopedic Signals

**TALOS** is an R/Shiny research application for automated analysis and visualization of wearable inertial measurement unit (IMU) gait recordings.

The application transforms annotated multi-sensor gait recordings into standardized bilateral temporal gait measurements, gait variability metrics, synchronized segmental-motion visualizations, cycle-level data, quality-control information, and reproducible trial reports.

TALOS was developed as an independent computational research project exploring how wearable sensor data can be transformed into interpretable measures of human gait across orthopedic and healthy populations.

---

## Overview

Wearable IMUs provide a relatively accessible method for studying human movement outside of traditional motion-capture laboratories. However, raw multi-sensor recordings can be difficult to interpret without substantial preprocessing and synchronization.

TALOS provides a streamlined workflow:

**Processed IMU Data + Trial Metadata**

↓

**Input Validation & Quality Control**

↓

**Gait-Event Reconstruction**

↓

**Straight Walking / U-Turn Segmentation**

↓

**Bilateral Gait-Cycle Analysis**

↓

**Temporal Metrics, Variability & Visualization**

↓

**Cycle-Level Data + Trial Report**

The application is designed to preserve transparency between the original gait-event annotations and the measurements presented to the researcher.

---

## Features

### Temporal Gait Analysis

TALOS calculates bilateral:

- Stride duration
- Swing duration
- Stance duration
- Stride-time asymmetry
- Swing-time asymmetry
- Stance-time asymmetry

Bilateral temporal asymmetry is calculated as:

`100 × |L − R| / ((|L| + |R|) / 2)`

---

### Gait Variability & Consistency

Within-trial temporal variability is summarized using the coefficient of variation:

`CV (%) = SD / Mean × 100`

TALOS reports:

- Stride-duration variability
- Swing-duration variability
- Stance-duration variability
- Number of complete analyzed strides

These measurements describe within-trial timing consistency and are not interpreted as diagnostic thresholds.

---

### Bilateral Gait-Event Visualization

TALOS reconstructs annotated left and right gait events and displays them alongside foot angular-velocity signals.

The visualization allows researchers to inspect the temporal relationship between:

- Left-foot movement
- Right-foot movement
- Swing intervals
- Stance intervals
- Straight walking
- U-turn periods

Straight walking before and after the U-turn is analyzed separately so that cycles spanning the turn are not included in straight-walking temporal measurements.

---

### Segmental IMU Motion

Synchronized movement signals can be visualized across four wearable sensor locations:

- Head
- Lower back
- Left foot
- Right foot

The synchronized display provides a trial-level view of segmental movement throughout the walking protocol.

---

### Cycle-Level Data

Researchers can inspect the individual gait cycles underlying the summary statistics.

Cycle-level outputs include:

- Side
- Walking segment
- Cycle number
- Swing duration
- Stance duration
- Stride duration

This allows summary measurements to be traced back to the individual annotated gait cycles from which they were calculated.

---

### Quality Control

TALOS performs automated checks on uploaded trial data, including metadata structure, gait-event annotations, sensor availability, sampling information, and processed signal compatibility.

The interface reports trial QC status before results are interpreted.

---

### Multipart Recording Support

Processed IMU recordings stored across multiple sequential files can be loaded together.

TALOS verifies compatible schemas and reconstructs a continuous recording before analysis.

---

### Trial Reports

Trial-level results can be exported as a structured HTML report containing the primary temporal gait measurements and analysis information.

---

## Clinical Gait Cohorts

TALOS is not restricted to a single orthopedic condition.

The software has been manually tested using multiple recordings representing:

- Anterior cruciate ligament (ACL) injury
- Knee osteoarthritis (KOA)
- Hip osteoarthritis (HOA)
- Healthy gait

The analysis pipeline operates on gait-event annotations and IMU signals rather than using the pathology label to determine the calculations performed.

---

## Software Verification

TALOS was manually tested across multiple orthopedic and healthy gait recordings to evaluate:

- Successful metadata and IMU-data ingestion
- Gait-event reconstruction
- Separation of straight walking and U-turn periods
- Bilateral stride, swing, and stance calculations
- Temporal asymmetry calculations
- Gait variability calculations
- Cycle-level output consistency
- IMU signal visualization
- Quality-control behavior
- Trial report generation

This testing represents **software verification rather than clinical validation**. TALOS has not been validated as a diagnostic or clinical decision-making system.

---

## Example

![TALOS interface](images/talos_interface.png)

*Example TALOS trial analysis showing temporal gait measurements, gait variability, bilateral gait events, and synchronized segmental IMU signals.*

---

## Dataset

TALOS was developed using the publicly available:

**A Dataset of Clinical Gait Signals with Wearable Sensors from Healthy, Neurological, and Orthopedic Cohorts**

The dataset contains clinical gait recordings acquired using wearable IMUs positioned at the:

- Head
- Lower back
- Left foot
- Right foot

Recordings were sampled at 100 Hz and include annotated gait events and U-turn boundaries.

**Dataset article:**  
Voisard et al., *Scientific Data* (2025)  
DOI: `10.1038/s41597-025-05959-w`

**Dataset:**  
Figshare  
DOI: `10.6084/m9.figshare.28806086`

**Original dataset repository:**  
`github.com/CyrilVoisard/dataset_gait_1`

TALOS is an independent project and is not affiliated with the creators of the original dataset.

---

## Built With

- R
- Shiny
- ggplot2
- Plotly
- dplyr
- tidyr
- jsonlite
- DT
- stringr

---

## Running TALOS

Clone or download this repository and open `app.R` in RStudio.

Install the required R packages if necessary:

    install.packages(c(
      "shiny",
      "jsonlite",
      "dplyr",
      "tidyr",
      "ggplot2",
      "plotly",
      "DT",
      "stringr"
    ))

Then run:

    shiny::runApp()

Within TALOS, select **Load trial** and provide:

1. Processed IMU data file(s)
2. Corresponding trial metadata (`.json`)

The application will validate the inputs and generate the trial analysis automatically.

---

## Methodological Notes

Gait-event annotations are used to identify swing intervals for each foot.

Stance duration is calculated between the end of one annotated swing interval and the beginning of the subsequent swing interval for the same side.

Stride duration is calculated between consecutive same-side swing initiations.

Cycles are calculated independently within each straight-walking segment. Cycles crossing the annotated U-turn are excluded from straight-walking temporal calculations.

Temporal measurements and variability statistics produced by TALOS are descriptive research outputs.

---

## Repository Structure

    TALOS-Gait/
    │
    ├── app.R
    ├── README.md
    │
    ├── images/
    │   └── talos_interface.png
    │
    ├── docs/
    │   └── methodology.md
    │
    ├── LICENSE
    └── .gitignore

---

## Limitations

TALOS currently relies on pre-existing gait-event annotations rather than independently detecting gait events from raw IMU signals.

The application has undergone manual software testing across multiple gait recordings but has not undergone formal clinical validation.

Measures of bilateral asymmetry and temporal variability should therefore be interpreted as descriptive research measurements rather than diagnostic indicators.

Future work may include automated gait-event detection, systematic cohort-level validation, and evaluation using independent wearable gait datasets.

---

## Author

**Nidhi Sree Perla**  
Cell & Molecular Biology  
University of South Florida

Independent computational research project, 2026.

---

## Disclaimer

TALOS is intended for research and educational use. It is not a medical device and should not be used for diagnosis, treatment decisions, or clinical decision-making.
