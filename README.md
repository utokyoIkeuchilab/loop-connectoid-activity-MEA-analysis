# MEA Analysis for Multi-Organoid Loop Cerebral Connectoids

MATLAB code for the Multi-Electrode Array (MEA) data analysis accompanying the paper:

> **Multi-Organoid Loop Cerebral Connectoids Exhibit Enhanced Neuronal Network Dynamics and Sequence-Specific Entrainment**
>
> Tomoya Duenki<sup>1,2,3,4</sup>, Yoshiho Ikeuchi<sup>1,2,3,4,\*</sup>
>
> <sup>1</sup> Institute of Industrial Science, The University of Tokyo, Meguro, Tokyo 153-8505, Japan  
> <sup>2</sup> Institute for AI and Beyond, The University of Tokyo, Bunkyo, Tokyo 113-8655, Japan  
> <sup>3</sup> Department of Chemistry and Biotechnology, The University of Tokyo, Bunkyo, Tokyo 113-8655, Japan  
> <sup>4</sup> LIMMS, CNRS–Institute of Industrial Science, IRL 2820, The University of Tokyo, Tokyo, Japan

---

## Repository Structure

```
.
├── Main_raw_trace_analysis.m                       # Raw trace analysis (MATLAB 2020a)
├── Main_raw_trace_analysis_Matlab2022a_and_newer.m # Raw trace analysis (MATLAB 2022a+)
├── Main_correlation_analysis.m                     # Correlation analysis
├── Main_burst_analysis.m                           # Burst analysis
├── helper_functions/
│   ├── filter_signal.m
│   ├── spike_detection.m
│   ├── wavelet_coherence.m
│   └── wavelet_transformation.m
├── example_data/
│   ├── example_signal_trace.mat   # 10-min raw signals from 4 organoids (12-week DIV)
│   └── example_spike_data.mat     # Extracted spike activity from 64 electrodes
└── README.md
```

---

## Requirements

- **MATLAB 2020a** (or 2022a+ — use `Main_raw_trace_analysis_Matlab2022a_and_newer.m` for newer versions)
- Required Toolboxes:
  - Signal Processing Toolbox
  - Wavelet Toolbox
  - Statistics and Machine Learning Toolbox
  - Parallel Computing Toolbox
  - Image Processing Toolbox

---

## Getting Started

1. Clone or download this repository.
2. Open MATLAB and add the repository root and `helper_functions/` to the path:
   ```matlab
   addpath(genpath('path/to/this/repo'))
   ```
3. Run the desired main script (see below). Example data in `example_data/` is loaded automatically.

---

## Analysis Scripts

### 1. Raw Trace Analysis — `Main_raw_trace_analysis.m`
Processes raw electrode signals and extracts signal features.

- Preprocessing and filtering
- Spike detection
- Signal plotting
- Time–frequency analysis
- Wavelet coherence
- Frequency separation

> For MATLAB 2022a or newer, use `Main_raw_trace_analysis_Matlab2022a_and_newer.m` instead.

### 2. Correlation Analysis — `Main_correlation_analysis.m`
Quantifies synchrony within and between organoids.

- Summarizing spiking data within each organoid
- Inter-organoid correlation (direct & indirect connections)
- Intra-organoid correlation
- Spike activity plots and correlation matrices

### 3. Burst Analysis — `Main_burst_analysis.m`
Characterizes network burst dynamics.

- Burst detection
- Burst frequency
- Interburst interval (IBI) and IBI coefficient of variation
- Burst peak size and burst size coefficient of variation
- Consecutive burst detection
- Spontaneous Activity Transient (SAT) duration

---

## Example Data

`example_data/` contains two `.mat` files from a 10-minute recording of 4 organoids connected and differentiated for 12 weeks:

| File | Contents |
|------|----------|
| `example_signal_trace.mat` | Raw signal traces (4 organoids) |
| `example_spike_data.mat` | Extracted spike activity (64 electrodes) |

---

## Citation

If you use this code in your work, please cite:

```
Duenki T, Ikeuchi Y. Multi-Organoid Loop Cerebral Connectoids Exhibit Enhanced Neuronal
Network Dynamics and Sequence-Specific Entrainment. [Journal, Year]. DOI: [add DOI]
```

---

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.  
Copyright (c) 2025 Univ Tokyo IIS Ikeuchi lab.
