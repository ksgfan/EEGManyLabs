# EEGManyLabs

This folder contains experimental code for the replication of Vogel, E. K., & Machizawa, M. G. (2004). Neural activity predicts individual differences in visual working memory capacity. Nature, 428(6984), 748–751.

Work is still in progress.

---

## Project Overview

We collected raw EEG data from **10 different laboratories**, each using different acquisition systems and file formats.

The processing workflow consists of:

1. **Loading and visually screening raw data**
2. **Loading the corresponding channel-location files**
3. **Merging CDA blocks (EEG and behavioral data)**
4. **Saving unified MATLAB data structures**

These steps produce the following files:

- `id_CDA.mat` – CDA task EEG + behavior  
- `id_Eye.mat` – Eye-tracking or EOG calibration data (if available)  
- `id_Resting.mat` – Resting-state EEG (if available)

All scripts for this stage are located in: `scripts/a_prepare_data/`


---

## BIDS Conversion

After harmonizing the raw data, all datasets were converted to **BIDS** using: `scripts/a1_convert_to_bids/`




