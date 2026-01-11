# Optoretinography Classification Method

This repository contains code and data for demonstrating our prolonged and multilayered optoretinography (ORG) method, which reveals light-evoked deformations in rod photoreceptors, pigment epithelium, and subretinal space.

## System Requirements

The code has been tested under the following configurations:

- **Red Hat Enterprise Linux 8.9**: MATLAB R2022b
- **macOS Ventura**: MATLAB R2021a

No non-standard hardware or customized software is required.

## File Descriptions

1. **example1.m**: Classify ORG signals from an example 5-second recording at 200 B-scans/second.
2. **example2.m**: Classify ORG signals from an example prolonged recording at 25 B-scans/second.
3. **example1-data** (will be downloaded when running example1.m):
    - `I_reg.mat`: Registered complex-valued OCT data (460 * 1000 * 1000, z * x * t).
    - `seg_mask_v2.mat`: Segmentation results.
4. **example2-data** (will be downloaded when running example2.m):
    - `I_reg.mat`: Registered complex-valued OCT data (480 * 1000 * 650, z * x * t).
    - `seg_mask_v2.mat`: Segmentation results.
5. **functions folder**: Functions for OCT data extraction and analysis.
6. **LICENSE**: License details.
7. **README**: This document.
8. **SVM folder**: Pre-defined parameters and pre-trained SVM model.
    - `PCA_coeff.mat`: Principal component coefficients and pre-defined parameters for constructing the spatiotemporal feature space.
    - `SVM_Mdl.mat`: Pre-trained SVM model.

## Instructions for Use

1. Run `example_1.m` or `example_2.m` in MATLAB.
2. Performance:
   - On macOS with an Apple M1 chip, `example_1.m` ~ 80 seconds, `example_2.m` ~ 50 seconds (excluding data download time).
3. Output:
   - Analysis results saved in "analysis" and "figures" folders within the dataset directory.

## Citation

Tan, B., Li, H., Zhuo, Y., Han, L., Mupparapu, R., Nanni, D., Barathi, V. A., Palanker, D., Schmetterer, L. & Ling, T. Light-evoked deformations in rod photoreceptors, pigment epithelium and subretinal space revealed by prolonged and multilayered optoretinography. Nat. Commun. 15, 5156 (2024). [https://doi.org/10.1038/s41467-024-49014-5](https://doi.org/10.1038/s41467-024-49014-5)


## Contact Information

**Huakun Li**
- Email: [HUAKUN001@e.ntu.edu.sg](mailto:HUAKUN001@e.ntu.edu.sg)

**Tong Ling**
- Email: [tong.ling@ntu.edu.sg](mailto:tong.ling@ntu.edu.sg)
