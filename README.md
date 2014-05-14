Lung-Fissure-Segmentation
=========================

Uses a variation of the watershed transform to segment lung fissures in low dose CT images

Data set is provided by Anthony Reeves and the Cornell Vision Image Analysis Group. The code is run on the VisionX system. Process involves:
1. Perform gradient operator on CT image
2. Opening then closing on CT image, then find regional maxima
3. Threshold to find regional minima
4. Run watershed transform using 2,3 as seed points and catchment basins
5. Post process with lung masks and removal of noise

Worked with Caleb Woo and Ryan Ashley
