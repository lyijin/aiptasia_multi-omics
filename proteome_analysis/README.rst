=================================
Excel sheet for proteome analysis
=================================
The Excel sheet is absolutely massive because of the formulas stored in the cells--we kept them in so that people could check the way we calculated stuff. It does take a while to load, so be patient!

The meaning of individual sheets are explained in the Excel sheet itself, but reproduced here for convenience:

Terminology
-----------
NOTE: names are mainly from Herbich et al., JPR, 2013. "Statistical Inference from Multiple iTRAQ Experiments without Using Common Reference Standards"

Experiment: that thing that we did 4 separate times, each containing 8 biological samples.

Channel: 113, 114, 115, …, 121.

Spectrum: every experiment has three spectral values per protein - think of this as technical replicates of each other.

Start
-----
0. Raw data: MaxQuant-corrected values

1. Reorganise raw data
   
   Does what it says on the tin really.
   
   Experiments are given distinct colours, spectras within experiments have a different hue.

2. Remove failed samples
   
   Each experiment had one channel that completely failed (115, 118, 118 and 121 respectively).

   2a. Get mean spectral count
       
       Do some prep work for 2b--calculate the per-channel, per-gene mean spectral count if the count was above 0.
       
       e.g. mean(100, 200, 300) = 100; mean(0, 0, 200) = 200; mean(100, 0, 300) = 200.
   
   2b. Compute SCPM
       
       SCPM = spectral counts per million. This is done so that we could compare transcript vs. protein expression on a per-sample basis.
       
       (These plots are all in /compare_abs_values.)

3. log2 transformation
   
   Changed raw intensities to log2-transformed values. All further mathematical stuff operates on log2 values.
   
   In preparation for (4), calculate per-channel median value, which is stored in last row of the sheet.

4. Remove loading effect
   
   Different channels might've had different amounts of protein loaded, i.e. loading effect.
   
   To account for this, each value had had its per-channel median subtracted. As a result, per-channel median is now 0.
   
   To prepare for (5), calculate per-spectra sample medians. This is stored on the RHS of the sheet, coloured grey.

5. Median-polished log2
   
   Every cell had had its per-spectra, per-protein median value subtracted from its original value.
   
   The new value - median-polished log2 - represents the relative abundance of the protein in that sample.
   
   From what I've digested from the paper, this value is calculated on a per-experiment basis, but it CAN BE used for comparisons across experiments.
   
   NOTE: I tried swapping (4) and (5) as the paper seemingly suggests that loading effect should be removed last, but it produced far worse results (much higher standard deviations) compared to this order.

(from this point onwards,we stop relying on the paper)

6. Sort by sample
   
   Does what it says on the tin - breaks the existing per-experiment grouping into per-Aiptasia strain/temperature combination.
   
   The colour coding hints at which experiment they belonged to.

7. Stats per replicate
   
   FILTER ALERT: ONLY SAMPLES WITH TWO OR MORE SPECTRAL VALUES PROCEED PAST THIS STEP. Reason: cannot calculate sample stdev from 1 sample (there's a 1/(n-1) term somewhere…).
   
   Having said that, most samples do still have 3 median-polished log2 from 3 spectra.
   
   Calculates the average/stdev/stderr [where stderr = stdev / sqrt(n)].
   
   The fancy formulaes are basically to prettify output - stuff that gets filtered out become blanks, instead of values that look like #ERROR! in those cells.

   7a. t-test with correction
       
       Based on the values in (7), perform a t-test to see which genes show temperature-dependent shifts in protein expression.
       
       TL;DR: only one gene in one condition had a corrected p < 0.05.
   
   7b. GLM with correction
       
       Based on the values in (7), carry out a GLM to model protein expression as a function of strain and temperature, or mix of both.
       
   	   TL;DR: again, only one gene in a dependent variable had a corrected p < 0.05.

8. Stats across replicates

   FILTER ALERT: ONLY SAMPLES WITH ONE OR MORE BIOLOGICAL REPLICATE HAVING A VALID VALUE PASSES THIS STEP.
   
   Means are calculated by adding the median values up and dividing them by the count (bit strange to take the mean of median values though).
   
   Standard errors are calculated by assigning equal weightage to each observation, and combined with stderr = sqrt[(stderr1 ^ 2 + stderr2 ^ 2 + …) / n] (i.e. "Pooled variance").

9. Calculate fold changes
   
   Computes a fold change value (bear in mind this is log2-ed) as (32 - 25) if both treatments have valid values.
