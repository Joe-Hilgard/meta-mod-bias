1_generate_data.R
- Uses runStudy() to generate simulated meta-analysis results, which are then saved to sim.RData
- I expect this is computationally expensive, especially with QRPs implemented
2_analysis.R
- Uses summarize_run() and plotCellMeans() and plotParamMeans() to investigate results across runs
- Generates summary stats of mean error, RMSE, and TypeI/(1-TypeII) error rates
- Saves results to final_output.csv
- Makes the plots for my old poster 

documentation.Rmd
- Seems to be a list of which parameters are taken by which functions
- Would be more typical to build with roxygen2 I guess
- TODO: won't knit, says it can't find theme_poster

Examples.R
- 