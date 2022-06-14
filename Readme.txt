# HIMA2
This is the main function to estimate and testing high-dimensional mediation effects in omic studies.

# simHIMA2
This includes the R codes to generate simulation data for high-dimensional mediation analysis.

# Sample
This is a sample to generate simulation data and run high-dimensional mediation analysis.
The results are shown as follows:

    M             alpha                          alpha_SE                        beta                             beta_SE                      alpha.beta                    p_val
1 M1  0.186202499530746  0.0233729246236753  0.190685300620237  0.0539252788249321  0.0355060795992599  0.000406077463864567
2 M2  0.23760843774173    0.0254982489444596  0.275968314139476  0.0498959717919173  0.0655723999888998  3.18639839765779e-08
3 M4  0.279243031888436  0.0237937624126595  0.330119236231144  0.0527676507773055  0.0921834964098794  3.94745935005543e-10
4 M5  0.314394446452975  0.0244683773118442  0.398917101220018  0.0503509385909565  0.125417321218693    2.32343976532013e-15

References
Perera C, Zhang H, Zheng Y, Hou L, Qu A, Zheng C, Xie K, Liu L. HIMA2: High-dimensional Mediation Analysis and its application in epigenome-wide DNA methylation data. BMC Bioinformatics. 2022. Accepted.

Depends R (>= 3.4) hdi, HDMT, MASS