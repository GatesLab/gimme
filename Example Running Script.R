install.packages("gimme", dependencies = T)
library(gimme)

dir.create("C:/Users/stlane/Desktop/t_120_n_25_v_5_out")

gimmeSEM(data     = "C:/Users/stlane/Desktop/t_120_n_25_v_5", 
         out      = "C:/Users/stlane/Desktop/t_120_n_25_v_5_out", 
         sep      = ",",
         header   = F,
         ar       = T, 
         subgroup = T)



         