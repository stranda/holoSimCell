###
###
pollen_pred <- readRDS("habitat_suitability_pollen_ICE_v1.RDS")
pextent <- readRDS("habitat_suitability_pollen_xyvals.RDS")
left <- min(pextent$xvals)
right <- max(pextent$xvals)
bot <- min(pextent$yvals)
top <- max(pextent$yvals)
corners <- matrix(c(left, top,
                    left, bot,
                    right, bot,
                    right, top), ncol = 2, byrow = T))
