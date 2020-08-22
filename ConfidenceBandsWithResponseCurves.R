#95$-Confidence Interval around predictions along an environment gradient for all three
#prediction types

#dataframe for the mean of the predictions for every observation

mean_pred <- data.frame(single = c(pred_cp_uv, p_cp_gj_un, p_cp_mvco),
                        type = c(rep("Univariate", length(pred_cp_uv)),
                                 rep("Unconditional Multivariate", length(p_cp_gj_un)),
                                 rep("Conditional Multivariate", length(p_cp_mvco))),
                        IA_500 = rep(test$IA_500, 3),
                        NDVIBEF_2000 = rep(test$NDVIBEF_2000, 3))

uv <- pred_cpuv %>% t %>% as.vector()
un <- pred_cp_gj_un$ychains[,1:140] %>% t %>% as.vector()
co <- pred_cp_mvco$ychains[,1:140] %>% t %>% as.vector()
d_gg_bpenv <- data.frame(pred = c(uv, un, co), IA_500 = rep(test$IA_500, nrow(pred_cpuv)),
                         NDVIBEF_2000 = rep(test$NDVIBEF_2000, nrow(pred_cpuv)),
                         type = c(rep("Univariate", length(uv)),
                                  rep("Unconditional Multivariate", length(un)),
                                  rep("Conditional Multivariate", length(co)))) 


#First Plot a line with the mean-predictions
ggplot(mean_pred, aes(IA_500, single, color = type))+
  geom_point()
#That doesnt make sense, maybe with the response curves? But, then it has nothing to do with 
#the test data