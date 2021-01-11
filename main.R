library(tidyverse)
library(moments)
library(lattice)
library(MASS)
library(lattice)
library(caret)
library(lmtest)

library(ggthemes)


old <- theme_set(theme_minimal() + theme(text = element_text(family = 'Arial', size = 20)))

# load and merge ----------------------------------------------------------------


df1 <- read_csv('/Users/han/Desktop/CASA/CASA0005 Geog Inform Sys & Sci/0005data/Crime LSOA.csv')
df2 <- read_csv('/Users/han/Desktop/CASA/CASA0005 Geog Inform Sys & Sci/0005data/Accessibility LSOA.csv')
df1$sum_2014 <- rowSums(df1[,(5:16)])
df_lsoa <- df1 %>% group_by(LSOA_Code) %>% 
  summarise(n = n(),
            LSOA_sums = sum(sum_2014)) %>%
  arrange(-LSOA_sums)

dfall <- merge(df_lsoa,df2,by.x = 'LSOA_Code',by.y = 'LSOA2011') %>% 
  dplyr::select(LSOA_Code,LSOA_sums,AvPTAI2015)

shapiro.test(dfall$LSOA_sums)
shapiro.test(dfall$AvPTAI2015)
df <- data.frame(dfall$LSOA_sums,dfall$AvPTAI2015,dfall$LSOA_Code)
par(mfrow=c(1,1))
histogram(dfall$LSOA_sums,xlab = 'LSOA crime')
qqnorm(dfall$LSOA_sums,xlab = 'LSOA crime')
histogram(dfall$AvPTAI2015,xlab = 'Average PTAL 2015')
qqnorm(dfall$AvPTAI2015,xlab = 'Average PTAL 2015')
# log transformer ------------------------------------------------------------
y <- log(dfall$LSOA_sums)
x <- log(dfall$AvPTAI2015)
z <- dfall$LSOA_Code


# symbox(~x, 
#        df_log, 
#        na.rm=T,
#        powers=seq(-3,3,by=.5))
# histogram((df_log$x)^1.5)




# plots ------------------------------------------------------------

df_log <- data.frame(x,y,z)
# OutVals = boxplot(df_log$x, plot=FALSE)$out
# df_log <- df_log %>% filter(!x %in% OutVals)
# OutVals1 = boxplot(df_log$y, plot=FALSE)$out
# df_log <- df_log %>% filter(!y %in% OutVals1)
par(mfrow=c(1,1))
histogram(df_log$x,xlab = 'Average PTAL 2015')
qqnorm(df_log$x,xlab = 'Average PTAL 2015')
histogram(df_log$y,xlab = 'LSOA crime')
qqnorm(df_log$y,xlab = 'LSOA crime')

# shapiro.test(df_log$x)
# b <- boxcox(LSOA_sums ~ AvPTAI2015, data=dfall)
# lambda <- b$x
# lik <- b$y
# bc <- cbind(lambda, lik)
# bc[order(-lik),]#结果λ=-0.1818时lik值最大，因此λ取值-0.3
# 



ggplot(df_log, aes(x = x)) + 
  geom_line(colour = "cadetblue3", stat = "density")+
  ggtitle("PTAI Probability Density Distribution")+
  theme(plot.title = element_text(hjust = 0.5))

# boxplot(df_log$x)

# outliers -------------------------------------------------------------------

mu_x <- mean(df_log$x)
sigma_x <- sd(df_log$x)
upper_line <-  mu_x + 3*sigma_x
lower_line <-  mu_x - 3*sigma_x

df_log_cleanned <- df_log %>% filter(x<upper_line ,x>lower_line) 

print(paste("Deleted abnormal data records nums：",length(df_log$x) - length(df_log_cleanned$x)))


# linear analysis ------------------------------------------------------------------

ggplot(df_log_cleanned) +
  geom_point(aes(x = x ,y = y),size = .1)+
  ggtitle("PTAI & Crimes Scatter plot")+
  theme(plot.title = element_text(hjust = 0.5))

cor_xy <- cor(df_log_cleanned$y, df_log_cleanned$x)
#cor.test(df_log_cleanned$y, df_log_cleanned$x, method = "spearman")
print(paste("pearson correlation coefficient：",round(cor_xy,2)))




# linear regression --------------------------------------------------------------------

# 初始模型
reg1 <- lm(y~x,data = df_log_cleanned)

ggplot(df_log_cleanned, aes(x = x ,y = y)) +
  geom_point(size = .1)+
  geom_smooth(method = lm)+
  ggtitle("Linear regression")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab('PTAI')+ylab('Crimes')
print(reg1)
print(summary(reg1))

print('forecast result')
print(df_log_cleanned %>% mutate(y_pred = reg1$fitted.values,
                           residuals = reg1$residuals))

library(sandwich)
library(stargazer)
library(lmtest)
#回归
reg<-lm(y~x,data=df_log_cleanned)
#求稳健标准误
rse1<-vcovHAC(reg)
rse2<-sqrt(diag(rse1))
#工作台显示
coeftest(reg,.vcov=rse1)
stargazer(reg,se=list(NULL,rse2),out="reg.doc",type="text")
# bootstrap -----------------------------------------------------------------
# 定义预测模型
train.control <- trainControl(method = "boot", number = 100)
# 训练模型
model <- train(y ~ x, data = df_log_cleanned, method = "lm",
               trControl = train.control)
# print(paste('自举重采样最终模型',model$finalModel))
print(model)
print('forecast result of Bootstrap resampling final model')
print(df_log_cleanned %>% mutate(y_pred = reg1$fitted.values,
                                 residuals = reg1$residuals))
model_coef <- function(data, index){
   coef(lm(y ~ x, data = data, subset = index))
 }
 model_coef(df_log_cleanned, 1:100)
 
 library(boot)
 boot(df_log_cleanned, model_coef, 500)




# residual&同方差性analysis -------------------------------------------------------------

print(reg1$residuals)
shapiro.test(reg1$residuals)
par(mfrow=c(2,2))
plot(reg1) # 4个图

bptest(reg1,studentize = FALSE)# p小于0.05说明具有异方差性


# spatial corr -----------------------------------------------------------------

library(sf)
library(tmap)
library(sp)
library(spdep)
par(mfrow=c(1,1))


# 读取地图数据
shp_name <- "/Users/han/Desktop/CASA/CASA0007 QM/QM_Assessment/statistical-gis-boundaries-london/ESRI/LSOA_2004_London_Low_Resolution.shp"
ncovr_sf <- st_read(shp_name)
current_style <- tmap_style("col_blind")
# plot(ncovr_sf) # 测试地图

# LSOA_CODE不能完全和ncovr_sf里面匹配上数据匹配上，后面想要做moran检验必选先在这里补全LSOA_CODE
df_log_cleanned_new <- left_join(data.frame(ncovr_sf$LSOA_CODE),df_log_cleanned,by = c("ncovr_sf.LSOA_CODE"="z"))
df_log_cleanned_new[is.na(df_log_cleanned_new)] <- 0
reg2<- lm(y~x,data = df_log_cleanned_new) 

# df_test 为增加了residuals和fitted value的dataframe
df_test <- df_log_cleanned_new 
df_test$res_reg2 <- residuals(reg2)
df_test$fitted_reg2 <- fitted(reg2)

# merge地图与原数据
ncovr_sf <- inner_join(ncovr_sf, df_test, by = c("LSOA_CODE"="ncovr_sf.LSOA_CODE"))
# 画图单独看犯罪数
tm_shape(ncovr_sf) +
  tm_polygons("y", style="quantile", title="crimes")

# 在一个图看犯罪数与可达性
# # 因为我们的区域太密集了，将可达性和犯罪数都花在一个图里面，几乎什么都看不清
# tm_shape(ncovr_sf) +
#   tm_polygons("x",
#               palette = "RdBu",title = 'PTAI') +
#   tm_shape(ncovr_sf) + tm_dots(size ='y',size.max = 200,title = 'crimes')
#   
ncovr_sf_plot <- ncovr_sf
ncovr_sf_plot$crimes <- ncovr_sf$y
tm_shape(ncovr_sf_plot) +
  tm_polygons("x",
              palette = "RdBu",title = 'PTAI') +
  tm_shape(ncovr_sf_plot) + tm_dots(col ='crimes',palette = "Set3")


# 看每个地区预测的残差
ncovr_sf$sd_breaks <- scale(ncovr_sf$res_reg2)[,1]
my_breaks <- c(-14,-3,-2,-1,1,2,3,14)
# 蓝色与红色越深都代表预测误差越大
tm_shape(ncovr_sf) + 
  tm_fill("sd_breaks", title = "Residuals", style = "fixed", breaks = my_breaks, palette = "-RdBu") +
  tm_borders(alpha = 0.1) +
  tm_layout(main.title = "Residuals", main.title.size = 0.7 ,
            legend.position = c("right", "bottom"), legend.title.size = 0.8)



# moran I

ncovr_sp <- as(ncovr_sf, "Spatial")
ncovr_sp$LSOA_CODE <- as.character(ncovr_sp$LSOA_CODE)
w <- poly2nb(ncovr_sp, row.names=ncovr_sp$LSOA_CODE)
wm <- nb2mat(neighbours = w, style='B')
rwm <- mat2listw(wm, style='W')

# moran检验，p<0.05代表存在空间自相关
lm.morantest(reg2, rwm, alternative="two.sided")

# moran值
moran_value <- moran(ncovr_sp$y, rwm, n=length(rwm$neighbours), S0=Szero(rwm))
print(moran_value$I) # 打印moran值
localmoran_value <- localmoran(ncovr_sp$y, rwm)
print(localmoran_value)
# spatial clustering
library(fpc)
mydata=data.frame(ncovr_sf$x,ncovr_sf$y)
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata,centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")


# 用fpc包中的dbscan函数进行密度聚类
dbscan_data_res <- dbscan(data.frame(ncovr_sf$x,ncovr_sf$y), eps=0.20, MinPts=6)
kmeans_data_res<-kmeans(data.frame(ncovr_sf$x,ncovr_sf$y), centers=4)
# 画出空间聚类结果
tm_shape(ncovr_sf %>% mutate(dbscan_cluster = as.character(dbscan_data_res$cluster))
) +
  tm_polygons("dbscan_cluster", style="quantile", title="dbscan")

tm_shape(ncovr_sf %>% mutate(kmeans_cluster = as.character(kmeans_data_res$cluster))
) +
  tm_polygons("kmeans_cluster", style="quantile", title="kmeans")

# 拉格朗日乘数选取替代模型
lm.LMtests(reg2, rwm, test = c("LMerr","LMlag","RLMerr","RLMlag","SARMA"))
#观测RLMerr和RLMlag统计量， 误差>滞后，所以选取空间误差模型


# 误差模型 --------------------------------------------------------------------

# 生成误差模型
reg2_err <- errorsarlm(y ~ x, data=ncovr_sp, rwm)
print(cor(reg2_err$y,reg2_err$fitted.values)^2) # 误差模型的R方
print(sum(reg2_err$y-reg2_err$fitted.values)^2)
summary(reg2_err)
ncovr_sf_err <- ncovr_sf %>%  
  mutate(reg2_err_residuals = residuals(reg2_err))
# 误差模型的残差在地图表示
qtm(ncovr_sf_err, fill = "reg2_err_residuals")
print('Error model forecast result')
print(reg2_err$fitted.values)
# 杜宾模型 --------------------------------------------------------------------
library(splm)
library('rgdal')
library('spatialreg')
library('texreg')
# 杜宾模型
sdm <- lagsarlm(y ~ x, ncovr_sp, rwm, Durbin=T)
print(cor(sdm$y,sdm$fitted.values)^2) # 杜宾模型的R方
ncovr_sf_Durbin <- ncovr_sf %>%  
  mutate(reg2_Durbin_residuals = residuals(sdm))
# 杜宾模型的残差在地图表示
qtm(ncovr_sf_Durbin, fill = "reg2_Durbin_residuals")
print(summary(sdm))
# 几个模型的数据比较
screenreg(list(reg2, reg2_err, sdm),
          custom.model.names=c("OLS","SEM","SDM"))
print(sum(sdm$y-sdm$fitted.values)^2)
# GWR ----------------------------------------------------------------

library(spgwr)

GWRbandwidth <- gwr.sel(y ~ x, 
                        data = ncovr_sp, 
                      adapt=T)
# GWR模型训练
gwr.model <-  gwr(y ~ x, 
                data = ncovr_sp, 
                adapt=GWRbandwidth, 
                hatmatrix=TRUE, 
                se.fit=TRUE)

print(gwr.model$results) # 打印GWR模型概要
print(gwr.model$SDF$localR2)
print(cor(gwr.model$lm$y,gwr.model$SDF$pred)^2) # 打印GWR模型的R方 ,R方显著提升，大于了0.7
print(summary(gwr.model))
print(sum(sdm$y-sdm$fitted.values)^2)
print(sum(gwr.model$lm$y-gwr.model$SDF$pred)^2)
# gwr.model的残差分布图
plot(x = gwr.model$lm$fitted.values,y = gwr.model$lm$residuals)
results <- as.data.frame(gwr.model$SDF)
names(results)
ncovr_sf2 <- ncovr_sf %>%
  mutate(coefx = results$x)
# 画出每个区域的回归系数
tm_shape(ncovr_sf2) +
  tm_polygons(col = "coefx", 
              palette = "RdBu", 
              alpha = 0.5)
print('GWR model forecast result')
print(gwr.model$SDF$pred)

#library(GWmodel)
#bw=bw.gwr(y ~ x, 
#          data = ncovr_sp, approach="CV",kernel="gaussian")
#gwr_model2<-gwr.basic(y ~ x, 
#                      data = ncovr_sp,  bw=bw, kernel='gaussian')
#print(gwr_model2)
# #run the significance test
# sigTest = abs(gwr.model$SDF$"x")-  2 * gwr.model$SDF$"x"
# #store significance results
# ncovr_sf2 <- ncovr_sf2 %>%
#   mutate(GWRxSig = sigTest)
# # 大于0（蓝色）意味着区域变量x（可达性）是显著的
# tm_shape(ncovr_sf2) +
#   tm_polygons(col = "GWRxSig", 
#               palette = "RdYlBu")




