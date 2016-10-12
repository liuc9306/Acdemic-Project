############################################### Set up ####################################################
data <- read.table("C:/Users/Administrator/Desktop/Computational Statistics_final paper/student-mat.csv", sep = ";", header = TRUE)
head(data)
sum(is.na(data)) # no ommmited value

# get rid of Mjob, Fjob
data <- data[ ,-c(9, 10)]
dim(data) # 395 observations, 31 variables

# regression with all variables without adjustment
reg.unadj <- lm(data = data, G3 ~ .)
summary(reg.unadj)

# look at the variables
summary(data)

hist(data$famrel)

hist(data$freetime)
data$freetime <- ifelse(data$freetime < 3, "low", "high")
data$freetime <- as.factor(data$freetime)
table(data$freetime)

hist(data$goout)

hist(data$failures)
data$failures <- ifelse(data$failures < 1, "no", "yes")
data$failures <- as.factor(data$failures)
table(data$failures)

hist(data$Dalc)
data$Dalc <- ifelse(data$Dalc > 1, "yes", "no")
data$Dalc <- as.factor(data$Dalc)
table(data$Dalc)

hist(data$Walc)
data$Walc <- ifelse(data$Walc > 1, "yes", "no")
data$Walc <- as.factor(data$Walc)
table(data$Walc)

hist(data$health)
data$health <- ifelse(data$health < 3, "poor", "good")
data$health <- as.factor(data$health)
table(data$health)

hist(data$absences)
median(data$absences) # 4
data$absences <- ifelse(data$absences < 5, "no", "yes")
data$absences <- as.factor(data$absences)
table(data$absences)

table(data$higher) # extremely unequivalent

hist(data$G3)
#data$G3 <- ifelse(data$G3 < 10, 0, 1)
#table(data$G3)

cor(data$G1 + data$G2, data$G3)

# regression after adjustment the variables
reg.adj <- lm(data = data, G3 ~ .)
summary(reg.adj)
# it turns out that after adjustment, the model fits the data better. So the adjustment is reasonable. 

############################################ Goal 1 #######################################################
# best subset selection
library("leaps")
reg.full <- regsubsets(data = data, G3 ~ ., nvmax = 33)
reg.summary <- summary(reg.full)
which.max(reg.summary$adjr2) # 13
which.min(reg.summary$bic) # 4
which.min(reg.summary$cp) # 10

par(mfrow = c(1, 3))
plot(reg.summary$adjr2, xlab = "Number of Variables", ylab = "adj-r2", type = "b", ylim = c(0.841, 0.844))
plot(reg.summary$bic, xlab = "Number of Variables", ylab = "BIC", type = "b", ylim = c(-690, -650))
plot(reg.summary$cp, xlab = "Number of Variables", ylab = "cp", type = "b", ylim = c(0, 20))
par(mfrow = c(1, 1))

reg.best <- lm(data = data, G3 ~ G2 + absences + famrel + failures)
summary(reg.best)



# foward stepwise
reg.fwd <- regsubsets(data = data, G3 ~ ., nvmax = 33, method = "forward")
regfwd.summary <- summary(reg.fwd)
which.min(regfwd.summary$bic) # 4
coef(reg.fwd, 4) # same result with best subset



# backward stepwise
reg.bwd <- regsubsets(data = data, G3 ~ ., nvmax = 33, method = "backward")
regbwd.summary <- summary(reg.bwd)
which.min(regbwd.summary$bic) # 4
coef(reg.bwd, 4) # same result with best subset


# OLS assumptions verify
# 1. Normality
library("car")
qqPlot(reg.best)
shapiro.test(reg.best$residuals) # significant, so not normal
# not normal, but since big data and CLT kicks in, the tests are still reliable.

# 2. Homogeneity
residualPlots(reg.best)
# The residuals do not have a clear patter, so homogeniety

# 3. Linearity
crPlots(reg.best) 
# The green line is very close to red line, which implies the best fit.




############################################ Goal 2 #######################################################
########################################
# subset selection by cross-validation #
########################################
library("boot")

reg.full <- regsubsets(data = data, G3 ~ ., nvmax = 33)
reg.summary <- summary(reg.full)

# create a function to predict regsubsets models
predict.regsubsets = function(object, newdata, id) {
  form  <-  as.formula(~.)
  mat  <-  model.matrix(form, newdata)
  coefi  <-  coef(object, id)
  xvars  <-  names(coefi)
  mat[, xvars] %*% coefi
}

# calculate the cross-validation value for each fold with different number of variables
set.seed(12374)
k <- 10
p <- 31
folds <- sample(1:k, nrow(data), replace = TRUE)
cv.error <- matrix(0, k, p)

for (j in 1:k) {
  mod <- regsubsets(data = data[folds != j, ], G3 ~ ., nvmax = p)
  for (i in 1:p) {
    pred <- predict.regsubsets(mod, data[folds == j, ], id = i)
    cv.error[j, i] <- mean((data$G3[folds == j] - pred) ^ 2)
  }
}

cv.error1 <- apply(cv.error, 2, mean)
min <- which.min(cv.error1) # 4
error1 <- cv.error1[min]
plot(cv.error1, type = "b", xlab = "Number of variables", ylab = "CV error")
points(x = min, y = cv.error1[min], col = 2, pch = 19)

# best model picked by cross-validation
coef(reg.full, id = 4)






#############
#   Lasso   #
#############
library("glmnet")
folds # generated before

# a single loop demonstration
y <- data$G3
x <- model.matrix(G3 ~ . - 1, data = data)
lasso <- cv.glmnet(x[folds != 1, ], y[folds != 1], alpha = 1)
lambda.best <- lasso$lambda.min
lasso.pred <- predict(lasso, s = lambda.best, newx = x[folds == 1, ])
error <- mean((lasso.pred - y[folds == 1]) ^ 2)

# calculate the CV error for each fold
cv.error2 <- vector(length = k)
for (i in 1:k) {
  lasso <- cv.glmnet(x[folds != i, ], y[folds != i], alpha = 1)
  lambda.best <- lasso$lambda.min
  lasso.pred <- predict(lasso, s = lambda.best, newx = x[folds == i, ])
  cv.error2[i] <- mean((lasso.pred - y[folds == i]) ^ 2)
}

error2 <- mean(cv.error2)

# model picked by Lasso
lasso <- cv.glmnet(x, y, alpha = 1)
lambda.best <- lasso$lambda.min
coef(lasso, s = lambda.best)

plot(lasso$cvm[25:74] ~ lasso$lambda[25:74], main = "lasso", xlab = "lambda", ylab = "cv error", type = "b")
ind <- which.min(lasso$cvm)
points(x = lasso$lambda[ind], y = lasso$cvm[ind], col = "red", pch = 19)


#################
# Random Forest #
#################
library("randomForest")
set.seed(987123)
rf <- randomForest(G3 ~ ., data = data, mtry = 10, importance = TRUE)
mean(rf$mse)

error3 <- vector(length = 12)
for (i in 7:18) {
  rf <- randomForest(G3 ~ ., data = data, mtry = i)
  pred <- rf$predicted
  error3[i - 6] <- mean((pred - data$G3) ^ 2)
}
which.min(error3) # 11
plot(y = error3, x = c(7:18), xlab = "m valule", ylab = "cv error", type = "b")
points(x = 17, y = error3[11], col = "red", pch = 19)

rf <- randomForest(G3 ~ ., data = data, mtry = 17, importance = TRUE)
rf$importance
varImpPlot(rf)
rf

error3 <- min(error3)


