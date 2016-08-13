#####################################################################################
############################# Univariate Optimization ###############################
#####################################################################################

################################################ Task 1 ########################################################
# An example of a unimodal function (on [0,2])
f <- function(x)
{
  1*x^4 - 14*x^3 + 60*x^2 - 70*x
}

### The derivative of f
f_prime <- function(x)
{
  4*x^3 - 42*x^2 + 120*x - 70
}

### The 2nd derivative of f
f_dbl <- function(x)
{
  12*x^2 - 84*x + 120
}

####################################################
####### FUNCTION FOR GOLDEN SECTION SEARCH #########
####################################################
golden <- function(f, int, precision = 1e-6, max_iter = 100)
{
  # ::: This function implements the golden section search for a 
  # ::: *minimum* for the function 'f' on the range [int]
  # ::: with precision no greater than 'precision'.
  # ::: Note: 'int' is an interval such as c(2,3).
  # ::: If you want to *maximize*, multiply your function by -1.
  
  rho <- (3-sqrt(5))/2 # ::: Golden ratio
  # ::: Work out first iteration here
  f_a <- f(int[1] + rho*(diff(int)))
  f_b <- f(int[2] - rho*(diff(int)))
  ### How many iterations will we need to reach the desired precision?
  N <- ceiling(log(precision/(diff(int)))/log(1-rho))
  a <- int[1] # interval starts from a
  b <- int[2] # interval ends at b
  for (i in 1:(N)) {                   # index the number of iterations
    if (i > max_iter) stop("Exceed maximum number of iterations.")
    if (f_a < f_b) {
      b <- b - rho*(diff(c(a, b)))
      f_b <- f_a
      f_a <- f( a + rho*(diff(c(a, b))) )
    } else {
      a <- a + rho*(diff(c(a, b)))
      f_a <- f_b
      f_b <- f( b - rho*(diff(c(a, b))) )
    }
    print(paste0("Iter", i, "; [", a, ", ", b, "]."))
  }
  c(a, b)
}



##########################################################
####### FUNCTION FOR BISECTION METHOD ####################
##########################################################
bisection <- function(f_prime, int, precision = 1e-6, max_iter = 100)
{
  # ::: f_prime is the function for the first derivative
  # ::: of f, int is an interval such as c(0,1) which 
  # ::: denotes the domain
  
  N <- ceiling(log(precision/(diff(int)))/log(.5))
  f_prime_a <- f_prime(int[1] + diff(int)/2)
  a <- int[1]
  b <- int[2]
  for (i in 1:N) {
    if (i > max_iter) stop("Exceed maximum number of iterations.")
    if (f_prime_a < 0) {
      a <- a + diff(c(a, b))/2
      f_prime_a <- f_prime(a + diff(c(a, b))/2)
    } else if (f_prime_a > 0) {
      b <- b - diff(c(a, b))/2
      f_prime_a <- f_prime(a + diff(c(a, b))/2)
      } else if (f_prime_a == 0) {
        print(paste0("Iter", i, "; [", a, ", ", b, "].", " Find point x = ", (a + b)/2 ))
        break
        }
    print(paste0("Iter", i, "; [", a, ", ", b, "]."))  
  }
  c(a, b)
}


##########################################################
####### FUNCTION FOR NEWTON'S METHOD #####################
##########################################################
newton <- function(f_prime, f_dbl, precision = 1e-6, start, max_iter = 100) {
  # ::: f_prime is first derivative function
  # ::: f_dbl is second derivitive function
  # ::: start is starting 'guess'
  
  x_old <- start
  x_new <- x_old - f_prime(x_old)/f_dbl(x_old)
    
    i <- 1 # ::: use 'i' to print iteration number
  print(paste0("Iteration ", i, "; Estimate = ", x_new) )
  while (abs(x_new - x_old) > precision) {
    if (i + 1 > max_iter) stop("Exceed maximum number of iterations.")
    x_old <- x_new
    x_new <- x_old - f_prime(x_old)/f_dbl(x_old)
      # ::: redefine variables and calculate new estimate
      
      # ::: keep track of iteration history
      print(paste0("Iteration ", i + 1, "; Estimate = ", x_new) )
    i <- i + 1
  }
  x_new
}



########################################################
####### FUNCTION FOR SECANT METHOD #####################
########################################################
secant <- function(f_prime, start1, start2, precision = 1e-6, max_iter = 100)
{
  x_old1 <- start1
  x_old2 <- start2
  x_new <- x_old2 - (x_old2 - x_old1)/(f_prime(x_old2) - f_prime(x_old1)) * f_prime(x_old2)
  i <- 1 # ::: use 'i' to print iteration number
  print(paste0("Iteration ", i, "; Estimate = ", x_new) )
  while (abs(x_new - x_old2) > precision) {
    if (i + 1 > max_iter) stop("Exceed maximum number of iterations.")
    x_old1 <- x_old2
    x_old2 <- x_new
    x_new <- x_old2 - (x_old2 - x_old1)/(f_prime(x_old2) - f_prime(x_old1)) * f_prime(x_old2)
    # ::: redefine variables and calculate new estimate
    
    # ::: keep track of iteration history
    print(paste0("Iteration ", i + 1, "; Estimate = ", x_new) )
    i <- i + 1
  }
  x_new
}


################################################ Task 2 ########################################################

# An example of a function (on [-10,10])
f1 <- function(x) {
  4*x^4 - 15*x^3 + 75*x^2 - 125*x
}

f_prime1 <- function(x) {
  16*x^3 - 45*x^2 + 150*x -125
}

f_dbl1 <- function(x) {
  48*x^2 - 90*x +150
}


### plot the function from -11 to 11
curve(f1, from = -11, to = 11, n = 1001, lwd = 3) # plot the function from -10 to 10
### Create diagram to depict starting points and first picks
a0 <- -10; b0 <- 10;
### plot dashed lines from a0 and b0 to curve
text(a0 - 1,-30, "a0")
segments(a0, -100, a0, f1(a0), lty = 2, lwd = 2)
text(b0 - 1,-30, "b0")
segments(b0, -100, b0, f1(b0), lty = 2, lwd = 2)



# find the minimum of the f1 by four methods
golden(f1, c(-10, 10))
bisection(f_prime1, c(-10, 10))
newton(f_prime1, f_dbl1, start = -10)
secant(f_prime1, start1 = -10, start2 = -9.5)




################################################ Task 3 ########################################################

# An example of a function (on [0,20])
f2 <- function(x) {
  -log(x)/(1+x)
}

f_prime2 <- function(x) {
  ( x*log(x) - 1 - x )/( x*(1+x)^2 )
}

f_dbl2 <- function(x) {
  ( (1+x)^4 + 2*x^2 + 2*x - 2*x^2*x*log(x) )/( x^2*(1+x)^3 )
}


### plot the function from 0 to 20
curve(f2, from = -0.5, to = 21, n = 1001, lwd = 3) # plot the function from -1 to 21
### Create diagram to depict starting points and first picks
a0 <- 0.01; b0 <- 20;
### plot dashed lines from a0 and b0 to curve
text(a0 - 1,-30, "a0")
segments(a0, -100, a0, f2(a0), lty = 2, lwd = 2)
text(b0 - 1,-30, "b0")
segments(b0, -100, b0, f2(b0), lty = 2, lwd = 2)



# find the minimum of the f1 by four methods
golden(f2, c(0, 20))
bisection(f_prime2, c(0, 20))
newton(f_prime2, f_dbl2, start = 1)
secant(f_prime2, start1 = 1, start2 = 1.1)












#######################################################################################
############################# Multivariate Optimization ###############################
#######################################################################################

################################################ Task 4 ########################################################
# a multivaraite function within the interval of -4 <= x <= 4, -2 <= y <= 3.2
mulfun <- function(x, y) {
  10*(x^2)*y - 5*x^2 - 4*y^2 - x^4 -2*y^4
}

z <- function(xy) {
  x <- xy[ , 1]
  y <- xy[ , 2]
  z <- mulfun(x, y)
  z <- ifelse(z > -30, z, -30) # truncate the plot with z >= -30
  return(z)
}


# draw the function in 3 dimentions
x <- seq(-4, 4, 0.1)
y <- seq(-2, 3.2, 0.1)
xy <- expand.grid(x, y)
z <- matrix(z(xy), nrow = length(x), ncol = length(y))

library("graphics")
par(mfrow = c(1,3))
persp(x, y, z, phi = 45, theta = 45)
persp(x, y, z, phi = 0, theta = 90)
persp(x, y, z, phi = 90, theta = 45)



################################################ Task 5 ########################################################

# rewrite the fucntion with input of c(x, y)
f.mulfun <- function(xy) {
  x <- xy[1]
  y <- xy[2]
  10*(x^2)*y - 5*x^2 - 4*y^2 - x^4 -2*y^4
}

# derive the gradient function
grad.mulfun <- function(xy) {
  x <- xy[1]
  y <- xy[2]
  c(-4*x^3 - 10*x + 20*x*y, -8*y^3 - 8*y + 10*x^2)
}


# Gradient method when alpha is fixed (eg: 0.2)
gradient <- function(f,               # original function
                     gradf,           # gradient function
                     strtpt,          # starting point
                     maxiter = 100,   # maximum number iterations
                     alpha = .2,      # fixed step size
                     minimize = TRUE, # set to FALSE to maximize
                     epsilon = 1e-4,  # stopping criterion
                     iterhist = TRUE) # print iteration history
{
  p_old <- strtpt; error <- 1; i <- 1; check <- 1
  while(check > epsilon)
  {
    if(iterhist == TRUE) print(paste0("Iter ", i, "; f(x) = ", f(p_old)))
    if(i > maxiter) stop("Exceeded maximum number of iterations")
    p_new <- gradf(p_old)
    p_new <- p_new/sqrt(p_new %*% p_new) # Normalize the gradient vector
    ### Subtract for minimization; add for maximization
    ifelse(minimize, p_new <- p_old - alpha*p_new, p_new <- p_old + alpha*p_new)
    ### Calculate stopping criterion
    check <- sqrt((p_new - p_old) %*% (p_new - p_old)) / sqrt(p_old %*% p_old)
    ### Redefine the old point to the new point for the next iteration
    p_old <- p_new
    i <- i + 1
  } 
}


gradient(f = f.mulfun, gradf = grad.mulfun, strtpt = c(3, 2), minimize = FALSE)




################################################ Task 6 ########################################################

# Steepest method, where alpha is varying for each iterations
steepest <- function(f,               # original function
                     gradf,           # gradient function
                     strtpt,          # starting point
                     maxiter = 100,   # maximum number iterations
#                    alpha = .2,      # fixed step size
                     minimize = TRUE, # set to FALSE to maximize
                     epsilon = 1e-4,  # stopping criterion
                     iterhist = TRUE) # print iteration history
{
  p_old <- strtpt; check <- 1; i <- 1
  while(check > epsilon)
  {
    if(iterhist == TRUE) print(paste0("Iter ", i, "; f(x) = ", f(p_old)))
    if(i > maxiter) stop("Exceeded maximum number of iterations")
    p_new <- gradf(p_old)
    p_new <- p_new/sqrt(p_new %*% p_new) # Normalize the gradient vector
    ### Subtract for minimization; add for maximization
    if(minimize) { f_to_minimize <- function(alpha) { f(p_old - alpha*p_new) }
    alpha <- optimize(f_to_minimize, interval = c(0, 1), maximum = FALSE)$minimum
    p_new <- p_old - alpha*p_new}
    if(!minimize) { f_to_maximize <- function(alpha) { f(p_old + alpha*p_new) }
    alpha <- optimize(f_to_maximize, interval = c(0, 1), maximum = TRUE)$maximum 
    p_new <- p_old + alpha*p_new}
    ### Calculate stopping criterion
    check <- sqrt((p_new - p_old) %*% (p_new - p_old)) / sqrt(p_old %*% p_old)
    ### Redefine the old point to the new point for the next iteration
    p_old <- p_new
    i <- i + 1
  } 
}

steepest(f.mulfun, grad.mulfun, strtpt = c(3, 2), minimize = FALSE)





################################################ Task 7 ########################################################

# -f(x, y) function
f.mulfun <- function(xy) {
  x <- xy[1]
  y <- xy[2]
  - (10*(x^2)*y - 5*x^2 - 4*y^2 - x^4 -2*y^4)
}

# Gradient function of -f(x, y)
grad.mulfun <- function(xy) {
  x <- xy[1]; y <- xy[2]
  -1 * c(-4*x^3 - 10*x + 20*x*y, -8*y^3 - 8*y + 10*x^2)
}

# Hessian function of -f(x, y)
hess.mulfun <- function(xy) {
  x <- xy[1]; y <- xy[2]
  matrix(-1 * c(-12*x^2 - 10 + 20*y, 20*x, 20*x, -24*y^2 - 8), 2, 2)
}

# Newton's Method
LRNewton <- function(strtpt = c(0, 0), epsilon = 1e-4, f, grad, hess, 
                     maxiter = 100, iterhist = TRUE)
{
  p_old <- t(t(strtpt)); check <- 1; i <- 1
  out <- t(p_old)
  
  while (check > epsilon)
  {
    if (iterhist == TRUE) 
      print(paste0("Iter ", i, "; f(x, y) = ", round(-f(p_old), 6),
                   "; b0 = ", round(p_old[1], 2), "; b1 = ", round(p_old[2], 2)))
    if (i > maxiter) stop("Exceeded maximum number of iterations")
    p_new <- p_old - solve(hess(p_old)) %*% grad(p_old)
    check <- sqrt(t(p_new - p_old) %*% (p_new - p_old)) / sqrt(t(p_old) %*% p_old)
    ### Redefine the old point to the new point for the next iteration
    p_old <- p_new
    out <- rbind(out, t(p_new))
    i <- i + 1
  } 
  out
}

out <- LRNewton(strtpt = c(3, 2), f = f.mulfun, grad = grad.mulfun, hess = hess.mulfun)


