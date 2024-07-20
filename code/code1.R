# Load necessary libraries
library(spdep)
library(lme4)

# Simulate data (for demonstration purposes)
set.seed(123)
n <- 100  # number of districts
x <- matrix(rnorm(n * 3), ncol = 3)  # auxiliary variables
beta <- c(0.5, -0.3, 0.2)
u <- rnorm(n, sd = 1)  # random effects
epsilon <- rnorm(n, sd = 0.5)  # errors
y <- x %*% beta + u + epsilon  # direct estimates

# Create a spatial weights matrix (binary contiguity)
coords <- cbind(runif(n), runif(n))  # random coordinates for districts
nb <- knn2nb(knearneigh(coords, k = 4))  # nearest neighbors
w <- nb2mat(nb, style = "W", zero.policy = TRUE)  # spatial weights matrix

# Fit the spatial model using lme4
fit <- lmer(y ~ x + (1|coords), REML = FALSE)
summary(fit)

# Extract fixed and random effects
beta_hat <- fixef(fit)
u_hat <- ranef(fit)$coords

# Estimate the spatial autocorrelation parameter (rho)
rho_hat <- cor(u_hat)

# Print results
print(beta_hat)
print(rho_hat)
