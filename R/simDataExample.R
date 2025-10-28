set.seed(123)

par(mfrow = c(2, 2), mar = c(4, 4, 4, 1))

# Scenario 1: Low within variation & high separation between groups

data1 <- simulate_functional_data(
  within_var = 0.3, between_sep = 2.5,
  start1_mean = 2, start1_sd = 0.3,
  end1_mean = 0, end1_sd = 0.3,
  start2_mean = 0, start2_sd = 0.3,
  end2_mean = 4, end2_sd = 0.3
)
plot_functional_data(data1)


# Scenario 2: High within variation & high separation between groups

data2 <- simulate_functional_data(
  within_var = 2.0, between_sep = 2.5,
  start1_mean = 2, start1_sd = 0.5,
  end1_mean = 0, end1_sd = 0.5,
  start2_mean = 0, start2_sd = 0.5,
  end2_mean = 4, end2_sd = 0.5
)
plot_functional_data(data2)


# Scenario 3: Low within variation & low separation between groups

data3 <- simulate_functional_data(
  within_var = 0.3, between_sep = 1,
  start1_mean = 2, start1_sd = 0.3,
  end1_mean = 0, end1_sd = 0.3,
  start2_mean = 0, start2_sd = 0.3,
  end2_mean = 1, end2_sd = 0.3
)

plot_functional_data(data3)


# Scenario 4: High within variation & low separation between groups

data4 <- simulate_functional_data(
  within_var = 2.0, between_sep =1,
  start1_mean = 2, start1_sd = 0.5,
  end1_mean = 0, end1_sd = 0.5,
  start2_mean = 0, start2_sd = 0.5,
  end2_mean = 1, end2_sd = 0.5
)

plot_functional_data(data4)
