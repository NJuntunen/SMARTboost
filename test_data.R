pacman::p_load(tidymodels)

data <- diamonds %>%
  select(x,y,z) %>%
  mutate(date = seq(as.Date("1990-01-01"), by = "day", length.out = length(.$x)))


SMARTrecipe <- recipes::recipe(data) %>%
  update_role(x, new_role = "outcome") %>%
  update_role(c(y,z), new_role = "predictor") %>%
  update_role(date, new_role = "date")

x <- SMARTrecipe


test <- SMARTboost(x ~ z+y,data = diamonds)

predict(test, diamonds)
