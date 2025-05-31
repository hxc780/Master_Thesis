```{r}
data_XXX <-merged_data_NEW %>%
  filter(whale_index == "SMP_XXX") %>%
  select(whale_index, Sr_DIV_Ca)

data_XXX$row_number <- seq(1, by=1, length.out=nrow(data_XXX))

subset(data_XXX, Sr_DIV_Ca > 0.002)
subset(data_XXX, Sr_DIV_Ca < 0.00)
```

```{r}
data_XXX <- data_XXX[-,]
data_XXX <- data_XXX[-2,]
data_XXX <- data_XXX[-1,]

data_XXX$Sr_DIV_Ca[1]    <- (data_XXX$Sr_DIV_Ca[2] + data_XXX$Sr_DIV_Ca[3]) / 2
data_XXX$Sr_DIV_Ca[1]    <- data_XXX$Sr_DIV_Ca[3] - 2 * (data_XXX$Sr_DIV_Ca[4] - data_XXX$Sr_DIV_Ca[3]) / 4
data_XXX$Sr_DIV_Ca[2]    <- data_XXX$Sr_DIV_Ca[3] - 1 * (data_XXX$Sr_DIV_Ca[4] - data_XXX$Sr_DIV_Ca[3]) / 4
data_XXX$Sr_DIV_Ca[1424] <- data_XXX$Sr_DIV_Ca[1423] + 1 * (data_XXX$Sr_DIV_Ca[1430] - data_XXX$Sr_DIV_Ca[1423]) / 7
data_XXX$Sr_DIV_Ca[2724] <- (data_XXX$Sr_DIV_Ca[2723] + data_XXX$Sr_DIV_Ca[2725]) / 2

data_XXX$row_number <- seq(1, by=1, length.out=nrow(data_XXX))
```


```{r}
data_XXX$millimeters <- seq(0.1, by=nrow(data_XXX)/10000, length.out=nrow(data_XXX))

ggplot(data_XXX, aes(x = millimeters, y = Sr_DIV_Ca)) +
  geom_point() +
  labs(title = paste("Sr_DIV_Ca vs Artificial Millimeters for SMP_XXX"),
       x = "Millimeters",
       y = "Sr_DIV_Ca") +
  theme_bw()
```

```{r}
trend_XXX <- loess(data_XXX$Sr_DIV_Ca ~ data_XXX$millimeters)
trend1_XXX <- data.frame(trend_XXX = trend_XXX$fitted, millimeters = data_XXX$millimeters)

ggplot(data_XXX, aes(x = millimeters, y = Sr_DIV_Ca)) +
  geom_line() +
  geom_line(data = trend1_XXX, aes(y= trend_XXX), color = "blue") +
  labs(title = "Sr_DIV_Ca vs Artificial Millimeters for SMP_XXX",
       x = "Millimeters",
       y = "Sr_DIV_Ca") +
  theme_minimal()
```

```{r}
# subtracting the trend
data_XXX$no_trend <- data_XXX$Sr_DIV_Ca - trend1_XXX$trend_XXX

# dividing the new timeseries with the sd(timeserie)
data_XXX$nt_sd <- data_XXX$no_trend/sd(data_XXX$no_trend)

# see if they oscillate around 0
ggplot(data_XXX, aes(x = millimeters, y = nt_sd)) +
  geom_line() +
  labs(title = "Sr_DIV_Ca vs Artificial Millimeters for SMP_XXX",
       x = "Millimeters",
       y = "Sr_DIV_Ca") +
  theme_minimal()

# Exporting
write.csv(data_XXX, "data_XXX.csv", row.names = FALSE)
```