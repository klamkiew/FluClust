library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(tidyverse)

frame <- read.csv(paste0('../../Desktop/covid.csv'), header = TRUE, sep = ',')

require(dplyr)

mass <- frame %>% select(Country, Confirmed)

gender_mass2 <- mass %>%
  group_by(Country) %>%
  summarise(
    Confirmed = max(Confirmed, na.rm = T),
  ) %>%
  arrange(Country)


gender_mass2[gender_mass2$Country=='Afghanistan',]$MaxConfirmed
world2 <- left_join(world, gender_mass2, by = c("sovereignt"="Country"))

world2[world2$sovereignt=='United States of America', ]$Confirmed = 78039888


world <- ne_countries(scale = "medium", returnclass = "sf")

#my_countries <- c("Aruba","Afghanistan", "Morocco", "Canada")

world_modified <- world %>% 
  mutate(my_selection = ifelse(admin %in% my_countries['Country'], 1, NA))
ggplot(data = world2) +
  geom_sf(aes(fill=Confirmed)) +
  theme_classic()

ggsave(paste0('../../Desktop/test.pdf'), width = 50, height = 25, units = "cm", limitsize = FALSE)


