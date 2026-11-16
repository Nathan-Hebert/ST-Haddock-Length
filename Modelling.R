################################################################################
#
# A spatiotemporal 4X5Y haddock length-at-age model in sdmTMB
#
################################################################################

# Load necessary libraries
library(sdmTMB); library(sdmTMBextra); library(sf); library(tidyverse)
library(rnaturalearth); library(viridis); library(patchwork)

# Set font size for figures
font_size <- 16

# Set plus group and start age
start_age <- 1
plus_group <- 8

###############################Helper functions#################################

# Converts sf UTM units from m to km without losing CRS
sf_m_to_km <- function(obj)
{
  obj_km <- st_geometry(obj)/1000
  st_crs(obj_km) <- st_crs(obj)
  st_geometry(obj) <- obj_km 
  
  return(obj)
}

# Gets data ready for a qqplot for each age based on scaled DHARMa residuals
get_qq_dat <- function(model, data, nsims) {
  
  # Calculate scaled DHARMa residuals and append to dataset
  dh <- simulate(model, nsim = nsims, type = "mle-mvn") %>%
    dharma_residuals(model, return_DHARMa = TRUE)
  data$scaled_residuals <- dh$scaledResiduals
  
  # Generate data to return for plotting
  group_by(data, TAGE) %>%
    group_split() %>%
    purrr::map_dfr(
      ~{
        u <- .x$scaled_residuals
        n <- length(u)
        m <- seq_len(n) / (n + 1)
        qqplot(m, u, plot.it = FALSE)
      },
      .id = "TAGE"
    ) %>%
    mutate(age_label = factor(TAGE,
                              levels = start_age:plus_group,
                              labels = age_labels_full))
}

################################Model-fitting###################################

# Grab land data for the barrier and figures, convert to UTM Zone 20N in km
land <- ne_countries(scale = "large", returnclass = "sf", 
                     continent = "North America") %>%
  st_crop(c(xmin = -75, ymin = 40, xmax = -50, ymax = 50)) %>%
  st_transform(crs = 32620) %>% sf_m_to_km()

# Load in the raw survey data, convert to UTM Zone 20N in km, and then filter to
# the study area, summer months, correct ages, and Western II-A gear... calculate 
# important cols
survey_dat <- readRDS("Data/Survey/raw_survey_data.rds") %>% 
  st_as_sf(coords= c("X","Y"),  crs = 4326) %>% 
  st_transform(crs = 32620) %>% sf_m_to_km() %>%
  filter(STRAT %in% 470:495 & MONTH %in% 6:8 & YEAR %in% c(1983:2020, 2022:2023)) %>%
  filter(AGE > start_age-1) %>%
  filter(!(FLEN == 7 & AGE == 4)) %>% # Drop likely measurement/aging error
  mutate(COHORT = YEAR - AGE,
         DayOfYear = as.integer(format(SDATE, "%j")),
         DayOfYearScaled = scale(DayOfYear),
         BirthDecade = paste0(floor(COHORT / 10) * 10, "s"),
         SamplingDecade = paste0(floor(YEAR / 10) * 10, "s"),
         TAGE = pmin(AGE, plus_group), f_tage = factor(TAGE),
         weight_norm = weight / mean(weight))
survey_dat[c("x","y")] <- st_coordinates(survey_dat)

# Set up and plot the mesh... include a barrier
mesh_sdmTMB <- make_mesh(survey_dat[,,drop = T], 
                         xy_cols = c("x","y"), 
                         cutoff = 16) %>% 
  add_barrier_mesh(land, range_fraction = 0.1, proj_scaling = 1, plot = TRUE)

# Fit the base model
nyear <- length(unique(survey_dat$YEAR))
map_list <- list(ln_tau_Z = factor(rep(1, nyear)))
mod_len <- sdmTMB(FLEN ~ 1 +
                    s(YEAR, f_tage, bs = "fs", k = 30) +
                    DayOfYearScaled:f_tage, 
              spatial = "on", mesh = mesh_sdmTMB, 
              spatial_varying = ~ 0 + factor(YEAR),
              time = "TAGE", spatiotemporal = "rw",
              data = survey_dat[,,drop = T],
              family = Gamma(link = "log"),
              control = sdmTMBcontrol(map = map_list),
              weights = survey_dat$weight_norm,
              dispformula = ~ 0 + f_tage)
sanity(mod_len)

# Fit the model without age-varying dispersion, to show it has worse residuals 
# for young fish (sensitivity)
mod_len_common <- sdmTMB(FLEN ~ 1 +
                    s(YEAR, f_tage, bs = "fs", k = 30) +
                    DayOfYearScaled:f_tage, 
                  spatial = "on", mesh = mesh_sdmTMB, 
                  spatial_varying = ~ 0 + factor(YEAR),
                  time = "TAGE", spatiotemporal = "rw",
                  data = survey_dat[,,drop = T],
                  family = Gamma(link = "log"),
                  control = sdmTMBcontrol(map = map_list),
                  weights = survey_dat$weight_norm)
sanity(mod_len_common)

################################Raw data plots##################################

# Set up plot labels for age groups
if (plus_group <= 8) {
  age_breaks <- start_age:plus_group
} else {
  age_breaks <- c(seq(start_age, plus_group - 1, by = 2), plus_group)
}
age_labels <- ifelse(age_breaks == plus_group, 
                     paste0(age_breaks, "+"), age_breaks)
age_labels_full <- c(paste0("Age ", start_age:(plus_group - 1)),
                     paste0("Age ", plus_group, "+"))

# Plot the sample size by age and year
survey_dat %>%
  group_by(YEAR, AGE) %>%
  summarise(n = n(), .groups = "drop") %>%
  complete(YEAR = 1983:2023, AGE, fill = list(n = 0)) %>%
  ggplot(aes(x = AGE, y = YEAR, fill = n)) +
  geom_tile() +
  geom_text(aes(label = n), size = 2.5, colour = "white") +
  scale_fill_viridis_c(limits = c(0, NA)) +
  scale_y_reverse(breaks = seq(1983, 2023, by = 10)) +
  coord_cartesian(expand = F) +
  scale_x_continuous(breaks = seq(start_age, max(survey_dat$AGE), by = 2)) +
  labs(x = "Age", y = "Survey year", fill = "Sample\nsize (n)") +
  theme_classic() +
  theme(text = element_text(size = font_size))
ggsave("Figs/sample_size.jpeg", width = 8, height = 6, units = "in")

# Plot mean normalized weights by length with loess smoother on the mean
survey_dat %>%
  st_drop_geometry() %>%
  group_by(FLEN) %>%
  summarise(mean_weight = mean(weight_norm), n = n(), .groups = "drop") %>%
  ggplot(aes(FLEN, mean_weight)) +
  geom_point(aes(size = n), colour = "grey40", alpha = 0.7) +
  geom_smooth(aes(weight = n),
              method = "loess", span = 0.4,
              se = FALSE, linewidth = 1.5, colour = "#3366FF") +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5) +
  scale_size_area(max_size = 3) +
  guides(size = "none") +
  labs(x = "Length (cm)", y = "Mean normalized HT weight") +
  theme_classic() +
  theme(text = element_text(size = font_size)) +
  coord_cartesian(expand = FALSE, ylim = c(0,1.7))
ggsave("Figs/weights.jpeg", width = 8, height = 6, units = "in")

# Plot the raw survey time series... gaps for 2021 and when there's no data 
# coverage
survey_dat %>%
  group_by(YEAR, TAGE) %>%
  summarize(mean = mean(FLEN, na.rm = TRUE), .groups = "drop") %>%
  complete(YEAR = seq(min(YEAR), max(YEAR), by = 1),
           TAGE = start_age:plus_group) %>%
  mutate(age_label = factor(TAGE, levels = start_age:plus_group,
                            labels = age_labels_full)) %>%
  arrange(TAGE, YEAR) %>% 
  ggplot(aes(x = YEAR, y = mean)) +
  geom_point(size = 1) + geom_path(size = 1.5) +
  theme_classic() + theme(text = element_text(size = font_size)) +
  coord_cartesian(xlim = c(1983, 2026), expand = F) +
  scale_x_continuous(breaks = c(1983,2003,2023)) +
  labs(x = "Survey year", y = "Mean length (cm)") +
  facet_wrap(~age_label, scales = "free_y")
ggsave("Figs/raw_time_series_faceted.jpeg", width = 8, height = 6.5, units = "in")

# Grab groundfish survey strata, filter out 4VW, and then convert to UTM Zone 20N 
# in km... ditch strata... used to delineate the study area
survey.strata <- st_read(paste0("R:/Science/Population Ecology Division/Shared/Ha",
                                "ddock 4X5Y/2023/Shapefiles/MaritimesRegionEcosyst",
                                "emAssessmentStrata(2014-)NAD83.shp"))
survey.strata.cropped <- survey.strata %>% 
  filter(StrataID %in% 470:495) %>%
  st_transform(crs = 32620) %>% sf_m_to_km() %>% 
  st_union()

# Common theme and map extent for spatial plots
spatial_theme <- theme_classic() +
  theme(text = element_text(size = font_size),
        plot.title = element_text(size = font_size),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_rect(fill = "lightblue", colour = NA),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        legend.position = "right",
        legend.key.height = unit(1.5, "cm"),
        legend.key.width = unit(0.5, "cm"))
spatial_coord <- coord_sf(xlim = c(130, 465), ylim = c(4650, 5055))

# Yearly plots of the survey tows overlaid on the study area
ggplot() +
  geom_sf(data = land, fill = "darkgrey") +
  geom_sf(data = survey.strata.cropped, 
          fill = "white", col = "black") +
  geom_sf(data = survey_dat, size = 1) +
  spatial_theme + spatial_coord +
  facet_wrap(~YEAR, ncol = 8) +
  ggtitle("Survey tows")
ggsave("Figs/survey_tow_locations.jpeg", width = 7, height = 7, units = "in")

################################Residual plots##################################

# Generate a DHARMa qqplot for each age for the base model
get_qq_dat(mod_len, survey_dat, nsims = 1000) %>%
  ggplot(aes(x, y)) + 
  facet_wrap(vars(age_label)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  theme(text = element_text(size = font_size)) +
  labs(x = "Theoretical quantiles",
       y = "Sample quantiles") +
  scale_x_continuous(breaks = c(0.1, 0.5, 0.9)) +
  scale_y_continuous(breaks = c(0.1, 0.5, 0.9))
ggsave("Figs/DHARMa_by_age.jpeg", width = 7.7, height = 6.6, units = "in")

# Generate a DHARMa qqplot for each age for the model without age-varying 
# dispersion (sensitivity)
get_qq_dat(mod_len_common, survey_dat, nsims = 1000) %>%
  ggplot(aes(x, y)) +  
  facet_wrap(vars(age_label)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  theme(text = element_text(size = font_size)) +
  labs(x = "Theoretical quantiles",
       y = "Sample quantiles") +
  scale_x_continuous(breaks = c(0.1, 0.5, 0.9)) +
  scale_y_continuous(breaks = c(0.1, 0.5, 0.9))
ggsave("Figs/DHARMa_by_age_common.jpeg", width = 7.7, height = 6.6, units = "in")

# Plot response residuals by age, faceted by decade of survey... weight smoother 
# by model weights
survey_dat$residuals_len <- survey_dat$FLEN-predict(mod_len, type = "response")$est
ggplot(survey_dat, aes(x = TAGE, y = residuals_len)) + 
  geom_point(size = 0.7) + geom_smooth(se = F, size = 1, method = "loess", 
                                       aes(weight = weight_norm)) + 
  geom_hline(yintercept = 0, size = 0.5) +
  labs(y = "Response residual (cm)", x = "Age") +
  theme_classic() + theme(text = element_text(size = font_size)) +
  scale_x_continuous(breaks = age_breaks, labels = age_labels) +
  facet_wrap(~SamplingDecade, nrow = 3)
ggsave("Figs/residuals_age_sampling_decade.jpeg", width = 6, height = 7, 
       units = "in")

# Plot response residuals by age, faceted by birth decade... weight smoother by 
# model weights
ggplot(survey_dat, aes(x = TAGE, y = residuals_len)) + 
  geom_point(size = 0.7) + geom_smooth(se = F, size = 1, method = "loess",
                                       aes(weight = weight_norm)) + 
  geom_hline(yintercept = 0, size = 0.5) +
  labs(y = "Response residual (cm)", x = "Age") +
  theme_classic() + theme(text = element_text(size = font_size)) +
  scale_x_continuous(breaks = age_breaks, labels = age_labels) +
  facet_wrap(~BirthDecade, nrow = 3)
ggsave("Figs/residuals_age_birth_decade.jpeg", width = 6, height = 7, 
       units = "in")

############################Model dispersion plot###############################

# Get CIs for age-specific CVs using the delta method
s <- summary(mod_len$sd_report)
disp <- as.data.frame(s[grep("^b_disp_k", rownames(s)), , drop = FALSE])
disp$age <- c(start_age:(plus_group - 1), paste0(plus_group, "+"))
disp$cv <- exp(-disp$Estimate / 2)
disp$cv_se <- 0.5 * disp$cv * disp$`Std. Error`
disp$cv_lwr <- disp$cv - 1.96 * disp$cv_se
disp$cv_upr <- disp$cv + 1.96 * disp$cv_se

# Using the above, generate the plot
ggplot(disp, aes(age, cv*100)) +
  geom_ribbon(aes(ymin = cv_lwr*100, ymax = cv_upr*100, group = 1),
              alpha = 0.3) +
  geom_line(aes(group = 1), linewidth = 2) +
  geom_point(size = 3) +
  labs(x = "Age",
       y = "Coefficient of variation (%)") +
  theme_classic() +
  theme(text = element_text(size = font_size)) +
  coord_cartesian(xlim = c(start_age, plus_group+0.2),
                  ylim = c(7, 12),
                  expand = FALSE, clip = "off") +
  scale_y_continuous(breaks = seq(7, 12, by = 1))
ggsave("Figs/model_CVs.jpeg", width = 7.7, height = 6.6, units = "in")

##############################Model predictions#################################

##### Plot growth curves for cohorts... highlight a few cohorts #####

# Set up the prediction grid across cohorts and ages
pred_dat_cohort <- expand.grid(COHORT = seq(min(survey_dat$COHORT), 
                                            max(survey_dat$COHORT), by = 1), 
                               AGE = start_age:(plus_group - 1), 
                               DayOfYearScaled = 0) %>% 
  mutate(TAGE = AGE, YEAR = COHORT + AGE, f_tage = factor(TAGE)) %>% 
  filter(YEAR %in% unique(survey_dat$YEAR))

# Make the predictions... split cohort trajectories where ages/years are missing
predictions <- predict(mod_len, pred_dat_cohort, type = "response", re_form = NA) %>%
  group_by(COHORT) %>% arrange(AGE, .by_group = TRUE) %>%
  mutate(segment = cumsum(c(0, diff(AGE)) > 1)) %>% ungroup()

# Ready the highlighted cohorts for plotting
highlight_cohorts <- seq(1983, 2013, by = 15)
label_dat <- predictions %>%
  filter(COHORT %in% highlight_cohorts) %>%
  group_by(COHORT) %>%
  arrange(AGE) %>%
  slice_tail(n = 3) %>%
  slice(if (cur_group_id() == n_distinct(COHORT)) 2 else 3) %>%
  ungroup() %>%
  mutate(label = COHORT,
         below = COHORT == max(highlight_cohorts),
         label_y = est + ifelse(below, -0.8, 0.8),
         label_vjust = ifelse(below, 1.5, -0.3))

# Generate the plot
len_age_effect <- ggplot(predictions, aes(x = AGE, y = est, 
                                          group = interaction(COHORT, segment), 
                                          colour = COHORT)) +
  geom_line(linewidth = 0.6, alpha = 0.25) +
  geom_line(data = subset(predictions, COHORT %in% highlight_cohorts),
            linewidth = 2, alpha = 1) +
  geom_text(data = label_dat, 
            aes(x = AGE, y = label_y, label = label, vjust = label_vjust, 
                colour = COHORT),
            hjust = 1, size = 6, fontface = "bold", inherit.aes = FALSE) +
  scale_colour_viridis_c("Cohort", breaks = as.numeric(round(with(pred_dat_cohort, 
                                       quantile(COHORT, probs = c(0, 0.5, 1)))))) +
  theme_classic() +
  theme(text = element_text(size = font_size), 
        legend.position = c(0.98, 0.05), legend.justification = c(1, 0)) +
  coord_cartesian(xlim = c(start_age, plus_group-0.7), ylim = c(NA, 62),
                  expand = FALSE, clip = "off") +
  scale_x_continuous(breaks = start_age:(plus_group - 1)) +
  labs(x = "Age", y = "Predicted length (cm)")

##### Plot the conditional effect of year across ages #####

# Set up the prediction grid across years and ages
pred_dat_year <- expand.grid(YEAR = seq(1983, 2023, by = 1), 
                             AGE = start_age:plus_group, 
                             DayOfYearScaled = 0) %>% 
  mutate(TAGE = AGE, COHORT = YEAR-AGE, f_tage = factor(TAGE),     
         age_label = factor(AGE, levels = start_age:plus_group, 
                            labels = age_labels_full)) %>% 
  filter(YEAR%in%unique(survey_dat$YEAR))

# Make the predictions, ensuring a gap in 2021
predictions <- predict(mod_len, pred_dat_year, se_fit = T, re_form = NA)
predictions_2021gap <- bind_rows(predictions, predictions %>%
                                   filter(YEAR == 2020) %>%
                                   mutate(YEAR = 2021, est = NA, est_se = NA)) %>% 
  arrange(AGE, YEAR)

# Generate the unfaceted plot
len_year_effect <- ggplot(data = predictions_2021gap, 
                          aes(x = YEAR, y = exp(est), 
                              ymin = exp(est-1.96*est_se), 
                              ymax = exp(est+1.96*est_se), 
                              fill = AGE, col = AGE, group = AGE)) + 
  geom_ribbon(alpha = 0.3, col = NA) + geom_path(size = 2) +
  theme_classic() + theme(text = element_text(size = font_size)) + 
  coord_cartesian(xlim = c(1983, 2025), expand = F) +
  scale_x_continuous(breaks = c(1983,2003,2023)) +
  labs(x = "Survey year", y = "Predicted length (cm)", 
       fill = "Age", col = "Age") + 
  scale_fill_viridis_c(breaks = age_breaks, labels = age_labels) + 
  scale_colour_viridis_c(breaks = age_breaks, labels = age_labels)
ggsave("Figs/year_effect.jpeg", width = 9, height = 7, units = "in")

# Now do the same plot, faceted by age
len_year_effect_faceted <- ggplot(data = predictions_2021gap, 
                                  aes(x = YEAR, y = exp(est), 
                                      ymin = exp(est-1.96*est_se), 
                                      ymax = exp(est+1.96*est_se))) + 
  geom_ribbon(alpha = 0.3, col = NA) + geom_path(size = 2) +
  theme_classic() + theme(text = element_text(size = font_size)) + 
  coord_cartesian(xlim = c(1983, 2025), expand = F) +
  scale_x_continuous(breaks = c(1983,2003,2023)) +
  labs(x = "Survey year", y = "Predicted length (cm)") + 
  facet_wrap(~age_label, scales = "free_y")
ggsave("Figs/year_effect_faceted.jpeg", width = 9, height = 7, units = "in")

##### Combine the year and cohort growth plots #####

len_year_effect_faceted + len_age_effect + plot_layout(ncol = 1) +
  plot_annotation(tag_levels = "a", tag_suffix = ")")
ggsave("Figs/year_cohort_combined.jpeg", width = 8, height = 12, units = "in")

##### Quantify length change between the early time series and late time series #####

predictions_2021gap %>%
  filter(YEAR %in% c(1983:1987, 2019:2023)) %>%
  mutate(period = ifelse(YEAR <= 1987, "early", "late")) %>%
  group_by(AGE, period) %>%
  summarise(mean_est = mean(exp(est), na.rm = T), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = period, values_from = mean_est) %>%
  mutate(pct_change = 100 * (late - early) / early)

##### Plot growth (length change per age increment) trend across cohorts #####

# Simulate predictions for cohorts from the model
set.seed(123)
pred_sim_cohort <- predict(mod_len, pred_dat_cohort, 
                           type = "link", re_form = NA, nsim = 1000)

# Get predictions on the response scale and then sort
pred_sim_cohort_cm <- exp(pred_sim_cohort)
ord_cohort <- order(pred_dat_cohort$COHORT, pred_dat_cohort$AGE)
pred_ord_cohort <- pred_dat_cohort[ord_cohort, ]
sim_ord_cohort <- pred_sim_cohort_cm[ord_cohort, ]

# Calculate raw annual growth between consecutive ages within a cohort
age_diff <- pred_ord_cohort$AGE[-1] - pred_ord_cohort$AGE[-nrow(pred_ord_cohort)]
same_cohort_year <- pred_ord_cohort$COHORT[-1] == pred_ord_cohort$COHORT[-nrow(pred_ord_cohort)]
valid_growth <- same_cohort_year & age_diff == 1
growth_sim_cohort <- sim_ord_cohort[-1, , drop = FALSE] - 
  sim_ord_cohort[-nrow(sim_ord_cohort), , drop = FALSE]
growth_sim_cohort <- growth_sim_cohort[valid_growth, , drop = FALSE]

# Summarize predicted growth across cohorts using median and quantiles
growth_cohort_dat <- pred_ord_cohort[-1, ] %>%
  filter(valid_growth) %>%
  mutate(age_interval = paste0("Age ", AGE - 1, "–age ", AGE),
         growth = apply(growth_sim_cohort, 1, median),
         growth_low = apply(growth_sim_cohort, 1, quantile, probs = 0.025),
         growth_high = apply(growth_sim_cohort, 1, quantile, probs = 0.975))

# Convert simulations to long format
sim_df_cohort <- as.data.frame(sim_ord_cohort)
names(sim_df_cohort) <- paste0("sim_", seq_len(ncol(sim_df_cohort)))
pred_sim_long_cohort <- pred_ord_cohort %>%
  select(COHORT, AGE) %>% bind_cols(sim_df_cohort) %>%
  pivot_longer(cols = starts_with("sim_"), names_to = "sim", values_to = "length")

# Estimate cohort slope and intercept for each simulation and age interval
growth_slopes_sim <- pred_sim_long_cohort %>%
  arrange(sim, COHORT, AGE) %>%
  group_by(sim, COHORT) %>%
  mutate(age_diff = AGE - lag(AGE),
         growth = ifelse(age_diff == 1, length - lag(length), NA_real_),
         age_interval = ifelse(age_diff == 1, paste0("Age ", AGE - 1, "–age ", 
                                                     AGE), NA_character_)) %>%
  ungroup() %>% filter(!is.na(growth), AGE > 0) %>%
  group_by(sim, age_interval) %>%
  summarise(intercept = coef(lm(growth ~ COHORT))[["(Intercept)"]],
            slope = coef(lm(growth ~ COHORT))[["COHORT"]],
            .groups = "drop")

# Summarise cohort slopes and intercepts across simulations
growth_slopes <- growth_slopes_sim %>%
  group_by(age_interval) %>%
  summarise(intercept = median(intercept),
            estimate = median(slope),
            low = quantile(slope, 0.025),
            high = quantile(slope, 0.975),
            .groups = "drop") %>%
  mutate(age_interval = factor(age_interval,
                               levels = paste0("Age ", start_age:(plus_group - 1), 
                                               "–age ", (start_age+1):plus_group)))

# Set factor order
growth_cohort_dat <- growth_cohort_dat %>%
  mutate(age_interval = factor(age_interval, 
                               levels = levels(growth_slopes$age_interval)))

# Save the observed data before adding missing cohort years
growth_cohort_observed <- growth_cohort_dat

# Add missing cohort years within each age interval so gaps are plotted as gaps
growth_cohort_dat <- growth_cohort_observed %>%
  group_by(age_interval) %>%
  complete(COHORT = seq(min(COHORT, na.rm = TRUE), 
                        max(COHORT, na.rm = TRUE), by = 1)) %>%
  ungroup() %>% arrange(age_interval, COHORT)

# Add a short extension after the final observed cohort for each age interval
# so each panel has its own ending segment
growth_extension_start <- growth_cohort_observed %>%
  filter(!is.na(growth)) %>%
  group_by(age_interval) %>%
  slice_max(COHORT, n = 1, with_ties = FALSE) %>%
  ungroup()
growth_extension_end <- growth_extension_start %>% mutate(COHORT = COHORT + 0.5)

# Combine the final observed point and the extension point for each panel
growth_extension <- bind_rows(growth_extension_start, growth_extension_end) %>%
  arrange(age_interval, COHORT)

# Generate slope labels
slope_labels <- growth_slopes %>%
  mutate(label = sprintf("Slope = %.2f cm/yr<br>95%% CI: %.2f to %.2f cm/yr",
                         estimate, low, high))

# Calculate endpoints for the simulation-derived linear trends
slope_lines <- growth_cohort_observed %>%
  filter(!is.na(growth)) %>%
  group_by(age_interval) %>%
  summarise(x = min(COHORT), xend = max(COHORT), .groups = "drop") %>%
  left_join(growth_slopes, by = "age_interval") %>%
  mutate(y = intercept + estimate * x, yend = intercept + estimate * xend)

# Only show zero line when zero is within the plotted uncertainty range
hline_dat <- growth_cohort_dat %>%
  group_by(age_interval) %>%
  summarise(ymin = min(growth_low, na.rm = TRUE),
            ymax = max(growth_high, na.rm = TRUE),
            .groups = "drop") %>%
  filter(ymin <= 0, ymax >= 0)

# Generate the plot
growth_plot <- ggplot(growth_cohort_dat,
                      aes(x = COHORT, y = growth)) +
  geom_ribbon(aes(ymin = growth_low,
                  ymax = growth_high), alpha = 0.3, na.rm = FALSE) +
  geom_line(linewidth = 2, na.rm = FALSE) +
  geom_ribbon(data = growth_extension,
              aes(x = COHORT, ymin = growth_low, ymax = growth_high),
              alpha = 0.3, inherit.aes = FALSE) +
  geom_line(data = growth_extension,
            aes(x = COHORT, y = growth),
            linewidth = 2, inherit.aes = FALSE) +
  geom_segment(data = slope_lines,
               aes(x = x, xend = xend, y = y, yend = yend),
               colour = "#3366FF", linewidth = 3, inherit.aes = FALSE) +
  geom_hline(data = hline_dat, aes(yintercept = 0),
             linetype = "dashed", linewidth = 0.5, inherit.aes = FALSE) +
  ggtext::geom_textbox(data = slope_labels, aes(x = Inf, y = Inf, label = label),
                       hjust = 1, vjust = 1, halign = 1, lineheight = 1.5,    
                       fill = NA, box.colour = NA, width = unit(0.75, "npc"),
                       inherit.aes = FALSE, size = font_size / 3) +
  facet_wrap(~ age_interval, scales = "free_y", ncol = 2) +
  labs(x = "Cohort", y = "Change in predicted length (cm/yr)") +
  theme_classic() +
  theme(text = element_text(size = font_size), 
        strip.background = element_blank(),
        strip.text = element_text(size = font_size, hjust = 1)) +
  coord_cartesian(expand = FALSE, xlim = c(NA, 2023)) +
  scale_x_continuous(breaks = as.numeric(round(with(growth_cohort_dat, 
         quantile(COHORT[!is.na(growth)], probs = c(0, 0.5, 1))))))
ggsave("Figs/growth_cohort_trends.jpeg", width = 9.5, height = 9.5, units = "in")

##### Plot spatial predictions (and fields) #####

# Setup a spatial prediction grid for a subset of the time series... 10-km res
pred_dat_space <- expand.grid(x = seq(min(survey_dat$x), max(survey_dat$x), by = 10),
                              y = seq(min(survey_dat$y), max(survey_dat$y), by = 10),
                              AGE = start_age:plus_group, YEAR = seq(1983, 2023, by = 3), 
                              DayOfYearScaled = 0) %>% 
  mutate(TAGE = AGE, f_tage = factor(TAGE),
         age_label = factor(AGE, levels = start_age:plus_group,
                            labels = age_labels_full))

# Make the predictions
predictions <- predict(mod_len, pred_dat_space, type = "response") %>% 
  st_as_sf(coords= c("x","y"), crs = 32620, remove = FALSE) %>%
  st_intersection(st_make_valid(survey.strata.cropped))

# Load in the depth contours, get them into utm coordinates and crop to the 
# study area
dc_crop <- st_read(
  "S:/Population_Ecology/Georges Bank/Useful R-scripts/Mapping Data/ctour100.shp") %>%
  st_transform(32620) %>% sf_m_to_km() %>%
  st_intersection(survey.strata.cropped)

# Plot the spatial fields by age, add depth contours
ggplot(data = predictions) +
  geom_raster(aes(x = x, y = y, fill = exp(epsilon_st + omega_s))) +
  geom_sf(data = land, fill = "darkgrey") +
  geom_sf(data = survey.strata.cropped, fill = NA, col = "black", size = 0.5) +
  geom_sf(data = dc_crop, size = 0.17) +
  spatial_theme + spatial_coord +
  scale_fill_gradient2("", midpoint = 1) +
  facet_wrap(~age_label, ncol = 4) +
  ggtitle("Spatial random effects (multiplicative scale)")
ggsave("Figs/spatial_fields.jpeg", width = 7.1, height = 4.8, units = "in")

# Plot the spatiotemporal fields for selected years, add depth contours
ggplot(data = predictions) +
  geom_raster(aes(x = x, y = y, fill = exp(est_rf - epsilon_st - omega_s))) +
  geom_sf(data = land, fill = "darkgrey") +
  geom_sf(data = survey.strata.cropped, fill = NA, col = "black", size = 0.5) +
  geom_sf(data = dc_crop, size = 0.17) +
  spatial_theme + spatial_coord +
  scale_fill_gradient2("", midpoint = 1) +
  facet_wrap(~YEAR, ncol = 5) +
  ggtitle("Spatiotemporal random effects (multiplicative scale)")
ggsave("Figs/spatiotemporal_fields.jpeg", width = 7.5, height = 6, units = "in")

# Plot the plus-group predictions for selected years, add depth contours
ggplot(data = predictions %>% filter(TAGE == plus_group)) +
  geom_raster(aes(x = x, y = y, fill = est)) +
  geom_sf(data = land, fill = "darkgrey") +
  geom_sf(data = survey.strata.cropped, fill = NA, col = "black", size = 0.5) +
  geom_sf(data = dc_crop, size = 0.17) +
  spatial_theme + spatial_coord +
  scale_fill_viridis_c("") +
  facet_wrap(~YEAR, ncol = 5) +
  ggtitle("Predicted length (cm) of age 8+ fish")
ggsave("Figs/plus_group_predictions.jpeg", width = 7.5, height = 6, units = "in")