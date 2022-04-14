library(drc)
library(rstatix)
library(tidyverse)
library(readxl)
library(ggpubr)
library(ggrepel)
library(scales)
library(ggprism)


# this is now a project, so all directory structure relative to home project directory
data_sheet_fn <- file.path("data", "2022-04-13 SFLT.xlsx") # 2022-03-11 sflt elisa.xlsx
output_dir_proj <- file.path("output", 
                        str_extract(string = basename(data_sheet_fn), 
                                    pattern = "[0-9]+-[0-9]+-[0-9]+"))
dir.create(output_dir_proj, recursive = T, showWarnings = F)


run_elisa_analysis <- function(data_sheet, output_dir) {
  #' *DEBUG*
  #' data_sheet = data_sheet_fn; output_dir = output_dir_proj
  standard_dat <- read_excel(data_sheet, sheet = 1)  %>%
    # assuming that Concentration of Standard is the 1st column
    mutate(log2concentration = log2(.[[1]]),
           Absorbance = .[[2]])
  results <- read_excel(data_sheet, sheet = 2) %>%
    mutate(Samples = .[[1]],
           Absorbance = .[[2]]) %>%
    mutate(Samples = make.unique(Samples,sep = "__"))
  
  #' *INFORMATION ABOUT LOGISTIC REGRESSION MODEL FOR ELISAS*
  # logistic regression model; estimate slope, ED50, and upper and lower bounds using non-linear method
  # logic behind this model: where we model the ED50 in log space,  four-parameter log-logistic model 
  # may be preferred for dose-response analysis involving very small datasets (<15–20) where a normally
  # distributed parameter estimate may more reasonably be assumed on a logarithm-transformed dose scale 
  # than the original dose scale.
  # drc paper: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0146021#sec001
  
  fit <- drm(formula = Absorbance ~ log2concentration, 
             data = standard_dat, 
             fct=LL.4(names=c("Slope", "Lower", "Upper", "ED50")),
             logDose = 2 # log 2 doses
  )
  
  # print the fit summary
  print(summary(fit))
  message("\n")
  
  # this is the adjusted absorbances of each of the samples; we want to predict their concentration based off the standard curve
  response <- results$Absorbance
  
  # log2 estimate along with associated error 
  DOSEx_temp <- ED(fit, response, type="absolute", display=F) %>% as_tibble()
  
  DOSEx <- DOSEx_temp %>% as.data.frame() %>%
    mutate(up = Estimate + `Std. Error`,
           down = Estimate - `Std. Error`) %>%
    mutate(across(.cols = everything(), log2))
  
  
  
  combined_res <- bind_cols(results, DOSEx_temp, 
                            DOSEx %>% rename_with(.fn = ~ str_c(.x, "_log2"), .cols = everything())) %>%
    dplyr::rename(`fit_to_curve_concentration_estimate` = Estimate,
                  `log2_fit_to_curve_concentration_estimate` = Estimate_log2)
  
  print(combined_res %>% select(1:3), n = Inf)
  message("\n")
  
  print(DOSEx)
  message("\n")
  
  # generate fit vals for plotting
  fit_vals <- expand.grid(conc=exp(seq(log2(0.5), log2(100), length=100)))
  pm <- predict(fit, newdata=fit_vals, interval = "confidence")
  fit_vals$p <- pm[,1]
  fit_vals$pmin <- pm[,2]
  fit_vals$pmax <- pm[,3]
  
  min_x <- min(combined_res$log2_fit_to_curve_concentration_estimate, na.rm = TRUE) - 1
  max_x <- max(combined_res$log2_fit_to_curve_concentration_estimate, na.rm = TRUE) + 1
  min_y <- min(combined_res$Absorbance, na.rm = TRUE) - 0.1
  max_y <- max(combined_res$Absorbance, na.rm = TRUE) + 0.5

   #theme info
  rel_size <- 1
  my_theme <- theme_prism(base_size = 7) +
    theme(strip.text.x = element_text(size = rel(rel_size*3)),
          title = element_text(size = rel(rel_size*3)),
          plot.subtitle = element_text(size = rel(rel_size*0.75)),
          plot.caption = element_text(size = rel(rel_size*0.4)),
          legend.box.spacing = unit(1, "cm"),
          legend.text = element_text(size = rel(rel_size*1.5)),
          axis.text.y = element_text(size = rel(rel_size*3), angle = 0, vjust = 0.2),
          axis.text.x = element_text(size = rel(rel_size*3), angle = 45),
          panel.grid.major = element_line(color = "gray",
                                    size = 0.6,
                                    linetype = 2),
          panel.grid.minor = element_line(color = "gray",
                                          size = 0.3,
                                          linetype = 2),
          panel.spacing = unit(1, "lines")) 
  # axes info
  nmajor <- ceiling(max_y); nminor <- ceiling(max_y*4)

  gplot_res <- ggplot(standard_dat, aes(x = log2concentration, y = Absorbance)) +
    
    geom_ribbon(data=fit_vals, aes(x=conc, y=p, ymin=pmin, ymax=pmax), alpha=0.3) +
    geom_line(data=fit_vals, aes(x=conc, y=p), color = "dodgerblue", lwd = 1.25) +
    
    geom_point(data = combined_res, 
               mapping = aes(x = log2_fit_to_curve_concentration_estimate, y = Absorbance),
               size = rel_size*2.25) +
    # scale_shape_manual(values=1:nlevels(factor(combined_res$Samples))) +
    geom_label_repel(data =combined_res,
                     mapping = aes(x = log2_fit_to_curve_concentration_estimate, y = Absorbance,
                                   label = Samples),
                     max.overlaps = Inf, force = 50, nudge_x = -0.5, nudge_y = 0.25, min.segment.length = 0,
                     size = rel_size*5) +
    geom_errorbar(data = combined_res, 
                  mapping = aes(x = log2_fit_to_curve_concentration_estimate, 
                                y = Absorbance, 
                                xmin = down_log2, xmax = up_log2),
                  size = rel_size) +

    coord_trans(x="log2") +
    # scale_y_continuous(breaks = extended_breaks(n = nmajor),
    #                    minor_breaks = extended_breaks(n = nmajor * nminor)) +
    # scale_x_continuous(breaks = seq(floor(min_x), floor(max_x), 1)) +
    xlim(min_x, max_x) +
    ylim(min_y, max_y) +
    
    labs(x = "Log2 Concentration",
         y = "Absorbance",
         subtitle = "4-paramter Log-logisitic Regression",
         caption = "Each point and their estimated log2 concentration with upper and lower bounds (black horizontal bars)\nBlue line is the line of best fit by the 4-param Log-logisitic model\nGray shaded ribbon is bounds of model predicitions on simulated data") + 
    my_theme +
    ggtitle("ELISA Analysis"); gplot_res
    
  ggsave(plot = gplot_res, filename = file.path(output_dir, "results_plot_ggplot.png"))
  
  plot_fn <- file.path(output_dir, "results_plot.pdf")
  # save the plot
  #' reading in global variables, not very good practice..
  save_plot <- function() {
    pdf(file = plot_fn, width = 6, height = 6)
    plot(fit, 
         xlim = c(min(DOSEx[,1], na.rm = T), max(DOSEx[,1], na.rm = T)),
         ylim = c(min(response, na.rm = T), max(response, na.rm = T)),
         main = "Plot of absorbance x concentration\nFitted sample values overlayed") +
      points(y=response,x=DOSEx[,1],col="blue",pch=19,cex=0.75) +
      # With error bars
      arrows(DOSEx[,1], response, DOSEx$up, response, length=0.05, angle=90, lwd=1,col="blue") +
      arrows(DOSEx[,1], response, DOSEx$down, response, length=0.05, angle=90, lwd=1,col="blue")
    dev.off()
  }
  save_plot()
  
  message("writing to output directory:")
  message(output_dir)
  write_tsv(combined_res %>% select(Samples, Absorbance, fit_to_curve_concentration_estimate), 
            file = file.path(output_dir, "results.tsv"))
  
  message("\nFinal result:")
  print(combined_res %>% select(Samples, Absorbance, fit_to_curve_concentration_estimate), n = Inf)
  
  # options(warn=0)
}


run_elisa_analysis(data_sheet = data_sheet_fn, 
                   output_dir = output_dir_proj)
