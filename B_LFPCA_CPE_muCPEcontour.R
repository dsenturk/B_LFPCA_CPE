muCPEcontour <- function(post_samples,                            # Outputs from FDpostSumms_BLFPCA_postSumms.R
                         alpha_cont = seq(0.05, 0.95, by = 0.10), # Vector of alpha values for the CPE contours, a vector of real numbers
                         s,                                       # Vector of longitudinal time points to be plotted, a vector of real numbers
                         t                                        # Vector of functional time points to be plotted, a vector of real numbers
){
  #############################################################################
  ## Description: Function for generating the CPE contour plots of the mean function 
  ##              posterior estimates calculated for a given simulation run displayed in 
  ##              Figures S10 of "Central posterior envelopes for Bayesian longitudinal 
  ##              functional principal component analysis" by Boland et al. (2024).  
  ##              Users can select to display at different alpha-levels using the arguments 
  ##              alpha_cont.
  ## Args:        (see above)
  ## Returns:     Plots of CPE contours for the mean function
  ## muCPEcontour Outline:
  ##              1. Clean and format data
  ##              3. Obtain legend
  ##              3. Plot the visualization
  #############################################################################
  
  # Install missing packages
  list.of.packages <- c("tidyverse", "reshape2", "viridis", "cowplot",
                        "latex2exp")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if (length(new.packages)) install.packages(new.packages) 
  
  # Load packages
  library(tidyverse)
  library(reshape2)
  library(viridis)
  library(cowplot)
  library(latex2exp)
  
  #############################################################################
  # 1. Clean and format data
  #############################################################################
  
  M <- 6000
  
  # Obtain data needed for plotting 
  data <- lapply(1:length(s), function(X){
    post_samples$mu_sample[which(rep(s, each = length(t)) == s[X]), ]
  })
  ranking <- post_samples$ranks$mu_MVD
  Truth <- lapply(1:length(s), function(X){
    data.frame(Time = t, est = truth_funcs$mu[which(rep(s, each = length(t)) == s[X])])
  })
  
  # Order the posterior sample based on MBD or MVD
  ranking <- ranking[order(-ranking$rank),]
  ordering_index <- ranking$id
  data_ordered <- lapply(1:length(s), function(X){
    as.matrix(t(data[[X]])[ordering_index, ])
  })
  
  # Obtain point estimate data
  point_est <- lapply(1:length(s), function(X){
    data.frame(Time = t, est = data_ordered[[X]][1, ]) 
  })
  
  # Calculate the CPE regions 
  total_alpha <- length(alpha_cont)
  CPE_regions <- lapply(1:length(s), function(S){
    lapply(X = 1:total_alpha, function(X){
      
      # Calculate upper bounds
      upper <- apply(data_ordered[[S]][1:(floor(alpha_cont[X] * M)),], 2, max)
      
      # Calculate lower bounds
      lower <- apply(data_ordered[[S]][1:(floor(alpha_cont[X] * M)),], 2, min)
      
      return(list(upper = upper, lower = lower))
    })
  })

  # Turn the data from wide to long format
  data <- lapply(1:length(s), function(X){
    subset <- data.frame(data[[X]])
    subset$Time <- t
    subset %>%
      melt(id.vars = c("Time"))
  })
  
  # Generate datasets to plot CPE contours 
  contour_min <- lapply(1:length(s), function(X){
    data.frame(cbind(Time = t, upper = CPE_regions[[X]][[1]]$upper,
                     lower = CPE_regions[[X]][[1]]$lower))
  })
 
  contour_bands <- lapply(1:length(s), function(S){
    lapply(2:total_alpha, function(X){
      contour_upper <- data.frame(cbind(Time = t, upper = CPE_regions[[S]][[X]]$upper,
                                        lower = CPE_regions[[S]][[X - 1]]$upper))
      
      contour_lower <- data.frame(cbind(Time = t, upper = CPE_regions[[S]][[X - 1]]$lower,
                                        lower = CPE_regions[[S]][[X]]$lower))
      
      return(list(contour_upper = contour_upper, contour_lower = contour_lower))
    })
  })

  
  #############################################################################
  # 2. Obtain legend
  #############################################################################
  
  # Colors used in contour plot
  cols <- viridis::viridis(n = total_alpha)
  
  # Generate a dataset and plot to obtain legend
  df <- data.frame(x = rnorm(1000, sd = 100), y = rnorm(1000, sd = 100))
  grid <- ggplot(df, aes(x = x, y = y)) +
    geom_density_2d_filled(bins = 10) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title=element_text(size=12),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 13),
          legend.background = element_rect(colour = 'black', fill = 'white', 
                                           linetype='solid'),
          axis.text=element_text(size=12)) +
    scale_fill_manual(values =  cols, 
                      labels = c(as.character(alpha_cont))) +
    guides(fill=guide_legend(title=TeX(r'($\alpha$ level)')))
  legend <- get_legend(grid)
  
  #############################################################################
  # 3. Plot the visualization
  #############################################################################
  
  plot_ids <- c("(a)", "(b)", "(c)", "(d)", 
                "(e)", "(f)", "(g)", "(h)",
                "(i)", "(j)", "(k)", "(l)",
                "(m)", "(n)", "(o)", "(p)",
                "(q)", "(r)", "(s)", "(t)")
  
  contour_plots <- lapply(1:length(s), function(X){
    contour_plot <- ggplot() +
      geom_line(data = data[[X]], 
                mapping = aes(x = Time, y = value, color = variable), alpha = 0.3) +
      scale_color_manual(values = c(rep("gray70", M))) +
      scale_x_continuous(expand = c(0, 0)) +
      theme_bw() + 
      labs(title = paste0(plot_ids[X], " MVE-CPE (s = ", round(s[X], 3), ")"), y = "Mean Function", x = "Time (t)") +
      theme(legend.position = "none", 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title = element_text(size = 12)) +
      geom_ribbon(data = contour_min[[X]], mapping = aes(x = Time, ymin = lower,
                                                    ymax = upper),
                  fill = cols[1], alpha = 0.8)
    
    for(i in 2:total_alpha){
      contour_plot <- contour_plot + 
        geom_ribbon(data = contour_bands[[X]][[i - 1]]$contour_upper, 
                    mapping = aes(x = Time, ymin = lower, ymax = upper),
                    fill = cols[i], alpha = 0.8) +
        geom_ribbon(data = contour_bands[[X]][[i - 1]]$contour_lower, 
                    mapping = aes(x = Time, ymin = lower, ymax = upper),
                    fill = cols[i], alpha = 0.8) 
    }
    
    contour_plot <- contour_plot + 
      geom_line(data = Truth[[X]], mapping = aes(x = Time, y = est), color = "red",
                size = 0.9) +
      geom_line(data = point_est[[X]], mapping = aes(x = Time, y = est), color = "black",
                size = 0.9) 
    
    
    return(contour_plot)
  })
  
  ggdraw() +
    draw_plot(contour_plots[[1]], x = 0, y = 1 - .2, width = 0.225 - 0.01, height = 0.2) +
    draw_plot(contour_plots[[2]], x = 0.225, y = 1 - .2, width = 0.225 - 0.01, height = 0.2) + 
    draw_plot(contour_plots[[3]], x = 0.45, y = 1 - .2, width = 0.225 - 0.01, height = 0.2) + 
    draw_plot(contour_plots[[4]], x = 0.675, y = 1 - .2, width = 0.225 - 0.01, height = 0.2) + 
    draw_plot(contour_plots[[5]], x = 0, y = 1 - .4, width = 0.225 - 0.01, height = 0.2) +
    draw_plot(contour_plots[[6]], x = 0.225, y = 1 - .4, width = 0.225 - 0.01, height = 0.2) + 
    draw_plot(contour_plots[[7]], x = 0.45, y = 1 - .4, width = 0.225 - 0.01, height = 0.2) + 
    draw_plot(contour_plots[[8]], x = 0.675, y = 1 - .4, width = 0.225 - 0.01, height = 0.2) + 
    draw_plot(contour_plots[[9]], x = 0, y = 1 - .6, width = 0.225 - 0.01, height = 0.2) +
    draw_plot(contour_plots[[10]], x = 0.225, y = 1 - .6, width = 0.225 - 0.01, height = 0.2) + 
    draw_plot(contour_plots[[11]], x = 0.45, y = 1 - .6, width = 0.225 - 0.01, height = 0.2) +
    draw_plot(contour_plots[[12]], x = 0.675, y = 1 - .6, width = 0.225 - 0.01, height = 0.2) +
    draw_plot(contour_plots[[13]], x = 0, y = 1 - .8, width = 0.225 - 0.01, height = 0.2) +
    draw_plot(contour_plots[[14]], x = 0.225, y = 1 - .8, width = 0.225 - 0.01, height = 0.2) +
    draw_plot(contour_plots[[15]], x = 0.45, y = 1 - .8, width = 0.225 - 0.01, height = 0.2) +
    draw_plot(contour_plots[[16]], x = 0.675, y = 1 - .8, width = 0.225 - 0.01, height = 0.2) +
    draw_plot(contour_plots[[17]], x = 0, y = 0, width = 0.225 - 0.01, height = 0.2) +
    draw_plot(contour_plots[[18]], x = 0.225, y = 0, width = 0.225 - 0.01, height = 0.2) +
    draw_plot(contour_plots[[19]], x = 0.45, y = 0, width = 0.225 - 0.01, height = 0.2) +
    draw_plot(contour_plots[[20]], x = 0.675, y = 0, width = 0.225 - 0.01, height = 0.2) +
    draw_grob(legend, x = 0.9, y = 0.2, width = 0.1 - 0.01, height = .6) 

}

