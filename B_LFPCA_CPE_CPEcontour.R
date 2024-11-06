CPEcontour <- function(post_samples,            # Outputs from FDpostSumms_BLFPCA_postSumms.R
                       type,                    # Type of eigenfunctions, chosen from "long_eigen" and "func_eigen".
                       jork,                    # Index for the eigenfunctions taking values 1, 2, 3
                       CPE,                     # Type of CPE calculated and plotted, chosen from "MBD", "Marginal MVD", and "Kernel MVD"
                       alpha_cont = seq(0.05, 0.95, by = 0.10), # alpha values for the CPE contours, a vector of real numbers
                       sort,                    # To denote the longitudinal (s) or functional (t) dimension 
                       title,                   # Title given to the plot, charactors
                       xlab                     # X-axis label, characters
){
  #############################################################################
  ## Description: Function for generating the CPE contour plots of the functional 
  ##              posterior estimates calculated for a given simulation run displayed in 
  ##              "Central posterior envelopes for Bayesian longitudinal functional principal 
  ##              component analysis" by Boland et al. (2024). 
  ##              Users can select which CPE type (MBD, MVD or kernel MVD) and 
  ##              respective  MBD/MVD/kernel MVD median (solid black line) and CPEs (shaded areas) 
  ##              they wish to display at different alpha-levels using the arguments 
  ##              CPE and alpha_cont, respectively. 
  ## Args:        (see above)
  ## Returns:     CPE contour plot for longitudinal or functional eigenfunctions
  ## CPEcontour Outline:
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
  
  if(jork == 1 & type == "long_eigen"){
    ylab <- TeX(r"( $\psi_1(s)$ )")
  } else if (jork == 2 & type == "long_eigen"){
    ylab <-  TeX(r"( $\psi_2(s)$ )")
  } else if (jork == 1 & type == "func_eigen"){
    ylab <- TeX(r"( $\phi_1(t)$ )")
  } else if (jork == 2 & type == "func_eigen"){
    ylab <- TeX(r"( $\phi_2(t)$ )")
  } else if (jork == 3 & type == "long_eigen"){
    ylab <-  TeX(r"( $\psi_3(s)$ )")
  } else if (jork == 3 & type == "func_eigen"){
    ylab <- TeX(r"( $\phi_3(t)$ )")
  }
  
  
  # Obtain data needed for plotting 
  if(type == "long_eigen"){
    data <- post_samples$Psi_long[, jork, ]
    Truth <- data.frame(Time = sort, est = truth_funcs$psi[[jork]]) 
    if(CPE == "MBD"){
      ranking <- post_samples$ranks$MBD.psi_j[[jork]]
    } else if (CPE == "Marginal MVD"){
      ranking <- post_samples$ranks$MVD.marginal_cov_long
    } else if (CPE == "Kernel MVD"){
      ranking <- post_samples$ranks$MVD.cov_kernel
    }
  } else if (type == "func_eigen"){
    data <- post_samples$Phi_func[, jork, ]
    Truth <- data.frame(Time = sort, est = truth_funcs$phi[[jork]]) 
    if(CPE == "MBD"){
      ranking <- post_samples$ranks$MBD.phi_k[[jork]]
    } else if (CPE == "Marginal MVD"){
      ranking <- post_samples$ranks$MVD.marginal_cov_func
    } else if (CPE == "Kernel MVD"){
      ranking <- post_samples$ranks$MVD.cov_kernel
    }
  }
  
  # Order the posterior sample based on MBD or MVD
  ranking <- ranking[order(-ranking$rank),]
  ordering_index <- ranking$id
  data_ordered <- as.matrix(t(data)[ordering_index, ])
  
  # Obtain point estimate data
  point_est <- data.frame(Time = sort, est = data_ordered[1, ]) 
  
  # Calculate the CPE regions 
  total_alpha <- length(alpha_cont)
  CPE_regions <- lapply(X = 1:total_alpha, function(X){
    
    # Calculate upper bounds
    upper <- apply(data_ordered[1:(floor(alpha_cont[X] * nrow(data_ordered))),], 2, max)
    
    # Calculate lower bounds
    lower <- apply(data_ordered[1:(floor(alpha_cont[X] * nrow(data_ordered))),], 2, min)
    
    return(list(upper = upper, lower = lower))
  })
  
  # Turn the data from wide to long format
  data <- data.frame(data)
  data$Time <- sort
  data <- data %>%
    melt(id.vars = c("Time"))
  
  # Generate datasets to plot CPE contours 
  contour_min <- data.frame(cbind(Time = sort, upper = CPE_regions[[1]]$upper,
                                  lower = CPE_regions[[1]]$lower))
  
  contour_bands <- lapply(2:total_alpha, function(X){
    contour_upper <- data.frame(cbind(Time = sort, upper = CPE_regions[[X]]$upper,
                                      lower = CPE_regions[[X - 1]]$upper))
    
    contour_lower <- data.frame(cbind(Time = sort, upper = CPE_regions[[X - 1]]$lower,
                                      lower = CPE_regions[[X]]$lower))
    
    return(list(contour_upper = contour_upper, contour_lower = contour_lower))
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
          axis.title=element_text(size=10),
          legend.text = element_text(size = 9),
          legend.title = element_text(size = 11),
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
  
  contour_plot <- ggplot() +
    geom_line(data = data, 
              mapping = aes(x = Time, y = value, color = variable), alpha = 0.3) +
    scale_color_manual(values = c(rep("gray70", nrow(data_ordered)))) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_bw() + 
    labs(title = title, y = ylab, x = xlab) +
    theme(legend.position = "none", 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 12)) +
    geom_ribbon(data = contour_min, mapping = aes(x = Time, ymin = lower,
                                                  ymax = upper),
                fill = cols[1], alpha = 0.8)
  
  for(i in 2:total_alpha){
    contour_plot <- contour_plot + 
      geom_ribbon(data = contour_bands[[i - 1]]$contour_upper, 
                  mapping = aes(x = Time, ymin = lower, ymax = upper),
                  fill = cols[i], alpha = 0.8) +
      geom_ribbon(data = contour_bands[[i - 1]]$contour_lower, 
                  mapping = aes(x = Time, ymin = lower, ymax = upper),
                  fill = cols[i], alpha = 0.8) 
  }
  
  contour_plot <- contour_plot + 
    geom_line(data = Truth, mapping = aes(x = Time, y = est), color = "black",
              size = 0.9) 
  
  return(contour_plot)
}
