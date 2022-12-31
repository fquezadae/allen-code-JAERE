library(ggplot2)
library(ggpubr)

# read data with no private information
par_res <- readRDS(paste0(dirout, "MC_", "0", ".rds"))
true <- lapply(par_res, `[[`, 1)

binwidth <- 0.03

# grab coefficeients
truecoef <- data.frame(Alpha =
  unlist(lapply(lapply(true, `[[`, 2), `[`, 1, 1)))

trueplot <- ggplot(data = truecoef) +
  geom_histogram(aes(x = Alpha, y = binwidth * ..density..,
    fill = "No catch error"), binwidth = binwidth, color = "black") +
  labs(y = "Probability", fill = "Methodology") +
  theme_bw() +
  theme(legend.position = c(.8, .8),
  legend.title = element_text(size = 7),
  legend.text = element_text(size = 7),
  legend.key.size = unit(0.2, "cm")) +
  scale_fill_grey(start = 0.2, end = 0.8)

# read data when dev=3
par_res <- readRDS(paste0(dirout, "MC_", "3", ".rds"))

correction <- lapply(par_res, `[[`, 1)
uncorrected <- lapply(par_res, `[[`, 2)
twostage <- lapply(par_res, `[[`, 3)
true <- lapply(par_res, `[[`, 4)

# make data with uncorrected and true marginal utility of catch
plotdata <- rbind(data.frame(Alpha =
  unlist(lapply(lapply(uncorrected, `[[`, 2), `[`, 1, 1)),
  Method = "Uncorrected"),
  data.frame(Alpha =
  unlist(lapply(lapply(true, `[[`, 2), `[`, 1, 1)),
  Method = "True catch parameters"))

binwidth <- 0.2

# note plots two subsets by "Method"
uncorplot <- ggplot(data = plotdata) +
    geom_histogram(aes(x = Alpha, y = binwidth * ..density..,
        fill = Method), binwidth = binwidth, color = "black",
        position = "identity", alpha = 0.75) +
    labs(y = "Probability", fill = "Methodology") +
    theme_bw() +
    theme(legend.position = c(.8, .8),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.2, "cm")) +
    scale_fill_grey(start = 0.2, end = 0.8)

figure <- ggarrange(trueplot, uncorplot,
    labels = c("A", "B"),
    ncol = 2, nrow = 1)

ggsave(paste0(dirout, "Figure_1.eps"), figure,
    width = 7.5, height = 4.21875, dpi = 800, device = cairo_ps)

# same but pull from marginal utility of catch from two-stage and true
plotdata <- rbind(data.frame(Alpha =
    unlist(lapply(lapply(twostage, `[[`, 2), `[`, 1, 1)),
    Method = "Two-stage"),
    data.frame(Alpha =
    unlist(lapply(lapply(true, `[[`, 2), `[`, 1, 1)),
    Method = "True catch parameters"))

binwidth <- 0.25

twostageplot <- ggplot(data = plotdata) +
    geom_histogram(aes(x = Alpha, y = binwidth * ..density..,
        fill = Method), binwidth = binwidth, color = "black",
        position = "identity", alpha = 0.75) +
    labs(y = "Probability", fill = "Methodology") +
    theme_bw() +
    theme(legend.position = c(.8, .8),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.2, "cm")) +
    scale_fill_grey(start = 0.2, end = 0.8)

# finally with full information
plotdata <- rbind(data.frame(Alpha =
    unlist(lapply(lapply(correction, `[[`, 2), `[`, 1, 1)),
    Method = "FIML"),
    data.frame(Alpha =
    unlist(lapply(lapply(true, `[[`, 2), `[`, 1, 1)),
    Method = "True catch parameters"))

binwidth <- 0.175

corplot <- ggplot(data = plotdata) +
    geom_histogram(aes(x = Alpha, y = binwidth * ..density..,
        fill = Method), binwidth = binwidth, color = "black",
        position = "identity", alpha = 0.75) +
    labs(y = "Probability", fill = "Methodology") +
    theme_bw() +
    theme(legend.position = c(.8, .8),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.2, "cm")) +
    scale_fill_grey(start = 0.2, end = 0.8)

figure2 <- ggarrange(twostageplot, corplot,
    labels = c("A", "B"),
    ncol = 2, nrow = 1)

ggsave(paste0(dirout, "Figure_2.eps"), figure2,
    width = 7.5, height = 4.21875, dpi = 800, device = cairo_ps)

# plot bias with different values of catch deviation
biasplot <- list()
for (it in c(1:7)) {

dev <- c(1, 1.5, 2, 2.5, 3, 3.5, 4)

i <- dev[it]

par_res <- readRDS(paste0(dirout, "MC_", i, ".rds"))

correction <- lapply(par_res, `[[`, 1)

uncorrected <- lapply(par_res, `[[`, 2)

twostage <- lapply(par_res, `[[`, 3)

true <- lapply(par_res, `[[`, 4)

biasplot[[it]] <- cbind(rbind(data.frame(
    Alpha = median(unlist(lapply(lapply(true, `[[`, 2), `[`, 1, 1))),
    se = median(unlist(lapply(lapply(true, `[[`, 2), `[`, 1, 2))),
    Methodology = "True"),
    data.frame(
    Alpha = median(unlist(lapply(lapply(uncorrected, `[[`, 2), `[`, 1, 1))),
    se = median(unlist(lapply(lapply(uncorrected, `[[`, 2), `[`, 1, 2))),
    Methodology = "Uncorrected"),
    data.frame(
    Alpha = median(unlist(lapply(lapply(twostage, `[[`, 2), `[`, 1, 1))),
    se = median(unlist(lapply(lapply(twostage, `[[`, 2), `[`, 1, 2))),
    Methodology = "Two-stage"),
    data.frame(
    Alpha = median(unlist(lapply(lapply(correction, `[[`, 2), `[`, 1, 1))),
    se = median(unlist(lapply(lapply(correction, `[[`, 2), `[`, 1, 2))),
    Methodology = "FIML")),
    data.frame(cdev = rep(i, 4)))

}

forbiasplot <- do.call(rbind, biasplot)

plotout <- ggplot(data = forbiasplot, aes(x = cdev, y = Alpha)) +
    geom_point(aes(y = Alpha, x = cdev,
        shape = Methodology, color = Methodology), fill = "white", size = 3) +
    geom_line(aes(y = Alpha, x = cdev,
        color = Methodology), size = 1) +
    geom_errorbar(aes(ymin = Alpha - (se), ymax = Alpha + (se),
      color = Methodology), width = 0.05, size = 0.5) +
    labs(y = "Alpha", x = "Catch error standard deviation") +
    scale_shape_manual(values = c(22, 21, 25, 24)) +
    theme_bw() +
    theme(legend.position = c(.15, .7)) +
    scale_color_grey(start = 0.2, end = 0.8)

ggsave(filename = paste0(dirout, "Figure_3.eps"), plotout,
    width = 7.5, height = 4.21875, dpi = 800, device = cairo_ps)