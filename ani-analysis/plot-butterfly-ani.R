library(vroom); library(dplyr); library(tidyr); library(ggplot2); library(patchwork); library(diptest)
library(ggpubr); library(cowplot); library(scales); library(latex2exp); library(purrr)  # theme_minimal_hgrid/vgrid, percent, TeX, map_dfr

SYM     <- "min"
BULK    <- "(0.1,0.15]"
LR_MIN  <- 6.635
# "#9D4030"
mc = c("#D0D55C", "#809D6F", "#BD4030", "#C3A97E", "#769DA6", "#A69E33", "#734002", "#595622", "#8C873F",  "#474B71")
BINS    <- c("[0,0.05)", "[0,0.1)", "(0.0,0.15]")
mpwd <- vroom("merged-distances.tsv")
gdiff_cfg <- "k29-l333-n10000-h3.gdiff"
gdiff_cfg <- "k25-l1000-n10000-h3.gdiff"
mpwd %>%
  mutate(cov=((query_coverage+subject_coverage)/2)) %>%
  filter(interaction(config, method) %in% c("k21-s10000.mash", "default.skani", gdiff_cfg, "default.FastANI", "v4.dashing2")) %>%
  ggplot() +
  aes(x=ANIm, y=estimate, color=method) +
  # facet_wrap(~cov, scale="free") +
  geom_point()+
  stat_cor(aes(label = after_stat(r.label)), show.legend = F) +
  geom_text(
    data = function(d) d %>% 
      mutate(method=factor(method, levels = c("dashing2", "FastANI", "gdiff", "mash", "skani"))) %>%
      group_by(method) %>% 
      summarize(
        mape = mean(abs(estimate-ANIm) / (ANIm)) * 100, 
        .groups = "drop"
      ) %>% 
      mutate(
        x = 0.04,
        y = 0.154 - ((row_number())* 0.00675)
      ),
    aes(x=x, y=y, label=paste0("MAPE = ", round(mape, 2), "%")),
    hjust=0, vjust=0, show.legend = F, 
  ) +
  # stat_smooth(method = "lm", se=F) +
  geom_abline() +
  theme_half_open() +
  # theme(panel.grid.major = element_line(color = "lightgray")) +
  scale_color_manual(values=mc) +
  labs(color="Method", label="", y=TeX(r'($\hat{D}$)'))

mmp <- merge(mpwd %>% filter(method=="gdiff") %>% select(genome_a, genome_b, estimate) %>% rename(gdiff=estimate), mpwd)
mmp %>% filter(method != "gdiff" & method != "ANIm") %>%
  ggplot() +
  facet_wrap(~method) +
  aes(x=gdiff, y=estimate, color=method) +
  geom_point(size=0.5, alpha=0.75) +
  geom_abline() +
  theme_minimal() +
  scale_color_manual(values=c("#D0D55C", "#809D6F", "#C3A97E", "#769DA6")) +
  labs(y="Alternative", color="Method") +
  coord_cartesian(xlim = c(0, 0.15), ylim=c(0, 0.15))
ggsave2("./butterfly-comparison.pdf", width = 6.5, height = 5.25)

mpwd %>%
  mutate(cov=((query_coverage+subject_coverage)/2)) %>%
  filter(interaction(config, method) %in% c("k21-s10000.mash", "default.skani", gdiff_cfg, "default.FastANI", "v4.dashing2")) %>%
  ggplot() +
  aes(x=cov, y=abs(estimate-ANIm)/ANIm, color=method) +
  # facet_wrap(~cov, scale="free") +
  geom_point()+
  # stat_summary_bin(breaks = c(0, 0.33, 0.67, 1.0), alpha=0.85) +
  scale_y_continuous(labels=percent) +
  scale_color_manual(values=mc) +
  theme_half_open() +
  theme(panel.grid.major = element_line(color = "lightgray")) +
  labs(color="Method", label="", x="Alignment coverage", y="MAPE")

dfr <- data.frame(
  method = c("dashing2", "FastANI", "gdiff", "mash", "skani"),
  wall_seconds = c(8.02, 931521, 40.36, 15.33, 25.84),
  peak_kb = c(192584, 186383012, 14461864, 414652, 17561772)
)
r <-merge(dfr, mpwd %>%
        mutate(cov=((query_coverage+subject_coverage)/2)) %>%
        filter(interaction(config, method) %in% c("k21-s10000.mash", "default.skani", "k27-delta4.gdiff", "default.FastANI", "v4.dashing2")) %>%
        group_by(method) %>%
        summarise(mape=mean(abs(estimate-ANIm)/ANIm))
)
r %>% ggplot() +
  aes(x=peak_kb/1000/1000, y=mape, color=method) +
  geom_point(size = 5) +
  scale_x_continuous(transform = "log10") +
  scale_y_continuous(labels = percent) +
  scale_color_manual(values=mc) +
  theme_minimal_grid() +
  labs(color="Method", x="Peak memory (GB)", y="MAPE")

r %>% ggplot() +
  aes(x=wall_seconds, y=mape, color=method) +
  geom_point(size = 5) +
  scale_x_continuous(transform = "log10") +
  scale_y_continuous(labels = percent) +
  scale_color_manual(values=mc) +
  theme_minimal_grid() +
  labs(color="Method", x="Running time (sec)", y="MAPE")

# Distributions
dfsk <- vroom("samples/samples-gdiff-k25-l1000-n10000-h3.tsv.gz", comment = "#") %>% 
  mutate(genome_a=gsub(x=genome_a, pattern="\\..*", replacement="")) %>%
  mutate(genome_b=gsub(x=genome_b, pattern="\\..*", replacement=""))
reconcile_directional_c <- function(df, group_cols = c("config", "genome_a", "genome_b"), direction_col = "direction") {
  missing_cols <- setdiff(c(group_cols, direction_col, "d"), names(df))
  if (length(missing_cols) > 0L) {
    stop("reconcile_directional(): missing column(s): ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  
  if (direction_col != "direction") {
    df <- dplyr::rename(df, direction = dplyr::all_of(direction_col))
  }
  
  placeable <- !is.na(df$direction) & df$direction %in% c("ab", "ba")
  if (any(!placeable)) {
    warning(sprintf(
      "reconcile_directional(): dropping %d row(s) whose `%s` is not \"ab\"/\"ba\"",
      sum(!placeable), direction_col), call. = FALSE)
    df <- df[placeable, , drop = FALSE]
  }
  
  df %>%
    group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    group_modify(~ {
      grp <- .x
      
      ab <- grp %>% filter(direction == "ab") %>% arrange(d)
      ba <- grp %>% filter(direction == "ba") %>% arrange(d)
      
      n_ab <- nrow(ab)
      n_ba <- nrow(ba)
      n_max <- max(n_ab, n_ba, 1L)
      
      pad <- function(block, n) {
        if (nrow(block) < n) {
          block <- bind_rows(block, block[rep(NA_integer_, n - nrow(block)), , drop = FALSE])
        }
        block
      }
      ab <- pad(ab, n_max)
      ba <- pad(ba, n_max)
      
      chosen <- map_dfr(seq_len(n_max), function(i) {
        row_ab <- ab[i, , drop = FALSE]
        row_ba <- ba[i, , drop = FALSE]
        d_ab <- row_ab$d
        d_ba <- row_ba$d
        
        if (is.na(d_ab) && is.na(d_ba)) {
          out <- row_ab
          out[] <- NA
        } else if (is.na(d_ab)) {
          out <- row_ba
        } else if (is.na(d_ba)) {
          out <- row_ab
        } else if (d_ab <= d_ba) {
          out <- row_ab
        } else {
          out <- row_ba
        }
        out
      })
      
      chosen$N_ab <- n_ab
      chosen$N_ba <- n_ba
      chosen
    }) %>%
    ungroup() %>%
    relocate(all_of(group_cols), N_ab, N_ba)
}
dfskw <- reconcile_directional_c(dfsk)

dfskw %>% filter(lr_ub > 6.635) %>%
  filter(
    (genome_b == "Heliconius_erato_2" & genome_a=="Heliconius_himera_2") |
      (genome_a == "Heliconius_erato_2" & genome_b=="Heliconius_himera_2") |
      (genome_b == "Heliconius_erato_x_himera-F1" & genome_a=="Heliconius_himera_2") |
      (genome_a == "Heliconius_erato_x_himera-F1" & genome_b=="Heliconius_himera_2") |
      (genome_b == "Heliconius_himera_2" & genome_a=="Heliconius_erato_1") |
      (genome_a == "Heliconius_himera_2" & genome_b=="Heliconius_erato_1")
  ) %>%
  mutate(ptx=interaction(genome_a, genome_b, sep = " / ")) %>%
  mutate(ptx=case_when(
    ptx == "Heliconius_erato_2 / Heliconius_himera_2" ~ "H. himera / H. erato etylus (East / East Andes)",
    ptx == "Heliconius_erato_x_himera-F1 / Heliconius_himera_2" ~ "H. himera / H. erato x himera (F1 Hybrid)",
    ptx == "Heliconius_erato_1 / Heliconius_himera_2" ~ "H. himera / H. erato demophoon (East / West Andes)",
    TRUE ~ ptx)
  ) %>%
  ggplot() +
  aes(x=d, fill=ptx) +
  facet_wrap(~factor(ptx, levels=c("H. himera / H. erato demophoon (East / West Andes)", "H. himera / H. erato etylus (East / East Andes)", "H. himera / H. erato x himera (F1 Hybrid)")), ncol=1, scale="free_y") +
  stat_bin(aes(y = after_stat(count / sum(count)*3)), bins=50, position = position_identity(), show.legend=F, alpha=0.5, color="black") +
  theme_bw() +
  scale_fill_manual(values=c("#994F00","#006CD1","#006C21")) +
  labs(x=TeX(r'($\hat{D}$)'), y="% of windows", fill="Pair") +
  scale_x_continuous(transform = "log10", labels = percent) +
  scale_y_continuous(labels = percent) +
  theme(legend.position = "bottom", legend.direction = "vertical") +
  theme(strip.text = element_text(face = "italic"))

dfskw %>% filter(lr_ub > 6.635) %>%
  filter(
    (genome_b == "Heliconius_erato_2" & genome_a=="Heliconius_himera_2") |
      (genome_a == "Heliconius_erato_2" & genome_b=="Heliconius_himera_2") |
      (genome_b == "Heliconius_erato_x_himera-F1" & genome_a=="Heliconius_himera_2") |
      (genome_a == "Heliconius_erato_x_himera-F1" & genome_b=="Heliconius_himera_2") |
      (genome_b == "Heliconius_himera_2" & genome_a=="Heliconius_erato_1") |
      (genome_a == "Heliconius_himera_2" & genome_b=="Heliconius_erato_1")
  ) %>%
  mutate(ptx=interaction(genome_a, genome_b, sep = " / ")) %>%
  mutate(ptx=case_when(
    ptx == "Heliconius_erato_2 / Heliconius_himera_2" ~ "H. himera / H. erato etylus (East / East Andes)",
    ptx == "Heliconius_erato_x_himera-F1 / Heliconius_himera_2" ~ "H. himera / H. erato x himera (F1 Hybrid)",
    ptx == "Heliconius_erato_1 / Heliconius_himera_2" ~ "H. himera / H. erato demophoon (East / West Andes)",
    TRUE ~ ptx)
  ) %>%
  ggplot() +
  aes(x=d, fill=ptx) +
  facet_wrap(~factor(ptx, levels=c("H. himera / H. erato demophoon (East / West Andes)", "H. himera / H. erato etylus (East / East Andes)", "H. himera / H. erato x himera (F1 Hybrid)")), ncol=1, scale="free_y") +
  # stat_bin(aes(y = after_stat(count / sum(count))*3), breaks = c(0:20)/200, position = position_identity(), show.legend=F, alpha=0.5, color="black") +
  stat_density(aes(y = after_stat(count / sum(count))*40), bins=40, position = position_identity(), show.legend=F, alpha=0.25, color="black", bw = 0.005) +
  stat_bin(aes(y = after_stat(count / sum(count))*3), bins=40, position = position_identity(), show.legend=F, alpha=0.5, color="black") +
  theme_cowplot() +
  labs(x=TeX(r'($\hat{D}$)'), y="% of windows", fill="Pair") +
  # scale_x_continuous(labels = percent) +
  scale_y_continuous(labels = percent) +
  theme(legend.position = "bottom", legend.direction = "vertical") +
  scale_fill_manual(values=c("#994F00","#006CD1","#006C21")) +
  theme(strip.text = element_text(face = "italic", size=7.5)) +
  coord_cartesian(xlim = c(0, 0.15))
ggsave2("./butterfly-example-dist.pdf", width = 3.5, height = 4)

dfskw %>% # filter(grepl("erato", genome_a) | grepl("erato", genome_b)) %>%
  group_by(genome_a, genome_b) %>%
  summarize(dd=mean(d), p=dip(d[is.finite(d)], full.result = FALSE, min.is.0 = FALSE, debug = FALSE)) %>% 
  filter((genome_a=="Heliconius_himera_2") | (genome_b=="Heliconius_himera_2") ) %>% arrange(p) %>% 
  arrange(p)

dfg <- dfskw %>% filter(lr_ub > 6.635) %>%
  group_by(genome_a, genome_b) %>%
  summarize(p = dip(d[is.finite(d)], full.result = FALSE, min.is.0 = TRUE, debug = FALSE), .groups = "drop", d=mean(d, na.rm=T)) %>%
  mutate(s_a=gsub(x=gsub(x=genome_a, pattern="Heliconius_", replacement = ""), "_.*", "")) %>%
  mutate(s_b=gsub(x=gsub(x=genome_b, pattern="Heliconius_", replacement = ""), "_.*", "")) %>%
  mutate(F1a=grepl("F1", genome_a), F1b=grepl("F1", genome_b)) %>%
  mutate(s_a = if_else(F1a, "F1", s_a), s_b = if_else(F1b, "F1", s_b)) %>%
  mutate(ss=s_a==s_b)

dfg %>%
  # filter(!(F1b&!grepl("himera_2", genome_a))) %>%
  # filter(!(F1a&grepl("himera_2", genome_b))) %>%
  # filter(!(genome_a == "Heliconius_erato_1" & F1b)) %>%
  # filter(!(genome_b == "Heliconius_himera_2" & F1a)) %>%
  mutate(cc = case_when(
    ss ~ "Same species",
    (F1a & grepl("himera", genome_b)) | (F1b & grepl("himera", genome_a)) ~ "F1 Hybrid",
    (F1a & grepl("erato", genome_b)) | (F1b & grepl("erato", genome_a)) ~ "F1 Hybrid",
    !F1a & !F1b& (grepl("erato_2", genome_a) & grepl("himera", genome_b)) | (grepl("erato_2", genome_b) & grepl("himera", genome_a)) ~ "H. himera / erato etylus",
    !F1a & !F1b& (grepl("erato_1", genome_a) & grepl("himera", genome_b)) | (grepl("erato_1", genome_b) & grepl("himera", genome_a)) ~ "H. himera / erato demophoon",
    TRUE ~ "Other")
  ) %>%
  ggplot() +
  aes(y=reorder(cc, p), x=p) +
  geom_point(alpha=0.85, size=2.5, aes(color=d), position = position_jitter(width = 0, height = 0.1)) +
  stat_summary(position = position_nudge(y=-0.25),
    fun.min = min,  fun.max = max,  geom = "errorbar",  width = 0.2) +
  stat_summary(
    position = position_nudge(y=-0.25),
    fun = mean, # median,
    geom = "point", shape = 21, fill = "white", color = "gray20", size = 3, stroke = 1.5
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "italic")) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  scale_color_viridis_c()+
  labs(color=TeX(r'($\hat{D}$)'), y="Pairs", x="Hartigans' dip statistic")
ggsave2("./butterfly-dip-dist.pdf", width = 4.5, height = 3.6)
