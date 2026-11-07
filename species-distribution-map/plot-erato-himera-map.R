#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Publication-ready distribution map: Heliconius erato etylus, H. erato demophoon,
## H. himera. Occurrence localities only (no range polygons), on a land-only
## topographic base (SRTM hillshade + contours, ocean blank white), single
## continuous panel cropped to the isthmus + Andean corridor: eastern edge at the
## Andean crest (no flat Amazon band), western Pacific left white and carrying the
## world locator inset, scale bar and north arrow. No caption.
##
## Input  (same directory): localities-erato-himera.tsv
##                          (source_type: gbif = GBIF occurrence records;
##                           literature = curated literature / collection localities;
##                           see references-erato-himera.md)
## Output (same directory): erato-himera-distribution.pdf / .png
##
## Symbols: filled circle = GBIF occurrence record; open square = iNaturalist
## research-grade observation (literature localities are kept in the TSV/docs but
## not plotted). Fill colour = taxon.
## Colours: demophoon #994F00 and etylus #006CD1 from the himera-pair palette of
## ani-analysis/plot-butterfly-ani.R; himera #D0D55C (yellow) from the same
## script's `mc` palette.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2); library(sf); library(terra); library(elevatr)
  library(rnaturalearth); library(cowplot)
})

script_dir <- {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd()
}

COL <- c("Heliconius erato demophoon" = "#994F00",
         "Heliconius erato etylus"    = "#006CD1",
         "Heliconius himera"          = "#D0D55C")   # yellow, from same script's mc palette
LAB <- c("Heliconius erato demophoon" = expression(italic("H. erato demophoon")),
         "Heliconius erato etylus"    = expression(italic("H. erato etylus")),
         "Heliconius himera"          = expression(italic("H. himera")))
XL <- c(-90.7, -74.7)        # open Pacific (white; carries inset/scale bar) to upper Amazon foothills
YL <- c(-7.7, 10.8)          # continuous latitude, no axis break
FW <- 6.0                    # figure size (in)
FH <- (FW - 0.7) * diff(YL) / diff(XL) + 1.35

## ---- occurrence localities --------------------------------------------------
pts <- read.delim(file.path(script_dir, "localities-erato-himera.tsv"), stringsAsFactors = FALSE) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
pts$lon <- st_coordinates(pts)[, 1]      # st_as_sf drops the coord columns
pts$lat <- st_coordinates(pts)[, 2]
print(table(pts$taxon, pts$source_type))
## plot GBIF + iNaturalist only; literature localities stay in the TSV/docs, unplotted
pts <- pts[pts$source_type != "literature", ]

## ---- land-only topographic base (hillshade + contours; ocean = white) ------
tmpl <- terra::rast(ext = terra::ext(XL[1], XL[2], YL[1], YL[2]),
                    crs = "EPSG:4326", res = (XL[2] - XL[1]) / 700)
el   <- terra::rast(get_elev_raster(tmpl, z = 7, src = "aws"))
el   <- crop(el, ext(XL[1], XL[2], YL[1], YL[2]))
el[el <= 0] <- NA                          # drop bathymetry: ocean stays blank white
hs   <- terra::shade(terrain(el, "slope", unit = "radians"),
                     terrain(el, "aspect", unit = "radians"),
                     direction = 315, angle = 45)
m    <- as.matrix(hs, wide = TRUE)                       # grid matrix; row 1 = north
q    <- quantile(m, c(0.02, 0.98), na.rm = TRUE)         # quantile stretch for contrast
v    <- pmax(0, pmin(1, (m - q[1]) / diff(q)))
v[is.na(v)] <- 1
cols <- grey(v * 0.55 + 0.42)
cols[is.na(m)] <- "#FFFFFF"
hill <- grDevices::as.raster(matrix(cols, nrow = nrow(m)))
ct   <- st_as_sf(as.contour(el, levels = seq(0, 6000, by = 500)))

## ---- small world locator inset ----------------------------------------------
world <- ne_countries(scale = 110, returnclass = "sf")
pW <- ggplot() +
  geom_sf(data = world, fill = "grey82", colour = "grey55", linewidth = 0.15) +
  annotate("rect", xmin = XL[1], xmax = XL[2], ymin = YL[1], ymax = YL[2],
           fill = NA, colour = "grey25", linewidth = 0.5) +
  coord_sf(xlim = c(-180, 180), ylim = c(-58, 80), expand = FALSE, datum = NA) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = "grey35", linewidth = 0.4),
        plot.margin = margin(1, 1, 1, 1))
pW
iw <- 7.0                                  # inset width, degrees of longitude

## ---- map ---------------------------------------------------------------------
p <- ggplot() +
  annotation_raster(hill, XL[1], XL[2], YL[1], YL[2], interpolate = TRUE) +
  geom_sf(data = ct, inherit.aes = FALSE, colour = "grey40", linewidth = 0.22, alpha = 0.85) +
  annotation_custom(grob = ggplotGrob(pW),
                    xmin = XL[1] + 0.3, xmax = XL[1] + 0.3 + iw,
                    ymin = mean(YL) - iw / 2.55 / 2, ymax = mean(YL) + iw / 2.55 / 2) +
  geom_point(data = pts, aes(x = lon, y = lat, fill = taxon, shape = source_type),
             inherit.aes = FALSE, size = 3.5, colour = "grey20", stroke = 0.4, alpha = 0.95) +
  scale_fill_manual(name = NULL, values = COL, labels = LAB,
                    guide = guide_legend(order = 1,
                                         override.aes = list(shape = 21, size = 4.2,
                                                             colour = "grey20", stroke = 0.4))) +
  scale_shape_manual(name = NULL, values = c(gbif = 21, inaturalist = 22),
                     labels = c(gbif = "GBIF occurrence records",
                                inaturalist = "iNaturalist research-grade observations"),
                     guide = guide_legend(order = 2,
                                          override.aes = list(size = 4.2, fill = NA,
                                                              colour = "grey20", stroke = 0.5))) +
  coord_sf(xlim = XL, ylim = YL, expand = FALSE, datum = NA, crs = 4326) +
  theme_map() +
  theme(panel.border = element_rect(colour = "grey35", fill = NA, linewidth = 0.5),
        plot.background = element_rect(fill = "white", colour = NA),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.direction = "horizontal",
        legend.text = element_text(size = 12),
        legend.margin = margin(0, 0, 0, 0))
p
## scale bar + north arrow (over the blank Pacific)
kmdeg <- 111.32 * cos(-5 * pi / 180)
sb_x  <- XL[1] + 0.04 * diff(XL); sb_y <- YL[1] + 0.05 * diff(YL); sb_L <- 300 / kmdeg
p <- p +
  annotate("segment", x = sb_x, xend = sb_x + sb_L, y = sb_y, yend = sb_y,
           linewidth = 1.25, colour = "grey10") +
  annotate("text", x = sb_x + sb_L / 2, y = sb_y + 0.025 * diff(YL),
           label = "300 km", size = 4, colour = "grey10") +
  annotate("segment", x = XL[1] + 0.9, xend = XL[1] + 0.9,
           y = YL[2] - 0.9, yend = YL[2] - 0.2,
           linewidth = 2, colour = "grey10",
           arrow = arrow(length = unit(2.5, "mm"), type = "closed")) +
  annotate("text", x = XL[1] + 0.9, y = YL[2] - 1.35, label = "N", size = 4, colour = "grey10")
p
ggsave(file.path(script_dir, "erato-himera-distribution.pdf"), p, width = FW, height = FH)
ggsave(file.path(script_dir, "erato-himera-distribution.png"), p, width = FW, height = FH, dpi = 300)
cat(sprintf("saved erato-himera-distribution.pdf/.png (%.2f x %.2f in)\n", FW, FH))


library(vroom); library(dplyr); library(tidyr); library(ggplot2); library(patchwork); library(diptest)
library(ggpubr); library(cowplot); library(scales); library(latex2exp); library(purrr)  # theme_minimal_hgrid/vgrid, percent, TeX, map_dfr

SYM     <- "min"
BULK    <- "(0.1,0.15]"
LR_MIN  <- 6.635
# "#9D4030"
mc = c("#D0D55C", "#809D6F", "#BD4030", "#C3A97E", "#769DA6", "#A69E33", "#734002", "#595622", "#8C873F",  "#474B71")
BINS    <- c("[0,0.05)", "[0,0.1)", "(0.0,0.15]")
mpwd <- vroom("../ani-analysis/merged-distances.tsv")
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
# ggsave2("./butterfly-comparison.pdf", width = 6.5, height = 5.25)

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
dfsk <- vroom("../ani-analysis/samples/samples-gdiff-k25-l1000-n10000-h3.tsv", comment = "#") %>% 
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
    ptx == "Heliconius_erato_2 / Heliconius_himera_2" ~ "H. himera / H. erato etylus (spatial adjacency)",
    ptx == "Heliconius_erato_x_himera-F1 / Heliconius_himera_2" ~ "H. himera / H. erato x himera (F1 hybrid)",
    ptx == "Heliconius_erato_1 / Heliconius_himera_2" ~ "H. himera / H. erato demophoon (strictly allopatric)",
    TRUE ~ ptx)
  ) %>%
  ggplot() +
  aes(x=d, fill=ptx) +
  facet_wrap(~factor(ptx, levels=c("H. himera / H. erato demophoon (strictly allopatric)", "H. himera / H. erato etylus (spatial adjacency)", "H. himera / H. erato x himera (F1 hybrid)")), ncol=1, scale="free_y") +
  stat_bin(aes(y = after_stat(count / sum(count)*3)), bins=50, position = position_identity(), show.legend=F, alpha=0.5, color="black") +
  theme_bw() +
  scale_fill_manual(values=c("#994F00","#006CD1","#006C21")) +
  labs(x=TeX(r'($\hat{D}$)'), y="% of windows", fill="Pair") +
  scale_x_continuous(transform = "log10", labels = percent) +
  scale_y_continuous(labels = percent) +
  theme(legend.position = "bottom", legend.direction = "vertical") +
  theme(strip.text = element_text(face = "italic"))

p1 <- dfskw %>% filter(lr_ub > 6.635) %>%
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
    ptx == "Heliconius_erato_2 / Heliconius_himera_2" ~ "H. himera / H. erato etylus (spatial adjacency)",
    ptx == "Heliconius_erato_x_himera-F1 / Heliconius_himera_2" ~ "H. himera / H. erato x himera (F1 hybrid)",
    ptx == "Heliconius_erato_1 / Heliconius_himera_2" ~ "H. himera / H. erato demophoon (strictly allopatric)",
    TRUE ~ ptx)
  ) %>%
  ggplot() +
  aes(x=d, fill=ptx) +
  facet_wrap(~factor(ptx, levels=c("H. himera / H. erato demophoon (strictly allopatric)", "H. himera / H. erato etylus (spatial adjacency)", "H. himera / H. erato x himera (F1 hybrid)")), ncol=1, scale="free_y") +
  # stat_bin(aes(y = after_stat(count / sum(count))*3), breaks = c(0:20)/200, position = position_identity(), show.legend=F, alpha=0.5, color="black") +
  stat_density(aes(y = after_stat(count / sum(count))*40), bins=40, position = position_identity(), show.legend=F, alpha=0.25, color="black", bw = 0.005) +
  stat_bin(aes(y = after_stat(count / sum(count))*3), bins=40, position = position_identity(), show.legend=F, alpha=0.5, color="black") +
  theme_cowplot(font_size = 18) +
  labs(x=TeX(r'($\hat{D}$)'), y="% of windows", fill="Pair") +
  # scale_x_continuous(labels = percent) +
  scale_y_continuous(labels = percent) +
  theme(legend.position = "bottom", legend.direction = "vertical") +
  scale_fill_manual(values=c("#994F00","#006CD1","#006C21")) +
  theme(strip.text = element_text(face = "italic")) +
  coord_cartesian(xlim = c(0, 0.15))
# ggsave2("./butterfly-example-dist.pdf", width = 3.5, height = 4)

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

p3 <- dfg %>%
  # filter(!(F1b&!grepl("himera_2", genome_a))) %>%
  # filter(!(F1a&grepl("himera_2", genome_b))) %>%
  # filter(!(genome_a == "Heliconius_erato_1" & F1b)) %>%
  # filter(!(genome_b == "Heliconius_himera_2" & F1a)) %>%
  mutate(cc = case_when(
    ss ~ "Same species",
    (F1a & grepl("himera", genome_b)) | (F1b & grepl("himera", genome_a)) ~ "F1 hybrid",
    (F1a & grepl("erato", genome_b)) | (F1b & grepl("erato", genome_a)) ~ "F1 hybrid",
    !F1a & !F1b& (grepl("erato_2", genome_a) & grepl("himera", genome_b)) | (grepl("erato_2", genome_b) & grepl("himera", genome_a)) ~ "H. himera / erato etylus",
    !F1a & !F1b& (grepl("erato_1", genome_a) & grepl("himera", genome_b)) | (grepl("erato_1", genome_b) & grepl("himera", genome_a)) ~ "H. himera / erato demophoon",
    TRUE ~ "Other")
  ) %>%
  ggplot() +
  aes(y=reorder(cc, p), x=p) +
  geom_point(alpha=0.85, size=2.5, aes(color=d), position = position_jitter(width = 0, height = 0.1)) +
  stat_summary(position = position_nudge(y=-0.33),
               fun.min = min,  fun.max = max,  geom = "errorbar",  width = 0.2) +
  stat_summary(
    position = position_nudge(y=-0.33),
    fun = mean, # median,
    geom = "point", shape = 21, fill = "white", color = "gray20", size = 3, stroke = 1.5
  ) +
  theme_minimal(base_size = 18) +
  theme(axis.text.y = element_text(face = "italic")) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  scale_color_viridis_c(option = "mako")+
  labs(color=TeX(r'($\hat{D}$)'), y="Pairs", x="Hartigans' dip statistic")
# ggsave2("./butterfly-dip-dist.pdf", width = 4.5, height = 3.6)
p3

plot_grid(
  p1,
  p + ,
  p3 + theme(axis.title.y = element_blank()) +
    theme(axis.text.y = element_text(vjust = 1)) +
    theme(legend.position = "top",  legend.key.width = unit(1, "cm"), legend.text = element_text(angle=30, hjust=1)),
nrow=1)
ggsave("./butterfly-combined.pdf", width = FW*3, height = FH)
