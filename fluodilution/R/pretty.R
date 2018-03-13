# Copyright (c) 2015-2018 Hadrien Chauvin
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#####################

#' @keywords internal
#' @export
print.summary.fd_data <- function (x, ...) {
  data <- x
  print(summary.data.frame(data))

  cat("Inoculums:\n")

  print(data.frame("Log mean" = attr(data, "fmm")$m0,
                   "Log std.dev." = attr(data, "fmm")$sd0))

  cat("\n(Default) cutoff points:\n")

  print(attr(data, "cutoff")[
    as.character(unique(data$Timepoint[data$Weight == "hist"]))])

  cat("\nProliferation model:\n")
  pro <- attr(data, "proliferation")
  cat("  Categories:", paste(pro$categories, collapse=", "), "\n")
  cat("  Times:", paste(pro$times, collapse=", "), "\n")
  cat("  Maximum generation number:", pro$mgen, "\n")
}

#' @export
#' @keywords internal
summary.fd_data <- function (object, ...) {
  structure(object, class="summary.fd_data")
}

#' @rdname fd_data
#' @export
#' @method plot fd_data
plot.fd_data <- function (x, type="overview",
              main = NULL,
              ...) {
  if (is.null(main)) main <- deparse(substitute(x))

  data <- x

  if (is.null(data$Individual)) {
    if (inherits(data, "groupedData")) {
      data$Individual <- data[[deparse(getGroupsFormula(data)[[2L]])]]
    } else {
      data$Individual <- "indiv_1"
    }
  }

  type <- match.arg(type,
            c("overview", "all", "hist", "N", "range", "balancing",
            "coverage", "cutoff"),
            several.ok=TRUE)

  if (("overview" %in% type || "all" %in% type) && length(type) != 1) {
    stop("when 'type' is \"overview\" or \"all\", ",
       "it cannot specify other types")
  }

  p_hist <- p_N <- p_range <- p_balancing <- p_coverage <- p_cutoff <- NULL

  overview <- FALSE
  if ("overview" %in% type) {
    overview <- TRUE
    type <- c("hist", "N", "range", "balancing")
  } else if ("all" %in% type) {
    type <- c("hist", "N", "range", "coverage", "balancing", "cutoff")
  } else if (any(duplicated(type)))
    stop("type: duplicated elements")

  htrans_inverse <- attr(data, "fmm")$htrans
  if (is.null(htrans_inverse)) htrans_inverse <- "identity"
  if (is.character(htrans_inverse))
    htrans_inverse <- match.fun(paste0(htrans_inverse, "_trans"))()
  htrans_inverse <- htrans_inverse$inverse

  cctrans_inverse <- attr(data, "fmm")$cctrans
  if (is.null(cctrans_inverse)) cctrans_inverse <- "identity"
  if (is.character(cctrans_inverse))
    cctrans_inverse <- match.fun(paste0(cctrans_inverse, "_trans"))()
  cctrans_inverse <- cctrans_inverse$inverse

  if ("hist" %in% type) {
    # NOTE: use asinh transformation instead?
    df_hists <- fd_transform(data %>% subset(Weight == "hist"),
                 hist = "identity")
    if (NROW(df_hists) > 0) {
      if (length(unique(df_hists$Category)) == 1) {
        facet <- facet_wrap(
          ~ Time,
          labeller = function (labels) {
            list(Time = paste0(round(labels$Time, digits=2), "h"))
          })
      } else {
        facet <- facet_grid(
          Time ~ Category,
          labeller = function (labels) {
            if (names(labels) == "Time") {
              list(Time =
                   paste0(round(labels$Time, digits=2), "h"))
            } else
              list(Category = as.character(labels$Category))
          })
      }
      df_hists2 <- df_hists %>%
        transform(xmean = (a + b) / 2.0) %>%
        subset(xmean > 0)
      (p_hist <- ggplot(df_hists2)+
        stat_unitarea(aes(x= xmean,
                  yu=y,
                  group=Timepoint, color=Type))+
        facet +
        labs(x="Fluorescence", title="Histograms (unit area)")+
        theme_classic()+
        scale_color_brewer(palette="Set1")+
        theme(axis.title.y=element_blank(),
            legend.title=element_blank(),
            legend.margin=margin(0, 0, 0, 0))+
        theme(legend.position="bottom")+
        annotation_logticks(sides="b")+
        scale_x_log10())
    }
  }

  if ("N" %in% type) {
    ss_Ns <- transform(subset.data.frame(data, Type == "Ns"),
               y = cctrans_inverse(y))
    if (NROW(ss_Ns) > 0) {
      p2_aes <- list(x =~ as.numeric(as.character(Time)),
               y =~ y)
      legend.position <- "bottom"
      if (nlevels(ss_Ns$Individual) > 1)
        p2_aes <- append(p2_aes, list(group =~ Individual))
      if (nlevels(ss_Ns$Category) > 1)
        p2_aes <- append(p2_aes, list(color =~ Category))
      else if (nlevels(ss_Ns$Individual) > 1) {
        p2_aes <- append(p2_aes, list(color =~ Individual))
        legend.position <- "none"
      }
      (p_N <- ggplot(ss_Ns, do.call(aes_, p2_aes))+
        stat_summary(fun.y=mean, geom="line")+
          labs(title="Cell counts")+
          theme_classic()+
        geom_point()+
        stat_summary(geom="pointrange", fun.data=function (dat) {
          data.frame(ymin = min(dat), ymax = max(dat), y = mean(dat))
        })+
        scale_color_brewer(palette="Set1")+
        labs(x="Time (hours)")+
        theme(axis.title.y=element_blank())+
          scale_y_log10()+
        annotation_logticks(sides="l")+
        theme(legend.position=legend.position))
    }
  }

  if ("range" %in% type) {
    sd0 <- attr(data, "fmm")$sd0
    m0 <- attr(data, "fmm")$m0
    if (!is.null(sd0) && !is.null(sd0)) {
      df_moments <- data.frame(m0 = exp(m0),
                   sd0 = exp(sd0))
      df_moments$Inoculum <- factor(rownames(df_moments),
                      levels = levels(data$Inoculum))
      mean_m0 <- exp(mean(log(df_moments$m0)))
      if (is.null(attr(data, "cutoff"))) {
        min_cutoff <- mean_m0 / 2.0 ^ 7
      } else {
        min_cutoff <- min(attr(data, "cutoff"))
      }
      df_gen <- transform(data.frame(
        gen = 0:ceiling(log2(mean_m0 / min_cutoff))),
        y = mean_m0 / 2.0 ^ gen)
      (p_range <- ggplot(df_moments)+
        geom_hline(data=df_gen,
               aes(yintercept=y), linetype="dashed")+
        annotate(geom="label", y = df_gen$y, label=df_gen$gen,
             x=levels(df_moments$Inoculum)[
               length(levels(df_moments$Inoculum))])+
        geom_pointrange(aes(x = Inoculum, y = m0,
                  ymin = m0 / sd0, ymax = m0 * sd0))+
        theme_classic()+
        labs(x="Inoculum #")+
        theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank())+
        scale_y_log10()+
        annotation_logticks(sides="l")+
        labs(y="Fluorescence range",
           title="m0, sd0 (cutoff)"))
      if (!is.null(attr(data, "cutoff"))) {
        df_cutoff <-
          plyr::ddply(
            subset.data.frame(data, Weight == "hist"),
            "Timepoint",
            function (df) {
              data.frame(
                Inoculum = df$Inoculum[1],
                Category = df$Category[1],
                Cutoff = attr(data, "cutoff")[df$Timepoint[1]])
            })
        p_range <- p_range +
          geom_point(data=subset(df_cutoff, is.finite(Cutoff)),
                 aes(x = Inoculum, y = Cutoff, color=Category),
                 shape=4, size=3
            )+
          scale_color_brewer(name="Cutoff", palette="Set1")+
          theme(legend.position="bottom")
      }
      (p_range)
    } else {
      type <- setdiff(type, "range")
    }
  }

  if ("balancing" %in% type) {
    df_balancing <-
      plyr::ddply(data, c("Type", "Time", "Category", "Individual"),
            plyr::summarise,
            Replicates = length(unique(Timepoint)))
     (p_balancing <- ggplot(df_balancing)+
      geom_hline(aes(yintercept = as.numeric(as.character(Time))),
             color="gray")+
      geom_point(aes(x = Individual,
               y = as.numeric(as.character(Time))+
                 c(0.2, -0.2, 0.4)[Type],
               size = as.character(Replicates),
               color = Type), alpha=0.8)+
      facet_wrap(~ Category)+
      scale_size_manual(breaks=c("1", "2", "3"), values=c(1, 2, 3))+
      theme_classic()+
       labs(x="Individual #", size="# Replicates") + # nolint
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+
      scale_color_brewer(palette="Set1")+
      labs(y="Time (hours)", title="Balancing"))
  }

  if ("coverage" %in% type) {
    if (is.null(attr(data, "counts")))
      stop("cannot plot coverage if 'counts' is missing")
    df_coverage <- plyr::ddply(data, c("Weight", "Time", "Category"),
              plyr::summarise,
              Replicates = length(unique(Timepoint)),
              Count = sum(attr(data, "counts")[
                unique(Timepoint)])) %>%
      transform(Time = factor(Time, levels=sort(unique(Time))),
            Label = as.character(Replicates),
            Size = (Replicates - min(Replicates)) /
              (max(Replicates) - min(Replicates)),
            Type = ifelse(Weight == "N",
                  "# cell counts",
                  "# histograms"))
    df_coverage %<>% rbind(transform(
      subset(df_coverage, Count > 0),
      Type = "# events",
      Label = paste0("scriptstyle(10)^", round(log10(Count))),
      Size = (Count - min(Count)) / (max(Count) - min(Count))))
     (p_coverage <- ggplot(df_coverage,
                 aes(x = Category,
               y = Time))+
      geom_point(aes(color = Type, size = Size), alpha=1)+
       geom_text(aes(label=Label), parse=TRUE, color="white")+
       facet_grid(Type ~ .)+
      theme_classic()+
       theme(axis.title.x=element_blank(), legend.position="none")+
       labs(size="# replicates")+
      scale_color_brewer(palette="Set1")+
       scale_size_continuous(range=c(7, 15))+
      labs(y="Time (hours)", title="Coverage"))
  }

  if ("cutoff" %in% type) {
    df_cutoff <- plyr::ddply(
      transform(subset.data.frame(data, Weight == "hist"),
            y = htrans_inverse(y)),
      c("Time", "Category", "Timepoint"),
      function (df) {
        ans <-
          data.frame(
            Cutoff =
              sum(df$y[df$a <= attr(data, "cutoff")[
                    df$Timepoint[1]]]))
        if (!is.null(attr(data, "counts")))
          ans$Count <- attr(data, "counts")[df$Timepoint[1]]
        else ans$Count <- NA
        ans
      }
    )
    (p_cutoff <- ggplot(df_cutoff,
         aes(x=factor(Time), y=Cutoff))+
      facet_grid(Category ~ .)+
      geom_point(stat="identity", color="black")+
      theme_classic()+
      labs(x="Time (hours)")+
      coord_flip()+
      theme(legend.position="bottom")+
      theme(panel.grid.major=element_line(color="gray"),
          panel.grid.major.x=element_blank())+
      scale_y_continuous(labels=scales::percent_format())+
      ggtitle("Percentage below cut off"))
  }

  if (overview) {
    if (is.null(p_range)) g2 <- p_balancing
    else g2 <- gridExtra::arrangeGrob(p_range, p_balancing, ncol=2)
    if (is.null(p_N) && is.null(p_hist))
      gridExtra::grid.arrange(g2, top=grid::textGrob(main))
    else {
      if (is.null(p_N))
        g1 <- p_hist
      else if (is.null(p_hist))
        g1 <- p_N
      else g1 <- gridExtra::arrangeGrob(p_hist, p_N, ncol=2,
                        widths=c(2, 1))
      gridExtra::grid.arrange(g1, g2, heights=c(1.2, 1),
                  top=grid::textGrob(main))
    }
  } else {
    if (!is.null(p_hist)) print(p_hist)
    if (!is.null(p_N)) print(p_N)
    if (!is.null(p_range)) print(p_range)
    if (!is.null(p_balancing)) print(p_balancing)
    if (!is.null(p_coverage)) print(p_coverage)
    if (!is.null(p_cutoff)) print(p_cutoff)
  }
}
