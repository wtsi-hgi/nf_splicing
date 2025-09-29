t_col <- function(col, rate)
{
    newcol <- rgb(col2rgb(col)["red", ],
                  col2rgb(col)["green", ],
                  col2rgb(col)["blue", ],
                  as.integer(rate * 255),
                  maxColorValue = 255)
    return(newcol)
}

select_colorblind <- function(col_id)
{
    col8 <- c("#D55E00", "#56B4E9", "#E69F00",
              "#009E73", "#F0E442", "#0072B2",
              "#CC79A7", "#000000")

    col12 <- c("#88CCEE", "#CC6677", "#DDCC77",
               "#117733", "#332288", "#AA4499",
               "#44AA99", "#999933", "#882255",
               "#661100", "#6699CC", "#888888")

    col15 <- c("red",       "royalblue", "olivedrab",
               "purple",    "violet",    "maroon1",
               "seagreen1", "navy",      "pink",
               "coral",     "steelblue", "turquoise1",
               "red4",      "skyblue",   "yellowgreen")

    col21 <- c("#F60239", "#009503", "#FFDC3D",
               "#9900E6", "#009FFA", "#FF92FD",
               "#65019F", "#FF6E3A", "#005A01",
               "#00E5F8", "#DA00FD", "#AFFF2A",
               "#00F407", "#00489E", "#0079FA",
               "#560133", "#EF0096", "#000000",
               "#005745", "#00AF8E", "#00EBC1")

    if (col_id == "col8") {
        return(col8)
    } else if (col_id == "col12") {
        return(col12)
    } else if (col_id == "col15") {
        return(col15)
    } else if (col_id == "col21") {
        return(col21)
    } else {
        stop(paste0("====> Error: wrong col_id"))
    }
}

panel.hist <- function(x, ...)
{
    x <- na.omit(x)

    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5))
    his <- hist(x, plot = FALSE)
    breaks <- his$breaks
    nB <- length(breaks)
    y <- his$counts
    y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = t_col("yellowgreen", 0.8), ...)
    lines(density(x), col = 2, lwd = 2)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))

    Cor <- abs(cor(x, y, use = "complete.obs"))
    txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])

    if(missing(cex.cor))
    {
        cex.cor <- 0.4 / strwidth(txt)
    }

    # Generate heatmap color scale based on correlation value
    col.range = c("blue", "white", "red")
    heatmap_col <- colorRampPalette(col.range)(100)
    col_index <- round(Cor * 100)
    bg.col <- heatmap_col[col_index]

    # Draw the rectangle with heatmap color
    rect(0, 0, 1, 1, col = bg.col, border = NA)
    
    # Draw correlation text on top
    text(0.5, 0.5, txt, cex = 0.6 + cex.cor * Cor)
}

panel.smooth <- function(x, y, 
                         col.line = "red", 
                         plot.type = c("point", "smooth"),
                         method = c("loess", "lm"), span = 0.5, degree = 2, level = 0.95, 
                         pch = 21, col.pch = t_col("royalblue", 0.8), col.bg = t_col("royalblue", 0.4), cex.pch = 0.8, ...)
{
    plot.type <- match.arg(plot.type, c("point", "smooth"))

    data <- na.omit(data.frame(x, y))
    x <- data$x
    y <- data$y

    if(plot.type == "point")
    {
        points(x, y, pch = pch, col = col.pch, bg = col.bg, cex = cex.pch, ...)
    } else {
        smoothScatter(x, y, add = TRUE, nrpoints = 0)
    }

    method <- match.arg(method)
    data <- subset(data, x != 0 & y != 0)
    x <- data$x
    y <- data$y
    if(method == "loess")
    {
        loess_fit <- loess(y ~ x, span = span, degree = degree)
        pred <- predict(loess_fit, se = TRUE)
        lines(x, pred$fit, col = col.line, ...)
    } else if (method == "lm") {
        lm_fit <- lm(y ~ x)
        pred <- predict(lm_fit, interval = "confidence", level = level)
        lines(x, pred[, "fit"], col = col.line, ...)
    }
}

bar_chart_pos_neg <- function(label,
                              value,
                              max_value = 1,
                              height = "1rem",
                              pos_fill = t_col("tomato", 0.8),
                              neg_fill = t_col("yellowgreen", 0.8)) {
    neg_chart <- div(style = list(flex = "1 1 0"))
    pos_chart <- div(style = list(flex = "1 1 0"))
    width <- paste0(abs(value / max_value) * 100, "%")

    if (value < 0) {
        bar <- div(style = list(marginLeft = "0.5rem", background = neg_fill, width = width, height = height))
        chart <- div(style = list(display = "flex", alignItems = "center", justifyContent = "flex-end"),
                     label,
                     bar)
        neg_chart <- tagAppendChild(neg_chart, chart)
    } else {
        bar <- div(style = list(marginRight = "0.5rem", background = pos_fill, width = width, height = height))
        chart <- div(style = list(display = "flex", alignItems = "center"), bar, label)
        pos_chart <- tagAppendChild(pos_chart, chart)
    }

    div(style = list(display = "flex"), neg_chart, pos_chart)
}

sort_paths_by_filename <- function(full_paths)
{
    file_names <- basename(full_paths)
    ord <- mixedorder(file_names)
    return(full_paths[ord])
}
