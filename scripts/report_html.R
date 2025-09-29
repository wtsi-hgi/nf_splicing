create_html_render <- function(file_summary_reads, 
                               file_summary_pct,
                               plot_reads_pct,
                               plot_barcodes_venn,
                               file_barcodes_summary,
                               plot_junctions_venn,
                               plot_junctions_corr,
                               plot_junctions_category,
                               file_junctions_category,
                               list_files_junctions_diagram,
                               list_files_junctions_scatter,
                               plot_psi_corr,
                               file_psi,
                               out_render_context)
{
    rmd_render_context <- glue(r"(
---
title: "Splicing Report"
date: "`r format(Sys.time(), '%d %B %Y -- %A -- %X')`"
output:
    html_document:
        toc: true
        toc_depth: 4
        theme: united
        highlight: tango
---

```{{r setup, include = FALSE}}
knitr::opts_chunk$set(echo = TRUE, fig.align = "center")
library(reactable)
library(sparkline)
library(UpSetR)
```

```{{js, echo = FALSE}}
function formatNumber(num, precision = 1) {{
    const map = [
        {{ suffix: 'T', threshold: 1e12 }},
        {{ suffix: 'B', threshold: 1e9 }},
        {{ suffix: 'M', threshold: 1e6 }},
        {{ suffix: 'K', threshold: 1e3 }},
        {{ suffix: '', threshold: 1 }},
    ];
    const found = map.find((x) => Math.abs(num) >= x.threshold);
    if (found) {{
        const formatted = (num / found.threshold).toFixed(precision) + found.suffix;
        return formatted;
    }}
    return num;
}}

function rangeMore(column, state) {{
    let min = Infinity
    let max = 0
    state.data.forEach(function(row) {{
        const value = row[column.id]
        if (value < min) {{
            min = Math.floor(value)
        }}
        if (value > max) {{
            max = Math.ceil(value)
        }}
    }})

    const filterValue = column.filterValue || min
    const input = React.createElement('input', {{
        type: 'range',
        value: filterValue,
        min: min,
        max: max,
        onChange: function(event) {{
            column.setFilter(event.target.value || undefined)
        }},
        style: {{ width: '100%', marginRight: '8px' }},
        'aria-label': 'Filter ' + column.name
    }})

    return React.createElement(
        'div',
        {{ style: {{ display: 'flex', alignItems: 'center', height: '100%' }} }},
        [input, formatNumber(filterValue)]
    )
}}

function filterMinValue(rows, columnId, filterValue) {{
    return rows.filter(function(row) {{
        return row.values[columnId] >= filterValue
    }})
}}

function rangeLess(column, state) {{
    let min = Infinity
    let max = 0
    state.data.forEach(function(row) {{
        const value = row[column.id]
        if (value < min) {{
            min = Math.floor(value)
        }}
        if (value > max) {{
            max = Math.ceil(value)
        }}
    }})

    const filterValue = column.filterValue || max
    const input = React.createElement('input', {{
        type: 'range',
        value: filterValue,
        min: min,
        max: max,
        onChange: function(event) {{
            column.setFilter(event.target.value || undefined)
        }},
        style: {{ width: '100%', marginRight: '8px' }},
        'aria-label': 'Filter ' + column.name
    }})

    return React.createElement(
        'div',
        {{ style: {{ display: 'flex', alignItems: 'center', height: '100%' }} }},
        [input, formatNumber(filterValue)]
    )
}}

function filterMaxValue(rows, columnId, filterValue) {{
    return rows.filter(function(row) {{
        return row.values[columnId] <= filterValue
    }})
}}
```

---

## 1. Introduction
This pipeline is designed to quantify splicing events within minigene-based mutagenesis libraries. 
Minigene systems are widely used to study the regulatory mechanisms of RNA splicing, 
particularly in the context of variant interpretation and functional genomics. 
By introducing specific mutations into synthetic constructs (minigenes), 
researchers can assess how sequence changes affect splicing outcomes in a controlled cellular environment.

---

## 2. Read Processing
This section summarises the distribution of reads according to the alignments.

* **total_reads:** the number of raw reads aftering adaptor trimming and quality filtering.
* **merged_reads:** the number of reads supporting the splicing events.
* **unmerged_reads:** the number of reads supporting the whole minigene transcript without any splicing event.
* **inclusion_reads:** the number of reads supporting canonical exon inclusion.
* **skipping_reads:** the number of reads supporting canonical exon skipping.
* **map_reads:** the number of reads supporting novel splicing events.
* **unexplain_reads:** the number of reads which cannot be categorized into any event.

```{{r, echo = FALSE}}
df <- as.data.frame(read.table("{file_summary_reads}", header = TRUE, sep = "\t", check.names = FALSE))
min_row <- ifelse(nrow(df) > 10, 10, nrow(df))
reactable(df, highlight = TRUE, bordered = TRUE, striped = TRUE, compact = TRUE, wrap = TRUE,
          minRows = min_row, defaultColDef = colDef(minWidth = 150, align = "left"))

df <- as.data.frame(read.table("{file_summary_pct}", header = TRUE, sep = "\t", check.names = FALSE))
min_row <- ifelse(nrow(df) > 10, 10, nrow(df))
reactable(df, highlight = TRUE, bordered = TRUE, striped = TRUE, compact = TRUE, wrap = TRUE,
          minRows = min_row, defaultColDef = colDef(minWidth = 150, align = "left"))
```
<br>

```{{r, echo = FALSE, fig.show = "hold", fig.align = "center", out.height = "50%", out.width = "50%"}}
knitr::include_graphics("{plot_reads_pct}", rel_path = FALSE)
```
<br>

---

## 3. Barcode Processing
Pipeline include a method of barcode detection which may differ from the method of input barcode association.
This section shows the discrepancy between the pipeline barcodes and the input barcodes

```{{r, echo = FALSE, fig.show = "hold", fig.align = "center", out.height = "100%", out.width = "100%"}}
knitr::include_graphics("{plot_barcodes_venn}", rel_path = FALSE)
```
<br>

---

## 4. Junction Processing
This section summarises all the novel splicing events

### 4.1. Correlations between replicates

```{{r, echo = FALSE, fig.show = "hold", fig.align = "default", out.height = "48%", out.width = "48%"}}
knitr::include_graphics(c("{plot_junctions_venn}", "{plot_junctions_corr}"), rel_path = FALSE)
```
<br>

---

### 4.2. Distribution of splicing junctions (appearing in 2 replicates at least)
```{{r, echo = FALSE}}
jucntion_category <- as_tibble(vroom("{file_junctions_category}", delim = "\t", col_names = TRUE, show_col_types = FALSE))
min_row <- ifelse(nrow(jucntion_category) > 20, 20, nrow(jucntion_category))
reactable(jucntion_category, highlight = TRUE, bordered = TRUE, striped = TRUE, compact = TRUE, wrap = TRUE,
          filterable = TRUE, minRows = min_row, defaultPageSize = 20, defaultColDef = colDef(minWidth = 80, align = "left"),
          columns = list("annotation" = colDef(minWidth = 400),
                         "avg_cov" = colDef(filterMethod = JS("filterMinValue"),
                                            filterInput = JS("rangeMore"),
                                            format = colFormat(digits = 2),
                                            minWidth = 200)))


```
<br>

```{{r, echo = FALSE, results = "asis", fig.align = "center", out.height = "80%", out.width = "80%"}}
    )")

    for(i in seq_along(list_files_junctions_diagram)) {
        rmd_render_context <- paste0(rmd_render_context, glue(r"(

cat("<details>\n")
cat("<summary>{names(list_files_junctions_diagram)[i]}</summary>\n")
cat("<br>\n")
knitr::include_graphics("{list_files_junctions_diagram[i]}", rel_path = FALSE)
knitr::include_graphics("{list_files_junctions_scatter[i]}", rel_path = FALSE)
cat("</details>\n\n")

        )"))
    }

    rmd_render_context <- paste0(rmd_render_context, glue(r"(
```
<br>

### 4.3. Splicing events  (appearing in 2 replicates at least) by variants
```{{r, echo = FALSE}}
jucntion_category <- as_tibble(vroom("{file_junctions_category}", delim = "\t", col_names = TRUE, show_col_types = FALSE))
jucntion_category_summary <- jucntion_category %>%
                               separate_rows(annotation, sep = ";") %>%
                               group_by(annotation) %>%
                               filter(annotation != "no_annotation") %>%
                               summarise(log_avg_cov = list(log10(avg_cov + 1)),
                                         no_of_variant = n_distinct(var_id),
                                         .groups = "drop")
boxplot_min <- min(unlist(jucntion_category_summary$log_avg_cov), na.rm = TRUE)
boxplot_max <- max(unlist(jucntion_category_summary$log_avg_cov), na.rm = TRUE)
reactable(jucntion_category_summary, highlight = TRUE, bordered = TRUE, striped = TRUE, compact = TRUE, wrap = TRUE,
          filterable = FALSE, sortable = FALSE, minRows = 12, defaultPageSize = 12, defaultColDef = colDef(minWidth = 200, align = "left"), 
          columns = list(annotation = colDef(name = "Splicing Type"),
                         log_avg_cov = colDef(name = "log10(avg_cov+1) ",
                                              cell = function(value, index) {{
                                                        if (length(value) > 5) {{
                                                            sparkline(value, type = "box", width = 180,
                                                                      chartRangeMin = boxplot_min, chartRangeMax = boxplot_max) }} }}),
                         no_of_variant = colDef(name = "No. of Variants") ))

knitr::include_graphics("{plot_junctions_category}", rel_path = FALSE)
```
<br>

## 5. PSI Correlation
This section summarises the correlation of PSI values between replicates.
```{{r, echo = FALSE, fig.show = "hold", fig.align = "center", out.height = "80%", out.width = "80%"}}
dt_psi <- as.data.table(vroom("{file_psi}", delim = "\t", col_names = TRUE, show_col_types = FALSE))
cols_to_scale <- names(dt_psi)[-1]
dt_psi[, (cols_to_scale) := lapply(.SD, `*`, 100), .SDcols = cols_to_scale]
reactable(dt_psi, highlight = TRUE, bordered = TRUE, striped = TRUE, compact = TRUE, wrap = TRUE,
          filterable = TRUE, minRows = 20, defaultPageSize = 20, 
          defaultColDef = colDef(minWidth = 100, align = "left", format = colFormat(digits = 2),
                                 filterMethod = JS("filterMinValue"), filterInput = JS("rangeMore")),
          columns = list(var_id = colDef(filterMethod = NULL, filterInput = "input")))

knitr::include_graphics("{plot_psi_corr}", rel_path = FALSE)
```
<br>
    )"))

    writeLines(rmd_render_context, out_render_context)
}