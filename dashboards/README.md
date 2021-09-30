# Genomic Dashboards
This document is designed to identify commonly used dashboard elements to improve situational awareness and provide facilitate understanding of geospatial and temporal patterns derived from SARS-CoV-2 genomic data sources. This document is intended as an introduction to dashboard elements, with simple examples using open source technologies.


## Table of Contents
> - [Numeric Counts](#1-numeric-counts)
> - [Maps](#2-maps)
> - [Epidemiologic Curves](#3-epidemiologic-curves)
> - [Variants over Time](#4-variants-over-time)
> - [Tables](#5-tables)
> - [Interaction and Filtering](#6-interaction-and-filtering)


### 1. Numeric Counts

Numeric counts are often highlighted on dashboards and frequently populate the upper- or left-most portion of the dashboard. These counts can represent totals of tests administered, cases, hospitalizations, deaths, or vaccinations. To draw the user's attention to these counts, they are often large, colorful, or accompanied by relevant icons.

![This is an example of dashboard counts from the Utah Department of Health's COVID-19 Case Count Dashboard on April 29th, 2021](images/utahCounts.png)

> Above is an example of dashboard counts from the [Utah Department of Health's COVID-19 Case Count Dashboard](https://coronavirus.utah.gov/case-counts/) on April 29, 2021.

***

In the context of genomic surveillance, count metrics often report the total number of SARS-CoV-2 genomes sequenced, the proportion of cases in a geographic region that have been sequenced, or the total number of sequences reported after a specific point in time.

![This is an example of counts on a SARS-CoV-2 sequencing dashboard developed by TGEN for the state of Arizona](images/tgen_arizona_counts.png)
> Above is an example of counts on a [SARS-CoV-2 sequencing dashboard](https://pathogen.tgen.org/covidseq-tracker/) developed by TGEN for the state of Arizona

It may be helpful to place numeric counts within a header that appear at the top of multiple dashboard views, providing additional context and a consistent reference point from which to compare other data summarized in charts, maps, or figures.

***

Count metrics can be further augmented with indicators of recent trends, as shown below.
![This shows 7-day average counts of numerous variables, with percentage trends indicated by numbers and arrows.](images/UWMad_augmented_counts.png)
> Above shows 7-day average counts of numerous variables taken from data collected at the University of Wisconsin-Madison and available on the [COVID-19 Response-UW-Madison Dashboard](https://covidresponse.wisc.edu/dashboard-2021-jan-may/), with percentage trends indicated by numbers and arrows.


### 2. Maps

A common feature of SARS-CoV-2 dashboards are choropleth maps, or maps where regions are colored according to a variable of interest. For example, variables shown in the tooltip, a text box displaying information when users hover or click a specific area, can detail the underlying data such as case incidence rates, cumulative counts per rolling time window, the proportion of cases where genomes are sequenced, or the number of variants of concern. They can be used to show these values within existing jurisdictional boundaries (e.g., states or counties) or custom areas (e.g., health regions or a “tri-county” area).

![kelseyMaps](images/kelsey_maps.png)

![utahMaps](images/utahMaps.png)

> Above are examples of standard jurisdictions with case counts and sequencing percentages (top-left), custom jurisdictions showing variant counts by for each region (top-right), and incidence rates over a 14-day rolling window, expressed in cases per 100,000 population by county (above). [Utah Department of Health's COVID-19 Case Count Dashboard](https://coronavirus.utah.gov/case-counts/).

Maps can also display information about specific genetic variants in a variety of ways. These may include lists of Variants of Concern (VOC) and counts, as shown above (top-right). Another option is to overlay pie charts that indicate the estimated proportion of specific variants  within a defined geographic region, as shown below.

![](images/dataTracker_map.png)

> Above is a map available on the [CDC COVID Data Tracker](https://covid.cdc.gov/covid-data-tracker/#variant-proportions) showing the distribution of select SARS-CoV-2 PANGO lineages with overlaid pie charts. Here, the geographic areas are defined by HHS Regions to prevent overplotting.

Below is the same map, but filtered to highlight only HHS Region 1, enabling the user to focus on variant proportions in one region of interest.

![](images/mapInteractive.png)
> Above is the same map available on the [CDC COVID Data Tracker](https://covid.cdc.gov/covid-data-tracker/#variant-proportions) showing the distribution of select SARS-CoV-2 PANGO lineages with overlaid pie charts, only this map has been filtered to show HHS Region 1.

When visualizing public health data in real-time with dashboards, the most recently available data is often incomplete and thus unreliable for trend analysis due to small sample size. Indicating the time point between which analyses are supported by sufficient data and uncertainty due to low sampling. In the example below, sampling density is aligned to provide users with context when evaluating the highlighted time period.

![](images/integratedWidgetMap.png)

> The above maps are filtered according to two different sampling densities (>=3 on the left, >=500 on the right) for Variant of Concern: B.1.1.7 (Alpha). Increasing the sampling density threshold helps to minimize random noise in the underlying data. These maps are available on the [CDC COVID Data Tracker](https://covid.cdc.gov/covid-data-tracker/#variant-proportions).


### 3. Epidemiologic Curves

Epidemiologic curves, or simply 'epi curves,' are a staple of any public health investigation. These visualizations show count metrics over time, typically broken out into daily, weekly, or monthly bins. To overcome the variability introduced by periodic reporting cycles, rolling averages (eg. 7-day) are commonly overlaid, which will often appear offset from the underlying bar plot.

![Utah epi curve](images/utahEpiCurve.png)

> Above is an example of epidemiologic curves displayed on [Utah Department of Health's COVID-19 Case Count Dashboard](https://coronavirus.utah.gov/case-counts/).

Epidemiologic curve visualizations can be further customized with colors to differentiate categorical variables of interest, including either clinical or demographic.

![](images/utah_coloredEpiCurve.png)
> Above is an example of epidemiologic curves colored according to categorical variable indicating test results. This figure is displayed on [Utah Department of Health's COVID-19 Case Count Dashboard](https://coronavirus.utah.gov/case-counts/).

### 4. Variants over Time

As the SARS-CoV-2 pandemic continues, focus has shifted to tracking the prevalence of particular Variants of Interest (VOIs) or Variants of Concern (VOCs). The ability to track changes in frequency over time is fundamental to genomic surveillance, and visualization of those data can vary.

A common approach to visualizing variant proportions over time is with a stacked bar chart. Each bar represents a unit of time, broken down by the proportion of variants identified during that period, with the most prevalent variants (often >=5%) labeled for clarity. Multiple time points are presented for comparison as discrete, adjacent bars. While helpful, this view is particularly hard to compare with quantitative precision because bar slices do not align from one time point to the next. To aid interpretation, it can be helpful to accompany a stacked bar chart with specific metrics organized in an adjacent table, as shown below, or via an interactive tooltip.

![variant tracker](images/cdc_variant_tracker_updated.png)

> Above is an example display of variant tracking over time, available on the [Variant Proportions](https://covid.cdc.gov/covid-data-tracker/#variant-proportions) subsection of the [CDC COVID Data Tracker](https://covid.cdc.gov/covid-data-tracker/#datatracker-home).

In the figure above, there are a number of subtle cues and interactive layers that provide additional context and detail. For example:
* The most recent two-week span is:
    1. highlighted with a bold selection border to focus the attention
    2. marked with a double-asterisk to denote that it is subject to change due a to sample collection and processing delay
* Only the most prevalent variants are labeled with text to prevent overplotting
* The stacked bar chart is accompanied by a color-matched table providing additional detail, including estimated proportions and confidence intervals


Alternatively, visualization of variant proportions can be depicted with continuous time, rather than cutting the dataset up into discrete interval bins. One prominent example of this are the variants tracking visualization powered by [Nextstrain](https://nextstrain.org), shown below. This approach provides smoothed, visual tracking of variant frequency over time but at the cost of direct comparative analyses enabled by discrete time bins.

![](images/nextstrain_smoothed_variants_over_time.png)

> Above is an example of an integrated, smoothed view of variant proportions over time produced by [Nextstrain](https://nextstrain.org), as integrated with the [Connecticut COVID Tracker](https://covidtrackerct.com/variant-surveillance/ )

As with maps, it is important to note that the array of Variants of Interest (VOIs) and Variants of Concern (VOCs) under surveillance are subject to change. Therefore, stacked bar charts can rapidly become complex and difficult to track. One way to combat information overload is to build interactive features and filters that empower users to focus on the most critical values and patterns in the visualization.

### 5. Tables

Tables are typically used on dashboards as a secondary visual to provide details and context to a primary visual, as shown above in the [Variants over Time](#4-variants-over-time) section. Tables, like the one shown below depicting variant frequencies in Washington State, are a staple format for epidemiological reports about public health investigations.

![washington_covid_variant_table](images/washingtonTable1.png)

***

Tables can also be used to provide very granular information across multiple variables, such as the cross-tabulation of variants by county one shown below, also from Washington State.

![washington covid county-level crosstab of variants](images/washingtonTable2.png)
> Above is a county-level cross-tabulation of VOIs and VOCs, which can be helpful in digging through very granular data. Tables are an ideal companion to provide context to broad, summary visualizations (as identified on the [Variants over Time](#4-variants-over-time) section).


### 6. Interaction and Filtering

 Dashboards are often designed to provide near real-time views of data, and their performance can vary dramatically depending on the scale of the underlying dataset. For example, when too much information is shown, visualization can fall victim to 'overplotting'. Overplotting describes the situation where data or labels overlap, making it difficult to discern individual data points or interpret patterns. Overplotting typically occurs with a large number of data points or a small number of unique values.

![](images/overplotExample.png)

> Above are two examples of overplotting; a scatter plot with an overwhelming number of data points that is difficult to read (left) and a pie chart with wedge labels  obscuring the data (right).

Solutions to overplotting include:
1. reducing data point size
1. changing data point shape, jitter, or transparency
1. tiling or subsetting data
1. algorithms to aggregate, cluster, or prevent label overlapping

***
Case Study: [Johns Hopkins COVID-19 Map](https://coronavirus.jhu.edu/map.html)

As the COVID-19 pandemic became widespread, the Johns Hopkins Coronavirus Resource Center's cumulative dashboard suffered from overplotting. The original map below better represents population distributions and territorial boundaries than progress of the pandemic.

![JHU_overplot_example](images/JHU_covid_map.png)

Johns Hopkins responded by producing a detailed heatmap visualization that is released in video format once per day. These videos effectively act as a rapid but effective walkthrough of multiple dashboard elements in a short period of time.

![](images/JHU_improved.png)
> Above is a screen capture of the [Johns Hopkins University's Daily COVID-19 Video](https://coronavirus.jhu.edu/covid-19-daily-video) captured on April 29th, 2021.

Variants of Interest and Variants of Concern can vary dynamically over time, and so it may be ineffective to report the frequencies for all variant at all time points. Rather, it may be useful to only report variants according to defined criteria, such as:
1. only the most frequent (top 5, top 10, etc.)
2. all variants above a threshold (1%, 5%, etc.)

It may be helpful to provide additional filtering parameters or graphical context to describe sequence data availability for the period of time or region under evaluation. The example below enables report-level filtering based on selected time intervals.

![](images/finalInteractive_outbreaksInfo.png)
> Above is an example of variant proportion visualization in one of the *Regional Reports* by [Outbreaks.Info](https://outbreak.info).
