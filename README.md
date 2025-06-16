# Cut-through corridors identification (Simple version)


### Overview

Identifying cut-through routes, commonly known as rat-runs, is a crucial
step in transport model development, as recommended in the Department
for Transport’s guidance (TAG Unit M3.1). These routes represent
diversions from main roads onto local or minor roads, often to bypass
congestion or avoid traffic signals. While traditionally, the
identification of such routes relies heavily on local knowledge, this
may not always be readily available or consistently applied. This
reproducible workflow, implemented in R, offers a data-driven approach.
It leverages open data sources and network analysis techniques to
systematically and efficiently identify potential cut-through corridors,
even in the absence of detailed local insights. This method provides a
transparent and replicable process for transport planners and modellers.

## Set-up

To begin our analysis, we first need to load the necessary R packages.
These packages provide the functions required for spatial data
manipulation, network analysis, and data visualization.

For this demonstration, we will focus our analysis on the city of Leeds.
We use the `zonebuilder` package to define our study area. This package
helps in creating consistent and reproducible geographical zones for
analysis.

``` r
selected_zones <- zonebuilder::zb_zone("Leeds",n_circles = 4)

# We then create a Well-Known Text (WKT) representation of the convex hull of these zones.
# This WKT filter will be used later to select network data only within our area of interest.
zones_wkt <- selected_zones |>
  st_union() |>
  st_convex_hull() |>
  st_transform(27700) |> 
  st_as_text()
```

Next, we define assumed speed limits for different road types. In a
real-world scenario, this data might come from more detailed datasets
like Ordnance Survey MasterMap Highways Network or OpenStreetMap (OSM).
For this example, we’ll use a simplified set of speeds. We also define
speeds for congested conditions, assuming a reduction of 10 mph for
Motorways, A roads, and B roads. These speeds are converted to meters
per second, which is a common unit for network analysis calculations.

``` r
custom_wp <- tibble(
  road_function = c("Local Road", "Minor Road", "B Road", "A Road",
                    "Motorway"),
  speed_ff = round(c(20,20,30,40,70)*(0.44704),5),
  speed_cg = round(c(20,20,20,30,60)*(0.44704),5),
  )
```

The road network data used in this analysis is OS OpenRoads, a freely
available dataset from Ordnance Survey. We load the road links within
our defined study area (using the `zones_wkt` filter) and exclude links
classified as ‘access roads’ as these are typically not part of through
routes. We then join our custom speed information to the network data.

``` r
selected_network <- st_read(
  "00_data/oproad_gb.gpkg",
  wkt_filter = zones_wkt,
  query = "SELECT * FROM \"road_link\" WHERE road_function NOT LIKE '%access%'") |> 
  left_join(custom_wp,by = "road_function") |> 
    # With the speeds assigned, we can calculate travel times for each road segment.
    # We calculate travel times for both free-flow (speed-limit) and congested conditions.
  mutate(
    tra_time_ff = length/speed_ff, # Travel time in free-flow (seconds)
    tra_time_cg = length/speed_cg  # Travel time in congested conditions (seconds)
         )
```

    Reading query `SELECT * FROM "road_link" WHERE road_function NOT LIKE '%access%''
    from data source `C:\temp_jf\MW-rapid-cut-through\00_data\oproad_gb.gpkg' using driver `GPKG'
    Re-reading with feature count reset from 33542 to 29003
    Simple feature collection with 29003 features and 20 fields
    Geometry type: LINESTRING
    Dimension:     XY
    Bounding box:  xmin: 418542.3 ymin: 421364.6 xmax: 441664.4 ymax: 443596
    Projected CRS: OSGB36 / British National Grid

Let’s take a quick look at the road network we’ve prepared. The map
below shows the different types of roads in our Leeds study area.

    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    not found in Windows font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    not found in Windows font database

    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family not found in Windows font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family not found in Windows font database

![](README_files/figure-commonmark/unnamed-chunk-4-1.png)

To perform network analysis, we need to convert our spatial lines data
into a graph representation. The `sfnetworks` package is ideal for this,
as it integrates `sf` spatial objects with `tidygraph` network analysis
capabilities. We convert the `selected_network` into an `sfnetwork`
object.

A common step in network simplification is to remove ‘pseudo nodes’.
These are nodes that connect only two edges and don’t represent true
junctions or decision points in the network. Removing them simplifies
the graph and can make subsequent analyses more efficient and
meaningful. The `to_spatial_smooth` function helps achieve this by
merging edges connected by pseudo nodes, while also summarizing their
attributes (like summing lengths and travel times). We ensure that edges
are only merged if they have the same `road_function`.

``` r
net = as_sfnetwork(selected_network)

# Simplifying the network by removing pseudo nodes
net_simplified = convert(net,
                         to_spatial_smooth,
                         summarise_attributes = list(length = "sum",
                                                     tra_time_ff = "sum",
                                                     tra_time_cg = "sum",
                                                     "first" # Keep the first value for other attributes
                                                     ),
                         require_equal = "road_function") # Only merge edges with the same road function
```

The core of our analysis involves calculating ‘betweenness centrality’.
This metric quantifies how often a particular road segment (edge) lies
on the shortest paths between pairs of other locations in the network. A
higher betweenness centrality suggests a road is more critical for
connecting different parts of the network. We calculate this for two
scenarios: 1. **Baseline (Free-flow):** Using travel times based on
speed limits. 2. **Congested State:** Using travel times that reflect
congestion on major roads.

By comparing the centrality scores between these two states, we can
identify roads that become disproportionately more ‘central’ or
‘important’ when major roads are congested. These are potential
candidates for cut-through routes. We use a `cutoff` of 5 minutes (300
seconds) to limit the shortest path calculations, focusing on relatively
local diversions.

``` r
# Calculate edge betweenness centrality for free-flow conditions
E(net_simplified)$centrality_ff <- edge_betweenness(
  net_simplified,
  directed = F, # Undirected graph
  weights = E(net_simplified)$tra_time_ff, # Use free-flow travel times as weights
  cutoff = 5*60 # Consider paths up to 5 minutes
  )

# Calculate edge betweenness centrality for congested conditions
E(net_simplified)$centrality_cg <- edge_betweenness(
  net_simplified,
  directed = F, # Undirected graph
  weights = E(net_simplified)$tra_time_cg, # Use congested travel times as weights
  cutoff = 5*60 # Consider paths up to 5 minutes
  )

# Now, we calculate the difference in betweenness centrality.
# A positive difference means the road segment became more central during congestion.
# We also calculate a log-transformed difference to handle potential large variations and skewness.
net_centrality <- net_simplified |> 
  activate(edges) |> 
  mutate(diff_bc = centrality_cg - centrality_ff,
         log_diff = if_else(diff_bc!=0,
                            (diff_bc/abs(diff_bc))*log(abs(diff_bc)),0 # Log transform, preserving sign
                            ))
```

Let’s examine the distribution of these centrality changes, specifically
for local roads. The histogram below shows the log-transformed
difference in betweenness centrality. Values to the right of zero
indicate local roads that gained importance (centrality) under congested
conditions.

``` r
sf_edges <- st_as_sf(net_centrality, "edges") |> 
  mutate(road_function = factor(road_function,
                                levels = road_levels,
                                ordered = T))

sf_edges |>
    filter(road_function == "Local Road") |> 
    ggplot(aes(log_diff))+
    geom_histogram(col = "white")+
    theme_minimal()+
    labs(title = "Distribution of Centrality Change on Local Roads",
         x = "Log-transformed Difference in Betweenness Centrality (Congested - Free Flow)",
         y = "Frequency")
```

    `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](README_files/figure-commonmark/unnamed-chunk-7-1.png)

Under congested conditions, most local roads loose *importance* as the
major roads provide less connectivity. However, a portion of them become
more important as they provide quicker routes for drivers, potentially
indicating a cut-through corridor.

We can explore the location of such corridors. These are the links with
the highest values (97th percentile).

``` r
cut_map <- sf_edges |>
    # We are interested in the change on local roads
    mutate(log_diff = if_else(road_function %in% c("Local Road"),log_diff,NA_real_)) |>
    # Identify local roads in the top 3% of centrality increase
     mutate(log_diff_p95 = log_diff > quantile(log_diff,0.97,na.rm = T)) |> 
    ggplot()+
  geom_sf(aes(col = log_diff_p95,
              linewidth = road_function,
              alpha = log_diff_p95
              ))+
  theme_void()+
  labs(col = "Potential Cut-Through",
       title = "Potential Cut-Through Corridors on Local Roads in Leeds",
       caption = "Highlighted roads are local roads in the 97th percentile of increased betweenness centrality under congestion.")+
  scale_color_manual(values = c("gray60","dodgerblue"))+
  scale_linewidth_manual(values = 1/c(1.5,1.8,2.05,2.5,2.1),guide = 'none')+
  scale_alpha_manual(
    values = c(0.2,1),guide = 'none',na.value = 0.4
    )+
  guides(col=guide_legend(nrow=3,byrow=TRUE))+
  theme(text = element_text(family = "Roboto Condensed"), legend.position = "right")

cut_map
```

    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family not found in Windows font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family not found in Windows font database

    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    not found in Windows font database

    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family not found in Windows font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family not found in Windows font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family not found in Windows font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family not found in Windows font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family not found in Windows font database

![](README_files/figure-commonmark/unnamed-chunk-9-1.png)

``` r
ggsave(plot = cut_map,filename = "cut_map.png",dpi = 320,width = 24,height = 14,units = "cm")
```
