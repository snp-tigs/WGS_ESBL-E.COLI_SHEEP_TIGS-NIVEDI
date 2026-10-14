"""
==============================================================================
Title: Figure 1 - Study Area Map (Karnataka, Chitradurga, and Tumakuru)
Description: Generates a geographic map of the sampling locations across the 
             Chitradurga and Tumakuru districts in Karnataka. Projects coordinates 
             to Web Mercator (EPSG:3857) for accurate spatial scaling and 
             exports as a transparent vector SVG for publication layouts.
              Script accompanies the publication: Whole-genome analysis reveals 
              shared extended-spectrum β-lactamase-producing Escherichia coli 
              lineages between migratory sheep and surrounding farm environments
==============================================================================
"""

import os
import geopandas as gpd
import matplotlib.pyplot as plt
import contextily as ctx
from pyproj import Transformer

# ==========================================
# 1. Configuration & File Paths
# ==========================================
# Relative paths for GitHub reproducibility
KARNATAKA_GEOJSON   = "./data_main_figures/Figure_1/karnataka.geojson"
CHITRADURGA_GEOJSON = "./data_main_figures/Figure_1/chitradurga_only.geojson"
TUMAKURU_GEOJSON    = "./data_main_figures/Figure_1/TUMAKURU_only.geojson"

# Output directories and files
RESULTS_DIR = "./results/Figure_1"
OUTPUT_SVG  = os.path.join(RESULTS_DIR, "Fig_1_Map_Overlay.svg")

# Create results directory if it doesn't exist
os.makedirs(RESULTS_DIR, exist_ok=True)

# ==========================================
# 2. Load and Project Spatial Data
# ==========================================
print("Loading GeoJSON boundaries...")
# Load State and District boundaries
gdf_karnataka   = gpd.read_file(KARNATAKA_GEOJSON)
gdf_chitradurga = gpd.read_file(CHITRADURGA_GEOJSON)
gdf_tumakuru    = gpd.read_file(TUMAKURU_GEOJSON)

# Reproject all data to Web Mercator projection (EPSG:3857) for accurate map rendering
gdf_karnataka   = gdf_karnataka.to_crs(epsg=3857)
gdf_chitradurga = gdf_chitradurga.to_crs(epsg=3857)
gdf_tumakuru    = gdf_tumakuru.to_crs(epsg=3857)

# Initialize transformer to convert Sample Coordinates from EPSG:4326 (lat/lon) to EPSG:3857
transformer = Transformer.from_crs("epsg:4326", "epsg:3857", always_xy=True)

# ==========================================
# 3. Define Sampling Sites
# ==========================================
sample_sites = [
    # --- Chitradurga Sites ---
    {'latitude': 14.2178, 'longitude': 76.6454, 'label': 'Ganjigunte'},
    {'latitude': 14.2161, 'longitude': 76.6595, 'label': 'Sondekere'},
    {'latitude': 14.3339, 'longitude': 76.6919, 'label': 'Bharamasagara'},
    {'latitude': 14.2954, 'longitude': 76.7044, 'label': 'Dodderi'},
    {'latitude': 14.1743, 'longitude': 76.8447690854118, 'label': 'Doddabeeranahalli'},
    {'latitude': 14.5310, 'longitude': 76.7643778669348, 'label': 'Gowrasamudra Kaval'},
    {'latitude': 14.2588, 'longitude': 76.8788, 'label': 'Parashuramapura kaval'},
    {'latitude': 14.2867, 'longitude': 76.7981010905807, 'label': 'Purlahalli'},
    {'latitude': 14.2904, 'longitude': 76.54589545, 'label': 'Ramjogihally'},
    {'latitude': 14.3683, 'longitude': 76.8037, 'label': 'Kalavehalli'},
    {'latitude': 14.4292, 'longitude': 76.8948, 'label': 'Donehalli'},
    {'latitude': 14.4660, 'longitude': 76.9319, 'label': 'Obalapura'},
    {'latitude': 13.8983, 'longitude': 76.6766, 'label': 'Bagganandu Kaval'},
    {'latitude': 13.8860, 'longitude': 76.5590, 'label': 'Paremenahalli'},
    {'latitude': 14.0422, 'longitude': 76.7705, 'label': 'Chillahalli'},
    {'latitude': 14.0786, 'longitude': 76.7710, 'label': 'Gulya'},
    {'latitude': 13.9938, 'longitude': 76.5190, 'label': 'Sooragondanahalli'},
    {'latitude': 14.2019, 'longitude': 76.4220, 'label': 'Kadehude'},
    {'latitude': 14.2016, 'longitude': 76.4802, 'label': 'Muddapura'},
    {'latitude': 13.7198, 'longitude': 76.4427, 'label': 'Bukkasagara'},
    {'latitude': 13.8772, 'longitude': 76.4187, 'label': 'Mallapura'},
    {'latitude': 13.9306, 'longitude': 76.3549, 'label': 'Janakal'},
    
    # --- Tumakuru Sites ---
    {'latitude': 13.8780, 'longitude': 76.7575, 'label': 'Dandikere'},
    {'latitude': 13.7288, 'longitude': 76.7368, 'label': 'Nerlegudda'},
    {'latitude': 13.9095, 'longitude': 76.7598, 'label': 'Herur'},
    {'latitude': 13.8116, 'longitude': 76.9119, 'label': 'Maddakkanahalli'},
    {'latitude': 14.0960, 'longitude': 77.2758, 'label': 'Pavgada ward no 1'},
    {'latitude': 14.1014, 'longitude': 77.2938, 'label': 'Pavgada ward no 2'},
    {'latitude': 14.18553,'longitude': 77.22837,'label': 'Ponnasamudra'},
    {'latitude': 13.43584,'longitude': 77.37192,'label': 'Kilarahalli'},
    {'latitude': 13.2087, 'longitude': 76.6578, 'label': 'Lokkamanahalli'},
    {'latitude': 13.1314, 'longitude': 76.7062, 'label': 'Kallankere'},
    {'latitude': 13.1943, 'longitude': 76.6123, 'label': 'Anekere'},
    {'latitude': 13.1700, 'longitude': 76.7778, 'label': 'Talakere'},
    {'latitude': 13.1059, 'longitude': 76.7850, 'label': 'Doddamalligere'},
    {'latitude': 13.3414, 'longitude': 77.1054, 'label': 'Tumkur ward no-9'},
    {'latitude': 13.4475, 'longitude': 77.1283, 'label': 'Thirumalapalya'},
    {'latitude': 13.3373, 'longitude': 77.1782, 'label': 'Janapanahalli'},
    {'latitude': 13.2342, 'longitude': 76.5527, 'label': 'Shivapura'},
    {'latitude': 13.2650, 'longitude': 76.9421, 'label': 'Hosahalli G'},
    {'latitude': 13.2987, 'longitude': 76.6864, 'label': 'Kattigenahalli'},
    {'latitude': 13.6202, 'longitude': 77.1868, 'label': 'Kammanakote'},
    {'latitude': 13.8018, 'longitude': 77.3551, 'label': 'Hosa Itakaloti'},
    {'latitude': 13.7187, 'longitude': 77.3497, 'label': 'Kalenahalli'},
    {'latitude': 13.6411, 'longitude': 77.0457, 'label': 'Avinamadu'}
]

# ==========================================
# 4. Rendering the Map
# ==========================================
print("Rendering map...")
fig, ax = plt.subplots(figsize=(20, 25))

# Plot Karnataka Base Map
gdf_karnataka.plot(ax=ax, linewidth=1.5, edgecolor='black', color='#93A79D', alpha=1, zorder=1)

# Plot Highlighted Districts (Chitradurga and Tumakuru)
gdf_chitradurga.plot(ax=ax, linewidth=2, edgecolor='#2c3e50', color='#A8BDB1', alpha=0.9, zorder=2)
gdf_tumakuru.plot(ax=ax, linewidth=2, edgecolor='#2c3e50', color='#A8BDB1', alpha=0.9, zorder=2)

# Plot Sample Coordinates
for site in sample_sites:
    x, y = transformer.transform(site['longitude'], site['latitude'])
    ax.plot(x, y, 'o', color='black', markeredgewidth=0.8, markeredgecolor='blue', markersize=8, zorder=3)
    # Uncomment the line below to render text labels for the villages directly on the map
    # ax.text(x, y, site['label'], fontsize=6, ha='right', va='bottom', fontweight='bold', color='black', rotation=-45, zorder=4)

# Optional: Add a contextual basemap (OpenStreetMap / CartoDB) underneath if desired
# ctx.add_basemap(ax, crs=gdf_karnataka.crs.to_string(), source=ctx.providers.CartoDB.PositronNoLabels, alpha=0.5, zorder=0)

# Format Axes
ax.set_title('Study Area: Sample Locations in Chitradurga and Tumakuru', fontsize=18, fontweight='bold', pad=20)
ax.set_xlabel('Longitude (Web Mercator)', fontsize=14)
ax.set_ylabel('Latitude (Web Mercator)', fontsize=14)

# Remove axis spines for a cleaner aesthetic
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# ==========================================
# 5. Export
# ==========================================
plt.savefig(OUTPUT_SVG, dpi=600, transparent=True, format="svg", bbox_inches="tight")
print(f"✅ Vector map successfully saved to: {OUTPUT_SVG}")
