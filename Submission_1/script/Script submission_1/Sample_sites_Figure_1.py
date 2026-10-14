import geopandas as gpd
import matplotlib.pyplot as plt
import contextily as ctx
from pyproj import Transformer

#GeoJSON file
gdf = gpd.read_file('/path/to/folder/karnataka.geojson')

# Reprojecting the data to Web Mercator projection (EPSG:3857)
gdf = gdf.to_crs(epsg=3857)

# Initialize transformer to convert from EPSG:4326 (lat/lon) to EPSG:3857
transformer = Transformer.from_crs("epsg:4326", "epsg:3857", always_xy=True)


fig, ax = plt.subplots(figsize=(20, 25))

#--------------------chitradurga---------------------------------
sample_sites = [
   #chitradurge
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

    {'latitude': 13.8860, 'longitude': 76.559, 'label': 'Paremenahalli'},
    {'latitude': 14.0422, 'longitude': 76.7705, 'label': 'Chillahalli'},
    {'latitude': 14.0786, 'longitude': 76.771, 'label': 'Gulya'},
    {'latitude': 13.9938, 'longitude': 76.519, 'label': 'Sooragondanahalli'},
    {'latitude': 14.2019, 'longitude': 76.422, 'label': 'Kadehude'},
    {'latitude': 14.2016, 'longitude': 76.4802, 'label': 'Muddapura'},
    {'latitude': 13.7198, 'longitude': 76.4427, 'label': 'Bukkasagara'},
    {'latitude': 13.8772, 'longitude': 76.4187, 'label': 'Mallapura'},
    {'latitude': 13.9306, 'longitude': 76.3549, 'label': 'Janakal'},
    #Tumkur
    {'latitude': 13.8780, 'longitude': 76.7575, 'label': 'Dandikere'},
    {'latitude': 13.7288, 'longitude': 76.7368, 'label': 'Nerlegudda'},
    {'latitude': 13.9095, 'longitude': 76.7598, 'label': 'Herur'},
    {'latitude': 13.8116, 'longitude': 76.9119, 'label': 'Maddakkanahalli'},

    {'latitude': 14.096, 'longitude': 77.2758, 'label': 'Pavgada ward no 1'},
    {'latitude': 14.1014, 'longitude': 77.2938, 'label': 'Pavgada ward no 2'},
    {'latitude': 14.18553, 'longitude': 77.22837, 'label': 'Ponnasamudra'},
    {'latitude': 13.43584, 'longitude': 77.37192, 'label': 'Kilarahalli'},
    {'latitude': 13.2087, 'longitude': 76.6578, 'label': 'Lokkamanahalli'},
    {'latitude': 13.1314, 'longitude': 76.7062, 'label': 'Kallankere'},
    {'latitude': 13.1943, 'longitude': 76.6123, 'label': 'Anekere'},
    {'latitude': 13.17, 'longitude': 76.7778, 'label': 'Talakere'},
    {'latitude': 13.1059, 'longitude': 76.785, 'label': 'Doddamalligere'},
    {'latitude': 13.3414, 'longitude': 77.1054, 'label': 'Tumkur ward no-9'},
    {'latitude': 13.4475, 'longitude': 77.1283, 'label': 'Thirumalapalya'},
    {'latitude': 13.3373, 'longitude': 77.1782, 'label': 'Janapanahalli'},
    {'latitude': 13.2342, 'longitude': 76.5527, 'label': 'Shivapura'},

    {'latitude': 13.265, 'longitude': 76.9421, 'label': 'Hosahalli G'},
    {'latitude': 13.2987, 'longitude': 76.6864, 'label': 'Kattigenahalli'},
    {'latitude': 13.6202, 'longitude': 77.1868, 'label': 'Kammanakote'},
    {'latitude': 13.8018, 'longitude': 77.3551, 'label': 'Hosa Itakaloti'},
    {'latitude': 13.7187, 'longitude': 77.3497, 'label': 'Kalenahalli'},
    {'latitude': 13.6411, 'longitude': 77.0457, 'label': 'Avinamadu'},


]

# Plot the GeoDataFrame with custom styling
#gdf.plot(ax=ax, linewidth= 1, edgecolor='black', color='#93A79D', alpha=1)

#karnataka_color
gdf.plot(ax=ax, linewidth= 1, edgecolor='black', color='#93A79D', alpha=1)
"""
# Convert lat/lon to EPSG:3857 and plot
for site in sample_sites:
    x, y = transformer.transform(site['longitude'], site['latitude'])
    ax.plot(x, y, 'o', color='black', markeredgewidth=0.2, markeredgecolor='blue', markersize=7)
    #ax.text(x, y, site['label'], fontsize=5, ha='right', va='bottom', fontweight='bold', color='red', rotation=-45)

# Add a basemap without labels (CartoDB Positron without labels)
#ctx.add_basemap(ax, crs=gdf.crs.to_string(), source=ctx.providers.CartoDB.VoyagerNoLabels, alpha=1)
ctx.add_basemap(ax, crs=gdf.crs.to_string(), source=ctx.providers.OpenStreetMap.CH, alpha=0)
"""
# Set labels and title
ax.set_title(' ')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

# Show plot
#plt.show()

# Save plot as PNG
plt.savefig ('/path/to/folder/CH_TUMKUR_OVERLAY_.svg', dpi=600, transparent=True)
