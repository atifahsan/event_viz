import os
import re
import math
import datetime as dt
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
import matplotlib.colors as colors


STATS_PATH = './stats'
MI_FILE = './tab_files/LTEDB_qgis.tab'

FIXED_COLS = ['Date', 'Time', 'Datetime', 'eNodeB Name', 'Cell FDD TDD Indication', 'Cell Name', 'LocalCell Id', 
              'eNodeB Function Name', 'Integrity', 'Start Time', 'Period(min)', 'Period(min)', 'NE Name', 'Cell']




def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map
       https://gist.github.com/jakevdp/91077b0cae40f8f8244a
    """
    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)


def get_end_point(lon1,lat1,bearing,d):
    R = 6371                     #Radius of the Earth
    brng = math.radians(bearing) #convert degrees to radians
    d = d/1000                  #convert meters to km
    lat1 = math.radians(lat1)    #Current lat point converted to radians
    lon1 = math.radians(lon1)    #Current long point converted to radians
    lat2 = math.asin( math.sin(lat1)*math.cos(d/R) + math.cos(lat1)*math.sin(d/R)*math.cos(brng))
    lon2 = lon1 + math.atan2(math.sin(brng)*math.sin(d/R)*math.cos(lat1),math.cos(d/R)-math.sin(lat1)*math.sin(lat2))
    lat2 = math.degrees(lat2)
    lon2 = math.degrees(lon2)
    return (lon2,lat2)


def lte_sec_from_cellname(name):
    pattern = r'_[L,T]\d_[A-F,O,X,Y]'
    sec = re.findall(pattern, name)
    if len(sec) == 1:
        sector = sec[0][-1:]
    else:
        sector = 'O'
    return sector


def convertnametoid(nodebname):
    name = nodebname
    nodebid = re.findall(r'\d{4,5}', name)
    if len(nodebid) > 0:
        siteid = int(nodebid[0])
    else:
        siteid = '0000'
    return str(siteid)


def read_u2020_stats(path):
    df = pd.DataFrame()
    for fname in os.listdir(path):
        if fname.startswith('Export'):
            print(fname)
            temp = pd.read_csv(os.path.join(path, fname), parse_dates=['Start Time'], na_values='NIL', compression='gzip', skiprows=6)
            df = df.append(temp)
    df.rename(columns={'Start Time': 'Datetime', 'NE Name': 'eNodeB Name', 'Cell': 'Cell Name'}, inplace=True)
    df.columns = df.columns.str.replace(' \(None\)', '')

    agg_dict = {col: 'sum' for col in df.columns if col not in FIXED_COLS}

    return df, agg_dict


def read_tab_file(path):
    gdf = gpd.read_file(MI_FILE, driver='Mapinfo File')
    gdf['Sector'] = gdf['Sector'].apply(lambda x: str(x[:1]).upper())
    gdf['Site'] = gdf['cellname'].apply(lambda x: convertnametoid(str(x)))
    gdf['Site_Sec'] = gdf['Site'].astype(str) + "_" + gdf['Sector']
    gdf = gdf.drop_duplicates(subset=['Site_Sec', 'geometry'])
    gdf = gdf[['Site', 'Site_Sec', 'geometry', 'Sector', 'Azimuth', 'Latitude', 'Longitude']]
    return gdf


def get_arcgis_map(geo_df, ax):
    buffer = 0.01
    map_retreived = False
    lon_min,lat_min,lon_max,lat_max = geo_df.geometry.buffer(buffer).total_bounds
    ax.set(xlim=[lon_min, lon_max], ylim=[lat_min, lat_max]) 
    
    while not map_retreived:
        try:
            lon_min,lat_min,lon_max,lat_max = geo_df.geometry.buffer(buffer).total_bounds
            bg_map = Basemap(llcrnrlon=lon_min, llcrnrlat=lat_min, urcrnrlon=lon_max, urcrnrlat=lat_max, epsg=4326)
            bg_map.arcgisimage(service='Canvas/World_Dark_Gray_Base', xpixels=1500, verbose=True, ax=ax)
            map_retreived = True
        except:
            buffer * 2


def visualize_kpi():
    bounds = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=10)

    gdf = read_tab_file(MI_FILE)

    df, agg_dict = read_u2020_stats(STATS_PATH)    
    df['Datetime'] = df['Datetime'].dt.strftime('%Y-%m-%d %H:%M:%S')
    df['Site'] = df['eNodeB Name'].apply(lambda x: convertnametoid(str(x)))
    df['Sec'] = df['Cell Name'].apply(lambda x: lte_sec_from_cellname(str(x)))
    df['Site_Sec'] = df['Site'] + "_" +  df['Sec']
    df = df.groupby(['Datetime', 'Site_Sec']).agg(agg_dict).reset_index()
    df['DL PRB Util'] = 100*df['L.ChMeas.PRB.DL.Used.Avg']/(df['L.ChMeas.PRB.DL.Avail']+0.0001)

    geo_df = pd.merge(gdf, df, on='Site_Sec')
    geo_df = geo_df.dropna()
    geo_df.crs = {'init':'epsg:4326'}



    time_ = list(geo_df['Datetime'].unique())
    time_.sort()

    f, ax = plt.subplots(1,1, figsize=(8, 8), facecolor='#232227')
    f.suptitle('NYE2020\nDL PRB Utilization%', size=20, color='whitesmoke')
    ax.set_axis_off()
    plt.axis('equal')
    plt.subplots_adjust(left=0.0, right=.95, top=0.9, bottom=0.0)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.1)
    cax.tick_params(color='whitesmoke', labelcolor='whitesmoke')


    get_arcgis_map(geo_df, ax)


    def init():
        plot_df = geo_df[(geo_df['Datetime'] == time_[0])]
        plot_df.plot(ax=ax, column='DL PRB Util', cmap=discrete_cmap(10, 'RdYlGn_r'), norm=norm, edgecolor='black', linewidth=.2, vmin=0, vmax=100, alpha=1, legend=True, cax=cax)
        ax.set_title(f'Time: {time_[t]}, Total Users = {plot_df["L.Traffic.User.Max"].sum():.2f}', color='whitesmoke')
        plot_df.apply(lambda x: ax.annotate(s=x.Sector, xy=get_end_point(*x.geometry.centroid.coords[0],x.Azimuth,200), size=5, va="center", ha="center", color='whitesmoke', alpha=0.55), axis=1);
        site_df = plot_df.copy(deep=True)
        site_df.drop_duplicates(subset=['Site'], inplace=True)
        for idx, row in site_df.iterrows():
            ax.annotate(s=row['Site'], xy=(row['Longitude'], row['Latitude']), size=6, va="center", ha="center", color='lightgrey', alpha=0.9)


    t = 1
    def update(t):
        plot_df = geo_df[(geo_df['Datetime'] == time_[t])]
        plot_df.plot(ax=ax, column='DL PRB Util', cmap=discrete_cmap(10, 'RdYlGn_r'), norm=norm, edgecolor='black', linewidth=.2, vmin=0, vmax=100, alpha=1, legend=False)
        ax.set_title(f'Time: {time_[t]}, Total Users = {plot_df["L.Traffic.User.Max"].sum():.2f}')


    animation = FuncAnimation(f, func=update, init_func=init, frames=len(time_), interval=1000)
    animation.save("./animation/viz.mp4", dpi=240, fps=1, codec="libx264", bitrate=5000, extra_args=['-pix_fmt', 'yuv420p'], savefig_kwargs={'facecolor' : f.get_facecolor()})
    animation.save('./animation/viz.gif', writer='pillow', savefig_kwargs={'facecolor' : f.get_facecolor()})
    plt.show()



if __name__ == "__main__":
    visualize_kpi()