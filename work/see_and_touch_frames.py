# from utils import *
import os

import plotly.graph_objects as go
import plotly
from datetime import datetime, timedelta

# MY MODULES
from plot_elements import *
from utils import create_epoch

# %%
wispr_epoch = np.load('./inner_epoch.npy', allow_pickle=True)
sub_wispr_epoch = wispr_epoch[::10]

# %%
crid = 2243
surface_dict = {'color': 'purple', 'opacity': 0.5, 'param': None}
dt = datetime(2021, 4, 27)
# dt_epoch = create_epoch([datetime(2021, 4, 27), datetime(2021, 4, 30)], timedelta(hours=4))

dt_epoch = sub_wispr_epoch
tracks = []
tracks.append({'r0': 11.3643, 'carrlon0': 63.305, 'carrlat0': -1.1469, 'v0': 296.49,
               't0': datetime(2021, 4, 27, 13, 12, 20)})  # TRACK 2
tracks.append({'r0': 8.8387, 'carrlon0': 104.845, 'carrlat0': -4.8701, 'v0': 241.51,
               't0': datetime(2021, 4, 27, 17, 12, 20)})  # TRACK 3
tracks.append({'r0': 8.2000, 'carrlon0': 87.856, 'carrlat0': -5.1566, 'v0': 288.29,
               't0': datetime(2021, 4, 27, 3, 42, 20)})  # TRACK 1b
tracks.append({'r0': 9.1360, 'carrlon0': 78.567, 'carrlat0': -2.005, 'v0': 282.77,
               't0': datetime(2021, 4, 27, 3, 42, 20)})  # TRACK 1a

track_t0_lst = [datetime(2021, 4, 27, 3, 42, 20),
                datetime(2021, 4, 27, 3, 42, 20),
                datetime(2021, 4, 27, 13, 12, 20),
                datetime(2021, 4, 27, 17, 12, 20)]

EXPORT_PATH = 'export/work/ST_frames/'
os.makedirs(EXPORT_PATH,exist_ok=True)

for i in range(len(dt_epoch)):
    dt_tmp = dt_epoch[i]

    plot = go.Figure(data=[plot_Spacecraft_model('SPP', dt_tmp, scale=5, color='silver', lighting=dict(specular=2))]
                          + plot_FOV('SPP', 'SPP_WISPR_INNER', dt_tmp,
                                     fov_dist_Rs=25, plot_type='wire+surf', color='lime')
                          + [plot_object_orbit('SPP', dt_epoch, type='line', color='white', width=6),
                             plot_HCS(crid, 'corona', **surface_dict, lighting=dict(diffuse=0.9)),
                             plot_sun(crid, scale=1.5)],
                     layout=go.Layout(title=dict(text=dt_tmp.strftime('%Y%m%d %H:%M:%S'),
                                                 x=.5,y=.9,xanchor='center',yanchor='top'),
                                      autosize=False,
                                      width=800,height=500,margin=dict(l=20,r=20,b=20,t=20),
                                      scene=dict(xaxis_title='X (Rs)', yaxis_title='Y (Rs)', zaxis_title='Z (Rs)',
                                                 xaxis_range=[-30, 30], yaxis_range=[-30, 30], zaxis_range=[-30, 30],
                                                 aspectratio=dict(x=1, y=1, z=1),
                                                 camera=dict(
                                                     up=dict(x=0, y=0, z=1),
                                                     center=dict(x=0, y=0, z=0),
                                                     eye=dict(x=.5, y=.8, z=.3))),
                                      template='plotly_dark',
                                      showlegend=False, ))
    for track in tracks:
        r_tmp = track['r0'] + track['v0'] * (dt_tmp - track['t0'])/timedelta(seconds=1)/Rs2km
        print(r_tmp)
        # if r_tmp>1. and r_tmp<30.:

        plot.add_traces(plot_flux_tube(track['carrlon0'], track['carrlat0'],
                                       np.arange(track['r0'], 30., .5),
                                       np.zeros_like(np.arange(track['r0'], 30., .5)) + track['v0'],
                                       delta_deg=3.,
                                       deg_vs_r_vect=1. - 1. / (
                                               abs(np.arange(track['r0'], 30., .5) - r_tmp)
                                               * abs(np.arange(track['r0'], 30., .5) - r_tmp - 3.) + 1)))
        plot.add_traces(plot_flux_tube(track['carrlon0'], track['carrlat0'],
                                                np.arange(track['r0'], 1., -.5),
                                                np.zeros_like(np.arange(track['r0'], 1., -.5)) + track['v0'],
                                                delta_deg=3.,
                                                deg_vs_r_vect=1. - 1. / (
                                                        abs(np.arange(track['r0'], 1., -.5) - r_tmp)
                                                        * abs(np.arange(track['r0'], 1., -.5) - r_tmp - 3.) + 1)))
        # else:
        #     plot.add_traces(plot_flux_tube(track['carrlon0'], track['carrlat0'],
        #                                    np.arange(1., 30., .5),
        #                                    np.zeros_like(np.arange(1., 30., .5)) + track['v0'],
        #                                    delta_deg=3.,))
    plot.write_image(EXPORT_PATH+'ST_frame_'+str(i)+'.png',scale=3.5)
    # plotly.offline.plot(plot)

