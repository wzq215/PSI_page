# from utils import *
import numpy as np
import plotly.graph_objects as go
import plotly
from datetime import datetime

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
fov_dt = datetime(2021, 4, 27, 12)

tracks = []

tracks.append({'r0': 11.3643, 'carrlon0': 63.305, 'carrlat0': -1.1469, 'v0': 296.49})
tracks.append({'r0': 8.8387, 'carrlon0': 104.845, 'carrlat0': -4.8701, 'v0': 241.51})
tracks.append({'r0': 8.2000, 'carrlon0': 87.856, 'carrlat0': -5.1566, 'v0': 288.29})
tracks.append({'r0': 9.1360, 'carrlon0': 78.567, 'carrlat0': -2.005, 'v0': 282.77})

config = {
  'toImageButtonOptions': {
    'format': 'svg', # one of png, svg, jpeg, webp
    'filename': 'custom_image',
    # 'height': 500,
    # 'width': 700,
    'scale': 2 # Multiply title/legend/axis/canvas sizes by this factor
  }
}

plot = go.Figure(data=[plot_Spacecraft_model('SPP', fov_dt, scale=5,color='silver',lighting=dict(specular=2))]
                      + plot_FOV('SPP', 'SPP_WISPR_INNER', fov_dt,
                                 fov_dist_Rs=25, plot_type='wire+surf', color='lime')
                      + sum([plot_flux_tube(track['carrlon0'], track['carrlat0'],
                                            np.arange(track['r0'], 25., .5),
                                            np.zeros_like(np.arange(track['r0'], 25., .5))+track['v0'],
                                            delta_deg=3.,
                                            deg_vs_r_vect=1. - 1. / (
                                                    abs(np.arange(track['r0'], 25., .5) - track['r0'])
                                                    * abs(np.arange(track['r0'], 25., .5) - track['r0'] - 3.) + 1)
                                            # deg_vs_r_vect=np.log(
                                            #     abs(np.linspace(track['r0'], 25., 20) - track['r0']) * abs(
                                            #         np.linspace(track['r0'], 25., 20) - track['r0'] - 3.) + 1)
                                            )
                             for track in tracks], [])
                      + sum([plot_flux_tube(track['carrlon0'], track['carrlat0'],
                                            np.arange(track['r0'], 1., -.5),
                                            np.zeros_like(np.arange(track['r0'], 1., -.5))+track['v0'],
                                            delta_deg=3.,
                                            deg_vs_r_vect=1. - 1. / (
                                                    abs(np.arange(track['r0'], 1., -.5) - track['r0'])
                                                    * abs(np.arange(track['r0'], 1., -.5) - track['r0'] - 3.) + 1)

                                            # deg_vs_r_vect=np.log(
                                            #     abs(np.linspace(track['r0'], 1., 20) - track['r0']) * abs(
                                            #         np.linspace(track['r0'], 1., 20) - track['r0'] - 3.) + 1)
                                            )
                             for track in tracks], [])
                      + [plot_object_orbit('SPP', dt_epoch, type='line', color='white', width=6),
                         plot_HCS(crid, 'corona', **surface_dict,lighting=dict(diffuse=0.9)),
                         plot_sun(crid, scale=1.5)],
                 layout=go.Layout(scene=dict(xaxis_title='X (Rs)', yaxis_title='Y (Rs)', zaxis_title='Z (Rs)',
                                             xaxis_range=[-30, 30], yaxis_range=[-30, 30], zaxis_range=[-30, 30],
                                             aspectratio=dict(x=1, y=1, z=1)),
                                  template='plotly_dark',
                                  showlegend=False,)
                                  # updatemenus=[dict(type='buttons',
                                  #                   buttons=[dict(label='Play', method='animate',
                                  #                                 args=[None])])])
                 )
from PIL import Image

# img = Image.open('/Users/ephe/Desktop/IMG_8756.JPG')
# plot.add_trace(go.Image(z=img),1,2)
plotly.offline.plot(plot,config=config)

# import io
# frames = []
# for s, fr in enumerate(plot.frames):
#     # set main traces to appropriate traces within plotly frame
#     plot.update(data=fr.data)
#     # move slider to correct place
#     plot.layout.sliders[0].update(active=s)
#     # generate image of current state
#     frames.append(Image.open(io.BytesIO(plot.to_image(format="png"))))
#
# # append duplicated last image more times, to keep animation stop at last status
# for i in range(5):
#     frames.append(frames[-1])
#
# # create animated GIF
# frames[0].save(
#         "test.gif",
#         save_all=True,
#         append_images=frames[1:],
#         optimize=False,
#         duration=500,
#         loop=0,
#     )
