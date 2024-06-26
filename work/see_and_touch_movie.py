# from utils import *
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
dt_epoch = create_epoch([datetime(2021, 4, 27), datetime(2021, 4, 30)], timedelta(hours=4))

dt_epoch = sub_wispr_epoch
tracks = []
tracks.append({'r0': 11.3643, 'carrlon0': 63.305, 'carrlat0': -1.1469, 'v0': 296.49})
tracks.append({'r0': 8.8387, 'carrlon0': 104.845, 'carrlat0': -4.8701, 'v0': 241.51})
tracks.append({'r0': 8.2000, 'carrlon0': 87.856, 'carrlat0': -5.1566, 'v0': 288.29})
tracks.append({'r0': 9.1360, 'carrlon0': 78.567, 'carrlat0': -2.005, 'v0': 282.77})

plot = go.Figure(data=[plot_Spacecraft_model('SPP', dt_epoch[0], scale=5)]
                      + plot_FOV('SPP', 'SPP_WISPR_INNER', dt_epoch[0],
                                 fov_dist_Rs=25, plot_type='wire+surf', color='lime')
                      + sum([plot_flux_tube(track['carrlon0'], track['carrlat0'], np.linspace(track['r0'], 25., 20),
                                            np.linspace(track['v0'], track['v0'], 20), ) for track in tracks], [])
                      + sum([plot_flux_tube(track['carrlon0'], track['carrlat0'], np.linspace(track['r0'], 1., 20),
                                            np.linspace(track['v0'], track['v0'], 20), ) for track in tracks], [])
                      + [plot_object_orbit('SPP', dt_epoch, type='line', color='white', width=6),
                         plot_HCS(crid, 'corona', **surface_dict), ],
                 frames=[go.Frame(data=[plot_Spacecraft_model('SPP', dt_tmp, scale=5)]
                                       + plot_FOV('SPP', 'SPP_WISPR_INNER', dt_tmp,
                                                  fov_dist_Rs=25, plot_type='wire+surf', color='lime'))
                         for dt_tmp in dt_epoch],
                 layout=go.Layout(scene=dict(xaxis_title='X (Rs)', yaxis_title='Y (Rs)', zaxis_title='Z (Rs)',
                                             xaxis_range=[-30, 30], yaxis_range=[-30, 30], zaxis_range=[-30, 30],
                                             aspectratio=dict(x=1, y=1, z=1)),
                                  template='plotly_dark',
                                  updatemenus=[dict(type='buttons',
                                                    buttons=[dict(label='Play', method='animate',
                                                                  args=[None])])])).set_subplots(1, 2)
from PIL import Image

# img = Image.open('/Users/ephe/Desktop/IMG_8756.JPG')
# plot.add_trace(go.Image(z=img),1,2)
plotly.offline.plot(plot)

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
