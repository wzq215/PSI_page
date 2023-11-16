import plotly.graph_objects as go
from datetime import datetime, timedelta
import plotly

from plot_elements import *
from utils import create_epoch

crid = 2269

dt_epoch = create_epoch([datetime(2023,3,25),datetime(2023,4,1)],timedelta(hours=3))

plot = go.Figure(data = [plot_Spacecraft_model('SOLO', datetime(2023,3,28), scale=50,color='silver'),
                         plot_planet_model('earth',datetime(2023,3,28),radius=10)]
                        +sum([plot_flux_tube_from_object('SOLO', datetime(2023,3,28),color='white')],[])
                        +sum([plot_flux_tube_from_object('earth', datetime(2023,3,28),color='white')],[])
                        +[plot_sun(crid=crid,surface_data='AIA_193',scale=10),
                         plot_object_orbit('SOLO', dt_epoch, type='line', color='red', width=6),
                         # plot_HCS(crid, 'corona', color='purple', opacity=0.5, param=None),
                         plot_HCS(crid, 'helio', color='purple', opacity= 0.5, param= None),
                         plot_zslice(crid,'corona',param='vr',z_const=0.,colorscale='plasma',clim=[200,400]),
                         # plot_slice_sun_object(crid,'helio',param='vr',object_name='earth',object_dt=datetime(2023,3,28),z_offset=10.,colorscale='jet',clim=[200,400]),
                         plot_object_orbit('earth',dt_epoch,type='line',color='blue',width=6),
                         ],
                 frames = [go.Frame(data=[plot_Spacecraft_model('SOLO', dt_tmp, scale=50, color='silver'),
                                          plot_planet_model('earth',dt_tmp,radius=10),]
                                         +sum([plot_flux_tube_from_object('SOLO', dt_tmp,color='white')],[])
                        +sum([plot_flux_tube_from_object('earth', dt_tmp,color='white')],[]))
                           for dt_tmp in dt_epoch[::6]],
                 layout=go.Layout(scene=dict(xaxis_title='X (Rs)', yaxis_title='Y (Rs)', zaxis_title='Z (Rs)',
                                             # xaxis_range=[-215, 215], yaxis_range=[-215, 215], zaxis_range=[-215, 215],
                                             aspectratio=dict(x=1, y=1, z=1)),
                                  template='plotly_dark',
                                  updatemenus=[dict(type='buttons',
                                                    buttons=[dict(label='Play', method='animate', args=[None])])]))
plotly.offline.plot(plot)
