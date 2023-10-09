'''
To plot PSI 3d plots.
Elements:
    - HCS
        - SC/IH
        - Surface data
    - Sun
        - Surface data
    - Objects
        - Planet [EARTH/VENUS/MERCURY] or Spacecrafts [PSP/SO/WIND]
            - Orbit
            - Point
            - B trace
            - Surface data or STL
    - Scene
Coordinates
    - Carrington
Platform
    - Pyvista
    - Plotly
'''

from utils import  *
from plot_elements import *
from get_traces import *

#@st.cache_data(ttl=3600)
def plot_ensemble(input_dict, plot):
    if len(input_dict['HCS']) != 0:
        HCS_dict = input_dict['HCS']
        trace_hcs_list = get_HCS_traces(HCS_dict)
        for i in range(len(trace_hcs_list)):
            plot.add_trace(trace_hcs_list[i])

    if len(input_dict['Spacecraft']) != 0:
        SC_dict = input_dict['Spacecraft']
        for SC_str in SC_dict['SC_list']:
            tmp_SC_dict = SC_dict[SC_str]
            trace_SC_list = get_SC_traces(tmp_SC_dict)
            for trace in trace_SC_list:
                plot.add_trace(trace)
    if len(input_dict['Sun']) != 0:
        Sun_dict = input_dict['Sun']
        trace_sun = get_sun_trace(Sun_dict)
        plot.add_trace(trace_sun)
    if len(input_dict['Planet']) != 0:
        Planet_dict = input_dict['Planet']
        trace_planet_list = get_planet_traces(Planet_dict)
        for trace in trace_planet_list:
            plot.add_trace(trace)
    plot_layout(plot)



def plot_layout(plot):
    camera = dict(
        up=dict(x=0, y=0, z=1),
        center=dict(x=0, y=0, z=0),
        eye=dict(x=0, y=-1.5, z=1.5)
    )
    plot.update_layout(
        autosize=False,
        width=1500,
        height=1000,
        scene=dict(
            xaxis=dict(showgrid=False, zeroline=False,visible=False),
            yaxis=dict(showgrid=False, zeroline=False,visible=False),
            zaxis=dict(showgrid=False, zeroline=False,visible=False),
            bgcolor='black',
            # xaxis_title='X (Rs)',
            # yaxis_title='Y (Rs)',
            # zaxis_title='Z (Rs)',
            xaxis_range=[-250, 250],
            yaxis_range=[-250, 250],
            zaxis_range=[-250, 250],
        ),
        plot_bgcolor='rgb(255,255,255)',
        scene_camera=camera,
        # title=dict(
        #     text='Carringtong Coordinate',
        #     font_size=30,
        #     y=0.9, x=0.5,
        #     xanchor='center',
        #     yanchor='top'),
        template='plotly_dark',
        showlegend=False,
        scene_aspectratio=dict(x=1, y=1, z=1),
    )
    plot.update_scenes(xaxis_tickfont_size=15,
                       yaxis_tickfont_size=15,
                       zaxis_tickfont_size=15,
                       xaxis_title_font_size=25,
                       yaxis_title_font_size=25,
                       zaxis_title_font_size=25,

                       )


if __name__ == '__main__':
    # HCS

    '''parameters for plotting Helio Current Sheet'''
    is_plot_HCS = 1
    if is_plot_HCS == 1:
        HCS_surface_dict = {'param': 'vr',
                            'cmap': 'jet', 'clim': [200, 500], 'opacity': .6,
                            'color': 'purple', 'showscale': False}
        HCS_dict = {'crid': 2233, 'plot_sc': False, 'plot_ih': True, 'surface': HCS_surface_dict}
    else:
        HCS_dict = {}

    '''parameters for plotting planets'''
    Planet_list = ['EARTH']
    Planet_epoch = create_epoch([datetime(2020, 8, 2), datetime(2020, 8, 30)], timedelta(days=0.5))
    Planet_trace_epoch = create_epoch([datetime(2020, 8, 2), datetime(2020, 8, 2)], timedelta(days=3))
    Planet_dict = {'Planet_list': Planet_list,
                   'Planet_epoch': Planet_epoch,  # Planet_epoch,
                   'Planet_trace_epoch': [],  # Planet_trace_epoch,
                   'Planet_model_dt': datetime(2022, 2, 23)}

    # Spacecraft
    SC_list = ['SOLO']
    Spacecraft_dict = {'SC_list': SC_list}
    if 'SPP' in SC_list:
        psp_orbit_epoch_range = [datetime(2020, 2, 20), datetime(2022, 2, 25)]
        psp_orbit_epoch_step = timedelta(days=1)
        psp_orbit_epoch = create_epoch(psp_orbit_epoch_range, psp_orbit_epoch_step)
        psp_orbit_marker_epoch_range = [datetime(2022, 2, 20), datetime(2022, 2, 25)]
        psp_orbit_marker_epoch_step = timedelta(days=3)
        psp_orbit_marker_epoch = create_epoch(psp_orbit_marker_epoch_range, psp_orbit_marker_epoch_step)
        psp_model_position_dt = datetime(2022, 2, 23)
        psp_trace_epoch = [datetime(2022, 2, 23)]
        SPP_dict = {'SC_str': 'SPP',
                    'orbit_epoch': psp_orbit_epoch,
                    'orbit_marker_epoch': psp_orbit_marker_epoch,
                    'model_position_dt': psp_model_position_dt,
                    'trace_epoch': psp_trace_epoch,
                    }
    if 'SOLO' in SC_list:
        solo_orbit_epoch_range = [datetime(2020, 8, 1), datetime(2020, 8, 3)]
        solo_orbit_epoch_step = timedelta(days=1)
        solo_orbit_epoch = create_epoch(solo_orbit_epoch_range, solo_orbit_epoch_step)
        solo_orbit_marker_epoch_range = [datetime(2020, 8, 2, 0, 0, 0), datetime(2020, 8, 3, 0, 0 ,0)]
        solo_orbit_marker_epoch_step = timedelta(days=3)
        solo_orbit_marker_epoch = create_epoch(solo_orbit_marker_epoch_range, solo_orbit_marker_epoch_step)
        solo_model_position_dt = datetime(2020, 8, 2, 12, 0, 0)
        solo_trace_epoch = [datetime(2020, 8, 2, 12, 0, 0)]
        SOLO_dict = {'SC_str': 'SOLO',
                     'orbit_epoch': solo_orbit_epoch,
                     'orbit_marker_epoch': solo_orbit_marker_epoch,
                     'model_position_dt': solo_model_position_dt,
                     'trace_epoch': solo_trace_epoch, }
    for SC_str in SC_list:
        Spacecraft_dict[SC_str] = eval(SC_str + '_dict')

    # Sun
    Sun_dict = {'scale': 20, 'crid': 2233, 'surface_data': 'AIA_193'}

    input_dict = {'HCS': HCS_dict, 'Sun': Sun_dict, 'Planet': Planet_dict, 'Spacecraft': Spacecraft_dict}

    import plotly.offline as py

    plot = go.Figure()
    plot_ensemble(input_dict, plot)

    py.plot(plot)

    # st.plotly_chart(plot)
