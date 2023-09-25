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

# from utils import  *
from get_traces import *

#@st.cache_data(ttl=3600)
def plot_ensemble(input_dict, plot):
    if len(input_dict['HCS']) != 0:
        HCS_dict = input_dict['HCS']
        trace_hcs_list = get_HCS_traces(HCS_dict)
        for i in range(len(trace_hcs_list)):
            plot.add_trace(trace_hcs_list[i])
    if len(input_dict['HPS']) != 0:
        HPS_dict = input_dict['HPS']
        trace_hps_list = get_HPS_traces(HPS_dict)
        for i in range(len(trace_hps_list)):
            plot.add_trace(trace_hps_list[i])
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
    plot_layout(plot, input_dict['Layout'])






if __name__ == '__main__':
    # HCS
    HCS_surface_dict = {'param': 'vr',
                        'cmap': 'jet', 'clim': [200, 500], 'opacity': .8,
                        'color': 'purple', 'showscale': True,
                        'r_power':0}
    HCS_dict = {'crid': 2254, 'plot_sc': False, 'plot_ih': True, 'surface': HCS_surface_dict}
    
    # HPS
    HPS_dict = {'crid': 2254}

    # Planets
    Planet_list = ['EARTH']
    Planet_epoch = create_epoch([datetime(2022, 2, 20), datetime(2022, 2, 25)], timedelta(days=1))
    Planet_trace_epoch = create_epoch([datetime(2022, 2, 20), datetime(2022, 2, 25)], timedelta(days=3))
    Planet_dict = {'Planet_list': Planet_list,
                   'Planet_epoch': [],  # Planet_epoch,
                   'Planet_trace_epoch': [],  # Planet_trace_epoch,
                   'Planet_model_dt': datetime(2022, 2, 23)}

    # Spacecraft
    SC_list = ['SPP', 'SOLO']
    Spacecraft_dict = {'SC_list': SC_list}
    psp_orbit_epoch_range = [datetime(2022, 2, 20), datetime(2022, 2, 25)]
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

    solo_orbit_epoch_range = [datetime(2022, 2, 20), datetime(2022, 2, 25)]
    solo_orbit_epoch_step = timedelta(days=1)
    solo_orbit_epoch = create_epoch(solo_orbit_epoch_range, solo_orbit_epoch_step)
    solo_orbit_marker_epoch_range = [datetime(2022, 2, 20), datetime(2022, 2, 25)]
    solo_orbit_marker_epoch_step = timedelta(days=3)
    solo_orbit_marker_epoch = create_epoch(solo_orbit_marker_epoch_range, solo_orbit_marker_epoch_step)
    solo_model_position_dt = datetime(2022, 2, 23)
    solo_trace_epoch = [datetime(2022, 2, 23)]
    SOLO_dict = {'SC_str': 'SOLO',
                 'orbit_epoch': solo_orbit_epoch,
                 'orbit_marker_epoch': solo_orbit_marker_epoch,
                 'model_position_dt': solo_model_position_dt,
                 'trace_epoch': solo_trace_epoch, }
    for SC_str in SC_list:
        Spacecraft_dict[SC_str] = eval(SC_str + '_dict')

    # Sun
    Spacecraft_dict = {}
    HCS_dict = {}
    HPS_dict = {}
    Sun_dict = {'scale': 1, 'crid': 2254,
                'surface_data': 'AIA_193'}
    Sun_dict = {}

    
    input_dict = {'HCS': HCS_dict,'HPS': HPS_dict, 'Sun': Sun_dict, 'Planet': Planet_dict, 'Spacecraft': Spacecraft_dict}

    import plotly.offline as py

    plot = go.Figure()
    plot_ensemble(input_dict, plot)

    py.plot(plot)

    # st.plotly_chart(plot)
