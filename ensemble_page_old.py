import pandas as pd

from ensemble import plot_ensemble, create_epoch
import streamlit as st
from datetime import datetime, timedelta, time
import sunpy.coordinates.sun as sun


def ask_for_epoch(key, ini_beg, ini_end):
    col1, col2 = st.sidebar.columns(2)
    with col1:
        beg_date = st.date_input('Begin Date', ini_beg, key=key + 'beg_d')
        end_date = st.date_input('End Date', ini_end, key=key + 'end_d')
    with col2:
        beg_time = st.time_input('Begin Time', time(0), step=timedelta(minutes=30), key=key + 'beg_t')
        end_time = st.time_input('End Time', time(0), step=timedelta(minutes=30), key=key + 'end_t')
    time_step = st.sidebar.selectbox('Time Step',
                                     ('1 Day', '12 Hours', '6 Hours', '3 Hours', '1 Hours', '30 Minutes', '15 Minutes'),
                                     key=key + 'step')
    time_step_timedelta = convert_time_step(time_step)
    beg_datetime = datetime.combine(beg_date, beg_time)
    end_datetime = datetime.combine(end_date, end_time)
    epoch = create_epoch([beg_datetime, end_datetime], time_step_timedelta)
    return epoch


def ask_for_datetime(key, ini_beg):
    col1, col2 = st.sidebar.columns(2)
    with col1:
        input_date = st.date_input('Date', ini_beg, key=key + 'd')

    with col2:
        input_time = st.time_input('Time', time(0), step=timedelta(minutes=30), key=key + 't')

    input_datetime = datetime.combine(input_date, input_time)
    return input_datetime


def ask_for_sc_settings(key, ini_beg_dt, ini_end_dt):
    SC_dict = {'SC_str': key}

    if st.sidebar.checkbox('Plot Orbit', key=key + 'orbit'):
        SC_dict['orbit_epoch'] = ask_for_epoch(key + 'orbit', ini_beg_dt, ini_end_dt)
    else:
        SC_dict['orbit_epoch'] = []

    if st.sidebar.checkbox('Plot Markers', key=key + 'marker'):
        SC_dict['orbit_marker_epoch'] = ask_for_epoch(key + 'marker', ini_beg_dt, ini_end_dt)
    else:
        SC_dict['orbit_marker_epoch'] = []

    if st.sidebar.checkbox('3D Model', key=key + '3d'):
        SC_dict['model_position_dt'] = ask_for_datetime(key + '3d', ini_beg_dt)
    else:
        SC_dict['model_position_dt'] = None

    if st.sidebar.checkbox('Trace Back', key=key + 'trace'):
        SC_dict['trace_epoch'] = ask_for_epoch(key + 'trace', ini_beg_dt, ini_end_dt)
    else:
        SC_dict['trace_epoch'] = []
    return SC_dict


def convert_time_step(str):
    if str == '1 Day':
        return timedelta(days=1)
    elif str == '12 Hours':
        return timedelta(hours=1)
    elif str == '6 Hours':
        return timedelta(hours=6)
    elif str == '3 Hours':
        return timedelta(hours=3)
    elif str == '1 Hour':
        return timedelta(hours=1)
    elif str == '30 Minutes':
        return timedelta(minutes=30)
    elif str == '15 Minutes':
        return timedelta(minutes=15)


if __name__ == '__main__':
    st.markdown('PKU Heliosphere Group')
    st.title('PSI Data 3D Visualization')
    st.markdown('测试中')

    st.sidebar.header('Elements to plot: ')

    col1, col2 = st.sidebar.columns(2)
    with col1:
        Sun_options = st.selectbox('Sun', ('Magnetogram (PSI)', 'AIA_193'))
    with col2:
        Sun_scale = st.number_input('Size (Rs)', 1, 10, 1)
    HCS_options = st.sidebar.multiselect('Heliospheric Current Sheet (HCS)',
                                         ('SC', 'IH'), 'SC')
    PLANETS_options = st.sidebar.multiselect('Planets',
                                             ('EARTH', 'VENUS', 'MERCURY'))
    Spacecrafts_options = st.sidebar.multiselect('Spacecrafts',
                                                 ('Parker Solar Probe', 'Solar Orbiter'))

    st.sidebar.divider()
    st.sidebar.header('Carrington Number For Sun/HCS')
    Ini_crid = st.sidebar.number_input('Carrington Rotation', 1500, 3000, 2254)
    Ini_beg_dt = sun.carrington_rotation_time(Ini_crid).to_datetime()
    Ini_end_dt = sun.carrington_rotation_time(Ini_crid + 1).to_datetime()
    st.sidebar.write('From ' + Ini_beg_dt.strftime('%Y-%m-%d/%H:%M:%S'))
    st.sidebar.write('To ' + Ini_end_dt.strftime('%Y-%m-%d/%H:%M:%S'))

    Sun_dict = {'crid': Ini_crid}
    if Sun_options == 'Magnetogram (PSI)':
        Sun_dict['surface_data'] = 'br'
    elif Sun_options == 'AIA_193':
        Sun_dict['surface_data'] = 'AIA_193'
    Sun_dict['scale'] = Sun_scale

    HCS_dict = {}
    HCS_surface_dict = {}
    if HCS_options != []:
        st.sidebar.divider()
        HCS_dict['plot_sc'] = ('SC' in HCS_options)
        HCS_dict['plot_ih'] = ('IH' in HCS_options)
        st.sidebar.header('HCS Settings')
        HCS_dict['crid'] = Ini_crid
        HCS_surface_style = st.sidebar.radio('Surface palettes',
                                             ('Mono', 'Solar Wind Speed', 'Solar Wind Proton Density'))
        HCS_default_clim = {'Solar Wind Speed': [200, 400], 'Solar Wind Proton Density': [0, 1000]}
        if HCS_surface_style == 'Mono':
            HCS_surface_color = st.sidebar.selectbox('Color', (
            'purple', 'white', 'blue', 'yellow'))  # [FUTURE] CHANGE WITH COMMON PLOTLY COLORS OR COLOR PICKER
            HCS_surface_dict['color'] = HCS_surface_color
            HCS_surface_dict['param'] = None
        else:
            if HCS_surface_style == 'Solar Wind Speed':
                HCS_surface_dict['param'] = 'vr'
            elif HCS_surface_style == 'Solar Wind Proton Density':
                HCS_surface_dict['param'] = 'rho'
            HCS_surface_cmap = st.sidebar.selectbox('Color Map', ('jet',
                                                                  'plasma'))  # [FUTURE] CHANGE WITH COMMON PLOTLY COLOR SCALES # [FUTURE] ON CHANGE, JUST CHANGE THE COLOR INSTEAD OF RERUN THE WHOLE CODE
            HCS_surface_dict['cmap'] = HCS_surface_cmap
            HCS_col1, HCS_col2 = st.sidebar.columns(2)
            with HCS_col1:
                HCS_surface_cmin = st.number_input('Min', None, None, HCS_default_clim[HCS_surface_style][0], step=50)
            with HCS_col2:
                HCS_surface_cmax = st.number_input('Max', None, None, HCS_default_clim[HCS_surface_style][1])
            HCS_surface_dict['clim'] = [HCS_surface_cmin, HCS_surface_cmax]
            HCS_surface_showscale = st.sidebar.checkbox('Show Colorbar')
            HCS_surface_dict['showscale'] = HCS_surface_showscale
        HCS_surface_opacity = st.sidebar.number_input('Opacity', 0.0, 1.0, value=0.9, step=0.1)
        HCS_surface_dict['opacity'] = HCS_surface_opacity
        HCS_dict['surface'] = HCS_surface_dict

    Planet_dict = {}
    if PLANETS_options != []:
        st.sidebar.divider()
        st.sidebar.header('Planets Settings')
        Planet_list = PLANETS_options
        Planet_dict['Planet_list'] = Planet_list
        if st.sidebar.checkbox('Plot Orbit'):
            planet_orbit_epoch = ask_for_epoch('planet_orbit', Ini_beg_dt, Ini_end_dt)
            Planet_dict['Planet_epoch'] = planet_orbit_epoch
        else:
            Planet_dict['Planet_epoch']=[]
        if st.sidebar.checkbox('Trace Back'):
            planet_trace_epoch = ask_for_epoch('planet_trace', Ini_beg_dt, Ini_end_dt)
            Planet_dict['Planet_trace_epoch'] = planet_trace_epoch
        else:
            Planet_dict['Planet_trace_epoch'] = []
        if st.sidebar.checkbox('Plot 3D Model'):
            Planet_model_dt = ask_for_datetime('planet_model', Ini_beg_dt)
            Planet_dict['Planet_model_dt'] = Planet_model_dt
        else:
            Planet_dict['Planet_model_dt'] = None


    Spacecraft_dict = {}
    if Spacecrafts_options != []:
        st.sidebar.divider()
        st.sidebar.header('Spacecrafts Settings')
        SC_list = []
        if 'Parker Solar Probe' in Spacecrafts_options:
            SC_list.append('SPP')
            st.sidebar.subheader('Parker Solar Probe')
            SPP_dict = ask_for_sc_settings('SPP', Ini_beg_dt, Ini_end_dt)
            Spacecraft_dict['SPP'] = SPP_dict

        if 'Solar Orbiter' in Spacecrafts_options:
            SC_list.append('SOLO')
            st.sidebar.subheader('Solar Orbiter')
            SOLO_dict = ask_for_sc_settings('SOLO', Ini_beg_dt, Ini_end_dt)
            Spacecraft_dict['SOLO'] = SOLO_dict
        Spacecraft_dict['SC_list'] = SC_list

    input_dict = {'Sun': Sun_dict, 'HCS': HCS_dict, 'Planet': Planet_dict, 'Spacecraft': Spacecraft_dict}
    import plotly.graph_objects as go

    if st.button('Plot'):
        plot = go.Figure()
        plot_ensemble(input_dict, plot)
        st.plotly_chart(plot,use_container_width=True)
        import plotly.offline
        plotly.offline.plot(plot,filename='ensemble/Ensemble.html',auto_open=False,auto_play=False)
        time_id = datetime.now().strftime('%Y%m%dT%H%M')

        input_file = open('ensemble/Ensemble_input.txt','w')
        input_file.write(str(input_dict))
        input_file.close()
        import zipfile
        zip_file = zipfile.ZipFile('ensemble.zip','w')
        zip_file.write('ensemble/Ensemble.html',compress_type=zipfile.ZIP_DEFLATED)
        zip_file.write('ensemble/Ensemble_input.txt', compress_type=zipfile.ZIP_DEFLATED)
        zip_file.close()

        with open('ensemble.zip','rb') as plot_file:
            btn = st.download_button(label='Download Zip',data=plot_file,file_name='Ensemble_'+time_id+'.zip')






