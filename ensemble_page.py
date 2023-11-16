import pandas as pd
import numpy as np

import streamlit as st
from datetime import datetime, timedelta, time
import sunpy.coordinates.sun as sun

# MY MODULES
from ensemble import plot_ensemble
from utils import create_epoch


def ask_for_epoch(key_, ini_beg_, ini_end_):
    col1, col2 = st.sidebar.columns(2)
    with col1:
        beg_date = st.date_input('Begin Date', ini_beg_, key=key_ + 'beg_d')
        end_date = st.date_input('End Date', ini_end_, key=key_ + 'end_d')
    with col2:
        beg_time = st.time_input('Begin Time', time(0), step=timedelta(minutes=30), key=key_ + 'beg_t')
        end_time = st.time_input('End Time', time(0), step=timedelta(minutes=30), key=key_ + 'end_t')
    time_step = st.sidebar.selectbox('Time Step',
                                     ('1 Day', '12 Hours', '6 Hours', '3 Hours', '1 Hours', '30 Minutes', '15 Minutes'),
                                     key=key_ + 'step')
    time_step_timedelta = convert_time_step(time_step)
    beg_datetime = datetime.combine(beg_date, beg_time)
    end_datetime = datetime.combine(end_date, end_time)
    epoch = create_epoch([beg_datetime, end_datetime], time_step_timedelta)
    return epoch


def ask_for_datetime(key_, ini_beg_):
    col1, col2 = st.sidebar.columns(2)
    with col1:
        input_date = st.date_input('Date', ini_beg_, key=key_ + 'd')

    with col2:
        input_time = st.time_input('Time', time(0), step=timedelta(minutes=30), key=key_ + 't')

    input_datetime = datetime.combine(input_date, input_time)
    return input_datetime


def ask_for_sc_settings(key_, ini_beg_dt_, ini_end_dt_):
    SC_dict = {'SC_str': key_}

    if st.sidebar.checkbox('Plot Orbit', key=key_ + 'orbit'):
        SC_dict['orbit_epoch'] = ask_for_epoch(key_ + 'orbit', ini_beg_dt_, ini_end_dt_)
    else:
        SC_dict['orbit_epoch'] = []

    if st.sidebar.checkbox('Plot Markers', key=key_ + 'marker'):
        SC_dict['orbit_marker_epoch'] = ask_for_epoch(key_ + 'marker', ini_beg_dt_, ini_end_dt_)
    else:
        SC_dict['orbit_marker_epoch'] = []

    if st.sidebar.checkbox('3D Model', key=key_ + '3d'):
        SC_dict['model_position_dt'] = ask_for_datetime(key_ + '3d', ini_beg_dt_)
    else:
        SC_dict['model_position_dt'] = None

    if st.sidebar.checkbox('Trace Back', key=key_ + 'trace'):
        SC_dict['trace_epoch'] = ask_for_epoch(key_ + 'trace', ini_beg_dt_, ini_end_dt_)
    else:
        SC_dict['trace_epoch'] = []
    return SC_dict


def convert_time_step(str_):
    if str_ == '1 Day':
        return timedelta(days=1)
    elif str_ == '12 Hours':
        return timedelta(hours=1)
    elif str_ == '6 Hours':
        return timedelta(hours=6)
    elif str_ == '3 Hours':
        return timedelta(hours=3)
    elif str_ == '1 Hour':
        return timedelta(hours=1)
    elif str_ == '30 Minutes':
        return timedelta(minutes=30)
    elif str_ == '15 Minutes':
        return timedelta(minutes=15)


if __name__ == '__main__':
    st.markdown('PKU Heliosphere Group')
    st.title('PSI Data 3D Visualization')
    st.markdown('Revised @ 23/11/16 \n CAUTION!!! LONGITUDINAL & LATITUDINAL SLICES NOT WORKING !!!')

    st.sidebar.header('Elements to plot: ')

    col1, col2 = st.sidebar.columns(2)
    with col1:
        Sun_options = st.selectbox('Sun', ('Magnetogram (PSI)', 'AIA_193'))
    with col2:
        Sun_scale = st.number_input('Size (Rs)', 1, 10, 1)
    HCS_options = st.sidebar.multiselect('Heliospheric Current Sheet (HCS)',
                                         ('SC', 'IH', 'HPS'), 'SC')
    col1, col2 = st.sidebar.columns(2)
    with col1:
        slice_options = st.number_input('Number of Slices', 0, 3, 0)
    with col2:
        slice_region = st.multiselect('Region', ('SC', 'IH'), None,)

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
    HPS_dict = {}
    if HCS_options != []:
        st.sidebar.divider()
        HCS_dict['plot_sc'] = ('SC' in HCS_options)
        HCS_dict['plot_ih'] = ('IH' in HCS_options)

        st.sidebar.header('HCS Settings')
        HCS_dict['crid'] = Ini_crid
        HCS_surface_style = st.sidebar.radio('Surface palettes',
                                             ('Mono', 'Solar Wind Speed', 'Solar Wind Proton Density'))
        HCS_default_clim = {'Solar Wind Speed': [200, 400], 'Solar Wind Proton Density': [1e2, 1e5]}
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
                                                                  'plasma'))
            # TODO: CHANGE WITH COMMON PLOTLY COLOR SCALES
            # TODO: ON CHANGE, JUST CHANGE THE COLOR INSTEAD OF RERUN THE WHOLE CODE
            HCS_surface_dict['cmap'] = HCS_surface_cmap
            HCS_col1, HCS_col2 = st.sidebar.columns(2)
            HCS_surface_islog = False
            with HCS_col1:
                HCS_surface_showscale = st.checkbox('Show Colorbar')
                HCS_surface_islog = st.checkbox('Log10')
            with HCS_col2:
                HCS_surface_rpower = st.number_input('*r^n', 0, 4, 0, step=1)

            HCS_default_clim[HCS_surface_style] = HCS_default_clim[HCS_surface_style] * 100 ** HCS_surface_rpower
            if HCS_surface_islog:
                HCS_default_clim[HCS_surface_style] = np.log10(HCS_default_clim[HCS_surface_style])

            HCS_col3, HCS_col4 = st.sidebar.columns(2)
            with HCS_col3:
                HCS_surface_cmin = st.number_input('Min', None, None, HCS_default_clim[HCS_surface_style][0])
            with HCS_col4:
                HCS_surface_cmax = st.number_input('Max', None, None, HCS_default_clim[HCS_surface_style][1])

            HCS_surface_dict['clim'] = [HCS_surface_cmin, HCS_surface_cmax]
            HCS_surface_dict['showscale'] = HCS_surface_showscale
            HCS_surface_dict['r_power'] = HCS_surface_rpower
            HCS_surface_dict['islog'] = HCS_surface_islog
        HCS_surface_opacity = st.sidebar.number_input('Opacity', 0.0, 1.0, value=0.9, step=0.1)
        HCS_surface_dict['opacity'] = HCS_surface_opacity
        HCS_dict['surface'] = HCS_surface_dict
        if 'HPS' in HCS_options:
            st.sidebar.header('HPS Settings')
            HPS_dict['crid'] = Ini_crid
            HPS_dict['isos_value'] = st.sidebar.number_input('Isos Value (1e4)', None, None, value=0.9, step=0.1,
                                                             key='isos_hps') * 1e4
            HPS_dict['opacity'] = st.sidebar.number_input('HPS Opacity', 0.0, 1.0, value=0.9, step=0.1,
                                                          key='opacity_hps')
            HPS_dict['color'] = st.sidebar.selectbox('HPS Color', ('azure', 'white'))

    slice_dict = {}
    if slice_options:
        st.sidebar.divider()
        slice_dict['plot_sc'] = ('SC' in slice_region)
        slice_dict['plot_ih'] = ('IH' in slice_region)
        slice_dict['number'] = slice_options
        st.sidebar.header('Slice Settings')

        for i in range(slice_options):
            slice_str = 'slice' + str(i + 1)
            slice_dict[slice_str] = dict()
            st.sidebar.subheader(slice_str)
            slice_dict[slice_str]['slice_type'] = st.sidebar.selectbox('Type',
                                                                 ('Latitudinal', 'Longitudinal', 'Object Plane'),key=slice_str+'_type')
            if slice_dict[slice_str]['slice_type'] == 'Latitudinal':
                slice_dict[slice_str]['z_const'] = st.sidebar.number_input('Z(Rs) = ', None, None, value=0.0, step=1.,key=slice_str+'z_const')
            elif slice_dict[slice_str]['slice_type'] == 'Longitudinal':
                slice_dict[slice_str]['lon_const'] = st.sidebar.number_input('Longitude(deg) = ', 0., 360., value=0.,
                                                                             step=1.,key=slice_str+'lon_const')
            elif slice_dict[slice_str]['slice_type'] == 'Object Plane':
                slice_object_col1, slice_object_col2 = st.sidebar.columns(2)
                with slice_object_col1:
                    slice_dict[slice_str]['object_name'] = st.sidebar.selectbox('Object', (
                    'Earth', 'SPP', 'SO', 'Venus', 'Mercury'),key=slice_str+'object_name')
                with slice_object_col2:
                    slice_dict[slice_str]['z_offset'] = st.sidebar.number_input('Z offset (Rs)', None, None, value=0.,
                                                                                step=1.,key=slice_str+'z_offset')
                slice_dict[slice_str]['object_dt'] = ask_for_datetime(slice_str + '_object', Ini_beg_dt)

        slice_param = st.sidebar.radio('Slice Param',
                                       ('Solar Wind Speed', 'Solar Wind Proton Density'))
        slice_default_clim = {'Solar Wind Speed': [200, 400], 'Solar Wind Proton Density': [1e2, 1e5]}

        if slice_param == 'Solar Wind Speed':
            slice_dict['param'] = 'vr'
        elif slice_param == 'Solar Wind Proton Density':
            slice_dict['param'] = 'rho'
        slice_cmap = st.sidebar.selectbox('Color Map', ('jet',
                                                        'plasma'),
                                          key='slice_cmap')

        slice_dict['colorscale'] = slice_cmap
        slice_col1, slice_col2 = st.sidebar.columns(2)

        with slice_col1:
            slice_showscale = st.checkbox('Show Colorbar')
            slice_islog = st.checkbox('Log10')
        with slice_col2:
            slice_rpower = st.number_input('*r^n', 0, 4, 0, step=1)
            slice_isinterp = st.checkbox('Interp')

        slice_default_clim[slice_param] = slice_default_clim[slice_param] * 100 ** slice_rpower
        if slice_islog:
            slice_default_clim[slice_param] = np.log10(slice_default_clim[slice_param])

        slice_col3, slice_col4 = st.sidebar.columns(2)
        with slice_col3:
            slice_cmin = st.number_input('Min', None, None, slice_default_clim[slice_param][0])
        with slice_col4:
            slice_cmax = st.number_input('Max', None, None, slice_default_clim[slice_param][1])
        slice_dict['clim'] = [slice_cmin, slice_cmax]
        slice_dict['showscale'] = slice_showscale
        slice_dict['r_power'] = slice_rpower
        slice_dict['islog'] = slice_islog
        slice_dict['interp'] = slice_isinterp
        slice_dict['crid'] = Ini_crid
        slice_dict['opacity'] = st.sidebar.number_input('Opacity', 0.0, 1.0, value=0.9, step=0.1, key='opacity_slc')

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
            Planet_dict['Planet_epoch'] = []
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

    layout_dict = {}
    st.sidebar.divider()
    st.sidebar.header('Layout Setting')
    layout_dict['theme'] = 'plotly'
    layout_dict['theme'] = st.sidebar.selectbox('Theme', (
        'plotly', 'plotly_white', 'plotly_dark', 'presentation', 'ggplot2', 'seaborn', 'simple_white', 'none'))
    layout_col1, layout_col2 = st.sidebar.columns(2)

    limit = 35
    try:
        if bool(slice_region) & ('IH' in slice_region):
            limit = 250
        if HCS_options & 'IH' in HCS_options:
            limit= 250
    except:
        pass
    with layout_col1:
        layout_xrng_min = st.number_input('Xmin', None, None, -limit)
        layout_yrng_min = st.number_input('Ymin', None, None, -limit)
        layout_zrng_min = st.number_input('Zmin', None, None, -limit)
    with layout_col2:
        layout_xrng_max = st.number_input('Xmax', None, None, limit)
        layout_yrng_max = st.number_input('Ymax', None, None, limit)
        layout_zrng_max = st.number_input('Zmax', None, None, limit)
    layout_dict['xrng'] = [layout_xrng_min, layout_xrng_max]
    layout_dict['yrng'] = [layout_yrng_min, layout_yrng_max]
    layout_dict['zrng'] = [layout_zrng_min, layout_zrng_max]

    input_dict = {'Sun': Sun_dict, 'HCS': HCS_dict, 'slice': slice_dict, 'HPS': HPS_dict,
                  'Planet': Planet_dict, 'Spacecraft': Spacecraft_dict, 'Layout': layout_dict}
    import plotly.graph_objects as go

    if st.button('Plot'):
        print('===============================Submitted at' + str(
            datetime.now()) + '====================================')
        print(str(input_dict))
        plot = go.Figure()
        plot_ensemble(input_dict, plot)
        st.plotly_chart(plot, use_container_width=True)
        import plotly.offline

        plotly.offline.plot(plot, filename='ensemble/Ensemble.html', auto_open=False, auto_play=False)
        time_id = datetime.now().strftime('%Y%m%dT%H%M')

        input_file = open('ensemble/Ensemble_input.txt', 'w')
        input_file.write(str(input_dict))
        input_file.close()
        import zipfile

        zip_file = zipfile.ZipFile('ensemble.zip', 'w')
        zip_file.write('ensemble/Ensemble.html', compress_type=zipfile.ZIP_DEFLATED)
        zip_file.write('ensemble/Ensemble_input.txt', compress_type=zipfile.ZIP_DEFLATED)
        zip_file.close()

        with open('ensemble.zip', 'rb') as plot_file:
            btn = st.download_button(label='Download Zip', data=plot_file, file_name='Ensemble_' + time_id + '.zip')
        print('=========================END at' + str(datetime.now()) + '============================')
