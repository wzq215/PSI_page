import os
import requests

import numpy as np
from datetime import datetime, timedelta
import pyvista

import sunpy.map
import sunpy.coordinates.sun as sun

import plotly.graph_objects as go
import spiceypy as spice

from ps_read_hdf_3d import ps_read_hdf_3d
from Trace_PSI_data import PSI_trace
import furnsh_kernels

Rs2km = 696300.
AU2km = 1.49597871e8
Rs2AU = Rs2km/AU2km
AU2Rs = AU2km/Rs2km
km2Rs = 1./Rs2km
km2AU = 1./AU2km
unit_lst = {'l': 696300,  # km
            'vr': 481.3711, 'vt': 481.3711, 'vp': 481.3711,  # km/s
            'ne': 1e14, 'rho': 1e8, 'p': 3.875717e-2, 'T': 2.807067e7,
            'br': 2.2068908e5, 'bt': 2.2068908e5, 'bp': 2.2068908e5,  # nT
            'j': 2.267e4}

from plot_elements import *

def get_HCS_traces(HCS_dict):
    HCS_trace_list = []
    crid = HCS_dict['crid']
    surface_dict = HCS_dict['surface']

    if HCS_dict['plot_sc']:
        HCS_trace_list.append(plot_HCS(crid, 'corona', surface_dict))
    if HCS_dict['plot_ih']:
        HCS_trace_list.append(plot_HCS(crid, 'helio', surface_dict))

    return HCS_trace_list

def get_HPS_traces(HPS_dict):
    HPS_trace_list = []
    crid = HPS_dict['crid']
    isos_value = HPS_dict['isos_value']
    opacity = HPS_dict['opacity']
    color = HPS_dict['color']
    HPS_trace_list.append(plot_HPS(crid,isos_value=isos_value,opacity=opacity,color=color))
    return HPS_trace_list

def get_SC_traces(SC_dict):
    SC_trace_list = []
    SC_str = SC_dict['SC_str']
    if len(SC_dict['orbit_epoch']) != 0:
        traces_orbit = plot_object_orbit(SC_str, SC_dict['orbit_epoch'], 'line')
        SC_trace_list.append(traces_orbit)

    if len(SC_dict['orbit_marker_epoch']) != 0:
        traces_marker = plot_object_orbit(SC_str, SC_dict['orbit_marker_epoch'], 'marker')
        SC_trace_list.append(traces_marker)
    if SC_dict['model_position_dt'] != None:
        trace_model = plot_Spacecraft_model(SC_str, SC_dict['model_position_dt'])
        SC_trace_list.append(trace_model)
    if len(SC_dict['trace_epoch']) != 0:
        traces_trace_list = plot_obj_magtrace(SC_str, SC_dict['trace_epoch'])
        for trace in traces_trace_list:
            SC_trace_list.append(trace)
    return SC_trace_list

def get_sun_trace(Sun_dict):
    return plot_sun(Sun_dict)


def get_planet_traces(Planet_dict):
    Planet_list = Planet_dict['Planet_list']
    Planet_orbit_epoch = Planet_dict['Planet_epoch']
    Planet_trace_epoch = Planet_dict['Planet_trace_epoch']
    Planet_model_dt = Planet_dict['Planet_model_dt']
    Planet_trace_list = []
    for Planet_name in Planet_list:
        if Planet_orbit_epoch != []:
            orbit_trace = plot_object_orbit(Planet_name, Planet_orbit_epoch, type='line')
            Planet_trace_list.append(orbit_trace)
        if Planet_trace_epoch != []:
            magtrace_trace_list = plot_obj_magtrace(Planet_name, Planet_trace_epoch)
            for magtrace_trace in magtrace_trace_list:
                Planet_trace_list.append(magtrace_trace)
        if Planet_model_dt != None:
            Planet_model_trace = plot_planet_model(Planet_name, Planet_model_dt)
            Planet_trace_list.append(Planet_model_trace)
    return Planet_trace_list

