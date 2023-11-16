"""
This script contain functions that get traces from plot_elements.py for ensemble.py.

parts:
the Sun (A sphere. go.Surface()) <- get_sun_trace
the HCS (SC surf and/or IH surf [go.Mesh3d,(go.Mesh3d)]) <- get_HCS_traces
the Planets (Spheres with surface images.) <- get_planet_traces
"""
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
Rs2AU = Rs2km / AU2km
AU2Rs = AU2km / Rs2km
km2Rs = 1. / Rs2km
km2AU = 1. / AU2km
unit_lst = {'l': 696300,  # km
            'vr': 481.3711, 'vt': 481.3711, 'vp': 481.3711,  # km/s
            'ne': 1e14, 'rho': 1e8, 'p': 3.875717e-2, 'T': 2.807067e7,
            'br': 2.2068908e5, 'bt': 2.2068908e5, 'bp': 2.2068908e5,  # nT
            'j': 2.267e4}

from plot_elements import *


def get_sun_trace(**Sun_dict):
    return plot_sun(**Sun_dict)


def get_HCS_traces(HCS_dict):
    HCS_trace_list = []
    crid = HCS_dict['crid']
    surface_dict = HCS_dict['surface']

    if HCS_dict['plot_sc']:
        HCS_trace_list.append(plot_HCS(crid, 'corona', **surface_dict))
    if HCS_dict['plot_ih']:
        HCS_trace_list.append(plot_HCS(crid, 'helio', **surface_dict))

    return HCS_trace_list


def get_HPS_traces(HPS_dict):
    HPS_trace_list = []
    crid = HPS_dict['crid']
    isos_value = HPS_dict['isos_value']
    opacity = HPS_dict['opacity']
    color = HPS_dict['color']
    HPS_trace_list.append(plot_HPS(crid, isos_value=isos_value, opacity=opacity, color=color))
    return HPS_trace_list


def get_slice_traces(slice_dict):
    slice_trace_list = []
    if slice_dict['plot_sc']:
        for i in range(slice_dict['number']):
            slice_str = 'slice' + str(i + 1)
            if slice_dict[slice_str]['slice_type'] == 'Latitudinal':
                slice_trace_list.append(plot_zslice(region='corona', **slice_dict[slice_str], **slice_dict))
            elif slice_dict[slice_str]['slice_type'] == 'Longitudinal':
                slice_trace_list.append(plot_lonslice(region='corona', **slice_dict[slice_str], **slice_dict))
            elif slice_dict[slice_str]['slice_type'] == 'Object Plane':
                slice_trace_list.append(plot_slice_sun_object(region='corona', **slice_dict[slice_str], **slice_dict))
    if slice_dict['plot_ih']:
        for i in range(slice_dict['number']):
            slice_str = 'slice' + str(i + 1)
            if slice_dict[slice_str]['slice_type'] == 'Latitudinal':
                slice_trace_list.append(plot_zslice(region='helio', **slice_dict[slice_str], **slice_dict))
            elif slice_dict[slice_str]['slice_type'] == 'Longitudinal':
                slice_trace_list.append(plot_lonslice(region='helio', **slice_dict[slice_str], **slice_dict))
            elif slice_dict[slice_str]['slice_type'] == 'Object Plane':
                slice_trace_list.append(plot_slice_sun_object(region='helio', **slice_dict[slice_str], **slice_dict))
    return slice_trace_list


def get_planet_traces(planet_list=[], planet_orbit_epoch=[], planet_trace_epoch=[], planet_model_dt=None):
    Planet_trace_list = []
    for Planet_name in planet_list:
        if planet_orbit_epoch != []:
            orbit_trace = plot_object_orbit(Planet_name, planet_orbit_epoch, type='line')
            Planet_trace_list.append(orbit_trace)
        if planet_trace_epoch != []:
            magtrace_trace_list = plot_obj_magtrace(Planet_name, planet_trace_epoch)
            for magtrace_trace in magtrace_trace_list:
                Planet_trace_list.append(magtrace_trace)
        if planet_model_dt != None:
            Planet_model_trace = plot_planet_model(Planet_name, planet_model_dt)
            Planet_trace_list.append(Planet_model_trace)
    return Planet_trace_list


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
