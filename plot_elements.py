import os
import requests

import numpy as np
from datetime import datetime, timedelta
import pyvista

import sunpy.coordinates.sun as sun

import plotly.graph_objects as go
import spiceypy as spice

# MY MODULES
from ps_read_hdf_3d import ps_read_hdf_3d
from Trace_PSI_data import PSI_trace
import furnsh_kernels
from utils import *

# TODO：提供从PFSS+ParkerSpiral生成电流片的选项
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


def plot_sun(Sun_dict):
    crid = Sun_dict['crid']
    if Sun_dict['surface_data'] == 'br':
        data_sc = ps_read_hdf_3d(crid, 'corona', 'br002', periodicDim=3)
        t_br = np.array(data_sc['scales2'])  # 150 in rad, latitude
        p_br = np.array(data_sc['scales3'])  # 256 in rad, Carrington longitude
        br = np.array(data_sc['datas'])  # 1CU = 2.205G = 2.205e-4T = 2.205e5nT
        br = br * 2.205e5  # nT
        sundata = np.squeeze(br[:, :, 0]).T
        x0, y0, z0 = get_sphere(Sun_dict['scale'], p_br, t_br, for_PSI=True)

        return go.Surface(x=x0, y=y0, z=z0, surfacecolor=sundata, colorscale='rdbu', opacity=1, showscale=False,
                          cmax=3e5,
                          cmin=-3e5)
    elif Sun_dict['surface_data'] == 'AIA_193':
        lon, lat, sundata, map_cmap = load_AIA_data(crid)
        x0, y0, z0 = get_sphere(Sun_dict['scale'], lon, lat)
        return go.Surface(x=x0, y=y0, z=z0, surfacecolor=np.squeeze(sundata), colorscale='turbid_r', opacity=1,
                          showscale=False)


def plot_planet_model(planet_str, dt):
    et = spice.datetime2et(dt)
    planet_pos, _ = spice.spkpos(planet_str, et, 'IAU_SUN', 'NONE', 'SUN')
    planet_pos = np.array(planet_pos).T / Rs2km

    import matplotlib.pyplot as plt
    img = plt.imread('Planet_model/' + planet_str + '.jpeg')
    from PIL import Image
    eight_bit_img = Image.fromarray(img).convert('P', palette='WEB', dither=None)
    idx_to_color = np.array(eight_bit_img.getpalette()).reshape(-1, 3)
    colorscale = [[i / 255.0, 'rgb({},{},{})'.format(*rgb)] for i, rgb in enumerate(idx_to_color)]
    lon = np.linspace(0, 2 * np.pi, img.shape[1])
    lat = np.linspace(np.pi / 2, -np.pi / 2, img.shape[0])
    x0, y0, z0 = get_sphere(1, lon, lat)
    print('Plotting ' + planet_str)
    trace = go.Surface(x=x0 + planet_pos[0], y=y0 + planet_pos[1], z=z0 + planet_pos[2],
                       surfacecolor=np.array(eight_bit_img), cmin=0, cmax=255, colorscale=colorscale, showscale=False)
    return trace


def plot_HPS(crid, isos_value=3e5, opacity=0.8, color='azure'):
    data = ps_read_hdf_3d(crid, 'corona', 'rho002', periodicDim=3)
    r = np.array(data['scales1'])  # 201 in Rs, distance from sun
    t = np.array(data['scales2'])  # 150 in rad, latitude
    p = np.array(data['scales3'])  # 256 in rad, Carrington longitude
    rho = np.array(data['datas'])  # 1CU = 2.205G = 2.205e-4T = 2.205e5nT
    rho = rho * 1e8  # cm^-3
    tv, pv, rv = np.meshgrid(t, p, r, indexing='xy')

    xv = rv * np.cos(pv) * np.sin(tv)
    yv = rv * np.sin(pv) * np.sin(tv)
    zv = rv * np.cos(tv)

    mesh = pyvista.StructuredGrid(xv, yv, zv)
    mesh.point_data['values'] = rho.ravel(order='F')  # also the active scalars
    isos = mesh.contour(isosurfaces=1, rng=[isos_value, isos_value])

    vertices = isos.points
    triangles = isos.faces.reshape(-1, 4)
    trace = go.Mesh3d(x=vertices[:, 0], y=vertices[:, 1], z=vertices[:, 2],
                      opacity=opacity,
                      i=triangles[:, 1], j=triangles[:, 2], k=triangles[:, 3],
                      color=color, )
    return trace


def plot_HCS(crid, region, surface_dict):
    data = ps_read_hdf_3d(crid, region, 'br002', periodicDim=3)
    r = np.array(data['scales1'])  # 201 in Rs, distance from sun
    t = np.array(data['scales2'])  # 150 in rad, latitude
    p = np.array(data['scales3'])  # 256 in rad, Carrington longitude
    br = np.array(data['datas'])  # 1CU = 2.205G = 2.205e-4T = 2.205e5nT
    br = br * 2.205e5  # nT

    tv, pv, rv = np.meshgrid(t, p, r, indexing='xy')

    xv = rv * np.cos(pv) * np.sin(tv)
    yv = rv * np.sin(pv) * np.sin(tv)
    zv = rv * np.cos(tv)

    mesh = pyvista.StructuredGrid(xv, yv, zv)
    mesh.point_data['values'] = br.ravel(order='F')  # also the active scalars
    isos = mesh.contour(isosurfaces=1, rng=[0, 0])

    vertices = isos.points
    triangles = isos.faces.reshape(-1, 4)

    if surface_dict['param'] != None:
        param = surface_dict['param']
        param_points = vertices[:, 0] * 0.

        data = ps_read_hdf_3d(crid, region, param + '002', periodicDim=3)
        r_param = np.array(data['scales1'])
        t_param = np.array(data['scales2'])
        p_param = np.array(data['scales3'])
        param_data = np.array(data['datas'])
        param_data = param_data * unit_lst[param]

        for i in range(len(vertices)):
            point = np.array(vertices[i])
            r_point, p_point, t_point = xyz2rtp_in_Carrington(point, for_psi=True)

            r_ind = np.argmin(abs(r_point - r_param))
            p_ind = np.argmin(abs(p_point - p_param))
            t_ind = np.argmin(abs(t_point - t_param))

            param_points[i] = param_data[p_ind, t_ind, r_ind] * (r_param[r_ind] ** surface_dict['r_power'])

        if surface_dict['islog']:
            param_points = np.log10(param_points)
        intensity = np.array(param_points).reshape(-1, 1)
        if surface_dict['clim'] == None:
            surface_dict['clim'] = [np.nanmin(intensity), np.nanmax(intensity)]

        trace = go.Mesh3d(x=vertices[:, 0], y=vertices[:, 1], z=vertices[:, 2],
                          opacity=surface_dict['opacity'],
                          colorscale=surface_dict['cmap'],
                          cmin=surface_dict['clim'][0], cmax=surface_dict['clim'][1],
                          i=triangles[:, 1], j=triangles[:, 2], k=triangles[:, 3],
                          intensity=intensity,
                          showscale=surface_dict['showscale'], )
    else:
        trace = go.Mesh3d(x=vertices[:, 0], y=vertices[:, 1], z=vertices[:, 2],
                          opacity=surface_dict['opacity'],
                          i=triangles[:, 1], j=triangles[:, 2], k=triangles[:, 3],
                          color=surface_dict['color'], )
    return trace


def plot_object_orbit(obj_name, dt_epoch, type='line',**kwargs):
    et_epoch = spice.datetime2et(dt_epoch)
    pos, _ = spice.spkpos(obj_name, et_epoch, 'IAU_SUN', 'NONE', 'SUN')
    pos = np.array(pos.T / Rs2km)
    pos = {'x': pos[0], 'y': pos[1], 'z': pos[2]}

    trace = []
    if type == 'line':
        trace = go.Scatter3d(x=pos['x'], y=pos['y'], z=pos['z'],
                             mode='lines',
                             line=dict(width=10,**kwargs))
    elif type == 'marker':
        trace = go.Scatter3d(x=pos['x'], y=pos['y'], z=pos['z'],
                             mode='markers',
                             marker=dict(size=3, **kwargs))
    return trace


def plot_Spacecraft_model(SC_str, dt, scale=10):
    et = spice.datetime2et(dt)
    SC_pos, _ = spice.spkpos(SC_str, et, 'IAU_SUN', 'NONE', 'SUN')
    SC_pos = np.array(SC_pos).T / Rs2km
    trans_arr = np.eye(3)
    if SC_str == 'SOLO':
        trans_arr = spice.sxform('SOLO_SRF', 'IAU_SUN', et)
        trans_arr = np.dot(trans_arr[0:3, 0:3], np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]]))
    elif SC_str == 'SPP':
        trans_arr = spice.sxform('SPP_SPACECRAFT', 'IAU_SUN', et)
        trans_arr = np.dot(trans_arr[0:3, 0:3], np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]]))

    trans_arr_4 = np.eye(4)
    trans_arr_4[0:3, 0:3] = trans_arr

    mesh = pyvista.read('SC_model/' + SC_str + '.stl')
    mesh.points = scale * mesh.points / np.max(mesh.points)

    axes = pyvista.Axes(show_actor=True, actor_scale=100.0, line_width=10)
    axes.origin = (0, 0, 0)

    trans = mesh.transform(trans_arr_4)

    # p = pyvista.Plotter()
    # p.add_actor(axes.actor, pickable=True)
    # p.add_mesh(mesh)
    # # p.add_mesh(trans)
    # p.show_grid()
    # p.show_axes()
    # p.show()

    vertices = trans.points
    triangles = trans.faces.reshape(-1, 4)
    trace = go.Mesh3d(x=vertices[:, 0] + SC_pos[0], y=vertices[:, 1] + SC_pos[1], z=vertices[:, 2] + SC_pos[2],
                      opacity=1,
                      color='gold',
                      i=triangles[:, 1], j=triangles[:, 2], k=triangles[:, 3],
                      showscale=False,
                      )
    return trace


def plot_FOV(SC_str, inst_str, dt, fov_dist_Rs=20., color='orange', plot_type='wire'):
    et = spice.datetime2et(dt)
    SC_pos, _ = spice.spkpos(SC_str, et, 'IAU_SUN', 'NONE', 'SUN')
    SC_pos = np.array(SC_pos).T / Rs2km
    if SC_str == 'SOLO':
        print('NOT SUPPORTED YET')
    elif SC_str == 'SPP':
        trans_arr = spice.sxform(inst_str, 'IAU_SUN', et)
        trans_arr = trans_arr[0:3, 0:3]
        if inst_str == 'SPP_WISPR_INNER':
            fov_parameter = spice.getfov(-96100, 4)
        elif inst_str == 'SPP_WISPR_OUTER':
            fov_parameter = spice.getfov(-96120, 4)
        else:
            print('NOT SUPPORTED YET')
    fov_edges = fov_parameter[4]
    fov_edges = fov_edges.T / fov_edges[:, 2] * fov_dist_Rs
    fov_edges_sun = np.dot(trans_arr, fov_edges).T

    trace_lst = []
    if 'wire' in plot_type:
        for i in range(4):
            trace_lst.append(go.Scatter3d(x=[SC_pos[0], SC_pos[0] + fov_edges_sun[i][0]],
                                          y=[SC_pos[1], SC_pos[1] + fov_edges_sun[i][1]],
                                          z=[SC_pos[2], SC_pos[2] + fov_edges_sun[i][2]],
                                          mode='lines',
                                          line=dict(width=5, color=color)))
            trace_lst.append(go.Scatter3d(x=[SC_pos[0] + fov_edges_sun[i - 1][0], SC_pos[0] + fov_edges_sun[i][0]],
                                          y=[SC_pos[1] + fov_edges_sun[i - 1][1], SC_pos[1] + fov_edges_sun[i][1]],
                                          z=[SC_pos[2] + fov_edges_sun[i - 1][2], SC_pos[2] + fov_edges_sun[i][2]],
                                          mode='lines',
                                          line=dict(width=5, color=color)))
    if 'surf' in plot_type:
        trace_lst.append(go.Mesh3d(x=[SC_pos[0],
                                      SC_pos[0] + fov_edges_sun[0][0],
                                      SC_pos[0] + fov_edges_sun[1][0],
                                      SC_pos[0] + fov_edges_sun[2][0],
                                      SC_pos[0] + fov_edges_sun[3][0], ],
                                   y=[SC_pos[1],
                                      SC_pos[1] + fov_edges_sun[0][1],
                                      SC_pos[1] + fov_edges_sun[1][1],
                                      SC_pos[1] + fov_edges_sun[2][1],
                                      SC_pos[1] + fov_edges_sun[3][1], ],
                                   z=[SC_pos[2],
                                      SC_pos[2] + fov_edges_sun[0][2],
                                      SC_pos[2] + fov_edges_sun[1][2],
                                      SC_pos[2] + fov_edges_sun[2][2],
                                      SC_pos[2] + fov_edges_sun[3][2], ],
                                   i=[0, 0, 0, 0, ], j=[1, 2, 3, 4], k=[2, 3, 4, 1],
                                   color='azure', opacity=0.5))

    return trace_lst


def plot_flux_tube(start_lon_deg, start_lat_deg, r_vect_Rs, v_vect_kmps, delta_deg=1., ):
    lon_vect_deg, lat_vect_deg = parker_spiral(r_vect_Rs * Rs2AU, start_lat_deg, start_lon_deg, v_vect_kmps)
    x0_vect_Rs, y0_vect_Rs, z0_vect_Rs = rlonlat2cart(r_vect_Rs, np.deg2rad(lon_vect_deg), np.deg2rad(lat_vect_deg))
    N = 8
    tube_edge_lst = []
    tube_circ_lst = []
    circ_edge = []
    for i in range(N):
        delta_lon_deg = np.cos(i / N * np.pi * 2) * delta_deg
        delta_lat_deg = np.sin(i / N * np.pi * 2) * delta_deg
        x_vect_Rs, y_vect_Rs, z_vect_Rs = rlonlat2cart(r_vect_Rs, np.deg2rad(lon_vect_deg + delta_lon_deg),
                                                       np.deg2rad(lat_vect_deg + delta_lat_deg))
        tube_edge_lst.append(go.Scatter3d(x=x_vect_Rs, y=y_vect_Rs, z=z_vect_Rs,
                                          mode='lines', line=dict(width=2, color='yellow')))
        circ_edge.append(np.vstack([np.hstack([x_vect_Rs[::2], x_vect_Rs[-1]]),
                                    np.hstack([y_vect_Rs[::2], y_vect_Rs[-1]]),
                                    np.hstack([z_vect_Rs[::2], z_vect_Rs[-1]])]))
    circ_edge = np.array(circ_edge)
    for i in range(np.shape(circ_edge)[2]):
        tube_circ_lst.append(go.Scatter3d(x=np.hstack([circ_edge[:, 0, i], circ_edge[0, 0, i]]),
                                          y=np.hstack([circ_edge[:, 1, i], circ_edge[0, 1, i]]),
                                          z=np.hstack([circ_edge[:, 2, i], circ_edge[0, 2, i]]),
                                          mode='lines', line=dict(width=2, color='yellow')))
    trace_lst = np.hstack([tube_edge_lst, tube_circ_lst])
    return list(trace_lst)


def plot_obj_magtrace(obj_str, dt_epoch, stop_at_obj=False):
    et_epoch = spice.datetime2et(dt_epoch)
    obj_pos, _ = spice.spkpos(obj_str, et_epoch, 'IAU_SUN', 'NONE', 'SUN')
    obj_pos = np.array(obj_pos.T / Rs2km)
    obj_magtrace_list = []
    for i in range(len(dt_epoch)):
        crid = int(sun.carrington_rotation_number(dt_epoch[i]))
        tmp_mag_trace_list = PSI_trace(obj_pos[:, i], crid)
        trace_points = []
        for mag_trace in tmp_mag_trace_list:
            for j in range(len(mag_trace.points)):
                trace_points.append(mag_trace.points[j, :])
        trace_points = np.array(trace_points)
        r_trace_points = np.sqrt(trace_points[:, 0] ** 2 + trace_points[:, 1] ** 2 + trace_points[:, 2] ** 2)
        trace_points = trace_points[np.argsort(r_trace_points)]
        r_obj, _, _ = xyz2rtp_in_Carrington(obj_pos)
        ind_obj = np.argmin(abs(r_obj - r_trace_points))
        if stop_at_obj:
            trace_points = trace_points[:ind_obj]
        obj_magtrace_list.append(go.Scatter3d(x=trace_points[:, 0], y=trace_points[:, 1], z=trace_points[:, 2],
                                              mode='lines', line=dict(color='white', width=1.)))
    return obj_magtrace_list


def plot_layout(plot, layout_dict):
    camera = dict(
        up=dict(x=0, y=0, z=1),
        center=dict(x=0, y=0, z=0),
        eye=dict(x=0, y=-1.5, z=1.5)
    )
    xrng = layout_dict['xrng']
    yrng = layout_dict['yrng']
    zrng = layout_dict['zrng']
    plot.update_layout(
        autosize=False,
        width=1000,
        height=1000,
        scene=dict(
            xaxis_title='X (Rs)',
            yaxis_title='Y (Rs)',
            zaxis_title='Z (Rs)',
            xaxis_range=[xrng[0], xrng[1]],
            yaxis_range=[yrng[0], yrng[1]],
            zaxis_range=[zrng[0], zrng[1]],
        ),
        scene_camera=camera,
        title=dict(
            text='Carrington Coordinate',
            font_size=30,
            y=0.9, x=0.5,
            xanchor='center',
            yanchor='top'),
        template=layout_dict['theme'],
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
