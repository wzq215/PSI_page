from datetime import datetime, timedelta

import numpy as np
import pandas as pd

import pyvista
import spiceypy as spice
from ps_read_hdf_3d import ps_read_hdf_3d
# from plot_body_positions import xyz2rtp_in_Carrington

import furnsh_kernels

Rs = 696300  # km

def xyz2rtp_in_Carrington(xyz_carrington, for_psi=False):
    """
    Convert (x,y,z) to (r,t,p) in Carrington Coordination System.
        (x,y,z) follows the definition of SPP_HG in SPICE kernel.
        (r,lon,lat) is (x,y,z) converted to heliographic lon/lat, where lon \in [0,2pi], lat \in [-pi/2,pi/2] .
    :param xyz_carrington:
    :return:
    """
    r_carrington = np.linalg.norm(xyz_carrington[0:3], 2)

    lon_carrington = np.arcsin(xyz_carrington[1]/np.sqrt(xyz_carrington[0]**2+xyz_carrington[1]**2))
    if xyz_carrington[0] < 0:
        lon_carrington = np.pi - lon_carrington
    if lon_carrington < 0:
        lon_carrington += 2*np.pi

    lat_carrington = np.pi / 2 - np.arccos(xyz_carrington[2] / r_carrington)
    if for_psi:
        lat_carrington = np.pi/2-lat_carrington
    return r_carrington, lon_carrington, lat_carrington

def rtp2xyz_in_Carrington(rtp_carrington, for_psi=False):
    if for_psi:
        rtp_carrington[2] = np.pi / 2 - rtp_carrington[2]

    z_carrington = rtp_carrington[0] * np.cos(np.pi / 2 - rtp_carrington[2])
    y_carrington = rtp_carrington[0] * np.sin(np.pi / 2 - rtp_carrington[2]) * np.sin(rtp_carrington[1])
    x_carrington = rtp_carrington[0] * np.sin(np.pi / 2 - rtp_carrington[2]) * np.cos(rtp_carrington[1])
    return x_carrington, y_carrington, z_carrington


def obtain_Brtp_and_grid_coordinate(cr, region='corona', for_trace=True):
    if region == 'corona':
        dict_corona = ps_read_hdf_3d(cr, 'corona', 'br002', 3)
        br_p_corona = dict_corona['scales3']
        br_r_corona = dict_corona['scales1']
        Br_corona = dict_corona['datas']

        dict_corona = ps_read_hdf_3d(cr, 'corona', 'bt002', 3)
        bt_r_corona = dict_corona['scales1']
        Bt_corona = dict_corona['datas']

        dict_corona = ps_read_hdf_3d(cr, 'corona', 'bp002', 3)
        bp_t_corona = dict_corona['scales2']
        Bp_corona = dict_corona['datas']

        if for_trace:
            rg_corona = br_r_corona[1:]
        else:
            rg_corona = bt_r_corona
        tg_corona = bp_t_corona
        pg_corona = br_p_corona

        if for_trace:
            br_g_corona = (Br_corona[:, :, 1:] + Br_corona[:, :, 1:]) / 2.

        else:
            br_g_corona = (Br_corona[:, :, 1:] + Br_corona[:, :, :-1]) / 2.
        bt_g_corona = (Bt_corona[:, 1:, :] + Bt_corona[:, :-1, :]) / 2.
        bp_g_corona = (Bp_corona[1:, :, :] + Bp_corona[:-1, :, :]) / 2.

        return rg_corona, tg_corona, pg_corona, br_g_corona, bt_g_corona, bp_g_corona
    elif region == 'helio':

        dict_helio = ps_read_hdf_3d(cr, 'helio', 'br002', 3)
        br_p_helio = dict_helio['scales3']
        br_r_helio = dict_helio['scales1']
        Br_helio = dict_helio['datas']

        dict_helio = ps_read_hdf_3d(cr, 'helio', 'bt002', 3)
        bt_r_helio = dict_helio['scales1']
        Bt_helio = dict_helio['datas']

        dict_helio = ps_read_hdf_3d(cr, 'helio', 'bp002', 3)
        bp_t_helio = dict_helio['scales2']
        Bp_helio = dict_helio['datas']

        if for_trace:
            rg_helio = br_r_helio[:-1]
        else:
            rg_helio = bt_r_helio
        tg_helio = bp_t_helio
        pg_helio = br_p_helio

        if for_trace:
            br_g_helio = (Br_helio[:, :, :-1] + Br_helio[:, :, :-1]) / 2.
        else:
            br_g_helio = (Br_helio[:, :, 1:] + Br_helio[:, :, :-1]) / 2.
        bt_g_helio = (Bt_helio[:, 1:, :] + Bt_helio[:, :-1, :]) / 2.
        bp_g_helio = (Bp_helio[1:, :, :] + Bp_helio[:-1, :, :]) / 2.
        # bp_g_helio = Bp_helio
        return rg_helio, tg_helio, pg_helio, br_g_helio, bt_g_helio, bp_g_helio
    else:
        print('Wrong Region (corona/helio)')
    return 0


def PSI_HCS(crid, region='sc'):
    hcs = []
    if region == 'sc' or region == 'scih':
        rg_corona, tg_corona, pg_corona, br_g_corona, bt_g_corona, bp_g_corona \
            = obtain_Brtp_and_grid_coordinate(crid, region='corona', for_trace=True)

        pv, tv, rv = np.meshgrid(pg_corona, tg_corona, rg_corona, indexing='ij')
        br = br_g_corona
        bt = bt_g_corona
        bp = bp_g_corona

        xv = rv * np.cos(pv) * np.sin(tv)
        yv = rv * np.sin(pv) * np.sin(tv)
        zv = rv * np.cos(tv)

        mesh = pyvista.StructuredGrid(xv, yv, zv)
        mesh.point_data['values'] = br.ravel(order='F')  # also the active scalars
        isos_br = mesh.contour(isosurfaces=1, rng=[0., 0.])
        hcs.append(isos_br)
    if region == 'ih' or region == 'scih':
        rg_helio, tg_helio, pg_helio, br_g_helio, bt_g_helio, bp_g_helio \
            = obtain_Brtp_and_grid_coordinate(crid, region='helio', for_trace=True)

        pv, tv, rv = np.meshgrid(pg_helio, tg_helio, rg_helio, indexing='ij')
        br = br_g_helio
        bt = bt_g_helio
        bp = bp_g_helio

        xv = rv * np.cos(pv) * np.sin(tv)
        yv = rv * np.sin(pv) * np.sin(tv)
        zv = rv * np.cos(tv)

        mesh = pyvista.StructuredGrid(xv, yv, zv)
        mesh.point_data['values'] = br.ravel(order='F')  # also the active scalars
        isos_br = mesh.contour(isosurfaces=1, rng=[0., 0.])
        hcs.append(isos_br)
    return hcs


def PSI_trace(start_point, crid):
    print('------PSI_trace------')
    start_r = np.linalg.norm(start_point)
    print('Starting Point: ', start_point)
    print('Starting Radius: ', start_r)
    # stream_points = []
    streams = []
    if start_r >= 30.:
        print('Tracing in Helio')
        rg_helio, tg_helio, pg_helio, br_g_helio, bt_g_helio, bp_g_helio = obtain_Brtp_and_grid_coordinate(crid,
                                                                                                           region='helio',
                                                                                                           for_trace=True)
        pv, tv, rv = np.meshgrid(pg_helio, tg_helio, rg_helio, indexing='ij')
        br = br_g_helio
        bt = bt_g_helio
        bp = bp_g_helio
        print(br.shape)
        print(bt.shape)
        print(bp.shape)

        xv = rv * np.cos(pv) * np.sin(tv)
        yv = rv * np.sin(pv) * np.sin(tv)
        zv = rv * np.cos(tv)

        Bxg = br * np.sin(tv) * np.cos(pv) + bt * np.cos(tv) * np.cos(pv) - bp * np.sin(pv)
        Byg = br * np.sin(tv) * np.sin(pv) + bt * np.cos(tv) * np.sin(pv) + bp * np.cos(pv)
        Bzg = br * np.cos(tv) - bt * np.sin(tv)

        mesh = pyvista.StructuredGrid(xv, yv, zv)
        vectors = np.empty((mesh.n_points, 3))
        vectors[:, 0] = Bxg.ravel(order='F')
        vectors[:, 1] = Byg.ravel(order='F')
        vectors[:, 2] = Bzg.ravel(order='F')
        mesh['vectors'] = vectors
        stream = mesh.streamlines('vectors', progress_bar=True, start_position=start_point, return_source=False,
                                  integration_direction='both',
                                  max_time=500., max_error=1e-2)
        # stream.tube(radius=1.).plot(show_grid=True)
        stream_rs = np.linalg.norm(stream.points, axis=1)
        min_r_ind = np.nanargmin(abs(stream_rs))
        mid_point = stream.points[min_r_ind, :]
        print('Reach Mid Point: ', mid_point)
        print('Mid Point Radius: ', np.linalg.norm(mid_point))
        if np.linalg.norm(mid_point) > 30.419035:
            print('=TRICK=')
            mid_point = mid_point*30.4/np.linalg.norm(mid_point)
        start_point = mid_point
        streams.append(stream)

    print('Tracing in Corona')
    rg_corona, tg_corona, pg_corona, br_g_corona, bt_g_corona, bp_g_corona = obtain_Brtp_and_grid_coordinate(crid,
                                                                                                             region='corona',
                                                                                                             for_trace=True)

    print('Corona outer boudary: ', rg_corona[-1])

    pv, tv, rv = np.meshgrid(pg_corona, tg_corona, rg_corona, indexing='ij')
    br = br_g_corona
    bt = bt_g_corona
    bp = bp_g_corona

    xv = rv * np.cos(pv) * np.sin(tv)
    yv = rv * np.sin(pv) * np.sin(tv)
    zv = rv * np.cos(tv)

    Bxg = br * np.sin(tv) * np.cos(pv) + bt * np.cos(tv) * np.cos(pv) - bp * np.sin(pv)
    Byg = br * np.sin(tv) * np.sin(pv) + bt * np.cos(tv) * np.sin(pv) + bp * np.cos(pv)
    Bzg = br * np.cos(tv) - bt * np.sin(tv)

    mesh = pyvista.StructuredGrid(xv, yv, zv)
    vectors = np.empty((mesh.n_points, 3))
    vectors[:, 0] = Bxg.ravel(order='F')
    vectors[:, 1] = Byg.ravel(order='F')
    vectors[:, 2] = Bzg.ravel(order='F')
    mesh['vectors'] = vectors
    stream = mesh.streamlines('vectors', progress_bar=True, start_position=start_point, return_source=False,
                              integration_direction='both',
                              max_time=100.)
    # stream.tube(radius=0.1).plot(show_grid=True)
    if stream.n_points > 0:
        stream_rs = np.linalg.norm(stream.points, axis=1)
        min_r_ind = np.nanargmin(abs(stream_rs))
        end_point = stream.points[min_r_ind, :]
        print('Reach End Point: ', end_point)
        print('End Point Radius: ', np.linalg.norm(end_point))
    else:
        print('Trace fail in corona')
    streams.append(stream)
    return streams


def PSI_trace_overview(crid, region='sc', sc_np=100, sc_src=30., ih_np=100, ih_src=200.):
    p = pyvista.Plotter()
    if region == 'ih' or region == 'scih':
        print('Helio')
        rg_helio, tg_helio, pg_helio, br_g_helio, bt_g_helio, bp_g_helio = obtain_Brtp_and_grid_coordinate(crid,
                                                                                                           region='helio',
                                                                                                           for_trace=True)
        pv, tv, rv = np.meshgrid(pg_helio, tg_helio, rg_helio, indexing='ij')
        br = br_g_helio
        bt = bt_g_helio
        bp = bp_g_helio

        xv = rv * np.cos(pv) * np.sin(tv)
        yv = rv * np.sin(pv) * np.sin(tv)
        zv = rv * np.cos(tv)

        Bxg = br * np.sin(tv) * np.cos(pv) + bt * np.cos(tv) * np.cos(pv) - bp * np.sin(pv)
        Byg = br * np.sin(tv) * np.sin(pv) + bt * np.cos(tv) * np.sin(pv) + bp * np.cos(pv)
        Bzg = br * np.cos(tv) - bt * np.sin(tv)

        mesh = pyvista.StructuredGrid(xv, yv, zv)
        vectors = np.empty((mesh.n_points, 3))
        vectors[:, 0] = Bxg.ravel(order='F')
        vectors[:, 1] = Byg.ravel(order='F')
        vectors[:, 2] = Bzg.ravel(order='F')
        mesh['vectors'] = vectors
        stream, src = mesh.streamlines('vectors', return_source=True, source_radius=ih_src, n_points=ih_np,
                                       progress_bar=True,
                                       max_time=400., max_error=1e-20)
        p.add_mesh(stream.tube(radius=1), color='white')

        mesh.point_data['values'] = br.ravel(order='F')  # also the active scalars
        isos_br = mesh.contour(isosurfaces=1, rng=[0., 0.])
        # isos_br.plot(show_grid=True,opacity=0.5)
        p.add_mesh(isos_br, opacity=0.5)
    if region == 'sc' or region == 'scih':
        print('Tracing in Corona')
        rg_corona, tg_corona, pg_corona, br_g_corona, bt_g_corona, bp_g_corona = obtain_Brtp_and_grid_coordinate(crid,
                                                                                                                 region='corona',
                                                                                                                 for_trace=True)

        pv, tv, rv = np.meshgrid(pg_corona, tg_corona, rg_corona, indexing='ij')
        br = br_g_corona
        bt = bt_g_corona
        bp = bp_g_corona

        xv = rv * np.cos(pv) * np.sin(tv)
        yv = rv * np.sin(pv) * np.sin(tv)
        zv = rv * np.cos(tv)

        Bxg = br * np.sin(tv) * np.cos(pv) + bt * np.cos(tv) * np.cos(pv) - bp * np.sin(pv)
        Byg = br * np.sin(tv) * np.sin(pv) + bt * np.cos(tv) * np.sin(pv) + bp * np.cos(pv)
        Bzg = br * np.cos(tv) - bt * np.sin(tv)

        mesh = pyvista.StructuredGrid(xv, yv, zv)
        vectors = np.empty((mesh.n_points, 3))
        vectors[:, 0] = Bxg.ravel(order='F')
        vectors[:, 1] = Byg.ravel(order='F')
        vectors[:, 2] = Bzg.ravel(order='F')
        mesh['vectors'] = vectors
        stream, src = mesh.streamlines('vectors', return_source=True, source_radius=sc_src, n_points=sc_np,
                                       progress_bar=True,
                                       max_time=50.)
        p.add_mesh(stream.tube(radius=0.05), color='white')
        mesh.point_data['values'] = br.ravel(order='F')  # also the active scalars
        isos_br = mesh.contour(isosurfaces=1, rng=[0., 0.])
        p.add_mesh(isos_br, opacity=0.5)

    p.add_mesh(pyvista.Sphere(1))
    p.show_grid()
    # p.show()

    return p


def lines_from_points(points):
    """Given an array of points, make a line set"""
    poly = pyvista.PolyData()
    poly.points = points
    cells = np.full((len(points) - 1, 3), 2, dtype=np.int_)
    cells[:, 1] = np.arange(0, len(points) - 1, dtype=np.int_)
    cells[:, 2] = np.arange(1, len(points), dtype=np.int_)
    poly.lines = cells
    return poly


def trace_psp_in_PSI(crid,datetime_beg,datetime_end,timestep):

    steps = (datetime_end - datetime_beg) // timestep + 1
    epoch = np.array([x * timestep + datetime_beg for x in range(steps)])
    ets = spice.datetime2et(epoch)

    psp_pos, _ = spice.spkpos('SPP', ets, 'IAU_SUN', 'NONE', 'SUN')  # km
    psp_pos = psp_pos.T / Rs

    line = lines_from_points(psp_pos.T)
    footpoints = psp_pos*0.
    p = pyvista.Plotter()
    for i in range(len(ets)):
        print(epoch[i])
        pos = psp_pos[:, i]
        streams = PSI_trace(pos, crid)
        if len(streams) == 1:
            p.add_mesh(streams[0].tube(radius=0.1), color='white')
            stream_rs = np.linalg.norm(streams[0].points, axis=1)
            end_r_ind = np.nanargmin(abs(stream_rs))
            end_point = streams[0].points[end_r_ind, :]
            footpoints[:,i] = end_point
        elif len(streams) > 1:
            p.add_mesh(streams[0].tube(radius=1.), color='white')
            if streams[1].n_points > 0:
                p.add_mesh(streams[1].tube(radius=0.1), color='white')
                stream_rs = np.linalg.norm(streams[1].points, axis=1)
                end_r_ind = np.nanargmin(abs(stream_rs))
                end_point = streams[1].points[end_r_ind, :]
                footpoints[:,i] = end_point
            else:
                footpoints[:,i] = [np.nan,np.nan,np.nan]
    hcs = PSI_HCS(crid, region='scih')
    p.add_mesh(hcs[0], opacity=0.5,color='purple')
    p.add_mesh(hcs[1], opacity=0.5,color='purple')
    p.add_mesh(pyvista.Sphere(radius=1.))
    p.add_mesh(line.tube(radius=2.), color='blue')
    p.show_grid()
    # p.show()

    psp_pos_rlonlat = np.array([xyz2rtp_in_Carrington(psp_pos[:, i]) for i in range(len(psp_pos.T))])
    psp_pos_rlonlat = psp_pos_rlonlat.T

    footpoints_rlonlat = np.array([xyz2rtp_in_Carrington(footpoints[:, i]) for i in range(len(footpoints.T))])
    footpoints_rlonlat = footpoints_rlonlat.T

    timestr_beg = datetime_beg.strftime('%Y%m%dT%H%M%S')
    timestr_end = datetime_end.strftime('%Y%m%dT%H%M%S')

    df = pd.DataFrame()
    df['epoch'] = epoch
    df['psp_x'] = psp_pos[0]
    df['psp_y'] = psp_pos[1]
    df['psp_z'] = psp_pos[2]
    df['psp_r'] = psp_pos_rlonlat[0]
    df['psp_lon'] = np.rad2deg(psp_pos_rlonlat[1])
    df['psp_lat'] = np.rad2deg(psp_pos_rlonlat[2])
    df['fpts_x'] = footpoints[0]
    df['fpts_y'] = footpoints[1]
    df['fpts_z'] = footpoints[2]
    df['fpts_r'] = footpoints_rlonlat[0]
    df['fpts_lon'] = np.rad2deg(footpoints_rlonlat[1])
    df['fpts_lat'] = np.rad2deg(footpoints_rlonlat[2])
    df.to_csv('export/Trace_PSI_data/psi_trace/trace_psp_fpts_('+timestr_beg+'-'+timestr_end+'-'+str(timestep//timedelta(minutes=1))+'min).csv')




    return p

if __name__ == '__main__':
    # p = PSI_trace_overview(2254,region='scih',sc_src=30.,sc_np=360,ih_src=200.,ih_np=0)
    # p.show()
    # quit()

    datetime_beg = datetime(2022,9,15,0,0,0)
    datetime_end = datetime(2022,10,11,0,0,0)
    timestep = timedelta(days=1)
    crid = 2262
    #
    p=trace_psp_in_PSI(crid,datetime_beg,datetime_end,timestep)
    # p.export_vtkjs('export/sample')
    # p.export_html('export/sample')

    p.show()

