import os
import requests
import numpy as np
# from datetime import datetime, timedelta
import sunpy.map

# Rs2km = 696300.
# AU2km = 1.49597871e8
# Rs2AU = Rs2km / AU2km
# AU2Rs = AU2km / Rs2km
# km2Rs = 1. / Rs2km
# km2AU = 1. / AU2km

# TODO: 单位标注

def dphidr(r, phi_at_r, Vsw_at_r):
    period_sunrot = 27. * (24. * 60. * 60)  # unit: s
    omega_sunrot = 2 * np.pi / period_sunrot
    result = omega_sunrot / Vsw_at_r  # unit: rad/km
    return result


def parker_spiral(r_vect_au, lat_beg_deg, lon_beg_deg, Vsw_r_vect_kmps):
    from_au_to_km = 1.49597871e8  # unit: km
    from_deg_to_rad = np.pi / 180.
    from_rs_to_km = 6.96e5
    from_au_to_rs = from_au_to_km / from_rs_to_km
    r_vect_km = r_vect_au * from_au_to_km
    num_steps = len(r_vect_km) - 1
    phi_r_vect = np.zeros(num_steps + 1)
    for i_step in range(0, num_steps):
        if i_step == 0:
            phi_at_r_current = lon_beg_deg * from_deg_to_rad  # unit: rad
            phi_r_vect[0] = phi_at_r_current
        else:
            phi_at_r_current = phi_at_r_next
        r_current = r_vect_km[i_step]
        r_next = r_vect_km[i_step + 1]
        r_mid = (r_current + r_next) / 2
        dr = r_current - r_next
        Vsw_at_r_current = Vsw_r_vect_kmps[i_step - 1]
        Vsw_at_r_next = Vsw_r_vect_kmps[i_step]
        Vsw_at_r_mid = (Vsw_at_r_current + Vsw_at_r_next) / 2
        k1 = dr * dphidr(r_current, phi_at_r_current, Vsw_at_r_current)
        k2 = dr * dphidr(r_current + 0.5 * dr, phi_at_r_current + 0.5 * k1, Vsw_at_r_mid)
        k3 = dr * dphidr(r_current + 0.5 * dr, phi_at_r_current + 0.5 * k2, Vsw_at_r_mid)
        k4 = dr * dphidr(r_current + 1.0 * dr, phi_at_r_current + 1.0 * k3, Vsw_at_r_next)
        phi_at_r_next = phi_at_r_current + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        phi_r_vect[i_step + 1] = phi_at_r_next
    lon_r_vect_deg = phi_r_vect / from_deg_to_rad  # from [rad] to [degree]
    lat_r_vect_deg = np.zeros(num_steps + 1) + lat_beg_deg  # unit: [degree]
    # r_footpoint_on_SourceSurface_rs = r_vect_au[-1] * from_au_to_rs
    # lon_footpoint_on_SourceSurface_deg = lon_r_vect_deg[-1]
    # lat_footpoint_on_SourceSurface_deg = lat_r_vect_deg[-1]
    return lon_r_vect_deg, lat_r_vect_deg


def appendSpherical_np(xyz):
    ptsnew = np.hstack((xyz, np.zeros(xyz.shape)))
    xy = xyz[:, 0] ** 2 + xyz[:, 1] ** 2
    ptsnew[:, 3] = np.sqrt(xy + xyz[:, 2] ** 2)
    # ptsnew[:,4] = np.arctan2(np.sqrt(xy), xyz[:,2]) # for elevation angle defined from Z-axis down
    ptsnew[:, 4] = np.arctan2(xyz[:, 2], np.sqrt(xy))  # for elevation angle defined from XY-plane up
    ptsnew[:, 5] = np.arctan2(xyz[:, 1], xyz[:, 0])
    return ptsnew


def get_sphere(R0, lon, lat, for_PSI=False):
    if for_PSI:
        lat = np.pi / 2 - lat
    llon, llat = np.meshgrid(lon, lat)
    x0 = R0 * np.cos(llon) * np.cos(llat)
    y0 = R0 * np.sin(llon) * np.cos(llat)
    z0 = R0 * np.sin(llat)
    return x0, y0, z0


def xyz2rtp_in_Carrington(xyz_carrington, for_psi=False):
    """
    Convert (x,y,z) to (r,p,t) in Carrington Coordination System.
        (x,y,z) follows the definition of SPP_HG in SPICE kernel.
        (r,lon,lat) is (x,y,z) converted to heliographic lon/lat, where lon \in [0,2pi], lat \in [-pi/2,pi/2] .
    :param xyz_carrington:
    :return:
    """
    r_carrington = np.linalg.norm(xyz_carrington[0:3], 2)

    lon_carrington = np.arcsin(xyz_carrington[1] / np.sqrt(xyz_carrington[0] ** 2 + xyz_carrington[1] ** 2))
    if np.shape(xyz_carrington[0]) != ():
        lon_carrington[xyz_carrington[0] < 0] = np.pi - lon_carrington[xyz_carrington[0] < 0]
    else:
        if xyz_carrington[0] < 0:
            lon_carrington = np.pi - lon_carrington
    # if lon_carrington < 0:
    lon_carrington = lon_carrington % (2 * np.pi)

    lat_carrington = np.pi / 2 - np.arccos(xyz_carrington[2] / r_carrington)
    if for_psi:
        lat_carrington = np.pi / 2 - lat_carrington
    return r_carrington, lon_carrington, lat_carrington


def rlonlat2cart(r, lon_rad, lat_rad):
    x = r * np.cos(lat_rad) * np.cos(lon_rad)
    y = r * np.cos(lat_rad) * np.sin(lon_rad)
    z = r * np.sin(lat_rad)
    return x, y, z


AIA_193_dir = 'AIA_193_map/'


def load_AIA_data(crid):
    aia_dir = AIA_193_dir
    aia_filename = 'aia193_synmap_cr' + str(crid) + '.fits'
    aia_file = os.path.join(aia_dir, aia_filename)
    if not os.path.exists(aia_file):
        aia_url = 'https://sun.njit.edu/coronal_holes/data/' + aia_filename
        aia_r = requests.get(aia_url)
        with open(aia_file, 'wb') as f:
            f.write(aia_r.content)
            f.close()
    syn_map = sunpy.map.Map(aia_file)
    map_data = syn_map.data
    map_cmap = syn_map.cmap
    map_x = np.linspace(syn_map.bottom_left_coord[0].value, syn_map.top_right_coord[0].value,
                        int(syn_map.dimensions.x.value))
    map_y = np.linspace(syn_map.bottom_left_coord[1].value, syn_map.top_right_coord[1].value,
                        int(syn_map.dimensions.y.value))
    lon = np.deg2rad(map_x)
    lat = np.deg2rad(map_y)
    return lon, lat, map_data, map_cmap


def create_epoch(range_dt, step_td):
    beg_dt = range_dt[0]
    end_dt = range_dt[1]
    return [beg_dt + n * step_td for n in range((end_dt - beg_dt) // step_td)]


if __name__ == '__main__':
    print('Hello World')
