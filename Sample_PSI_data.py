import os

import numpy as np
import pandas as pd
import sunpy.coordinates.sun as sun
from datetime import datetime,timedelta

import spiceypy as spice

from utils import appendSpherical_np
from ps_read_hdf_3d import ps_read_hdf_3d
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

def sample_in_psi(epoch,position,param_lst,lat2psilat=True,interp=False,save2csv='./export/Sample_PSI_data/'):
    os.makedirs(save2csv,exist_ok=True)
    sample={param:np.zeros_like(epoch) for param in param_lst}
    if lat2psilat:
        lat_sample = np.pi/2-position['lat']
    r_sample = position['r']
    lon_sample = position['lon']

    crid_beg = int(sun.carrington_rotation_number(epoch[0]))
    crid_end = int(sun.carrington_rotation_number(epoch[-1]))

    for crid in range(crid_beg,crid_end+1):
        subepochbin = (epoch>=sun.carrington_rotation_time(crid)) & (epoch<sun.carrington_rotation_time(crid+1))
        subepochbin_scbin = subepochbin & (r_sample<=30)

        subepochbin_ihbin = subepochbin & (r_sample>30)

        if sum(subepochbin_ihbin) > 0:
            for param in param_lst:
                data_ih = ps_read_hdf_3d(crid, 'helio', param+'002', periodicDim=3)
                r = np.array(data_ih['scales1'])
                t = np.array(data_ih['scales2'])
                p = np.array(data_ih['scales3'])
                datas = np.array(data_ih['datas'])
                datas = datas * unit_lst[param]


                r_inds = list(map(lambda x: np.argmin(abs(x-r)),np.array(r_sample)))
                t_inds = list(map(lambda x: np.argmin(abs(x-t)),np.array(lat_sample)))
                p_inds = list(map(lambda x: np.argmin(abs(x-p)),np.array(lon_sample)))

                sample[param][subepochbin_ihbin] = datas[p_inds,t_inds,r_inds][subepochbin_ihbin]
        if sum(subepochbin_scbin)>0:
            for param in param_lst:
                data_sc = ps_read_hdf_3d(crid, 'corona', param+'002', periodicDim=3)
                r = np.array(data_sc['scales1'])  # 201 in Rs, distance from sun
                t = np.array(data_sc['scales2'])  # 150 in rad, latitude
                p = np.array(data_sc['scales3'])  # 256 in rad, Carrington longitude
                datas = np.array(data_sc['datas'])  # 1CU = 2.205G = 2.205e-4T = 2.205e5nT
                datas = datas * unit_lst[param]  #

                r_inds = list(map(lambda x: np.argmin(abs(x - r)), r_sample))
                t_inds = list(map(lambda x: np.argmin(abs(x - t)), lat_sample))
                p_inds = list(map(lambda x: np.argmin(abs(x - p)), lon_sample))

                sample[param][subepochbin_scbin] = datas[p_inds, t_inds, r_inds][subepochbin_scbin]
    sample_df = pd.DataFrame(sample)
    sample_df['epoch'] = epoch
    sample_df['x'] = position['x']
    sample_df['y']  = position['y']
    sample_df['epoch'] = epoch
    sample_df['z'] = position['z']
    sample_df['r'] = position['r']
    sample_df['lon'] = position['lon']
    sample_df['lat'] = position['lat']

    if save2csv:
        sample_df.to_csv(save2csv+'export.csv')
    return sample_df

def sample_in_psi_psp(epoch,param_lst,lat2psilat=True,save2csv='export/Sample_PSI_data/'):

    et_epoch = spice.datetime2et(epoch)
    pos, _ = spice.spkpos('SPP', et_epoch, 'IAU_SUN', 'NONE', 'SUN')
    pos = np.array(pos.T/ Rs2km)
    pos = appendSpherical_np(pos.T)
    print(pos)
    position = {'x': pos[:,0], 'y': pos[:,1], 'z': pos[:,2],'r':pos[:,3],'lat':pos[:,4],'lon':pos[:,5]}
    return sample_in_psi(epoch,position,param_lst,lat2psilat=True,save2csv='export/Sample_PSI_data/')

if __name__ == '__main__':
    from utils import create_epoch
    epoch = create_epoch([datetime(2020,1,21),datetime(2020,1,31)],timedelta(minutes=60))
    sample_df = sample_in_psi_psp(epoch,['br'])
    sample_df.plot()
