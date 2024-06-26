"""
Name:
    ps_read_hdf_3d.
Purpose:
    Reads a 3-dimensional HDF file.
Inputs:
    crid: ID of the Carrington Rotation
    region: 'corona' or 'helio'
    param: e.g. 'br002' (https://www.predsci.com/mhdweb/data_access.php)
Keywords:
    PERIODIC_DIM:
        (Heritage from idl routine)Dimension number to be made periodic (0 - 2PI). Valid values are:
            1 - The first dimension.
            2 - The second dimension.
            3 - The third dimension.
Outputs:
    Dict containing:
        scales1:
            Scale values of the first dimension.
        scales2:
            Scale values of the second dimension.
        scales3:
            Scale values of the third dimension.
        datas:
            Data values.

Name:
    ps_load_hdf.
Purpose:
    Load HDF files from PSI site
Input:
    crid:
    Range: 'corona' or 'helio'
    param: e.g. 'br002'
Output:
    HDF files downloaded
"""
import re

from matplotlib import pyplot as plt
from pyhdf.HDF import *
from pyhdf.SD import *
import numpy as np
import pandas as pd
import glob
import requests
import os, sys
from tqdm import tqdm


def ps_read_hdf_3d(crid, region, param, periodicDim=[],psi_datapath='/Users/ephe/PSI_Data/'):
    path = psi_datapath+'cr' + str(crid) + '/' + region + '/'
    # path = '../PSI_Data/cr' + str(crid) + '/' + region + '/'
    filename = param + '.hdf'
    # print(path + filename)
    # print(os.path.exists(path + filename))
    # print(path)
    if not os.path.exists(path + filename):
        ps_load_hdf(crid, region, param)

    sdId = SD(path + filename)
    sdkeys = sdId.datasets().keys()

    scales3 = sdId.select('fakeDim0').get()
    scales2 = sdId.select('fakeDim1').get()
    scales1 = sdId.select('fakeDim2').get()
    datas = sdId.select('Data-Set-2').get()

    sdId.end()

    scales3 = np.array(scales3)
    scales2 = np.array(scales2)
    scales1 = np.array(scales1)
    datas = np.array(datas)

    # Periodic is assumed to have a range from 0 to 2PI
    if type(periodicDim) == int:
        nx = len(scales3)
        ny = len(scales2)
        nz = len(scales1)
        tmpDatas = datas

        if periodicDim == 1:
            scales1 = np.append(scales1, scales1[0] + 2.0 * np.pi)
            datas = np.zeros((nx, ny, nz + 1))
            # print(datas.shape)
            datas[:, :, 0:-1] = tmpDatas
            datas[:, :, -1] = tmpDatas[:, :, 0]
        elif periodicDim == 2:
            scales2 = np.append(scales2, scales2[0] + 2.0 * np.pi)
            datas = np.zeros((nx, ny + 1, nz))
            datas[:, 0:- 1, :] = tmpDatas
            datas[:, -1, :] = tmpDatas[:, 0, :]
        elif periodicDim == 3:
            scales3 = np.append(scales3, scales3[0] + 2.0 * np.pi)
            datas = np.zeros((nx + 1, ny, nz))
            datas[0:-1, :, :] = tmpDatas
            datas[-1, :, :] = tmpDatas[0, :, :]

    dict = {'scales1': scales1, 'scales2': scales2, 'scales3': scales3, 'datas': datas}
    # print('Shape of Scales:',scales1.shape,scales2.shape,scales3.shape)
    # print('Shape of Datas',datas.shape)

    return dict


def ps_load_hdf(crid, region, param,psi_datapath='/Users/ephe/PSI_Data/'):
    # https://www.predsci.com/data/runs/cr2246-high/hmi_mast_mas_std_0201/corona/
    url = 'https://www.predsci.com/data/runs/'
    url += 'cr' + str(crid) + '-'
    # path = '../PSI_Data/cr' + str(crid)
    path = psi_datapath+'cr' + str(crid)
    try:
        os.mkdir(path)
    except:
        pass

    if str(requests.get(url + 'high/', verify=False)) != '<Response [404]>':
        url += 'high/'
    else:
        print('No High Data, try Medium')
        url += 'medium/'
        if str(requests.get(url, verify=False)) == '<Response [404]>':
            print('Nothing Found!')
            return

    r_links = requests.get(url,verify=False)
    res_url =  r"(?<=href=\").+?(?=\")|(?<=href=\').+?(?=\')"
    links = re.findall(res_url,str(r_links.content),re.I|re.S|re.M)
    l = [s for s in links if 'mas' in s]
    for i in range(len(l)):
        print(i,l[i])
    # s=int(input('choose one file: (e.g. 0)'))
    print('Downloading [1]')
    s=1

    url += l[s] + region + '/'
    path += '/' + region + '/'
    filename = param + '.hdf'
    r = requests.get(url + filename, verify=False, stream=True)
    total = int(r.headers.get('content-length', 0))
    if str(r) != '<Response [404]>':
        try:
            os.mkdir(path)
        except:
            pass
        f = open(path + filename, 'wb')
        bar = tqdm(desc=path + filename, total=total, unit='iB', unit_scale=True, unit_divisor=1024, )
        for data in r.iter_content(chunk_size=1024):
            size = f.write(data)
            bar.update(size)
        f.close


if __name__ == '__main__':
    dict = ps_read_hdf_3d(2239, 'corona', 'br002')
    br_r = np.array(dict['scales1'])
    br_t = np.array(dict['scales2'])
    br_p = np.array(dict['scales3'])
    br = np.array(dict['datas'])
    print('---Br(r,t,p)---')
    print(br_r.shape,br_t.shape,br_p.shape)
    print(br.shape)
    print(br_r[-11:-2])
    print((br_r[-11:-2]+br_r[-10:-1])/2.)
    dict = ps_read_hdf_3d(2239, 'corona', 'bt002')
    bt_r = np.array(dict['scales1'])
    bt_t = np.array(dict['scales2'])
    bt_p = np.array(dict['scales3'])
    bt = np.array(dict['datas'])
    print('---Bt(r,t,p)---')
    print(bt_r.shape,bt_t.shape,bt_p.shape)
    print(bt.shape)
    print(bt_r[-10:-1])
    dict = ps_read_hdf_3d(2239, 'corona', 'bp002')
    bp_r = np.array(dict['scales1'])
    bp_t = np.array(dict['scales2'])
    bp_p = np.array(dict['scales3'])
    bp = np.array(dict['datas'])
    print('---Bp(r,t,p)---')
    print(bp_r.shape,bp_t.shape,bp_p.shape)
    print(bp.shape)
    print(bp_r[:10])
    plt.subplot(3,2,1)
    plt.title('rc(n+-1)-rg(n)')
    plt.plot(br_r[1:]-bt_r)
    plt.plot(br_r[:-1]-bt_r)
    plt.subplot(3,2,2)
    plt.title('[rc(n-1)+rc(n+1)]/2-rg(n)')
    plt.plot((br_r[1:]+br_r[:-1])/2.-bt_r)
    plt.subplot(3,2,3)
    plt.title('tc(n+-1)-tg(n)')
    plt.plot(bt_t[1:]-bp_t)
    plt.plot(bt_t[:-1]-bp_t)
    plt.subplot(3,2,4)
    plt.title('[tc(n-1)+tc(n+1)]/2-tg(n)')
    plt.plot((bt_t[1:]+bt_t[:-1])/2.-bp_t)
    plt.subplot(3,2,5)
    plt.title('pc(n+-1)-pg(n)')
    plt.plot(bp_p[1:]-br_p)
    plt.plot(bp_p[:-1]-br_p)
    plt.subplot(3,2,6)
    plt.title('[pc(n-1)+pc(n+1)]/2-pg(n)')
    plt.plot((bp_p[1:]+bp_p[:-1])/2.-br_p)
    plt.show()

    rg = bt_r
    tg = bp_t
    pg = br_p

    br_g = (br[:,:,1:]+br[:,:,:-1])/2.
    bt_g = (bt[:,1:,:]+bt[:,:-1,:])/2.
    bp_g = (bp[1:,:,:]+bp[:-1,:,:])/2.
    # print(dict['scales1'])
