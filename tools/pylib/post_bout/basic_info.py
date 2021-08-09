from __future__ import print_function
from __future__ import division
from builtins import zip
from builtins import range
from past.utils import old_div

# basic_info return some statistical averages and harmonic info
import numpy as np
import math


def basic_info(data, meta, rescale=True, rotate=False, user_peak=0, nonlinear=None):

    print("in basic_info")
    # from . import read_grid,parse_inp,read_inp,show

    dims = data.shape
    ndims = len(dims)

    mxg = meta["MXG"]["v"]

    if ndims == 4:
        nt, nx, ny, nz = data.shape
        print(nt, nx, ny)
    else:
        print("something with dimesions")

    dc = (
        data.mean(1).mean(1).mean(1)
    )  # there MUST be a way to indicate all axis at once
    amp = abs(data).max(1).max(1).max(1)
    dt = meta["dt"]["v"]

    if rescale:
        amp_o = amp - dc
        fourDamp = np.repeat(amp_o, nx * ny * nz)
        fourDamp = fourDamp.reshape(nt, nx, ny, nz)
        dc_n = old_div(dc, amp_o)
        data_n = old_div(data, fourDamp)

        print(data.shape)
        dfdt = np.gradient(data)[0]
        dfdt = abs(dfdt).max(1).max(1).max(1)

        ave = {"amp": amp, "dc": dc, "amp_o": amp_o, "dfdt": dfdt}

    else:
        print("no rescaling")
        ave = {"amp": amp, "dc": dc}

    if nonlinear is not None:  # add nonlinear part if user provides
        nl = abs(nonlinear[:, mxg : -1.0 * mxg, :, :]).max(1).max(1).max(1)
        nl_norm = (old_div(nl, dfdt)) * dt

        ave["nl"] = nl
        ave["nl_norm"] = nl_norm

    if rotate:
        print("rotate stuff")
        # will need to provide some grid geometry to do this one
    else:
        print("or not")

    # let's identify the dominant modes, for now at every [t,x] slice
    # if the data set is too large we can average over x
    peaks_db = fft_info(
        data, user_peak, meta=meta
    )  # Nt X Nx X (# of loc. max) list of dict

    # print peaks[0]['gamma']
    return peaks_db, ave


def fft_info(
    data,
    user_peak,
    dimension=[3, 4],
    rescale=False,
    wavelet=False,
    show=False,
    meta=0,
    edgefix=False,
):
    import numpy as np
    import math

    print("in fft_inf0")

    dims = data.shape
    ndims = len(dims)

    if ndims == 4:
        nt, nx, ny, nz = data.shape
        print(data.shape)
    else:
        print("something with dimesions")

    # data2 = data
    # if edgefix:
    #     data2 = np.zeros((nt,nx,ny+1,nz+1))

    #     for t in range(nt):
    #         for x in range(nx):
    #             temp = np.append(data[t,x,:,:],[data[t,x,0,:]],0)
    #             data2[t,x,:,:] = np.append(temp,
    #                                        np.transpose([temp[:,0]]),1)

    # dt, k labels for the revelant dimensions

    dt = meta["dt"]["v"]
    dz = meta["dz"]
    # IC = meta['IC']
    ky_max = old_div(ny, 2)
    kz_max = old_div(nz, 2)
    amp = abs(data).max(2).max(2)  # nt x nx

    print("dt: ", dt)

    # print data[0,2,:,:]

    IC = amp[0, :].max()  # intial condition, set
    print(IC)

    fft_data = np.fft.fft2(data)[
        :, :, 0:ky_max, 0:kz_max
    ]  # by default the last 2 dimensions

    power = fft_data.conj() * fft_data
    # print power[0].max(), (IC*(ky_max)*(kz_max))**2

    cross_pow = old_div((fft_data * (np.roll(fft_data, 1, axis=0)).conj()), (ny * nz))

    if rescale:
        fft_data_n = np.fft.fft2(data_n)[:, :, 0:ky_max, 0:kz_max]
        pow_n = np.sqrt((fft_data_n.conj() * fft_data_n).real)

    peaks = [[[] for i in range(nx)] for j in range(nt)]  # a list of dictionaries
    peaks_db = []

    peak_hist = [[0 for i in range(kz_max)] for j in range(ky_max)]  # a 2d bin array

    # for now using a lame 2x loop method

    if user_peak != 0:
        for mem in user_peak:
            print(mem)
            peak_hist[int(mem[0])][int(mem[1])] = abs(
                power.mean(0).mean(0)[int(mem[0]), int(mem[1])]
            )

        # floor =  ((IC*(kz_max*ky_max))**2)/10000

    else:
        for t in range(nt):
            for x in range(nx):
                peaks[t][x] = local_maxima(
                    power[t, x, :, :], 0, floor=(IC * (kz_max * ky_max)) ** 2
                )
                for p in peaks[t][
                    x
                ]:  # looping over each returned peakset at some fixed t,x pair
                    peak_hist[p["y_i"]][
                        p["z_i"]
                    ] += 1  # average across t and x, at least exclude pad
        floor = 0

    # this array is usefull for determining what the dominant modes are
    # but we want to retain the option of observing how the amplitude
    # of any give harmonic varies in space

    peak_hist = np.array(peak_hist)

    # let's find the top N overall powerfull harmnonics
    net_peak = local_maxima(peak_hist, user_peak, bug=False)

    print("net_peak: ", net_peak, user_peak != 0)
    # dom_mode = [{'amp':[],'amp_n':[],'phase':[],'freq':[],'gamma':[]} for x in net_peak]
    dom_mode_db = []

    Bp = meta["Bpxy"]["v"][:, old_div(ny, 2)]
    B = meta["Bxy"]["v"][:, old_div(ny, 2)]
    Bt = meta["Btxy"]["v"][:, old_div(ny, 2)]

    rho_s = meta["rho_s"]["v"]

    L_z = old_div(meta["L_z"], rho_s)
    # L_z =
    L_y = meta["lpar"]  # already normalized earlier in read_inp.py
    L_norm = old_div(meta["lbNorm"], rho_s)

    hthe0_n = 1e2 * meta["hthe0"]["v"] / rho_s  # no x dep
    hthe0_n_x = old_div(L_y, (2 * np.pi))  # no x dep

    print("L_z,Ly: ", L_z, L_y)

    # if user provides the harmo    nic info overide the found peaks

    # thi is where all the good stuff is picked up

    # look at each mode annd pull out some usefull linear measures
    for i, p in enumerate(net_peak):
        # print i,p['y_i'],p['z_i'],fft_data.shape,fft_data[:,:,p['y_i'],p['z_i']].shape

        amp = (
            old_div(np.sqrt(power[:, :, p["y_i"], p["z_i"]]), (kz_max * ky_max))
        ).real

        # print (np.angle(fft_data[:,:,p['y_i'],p['z_i']],deg=False)).real

        phase = -np.array(
            np.gradient(
                np.squeeze(np.angle(fft_data[:, :, p["y_i"], p["z_i"]], deg=False))
            )[0].real
        )  # nt x nx

        gamma_instant = np.array(np.gradient(np.log(np.squeeze(amp)))[0])

        # loop over radaii
        phasenew = []
        gammanew = []
        from scipy.interpolate import interp2d, interp1d
        from scipy import interp

        gamma_t = np.transpose(gamma_instant)
        for i, phase_r in enumerate(np.transpose(phase)):
            gamma_r = gamma_t[i]
            jumps = np.where(abs(phase_r) > old_div(np.pi, 32))
            # print jumps
            if len(jumps[0]) != 0:

                all_pts = np.array(list(range(0, nt)))
                good_pts = (np.where(abs(phase_r) < old_div(np.pi, 3)))[0]
                # print good_pts,good_pts
                # f = interp1d(good_pts,phase_r[good_pts],fill_value=.001)
                # print max(all_pts), max(good_pts)
                # phasenew.append(f(all_pts))
                try:
                    phase_r = interp(all_pts, good_pts, phase_r[good_pts])
                    gamma_r = interp(all_pts, good_pts, gamma_r[good_pts])
                except:
                    "no phase smoothing"
            phasenew.append(phase_r)
            gammanew.append(gamma_r)

        phase = old_div(np.transpose(phasenew), dt)

        gamma_i = old_div(np.transpose(gammanew), dt)

        amp_n = (
            old_div(np.sqrt(power[:, :, p["y_i"], p["z_i"]]), (kz_max * ky_max * amp))
        ).real
        # amp_n = dom_mode[i]['amp_n'] #nt x nx

        # let just look over the nx range
        # lnamp = np.log(amp[nt/2:,2:-2])
        try:
            lnamp = np.log(amp[old_div(nt, 2) :, :])
        except:
            print("some log(0) stuff in basic_info")

        t = dt * np.array(list(range(nt)))  # dt matters obviouslyww
        r = np.polyfit(t[old_div(nt, 2) :], lnamp, 1, full=True)

        gamma_est = r[0][0]  # nx
        f0 = np.exp(r[0][1])  # nx
        res = r[1]
        pad = [0, 0]
        # gamma_est = np.concatenate([pad,gamma_est,pad])
        # f0 = np.concatenate([pad,f0,pad])
        # res = np.concatenate([pad,res,pad])

        # sig = res/np.sqrt((x['nt']-2))
        sig = np.sqrt(old_div(res, (nt - 2)))
        # sig0 = sig*np.sqrt(1/(x['nt'])+ ) # who cares
        sig1 = sig * np.sqrt(old_div(1.0, (nt * t.var())))
        nt = np.array(nt)
        print("shapes ", nt.shape, nt, lnamp.shape, res.shape, gamma_est)
        # print r
        res = 1 - old_div(res, (nt * lnamp.var(0)))  # nx
        res[0:2] = 0
        res[-2:] = 0

        gamma = [gamma_est, sig1, f0, res]

        # gamma_est2 = np.gradient(amp)[0]/(amp[:,:]*dt)
        # gamma_w =  np.gradient(gamma_est2)[0]

        # gamma_i = np.abs(gamma_w).argmin(0) #index of the minimum for any given run
        # for j in range(nx):
        #     gamma_w[0:max([gamma_i[j],nt/3]),j] = np.average(gamma_w)*100000.0

        freq = np.array(
            weighted_avg_and_std(phase[-10:, :], weights=np.ones(phase[-10:, :].shape))
        )

        # gamma = weighted_avg_and_std(
        #         gamma_est2[-5:,:],weights=np.ones(gamma_est2[-5:,:].shape))

        k = [
            [p["y_i"], p["z_i"]],
            [2 * math.pi * float(p["y_i"]) / L_y, 2 * math.pi * p["z_i"] / L_z],
        ]
        # L_y is normalized

        # simple k def, works in drift-instability fine
        # k = [[p['y_i'],p['z_i']],
        #      [(B/Bp)**-1*2*math.pi*float(p['y_i'])/(L_y),(B/Bp)*2*math.pi*p['z_i']/L_z]]

        # k_r = [[p['y_i'],p['z_i']],
        #      [(Bp/B)*2*math.pi*float(p['y_i'])/(L_y),
        #       (B/Bp)*2*math.pi*p['z_i']/L_z]]

        k_r = [
            [p["y_i"], p["z_i"]],
            [
                2 * math.pi * float(p["y_i"]) / (L_y),
                (old_div(B, Bp)) * 2 * math.pi * p["z_i"] / L_z
                + (old_div(Bt, B)) * 2 * math.pi * float(p["y_i"]) / (L_y),
            ],
        ]

        # revised
        # k_r = [[p['y_i'],p['z_i']],
        #      [2*math.pi*float(p['y_i'])/(L_y),
        #       (Bp/B)*2*math.pi*p['z_i']/L_z
        #       - (Bt/Bp)*2*math.pi*float(p['y_i'])/(L_y)]]
        # revised

        # what I think is the most general one, works in drift-instability again
        # seems to work for Bz only helimak, now trying Bp = Bt
        # k = [[p['y_i'],p['z_i']],
        #      [((Bp/B)*float(p['y_i'])/(hthe0_n)) +
        #       2*np.pi*p['z_i']*np.sqrt(1-(Bp/B)**2)/L_z,
        #       2*math.pi*p['z_i']/L_norm -
        #       float(p['y_i'])*np.sqrt(1-(Bp/B)**2)/(hthe0_n)]]
        # k = [[p['y_i'],p['z_i']],
        #      [((Bp/B)*float(p['y_i'])/(hthe0_n)),
        #       2*math.pi*p['z_i']/L_norm]]
        # BOTH SEEM TO PRODOCE SAME RESULTS

        # k = [[p['y_i'],p['z_i']],
        #      [(float(p['y_i'])/(hthe0_n_x)),
        #       2*math.pi*float(p['z_i'])/L_norm]]

        dom_mode_db.append(
            {
                "modeid": i,
                "k": k[1],
                "gamma": gamma,
                "freq": freq,
                "amp": amp,
                "amp_n": amp_n,
                "phase": phase,
                "mn": k[0],
                "nt": nt,
                "k_r": k_r[1],
                "gamma_i": gamma_i,
            }
        )

    return dom_mode_db


# return a 2d array fof boolean values, a very simple boolian filter
def local_maxima(array2d, user_peak, index=False, count=4, floor=0, bug=False):

    from operator import itemgetter, attrgetter

    if user_peak == 0:
        where = (
            (array2d >= np.roll(array2d, 1, 0))
            & (array2d >= np.roll(array2d, -1, 0))
            & (array2d >= np.roll(array2d, 0, 1))
            & (array2d >= np.roll(array2d, 0, -1))
            & (array2d >= old_div(array2d.max(), 5.0))
            & (array2d > floor * np.ones(array2d.shape))
            & (array2d >= array2d.mean())
        )
    else:  # some simpler filter if user indicated some modes
        where = array2d > floor

    # ignore the lesser local maxima, throw out anything with a ZERO
    if bug == True:
        print(array2d, array2d[where.nonzero()], where.nonzero()[0])

    peaks = list(zip(where.nonzero()[0], where.nonzero()[1], array2d[where.nonzero()]))

    peaks = sorted(peaks, key=itemgetter(2), reverse=True)

    if len(peaks) > count and user_peak == 0:
        peaks = peaks[0:count]

    keys = ["y_i", "z_i", "amp"]

    peaks = [dict(list(zip(keys, peaks[x]))) for x in range(len(peaks))]

    return peaks
    # return np.array(peak_dic)


def weighted_avg_and_std(values, weights):
    """
    Returns the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """

    if len(values.shape) == 2:
        average = np.average(values, 0)  # , weights=weights)
        variance = old_div(
            (
                np.inner(
                    weights.transpose(), ((values - average) ** 2).transpose()
                ).diagonal()
            ),
            weights.sum(0),
        )
    else:
        average = np.average(values, weights=weights)
        variance = old_div(
            np.dot(weights, (values - average) ** 2), weights.sum()
        )  # Fast and numerically precise

    return [average, variance]
