'''
Read in stuff, do some calculations, and make plots for comparison
between Matt Rayson's Bay model and the shelf model.
'''

from mpl_toolkits.basemap import pyproj
import tracpy
import init
import netCDF4 as netCDF
from matplotlib.mlab import *
import projects
from matplotlib.tri import Triangulation

def run():

    # Location of TXLA model output
    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'

    # Matt's
    d = netCDF.MFDataset('tracks/matt/2010-05-*.nc', aggdim='ntrac')
    p = pyproj.Proj(proj='utm', zone='15')
    lonm,latm = p(d.variables['xp'][:], d.variables['yp'][:], inverse=True)
    # want one simulation to have a single set of the starting points
    d1 = netCDF.MFDataset('tracks/matt/2010-05-23T00.nc', aggdim='ntrac')
    lonmstart, latmstart = p(d1.variables['xp'][:,0], d1.variables['yp'][:,0], inverse=True)

    # Kristen's
    db = netCDF.MFDataset('tracks/galv_b/2010-05-*.nc',aggdim='ntrac')
    lonb = np.fliplr(db.variables['lonp'][:]) # flip to be forward in time
    latb = np.fliplr(db.variables['latp'][:]) # flip to be forward in time
    df = netCDF.MFDataset('tracks/galv_f/2010-05-*.nc',aggdim='ntrac')
    lonf = df.variables['lonp'][:]
    latf = df.variables['latp'][:]
    # Combined
    lonk = np.concatenate((lonb,lonf),axis=1)
    latk = np.concatenate((latb,latf),axis=1)
    grid = tracpy.inout.readgrid(loc)

    ## Shelf tracks ##
    # Eliminate drifters in the bay and outside Matt's domain
    # lonvert = np.array([-94.44896376, -94.38868859, -94.38625664, -94.35380894,
    #        -94.35357294, -94.32114683, -94.32099761, -94.2839534 ,
    #        -94.28371841, -94.32073196, -94.32015033, -94.35247582,
    #        -94.35185901, -94.38873039, -94.3886522 , -94.41861373,
    #        -94.41854639, -94.45541454, -94.45530469, -94.52433946,
    #        -94.52422355, -94.59322436, -94.59308138, -94.65976418,
    #        -94.65956549, -94.76770875, -94.76763167, -94.96778311,
    #        -94.96789969, -95.07604071, -95.07642981, -95.17770396,
    #        -95.17810623, -95.2795571 , -95.2798657 , -95.31901876,
    #        -95.3199739 , -95.14404646, -95.14434603, -95.0427632 ,
    #        -95.04318363, -95.00854345, -95.00884843, -94.9349127 ,
    #        -94.935271  , -94.90060651, -94.90089718, -94.83385365,
    #        -94.83415664, -94.76476871, -94.76501578, -94.72799532,
    #        -94.7285101 , -94.69609675, -94.69634997, -94.62918694,
    #        -94.62940326, -94.59001667, -94.59019564, -94.45345669, -94.44896376])
    # latvert = np.array([ 29.51080435,  29.51099607,  29.48073304,  29.48082425,
    #         29.41424186,  29.41432502,  29.3679228 ,  29.36800802,
    #         29.28529434,  29.28520919,  29.10363677,  29.10355388,
    #         28.92804192,  28.92793764,  28.90776056,  28.90766817,
    #         28.89152842,  28.89140532,  28.86720038,  28.85483934,
    #         28.8326488 ,  28.81822849,  28.7940243 ,  28.78361653,
    #         28.75336181,  28.75277136,  28.74268359,  28.74135471,
    #         28.75345628,  28.75261073,  28.78892442,  28.78805142,
    #         28.82234232,  28.83349224,  28.85769558,  28.85730671,
    #         28.92992664,  29.06675435,  29.0929739 ,  29.09382103,
    #         29.13417612,  29.134447  ,  29.16470469,  29.16525227,
    #         29.20358522,  29.20382761,  29.23610029,  29.23654313,
    #         29.27286272,  29.273285  ,  29.30556746,  29.30577777,
    #         29.37638645,  29.37656202,  29.41287395,  29.41321234,
    #         29.44751927,  29.44770176,  29.47796222,  29.47850429,  29.51080435])
    # # xvert, yvert = grid['basemap'](lonvert, latvert)
    # # These are the x,y locations for the normal projection, but want to save in lon/lat
    # # xvert = np.array([422614, 428447, 428671, 431812, 431812, 434953, 434953, 438543,
    # #                 438543, 434953, 434953, 431812, 431812, 428223, 428223, 425306,
    # #                 425306, 421716, 421716, 414986, 414986, 408255, 408255, 401749,
    # #                 401749, 391204, 391204, 371685, 371685, 361140, 361140, 351268,
    # #                 351268, 341397, 341397, 337583, 337583, 354858, 354858, 364730,
    # #                 364730, 368095, 368095, 375275, 375275, 378640, 378640, 385146,
    # #                 385146, 391877, 391877, 395467, 395467, 398608, 398608, 405114,
    # #                 405114, 408928, 408928, 422165, 422614])
    # # yvert = np.array([772846, 772846, 769480, 769480, 762076, 762076, 756916, 756916,
    # #                 747718, 747718, 727526, 727526, 708007, 708007, 705763, 705763,
    # #                 703968, 703968, 701276, 699930, 697462, 695891, 693199, 692077,
    # #                 688712, 688712, 687590, 687590, 688936, 688936, 692975, 692975,
    # #                 696789, 698135, 700827, 700827, 708904, 723936, 726852, 726852,
    # #                 731340, 731340, 734705, 734705, 738968, 738968, 742557, 742557,
    # #                 746596, 746596, 750186, 750186, 758038, 758038, 762076, 762076,
    # #                 765891, 765891, 769256, 769256, 772846])
    # # xyverts = np.vstack((xvert,yvert))
    # verts = np.vstack((lonvert,latvert))
    # # Form path
    # path = Path(verts.T)
    # # Loop through the drifters as represented by their initial location for lonf,latf
    # # Establish a set of indices of the drifters that are outside the desired domain
    # ind = np.zeros(lonf.shape[0])
    # for i in xrange(lonf.shape[0]):
    #     if path.contains_point(np.vstack((lonf[i,0],latf[i,0]))):
    #         ind[i] = 1
    # # These only contain tracks that go through the starting points of Matt's runs
    # lonknew = lonk[ind.astype(bool),:]
    # latknew = latk[ind.astype(bool),:]

    # Do plot
    # Find initial points (since this is a backward run) to plot as starting points
    # moving forward
    lonki = []
    latki = []
    for idrift in xrange(lonk.shape[0]):
        # pdb.set_trace()
        # print idrift
        # Find last non-nan and make sure it is in the desired month start time
        ind3 = ~np.isnan(lonk[idrift,:])
        #pdb.set_trace()
        # only plot if last non-nan (back in time) is in 1 month period
        # in order to plot the tracks that "started" in the plotted month
        if np.sum(ind3) > 1: # don't want tracks that start on land
            # This is for if we care when the drifter stopped
            # if t[find(ind3)[-1]] >= datetime(year,startMonth,startDay,0) and \
            #   t[find(ind3)[-1]] <= datetime(year,startMonth+1,startDay,0):
            # ind2 = ~np.isnan(lonk[idrift,:])
            if np.sum(np.isnan(lonk[idrift,:])) > 0 and np.sum(np.isnan(lonk[idrift,:])) < lonk.shape[1]: # if there is a nan
                # ax.plot(lonk[idrift,find(ind2)[-1]].T,latk[idrift,find(ind2)[-1]].T,'o',color='orange',liidth=.5,label='_nolegend_')
                lonki.append(lonk[idrift,find(ind3)[0]])
                latki.append(latk[idrift,find(ind3)[0]])
            else:
                # ax.plot(lonk[idrift,0].T,latk[idrift,0].T,'o',color='orange',liidth=.5,label='_nolegend_')
                lonki.append(lonk[idrift,find(ind3)[0]])
                latki.append(latk[idrift,find(ind3)[0]])
    xki, yki = grid['basemap'](lonki,latki)

    # Overall shelf tracks
    # Backward to forward
    tracpy.plotting.tracks(lonk, latk, 'matt/shelf', grid=grid)
    plot(xki,yki,'go',alpha=.4)
    savefig('figures/matt/shelftracks.png',bbox_inches='tight')
    # Backward to Bay opening only
    tracpy.plotting.tracks(lonb, latb, 'matt/shelf_backonly', grid=grid)
    plot(xki,yki,'go',alpha=.4)
    savefig('figures/matt/shelf_backonlytracks.png',bbox_inches='tight')
    # Bay forward only
    tracpy.plotting.tracks(lonf, latf, 'matt/shelf_forwardonly', grid=grid)
    plot(lonf[:,0],latf[:,0],'go',alpha=.4)
    savefig('figures/matt/shelf_forwardonlytracks.png',bbox_inches='tight')

    # Overall shelf transport
    projects.transport.run(name='galv_b', Title='Transport to Galveston', dmax=1.5)
    projects.transport.run(name='galv_f', Title='Transport from Galveston', dmax=1.5)


    ## Galveston histograms ##
    # Smaller basemap parameters.
    llcrnrlon=-95.6; llcrnrlat=28.3; urcrnrlon=-93.8; urcrnrlat=29.8;
    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
    smallgrid = tracpy.inout.readgrid(loc, llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, 
                                    urcrnrlat=urcrnrlat, urcrnrlon=urcrnrlon)

    # Matt histogram
    tracpy.plotting.hist(lonm, latm, 'matt/matt', grid=smallgrid, which='hexbin', bins=(160,160))
    xmstart, ymstart = smallgrid['basemap'](lonmstart,latmstart)
    f = gcf()
    # last axis was for colorbar so grab inital one
    f.axes[0].plot(xmstart,ymstart,'go',alpha=.3) 
    # add grid
    points = loadtxt('points.dat')
    lonp,latp = p(points[:,0], points[:,1], inverse=True)
    triang = Triangulation(lonp,latp)
    f.axes[0].triplot(triang)
    savefig('figures/matt/matthexbin.png',bbox_inches='tight')

    # Matt tracks
    tracpy.plotting.tracks(lonm, latm, 'matt/matt', grid=smallgrid)

    # Shelf:
    # Need to stop drifters as his domain boundary BUT not from going into the Bay
    # x and y verts on smallgrid projection, but this projection might change over time
    # xvertbig = np.array([ 113263.66720431,  119096.66652423,  119320.66686604,
    #         122461.66678215,  122461.66685854,  125602.66643451,
    #         125602.66642388,  129192.66658589,  129192.66626941,
    #         125602.66643373,  125602.66631186,  122461.66719913,
    #         122461.66712714,  118872.66642492,  118872.66647544,
    #         115955.66697976,  115955.6668092 ,  112365.66628673,
    #         112365.66670544,  105635.66662346,  105635.66692937,
    #          98904.66691465,   98904.66683249,   92398.66693207,
    #          92398.66669606,   81853.66708999,   81853.66694612,
    #          62334.66648085,   62334.66659896,   51789.66706805,
    #          51789.66695047,   41917.66679772,   41917.66677504,
    #          32046.66715647,   32046.66652647,   28232.66691345,
    #          28232.66655248,   28232.66655248,   97180.        ,
    #         113263.66720431])
    # yvertbig = np.array([ 133653.2476132 ,  133653.24740824,  130287.24684561,
    #         130287.24661835,  122883.24689633,  122883.24719082,
    #         117723.24701227,  117723.247473  ,  108525.24747157,
    #         108525.24679182,   88333.24664775,   88333.24703864,
    #          68814.24678444,   68814.2469479 ,   66570.24736787,
    #          66570.24692043,   64775.24764515,   64775.2469847 ,
    #          62083.24742793,   60737.24763316,   58269.24688878,
    #          56698.24684586,   54006.24684154,   52884.24719558,
    #          49519.24754775,   49519.24728893,   48397.24693971,
    #          48397.24676136,   49743.24761656,   49743.24760865,
    #          53782.24766385,   53782.24756101,   57596.24732823,
    #          58942.24741941,   61634.24691447,   61634.24696074,
    #          69711.24711498,  136328.        ,  169068.        ,
    #         133653.2476132 ])
    # lonvertbig, latvertbig = smallgrid['basemap'](xvertbig, yvertbig, inverse=True)
    lonvertbig = np.array([-94.44896376, -94.38868859, -94.38625664, -94.35380894,
           -94.35357294, -94.32114683, -94.32099761, -94.2839534 ,
           -94.28371841, -94.32073196, -94.32015033, -94.35247582,
           -94.35185901, -94.38873039, -94.3886522 , -94.41861373,
           -94.41854639, -94.45541454, -94.45530469, -94.52433946,
           -94.52422355, -94.59322436, -94.59308138, -94.65976418,
           -94.65956549, -94.76770875, -94.76763167, -94.96778311,
           -94.96789969, -95.07604071, -95.07642981, -95.17770396,
           -95.17810623, -95.2795571 , -95.2798657 , -95.31901876,
           -95.3199739 , -95.32790465, -94.61713345, -94.44896376])
    latvertbig = np.array([ 29.51080435,  29.51099607,  29.48073304,  29.48082425,
            29.41424186,  29.41432502,  29.3679228 ,  29.36800802,
            29.28529434,  29.28520919,  29.10363677,  29.10355388,
            28.92804192,  28.92793764,  28.90776056,  28.90766817,
            28.89152842,  28.89140532,  28.86720038,  28.85483934,
            28.8326488 ,  28.81822849,  28.7940243 ,  28.78361653,
            28.75336181,  28.75277136,  28.74268359,  28.74135471,
            28.75345628,  28.75261073,  28.78892442,  28.78805142,
            28.82234232,  28.83349224,  28.85769558,  28.85730671,
            28.92992664,  29.52892904,  29.82861514,  29.51080435])
    vertsbig = np.vstack((lonvertbig,latvertbig))
    # Form path
    pathbig = Path(vertsbig.T)

    # # Only use drifters starting where Matt started them
    # lonfnew = lonf[ind.astype(bool),:]
    # latfnew = latf[ind.astype(bool),:]
    # update lonfnew to stop drifters that exit matt's domain (but they can enter bay)
    for i in xrange(lonf.shape[0]):
        for t in xrange(lonf.shape[1]):
            if not pathbig.contains_point(np.vstack((lonf[i,t],latf[i,t]))):
                lonf[i,t:] = np.nan
                latf[i,t:] = np.nan
                break

    # Shelf histogram
    tracpy.plotting.hist(lonf, latf, 'matt/kristen', grid=smallgrid, which='hexbin', bins=(160,160))
    xf, yf = smallgrid['basemap'](lonf,latf)
    f = gcf()
    # plot grid
    xr = np.ma.masked_where(smallgrid['mask']==0,smallgrid['xr'])
    yr = np.ma.masked_where(smallgrid['mask']==0,smallgrid['yr'])
    f.axes[0].plot(xr,yr,xr.T,yr.T,color='lightgrey',zorder=0)
    # last axis was for colorbar so grab inital one
    f.axes[0].plot(xf[0:384,0],yf[0:384,0],'go',alpha=.3) # just plot them once
    savefig('figures/matt/kristenhexbin.png',bbox_inches='tight')

    # Shelf tracks
    tracpy.plotting.tracks(lonf, latf, 'matt/kristen', grid=smallgrid)

if __name__ == "__main__":
    run()
