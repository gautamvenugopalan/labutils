#import cdsutils as cds
import scipy.io as sio
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interp
import sys
import yaml
from matplotlib.ticker import FormatStrFormatter
from scipy.spatial import ConvexHull

if 'gvELOG' in plt.style.available:
	plt.style.use('gvELOG')
else:
	plt.style.use('bmh')

def main():
	if len(sys.argv) < 2:
		print(plotHeatMap.__doc__)
	else:
		paramFile = str(sys.argv[1])
		
		with open(paramFile,'r') as f:
			reader = yaml.load_all(f)
			params = reader.next()
			reader.close()

		if params['downloadData']:
			chans = params['chans']
			tStart = params['tStart']
			tStop = params['tStop']
			CDSdata = cds.getdata(CHANS,tData,tStart)
			data = CDSdata[0].data
			for i in range(len(chans)-1):
				data = np.vstack(data,CDSdata[i+1].data)
			data = data.T
		else:
			dataFile = params['dataFileName']
			data = np.loadtxt(dataFile)
		
		if np.shape(data)[1]>7:
			z, points, hull, PIT, YAW, grid_x, grid_y, pitRange, yawRange = shapeData(data, params) 
			plotHeatMap(z, points, hull, PIT, YAW, grid_x, grid_y, pitRange, yawRange, params)
		else:
			data = data[:,1:7]
			plotHeatMap(data,params)
	return

def plotHeatMap(z, points, hull, PIT, YAW, grid_x, grid_y, pitRange, yawRange, params):
	'''
	Script to generate heat maps for CP walking.
	Heat maps are generated using 2D interpolation from scipy.interpolate
	Example usage:
		python CPheatMap.py params.yml

	where params.yml is a paramter file specifying channels, plotting parameters etc.

	Option to download data from CDS on execution, or use an existing data file
	
	Supports interpolation using scipy.interpolate.Rbf or scipy.interpolate.griddata.
	Judicious use of either function is required to guarantee physical results.

	GV Feb 2018
	'''
	fig, ax = plt.subplots(2,3,figsize=(16,10), sharex=True, sharey=True)
	BLRMSrange = params['legends']
	nChans = len(BLRMSrange)
	for kk in range(nChans):
		print('Processing {}'.format(BLRMSrange[kk]))
		blrms = z[kk,:-1]
		pts = np.vstack((PIT,YAW)).T
		if params['interpFunc'] == 'griddata':
			grid_z = interp.griddata(pts,blrms,(grid_x,grid_y),method=params['interpMethod'],rescale=True)
		elif params['interpFunc'] == 'Rbf':
			smoo = interp.Rbf(PIT, YAW, blrms,function=params['interpMethod'],smooth=5)
			grid_z = smoo(grid_x, grid_y)
		else:
			print('Interpolation function not recognized - use "Rbf" or "griddata".')
			return
		heatmap = ax[int(np.heaviside(kk-2,0))][kk%3].imshow(grid_z,extent=(pitRange[0],pitRange[1],yawRange[0],yawRange[1]),cmap='magma',aspect=0.3, rasterized=True)#,norm=LogNorm())
		c1=fig.colorbar(heatmap, ax=ax[int(np.heaviside(kk-2,0))][kk%3])
		c1.set_label('RLP OUT16',fontsize=16, labelpad=6)
		ax[int(np.heaviside(kk-2,0))][kk%3].axes.set_aspect('equal','datalim')
		ax[int(np.heaviside(kk-2,0))][kk%3].set_title(BLRMSrange[kk],fontsize=14)
		fig.subplots_adjust(hspace=0.1, wspace=0.33)
		ax[1][2].axis('off')
		if params['plotPoints']:
			ax[int(np.heaviside(kk-2,0))][kk%3].plot(PIT,YAW,'ko',markersize=8,fillstyle='none')
		if params['plotHull']:
			for simplex in hull.simplices:
				ax[int(np.heaviside(kk-2,0))][kk%3].plot(points[simplex, 0], points[simplex, 1], 'r--')
	ax[1,0].set_xlabel(params['optic'] + ' PIT slider offset [urad]',fontsize=14)
	ax[1,1].set_xlabel(params['optic'] + ' PIT slider offset [urad]',fontsize=14)
	ax[0,0].set_ylabel(params['optic'] + ' YAW slider offset [urad]',fontsize=14)
	ax[1,0].set_ylabel(params['optic'] + ' YAW slider offset [urad]',fontsize=14)
	plt.suptitle(params['optic'] + ' CP walk interpolated')

	#Some formatting
	for axx in ax.flat:
		axx.yaxis.set_major_formatter(FormatStrFormatter('%3d'))
		axx.xaxis.set_major_formatter(FormatStrFormatter('%3d'))
		axx.grid('on',which='both',linestyle='--',alpha=0.4)

	if params['saveFig']:
		fig.savefig(params['optic'] + '_CP_BLRMS_interpolated.pdf')

	return

def findUniques(dataIn):
	'''
	Subroutine to find unique instances of PIT and YAW offsets to the CP in the timespan of interest
	'''
	uP,ind = np.unique(dataIn,return_inverse=True)
	ind = np.nonzero(np.diff(ind))
	#Append first and last indices before returning
	ind = np.append([0],ind)
	ind = np.append(ind,[np.shape(dataIn)[0]])
	return ind

def normalize(arr):
	maxx = np.max(arr)
	minn = np.min(arr)
	arrNorm = np.ones(len(arr))
	for ii, aa in enumerate(arr):
		arrNorm[ii] = (aa - minn) / (maxx - minn)
	return arrNorm

def shapeData(data, params):
	'''
	Subroutine to format input data into the format expected by plotHeatMap()
	'''
	if np.shape(data)[1]>7:
		if params['optic'] == 'ITMX':
			data = np.vstack((data[:,1],data[:,2],data[:,5],data[:,6],data[:,7],data[:,8],data[:,9])).T
		elif params['optic'] == 'ITMY':
			data = np.vstack((data[:,3],data[:,4],data[:,5],data[:,6],data[:,7],data[:,8],data[:,9])).T
	else:
		data = data[:,1:7]
	
	#Correctly formatted data should be 7 columns - PIT, YAW, and 5 BLRMS.
	#Some data processing - arrange everything into nice arrays
	Pindex = findUniques(data[:,0])
	Yindex = findUniques(data[:,1])

	#make some arrays for storing reprocessed data
	nPts = np.shape(Pindex)[0] - 1
	nChans = np.shape(params['chans'])[0] - 2

	z = np.zeros([nChans, np.shape(Pindex)[0]-1])
	PIT = np.zeros([nChans, np.shape(Pindex)[0]-1])
	YAW = np.zeros([nChans, np.shape(Yindex)[0]])
	#YAW = np.zeros([nChans, np.shape(Yindex)[0]-1])
	
	#Seems like some DQ check was required to avoid negative BLRMS values, not sure if this is still needed
	for kk in range(nChans):
		for jj, ii in enumerate(Pindex[:-1]):
			if(np.min(data[ii:Pindex[jj+1], kk+2]) and data[ii,0]!=0 and data[ii,1]!=0) > 0:
				z[kk, jj] = np.median(data[ii:Pindex[jj+1], kk+2]) #Median averaging for glitch robustness
				PIT[kk, jj] = data[ii,0]
				YAW[kk, jj] = data[ii,1]
	
	#Now define the grid of points over which to interpolate
	pitRange = [params['pitMin'],params['pitMax']]
	yawRange = [params['yawMin'],params['yawMax']]
	nPts_interpx = params['interpPts'] * 1j
	grid_x,grid_y = np.mgrid[pitRange[0]:pitRange[1]:nPts_interpx, yawRange[0]:yawRange[1]:nPts_interpx]

	#Since all the BLRMS have the same steps, we can simplify the PIT and YAW arrays
	PIT = PIT[0,:-1].T
	YAW = YAW[0,:-1].T
	
	points = np.vstack((PIT,YAW)).T
	hull = ConvexHull(points)
	return z, points, hull, PIT, YAW, grid_x, grid_y, pitRange, yawRange

if __name__ == "__main__":
	main()
