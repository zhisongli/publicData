'''
exec(open('singleSlopeRunup.py').read())
'''
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
from scipy.interpolate import LinearNDInterpolator

plt.rc("font", family = "Times New Roman", size = 10.0)
# plt.rcParams['text.usetex'] = True
# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "Times New Roman"
# })


class myClassSolitaryWaveInfo:
    def __init__(self, waveHeight=0.1, waterDepth=1.0):
        self.waterDepth = waterDepth
        self.waveHeight = waveHeight

class myInterfaceHeight(myClassSolitaryWaveInfo):
    def __init__(self, inDir, inFile, waveHeight=0.1, waterDepth=1.0):
        myClassSolitaryWaveInfo.__init__(self,waveHeight=waveHeight, waterDepth=waterDepth)
        if os.path.isdir(inDir):
            pass
        else:
            inDir = inDir.replace('postProcessing', 'myPostProcessing')

        thisFile = os.path.join(inDir, inFile)
        # print(thisFile)
        tmpData = np.loadtxt(
                        open(thisFile, "r"),
                        # delimiter=",",
                        # skiprows=1,
                        # usecols=[0, 2]
                    )
        index = -1
        index += 1; self.time = tmpData[:, index]
        index += 1; self.eta  = tmpData[:, index+1::2]
        self.nSta = np.shape(self.eta)[1]
        self.effectiveH = np.zeros((self.nSta, 2), dtype=float)

        if '00' in inDir:
            pass
        elif '01' in inDir:
            tmpStaGroupIndex = 0
        elif '12' in inDir:
            tmpStaGroupIndex = 1
        elif '23' in inDir:
            tmpStaGroupIndex = 2
        elif '34' in inDir:
            tmpStaGroupIndex = 3
        elif '45' in inDir:
            tmpStaGroupIndex = 4
        
        for i in range(self.nSta):
            self.effectiveH[i, 0] = tmpStaGroupIndex * 10 + i + 0.01
            self.effectiveH[i, 1] = np.max(self.eta[:, i])
        
        self.meanH = np.mean(self.effectiveH[1:-1,1])
        
class dataPlot:
    def __init__(self, inRoot, waterDepth=1.0, waveHeight=0.5, inLabel='none'):
        self.case = inLabel

        inFile = 'height.dat'
        dirIndex = -1

        dirIndex += 1
        inDir = f'{inRoot}/postProcessing/interfaceHeight{dirIndex}{dirIndex+1}/0/'
        self.tmpData01 = myInterfaceHeight(inDir, inFile, waterDepth=waterDepth, waveHeight=waveHeight)

        dirIndex += 1
        inDir = f'{inRoot}/postProcessing/interfaceHeight{dirIndex}{dirIndex+1}/0/'
        self.tmpData12 = myInterfaceHeight(inDir, inFile, waterDepth=waterDepth, waveHeight=waveHeight)

        dirIndex += 1
        inDir = f'{inRoot}/postProcessing/interfaceHeight{dirIndex}{dirIndex+1}/0/'
        self.tmpData23 = myInterfaceHeight(inDir, inFile, waterDepth=waterDepth, waveHeight=waveHeight)

        dirIndex += 1
        inDir = f'{inRoot}/postProcessing/interfaceHeight{dirIndex}{dirIndex+1}/0/'
        self.tmpData34 = myInterfaceHeight(inDir, inFile, waterDepth=waterDepth, waveHeight=waveHeight)

        dirIndex += 1
        inDir = f'{inRoot}/postProcessing/interfaceHeight{dirIndex}{dirIndex+1}/0/'
        self.tmpData45 = myInterfaceHeight(inDir, inFile, waterDepth=waterDepth, waveHeight=waveHeight)
        # eH = np.mean(np.array([tmpData12.meanH, tmpData23.meanH, tmpData34.meanH])) / waveHeight
        
waterDepth = 0.5
waveHeight = 0.30 * waterDepth

inRoot = 'demo.ee=0.30.new'
staNew = dataPlot(inRoot, waterDepth=waterDepth, waveHeight=waveHeight, inLabel='new')

inRoot = 'demo.ee=0.30.newB'
staNewB = dataPlot(inRoot, waterDepth=waterDepth, waveHeight=waveHeight, inLabel='new')


inRoot = 'demo.ee=0.30.old'
staOld = dataPlot(inRoot, waterDepth=waterDepth, waveHeight=waveHeight, inLabel='old')



dirOut = './'
if True:

    FigureName = os.path.join(dirOut, 'sta')


    FigureNx, FigureNy = 1, 1
    FigureSize = (9, 4)

    FigureNx, FigureNy = 3, 1
    FigureSize = (8, 8)

    FigureNx, FigureNy = 3, 3
    FigureSize = (8, 7)

    # FigureNx, FigureNy = 1, 1
    # FigureSize = (8, 4)

    FigureNx, FigureNy = 2, 2
    FigureSize = (9, 6)

    # FigureNx, FigureNy = 3, 2
    # FigureSize = (9, 8)

    fig, axs = plt.subplots(ncols=FigureNy, nrows=FigureNx, 
                            figsize=FigureSize,
                            # layout="constrained",
                            sharex=False, sharey=False,
                            )
    
    # 
    for FigureIx in range(FigureNx):
        for FigureIy in range(FigureNy):

            FigureIndex = (FigureIx -1 + 1 ) * FigureNy + FigureIy+1
            print(f'ix={FigureIx}, iy={FigureIy}, index={FigureIndex}')             

            if FigureNx == FigureNy == 1:
                myAx = axs
            elif FigureNx == 1:
                myAx = axs[FigureIy]
            elif FigureNy == 1:
                myAx = axs[FigureIx]
            else:
                myAx = axs[FigureIx, FigureIy]
            
            myAx.grid()

            textLabelY = "$\eta$/h (-)"
            textLabelX = "Time (s)"  
            # data involking
            if FigureIndex == 0:
                pass

            elif FigureIndex in range(1, 10, 1):
                textTitle=''
                if FigureIndex == 0:
                    pass
                elif FigureIndex == 1:
                    textTitle = 'x=10'

                    tmpDD = staNew
                    tmpData = tmpDD.tmpData12
                    tmp_x = tmpData.time
                    tmp_y = tmpData.eta[:, 0] / tmpData.waterDepth
                    maxTime = tmp_x[np.argmax(tmp_y)]
                    myAx.plot(tmp_x, tmp_y, 'k-',
                            label = tmpDD.case + ', maxEta=%.3f'%(np.nanmax(tmp_y)),
                            linewidth=2.0, markersize=12)
                    
                    tmpDD = staNewB
                    tmpData = tmpDD.tmpData12
                    tmp_x = tmpData.time
                    tmp_y = tmpData.eta[:, 0] / tmpData.waterDepth
                    maxTime = tmp_x[np.argmax(tmp_y)]
                    myAx.plot(tmp_x, tmp_y, 'b-',
                            label = tmpDD.case + ', maxEta=%.3f'%(np.nanmax(tmp_y)),
                            linewidth=2.0, markersize=12)
                    
                    tmpDD = staOld
                    tmpData = tmpDD.tmpData12
                    tmp_x = tmpData.time
                    tmp_y = tmpData.eta[:, 0] / tmpData.waterDepth
                    maxTimeNew = tmp_x[np.argmax(tmp_y)]
                    tmp_x = tmp_x - maxTimeNew + maxTime
                    myAx.plot(tmp_x, tmp_y, 'r--',
                            label = tmpDD.case + ', maxEta=%.3f'%(np.nanmax(tmp_y)),
                            linewidth=2.0, markersize=12)

                if FigureIndex == 2:
                    textTitle = 'x=20'

                    tmpDD = staNew
                    tmpData = tmpDD.tmpData23
                    tmp_x = tmpData.time
                    tmp_y = tmpData.eta[:, 0] / tmpData.waterDepth
                    maxTime = tmp_x[np.argmax(tmp_y)]
                    myAx.plot(tmp_x, tmp_y, 'k-',
                            label = tmpDD.case + ', maxEta=%.3f'%(np.nanmax(tmp_y)),
                            linewidth=2.0, markersize=12)
                    
                    tmpDD = staNewB
                    tmpData = tmpDD.tmpData23
                    tmp_x = tmpData.time
                    tmp_y = tmpData.eta[:, 0] / tmpData.waterDepth
                    maxTime = tmp_x[np.argmax(tmp_y)]
                    myAx.plot(tmp_x, tmp_y, 'b-',
                            label = tmpDD.case + ', maxEta=%.3f'%(np.nanmax(tmp_y)),
                            linewidth=2.0, markersize=12)

                    tmpDD = staOld
                    tmpData = tmpDD.tmpData23
                    tmp_x = tmpData.time
                    tmp_y = tmpData.eta[:, 0] / tmpData.waterDepth
                    maxTimeNew = tmp_x[np.argmax(tmp_y)]
                    tmp_x = tmp_x - maxTimeNew + maxTime
                    myAx.plot(tmp_x, tmp_y, 'r--',
                            label = tmpDD.case + ', maxEta=%.3f'%(np.nanmax(tmp_y)),
                            linewidth=2.0, markersize=12)
                    

                if FigureIndex == 3:
                    textTitle = 'x=30'

                    tmpDD = staNew
                    tmpData = tmpDD.tmpData34
                    tmp_x = tmpData.time
                    tmp_y = tmpData.eta[:, 0] / tmpData.waterDepth
                    maxTime = tmp_x[np.argmax(tmp_y)]
                    myAx.plot(tmp_x, tmp_y, 'k-',
                            label = tmpDD.case + ', maxEta=%.3f'%(np.nanmax(tmp_y)),
                            linewidth=2.0, markersize=12)
                    
                    tmpDD = staNewB
                    tmpData = tmpDD.tmpData34
                    tmp_x = tmpData.time
                    tmp_y = tmpData.eta[:, 0] / tmpData.waterDepth
                    maxTime = tmp_x[np.argmax(tmp_y)]
                    myAx.plot(tmp_x, tmp_y, 'b-',
                            label = tmpDD.case + ', maxEta=%.3f'%(np.nanmax(tmp_y)),
                            linewidth=2.0, markersize=12)
                    
                    tmpDD = staOld
                    tmpData = tmpDD.tmpData34
                    tmp_x = tmpData.time
                    tmp_y = tmpData.eta[:, 0] / tmpData.waterDepth
                    maxTimeNew = tmp_x[np.argmax(tmp_y)]
                    tmp_x = tmp_x - maxTimeNew + maxTime
                    myAx.plot(tmp_x, tmp_y, 'r--',
                            label = tmpDD.case + ', maxEta=%.3f'%(np.nanmax(tmp_y)),
                            linewidth=2.0, markersize=12)
                    

                if FigureIndex == 4:
                    textTitle = 'x=40'

                    tmpDD = staNew
                    tmpData = tmpDD.tmpData45
                    tmp_x = tmpData.time
                    tmp_y = tmpData.eta[:, 0] / tmpData.waterDepth
                    maxTime = tmp_x[np.argmax(tmp_y)]
                    myAx.plot(tmp_x, tmp_y, 'k-',
                            label = tmpDD.case + ', maxEta=%.3f'%(np.nanmax(tmp_y)),
                            linewidth=2.0, markersize=12)

                    tmpDD = staNewB
                    tmpData = tmpDD.tmpData45
                    tmp_x = tmpData.time
                    tmp_y = tmpData.eta[:, 0] / tmpData.waterDepth
                    maxTime = tmp_x[np.argmax(tmp_y)]
                    myAx.plot(tmp_x, tmp_y, 'b-',
                            label = tmpDD.case + ', maxEta=%.3f'%(np.nanmax(tmp_y)),
                            linewidth=2.0, markersize=12)
                    
                    tmpDD = staOld
                    tmpData = tmpDD.tmpData45
                    tmp_x = tmpData.time
                    tmp_y = tmpData.eta[:, 0] / tmpData.waterDepth
                    maxTimeNew = tmp_x[np.argmax(tmp_y)]
                    tmp_x = tmp_x - maxTimeNew + maxTime
                    myAx.plot(tmp_x, tmp_y, 'r--',
                            label = tmpDD.case + ', maxEta=%.3f'%(np.nanmax(tmp_y)),
                            linewidth=2.0, markersize=12)
                    
                
                if False:
                    FigureXmin, FigureXmax = np.min(tmp_x), np.max(tmp_x)
                    FigureYmin, FigureYmax = np.min(tmp_y), np.max(tmp_y)
                    FigureTickXN = 10
                    FigureTickYN = 10
                    FigureDx = (FigureXmax - FigureXmin) / FigureTickXN
                    FigureDy = (FigureYmax - FigureYmin) / FigureTickYN
                else:
                    FigureXmin, FigureXmax = 0, 20
                    # FigureXmin, FigureXmax = 2, 18
                    # FigureXmin, FigureXmax = 0, 16
                    FigureYmin, FigureYmax = -0.2, 0.8
                    FigureYmin, FigureYmax = -0.1, 0.4
                    # FigureYmin, FigureYmax = -0.05, 0.25

                    FigureTickXN = 10
                    FigureTickYN = 10
                    FigureDx = (FigureXmax - FigureXmin) / FigureTickXN
                    FigureDy = (FigureYmax - FigureYmin) / FigureTickYN

            # end

            # * * * * + * * * * + * * * * + * * * * + * * * * + * * * * + * * * * +
            if FigureIndex >= 1:
                myAx.legend(
                    markerscale = 1.0,
                    loc = "lower center",
                    fontsize = "medium",
                    edgecolor = "black",
                    facecolor = "white",
                    ncol = 2,
                    # title = "Velocity CDF",
                )

            if FigureNx * FigureNy > 1:
                    textTitle = f'({FigureIndex}) ' + textTitle
            myAx.set_title(textTitle, loc = "left", fontsize = 'small')
            '''
                WARN: following lines can not be modified, if YOU know how to do them.
            '''
            # - - - - + - - - - + - - - - + - - - - + - - - - + - - - - + - - - - +
            myAx.set_xlim(FigureXmin, FigureXmax)
            FigureTickX = np.linspace(FigureXmin, FigureXmax, num=FigureTickXN+1)
            # print(f'FigureTickX = {FigureTickX}')
            myAx.set_xticks(FigureTickX)
            FigureTickXLabel = ['%04d'%i for i in FigureTickX] if np.mod(FigureDx,1) == 0 else ['%.1f'%i for i in FigureTickX]
            # myAx.set_xticklabels(FigureTickXLabel)
            if FigureIy == 1 - 1:
                myAx.set_ylabel(textLabelY)
            # end

            myAx.set_ylim(FigureYmin, FigureYmax)
            FigureTickY = np.linspace(FigureYmin, FigureYmax, num=FigureTickYN+1)
            myAx.set_yticks(FigureTickY)
            FigureTickYLabel = ['%04d'%i for i in FigureTickY] if np.mod(FigureDy,1) == 0 else ['%.1f'%i for i in FigureTickY]
            # myAx.set_yticklabels(FigureTickYLabel)
            if FigureIx == FigureNx-1:
                myAx.set_xlabel(textLabelX)
            # end
            # - - - - + - - - - + - - - - + - - - - + - - - - + - - - - + - - - - +
            
            
    plt.tight_layout()
    # plt.savefig(f'{FigureName}200.png', bbox_inches='tight', dpi=200)
    # plt.savefig(f'{FigureName}400.png', bbox_inches='tight', dpi=400)
    # plt.savefig(f'{FigureName}600.png', bbox_inches='tight', dpi=600)
    plt.savefig(f'{FigureName}.png', bbox_inches='tight', dpi=600)
    plt.savefig(f'{FigureName}.pdf', bbox_inches='tight')
    plt.show()
# end


# 
